import json
from functools import lru_cache
from pathlib import Path

import pytest

from softarm.config import ConfigError, load_config
from softarm.derive import derive
from softarm.pipeline import derive_system

ROOT = Path(__file__).parents[1]


@lru_cache
def _reference_system(path_text: str):
    config = load_config(Path(path_text))
    return config, derive_system(config)


def test_reference_configs_are_valid():
    for path in (ROOT / "examples" / "config").glob("*.toml"):
        config, system = _reference_system(str(path))
        assert system.plant.config == config
        assert system.plant.config.segments >= 1


def test_readme_example_index_covers_every_reference_config():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    example_index = readme.split("## 8. 示例索引", maxsplit=1)[1]
    for path in (ROOT / "examples" / "config").glob("*.toml"):
        row = next(
            line for line in example_index.splitlines()
            if f"`{path.stem}`" in line
        )
        config = load_config(path)
        assert f"`{config.dynamics.formulation}`" in row


def test_reference_bundle_manifests_match_configs():
    for path in (ROOT / "examples" / "config").glob("*.toml"):
        config, system = _reference_system(str(path))
        manifest_path = ROOT / "examples" / "generated" / path.stem / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        model = manifest["model"]
        assert model["rod"] == config.rod
        assert model["parameterization"] == config.parameterization
        assert model["segments"] == config.segments
        assert model["base_mode"] == config.base.mode
        assert model["dynamics_formulation"] == config.dynamics.formulation
        assert model["mount_xyz"] == list(config.base.mount_xyz)
        assert model["mount_rpy"] == list(config.base.mount_rpy)
        expected_parameters = list(system.plant.parameters)
        if system.actuation is not None:
            expected_parameters.extend(system.actuation.parameters)
        if system.constraint is not None:
            expected_parameters.extend(system.constraint.parameters)
        assert [item["name"] for item in manifest["parameters"]] == [
            item.name for item in expected_parameters
        ]
        assert [item["default"] for item in manifest["parameters"]] == [
            item.default for item in expected_parameters
        ]
        expected_actuation = None if system.actuation is None else {
            "family": system.actuation.family,
            "acceleration": system.actuation.acceleration,
            "channels": [
                {"name": name, "kind": kind}
                for name, kind in zip(
                    system.actuation.channel_names,
                    system.actuation.channel_kinds,
                    strict=True,
                )
            ],
        }
        assert manifest["actuation"] == expected_actuation
        expected_constraint = None if system.constraint is None else {
            "family": system.constraint.family,
            "channels": [
                {"name": name, "kind": kind}
                for name, kind in zip(
                    system.constraint.channel_names,
                    system.constraint.channel_kinds,
                    strict=True,
                )
            ],
        }
        assert manifest["constraint"] == expected_constraint


def test_base_and_constraint_validation(tmp_path):
    path = tmp_path / "floating.toml"
    path.write_text(
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[base]\nmode="floating_rpy"\nmount_xyz=[0,0,0]\nmount_rpy=[0,0,0]\n'
        '[constraint]\nfamily="plane_point_contact"\nplane_normal=[0,0,0]\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="plane_normal.*nonzero"):
        load_config(path)


def test_euler_bernoulli_requires_explicit_normalized_ritz(tmp_path):
    path = tmp_path / "bad.toml"
    path.write_text('[model]\nrod="euler_bernoulli"\nparameterization="ritz"\nsegments=1\n[dynamics]\nformulation="symbolic_lagrange"\n[ritz]\nx=[0,0,1]\ny=[0,0,2]\n')
    with pytest.raises(ConfigError, match="normalized"):
        derive(load_config(path))


@pytest.mark.parametrize("actuation,match", [
    ('family="tendon"\n', "acceleration"),
    ('family="tendon"\nacceleration="auto"\n', "strict.*none"),
    ('family="tendon"\nacceleration="none"\n', "at least one channel"),
])
def test_tendon_actuation_requires_explicit_valid_structure(tmp_path, actuation, match):
    path = tmp_path / "bad_actuation.toml"
    path.write_text(
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
        '[actuation]\n' + actuation
        + '\n[dynamics]\nformulation="symbolic_lagrange"\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match=match):
        load_config(path)


def test_tendon_config_rejects_duplicate_span(tmp_path):
    path = tmp_path / "duplicate.toml"
    path.write_text(
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[actuation]\nfamily="tendon"\nacceleration="none"\n'
        '[[actuation.channels]]\nname="t1"\nkind="unilateral"\n'
        '[[actuation.channels.spans]]\nsection=1\nradius=0.02\nangle=0\n'
        '[[actuation.channels.spans]]\nsection=1\nradius=0.03\nangle=1\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="repeats section"):
        load_config(path)


def test_registered_cosserat_pcs_validator_checks_inertia_integration(tmp_path):
    distributed = tmp_path / "distributed.toml"
    distributed.write_text(
        '[model]\nrod="cosserat"\nparameterization="pcs"\nsegments=1\ninertia="distributed"\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[integration]\nmethod="analytic"\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="requires integration.method='gauss'"):
        derive(load_config(distributed))

    order_one = tmp_path / "order_one.toml"
    order_one.write_text(
        '[model]\nrod="cosserat"\nparameterization="pcs"\nsegments=1\ninertia="distributed"\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[integration]\nmethod="gauss"\norder=1\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="order at least 2"):
        derive(load_config(order_one))

    lumped = tmp_path / "lumped.toml"
    lumped.write_text(
        '[model]\nrod="cosserat"\nparameterization="pcs"\nsegments=1\ninertia="lumped"\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        '[integration]\nmethod="gauss"\norder=2\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="not applicable"):
        derive(load_config(lumped))


def test_unregistered_combination_is_discovered_from_registry(tmp_path):
    path = tmp_path / "unsupported.toml"
    path.write_text(
        '[model]\nrod="cosserat"\nparameterization="ritz"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n',
        encoding="utf-8",
    )
    config = load_config(path)
    with pytest.raises(ValueError, match="unregistered model combination"):
        derive(config)


def test_pac_validation_and_explicit_cosserat_rejection(tmp_path):
    low_order = tmp_path / "low_order.toml"
    low_order.write_text(
        '[model]\nrod="euler_bernoulli"\nparameterization="pac"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n'
        'inertia="distributed"\n[integration]\nmethod="gauss"\norder=3\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="order at least 4"):
        derive(load_config(low_order))

    unsupported = tmp_path / "cosserat_pac.toml"
    unsupported.write_text(
        '[model]\nrod="cosserat"\nparameterization="pac"\nsegments=1\n'
        '[dynamics]\nformulation="symbolic_lagrange"\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unregistered model combination"):
        derive(load_config(unsupported))


def test_legacy_model_family_is_rejected(tmp_path):
    path = tmp_path / "legacy.toml"
    path.write_text('[model]\nfamily="legacy"\nsegments=1\n', encoding="utf-8")
    with pytest.raises(ConfigError, match="model.rod"):
        load_config(path)


@pytest.mark.parametrize(
    ("dynamics", "match"),
    [
        ("", "dynamics must be a TOML table"),
        ("[dynamics]\n", "dynamics.formulation"),
        ("[dynamics]\nformulation='symbolic_el'\n", "symbolic_lagrange.*recursive"),
    ],
)
def test_dynamics_formulation_is_required_and_explicit(tmp_path, dynamics, match):
    path = tmp_path / "dynamics.toml"
    path.write_text(
        '[model]\nrod="euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
        + dynamics,
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match=match):
        load_config(path)
