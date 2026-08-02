from pathlib import Path

import pytest

from softarm.config import ConfigError, load_config
from softarm.derive import derive

ROOT = Path(__file__).parents[1]


def test_reference_configs_are_valid():
    for path in (ROOT / "examples" / "config").glob("*.toml"):
        config = load_config(path)
        assert derive(config).config.segments >= 1


def test_base_and_constraint_validation(tmp_path):
    path = tmp_path / "floating.toml"
    path.write_text(
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
        '[base]\nmode="floating_rpy"\nmount_xyz=[0,0,0]\nmount_rpy=[0,0,0]\n'
        '[constraint]\nfamily="plane_point_contact"\nplane_normal=[0,0,0]\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="plane_normal.*nonzero"):
        load_config(path)


def test_euler_bernoulli_requires_explicit_normalized_ritz(tmp_path):
    path = tmp_path / "bad.toml"
    path.write_text('[model]\nrod="euler_bernoulli"\nparameterization="ritz"\nsegments=1\n[ritz]\nx=[0,0,1]\ny=[0,0,2]\n')
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
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n[actuation]\n' + actuation,
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match=match):
        load_config(path)


def test_tendon_config_rejects_duplicate_span(tmp_path):
    path = tmp_path / "duplicate.toml"
    path.write_text(
        '[model]\nrod="extensible_euler_bernoulli"\nparameterization="pcs"\nsegments=1\n'
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
        '[integration]\nmethod="analytic"\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="requires integration.method='gauss'"):
        derive(load_config(distributed))

    order_one = tmp_path / "order_one.toml"
    order_one.write_text(
        '[model]\nrod="cosserat"\nparameterization="pcs"\nsegments=1\ninertia="distributed"\n'
        '[integration]\nmethod="gauss"\norder=1\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="order at least 2"):
        derive(load_config(order_one))

    lumped = tmp_path / "lumped.toml"
    lumped.write_text(
        '[model]\nrod="cosserat"\nparameterization="pcs"\nsegments=1\ninertia="lumped"\n'
        '[integration]\nmethod="gauss"\norder=2\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="not applicable"):
        derive(load_config(lumped))


def test_unregistered_combination_is_discovered_from_registry(tmp_path):
    path = tmp_path / "unsupported.toml"
    path.write_text(
        '[model]\nrod="cosserat"\nparameterization="ritz"\nsegments=1\n',
        encoding="utf-8",
    )
    config = load_config(path)
    with pytest.raises(ValueError, match="unregistered model combination"):
        derive(config)


def test_pac_validation_and_explicit_cosserat_rejection(tmp_path):
    low_order = tmp_path / "low_order.toml"
    low_order.write_text(
        '[model]\nrod="euler_bernoulli"\nparameterization="pac"\nsegments=1\n'
        'inertia="distributed"\n[integration]\nmethod="gauss"\norder=3\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="order at least 4"):
        derive(load_config(low_order))

    unsupported = tmp_path / "cosserat_pac.toml"
    unsupported.write_text(
        '[model]\nrod="cosserat"\nparameterization="pac"\nsegments=1\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unregistered model combination"):
        derive(load_config(unsupported))


def test_legacy_model_family_is_rejected(tmp_path):
    path = tmp_path / "legacy.toml"
    path.write_text('[model]\nfamily="legacy"\nsegments=1\n', encoding="utf-8")
    with pytest.raises(ConfigError, match="model.rod"):
        load_config(path)
