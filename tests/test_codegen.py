import json
from dataclasses import replace
from pathlib import Path

import sympy as sp

from softarm import derive_actuation, load_config
from softarm.backends.protocol import decode_dag, encode_dag
from softarm.codegen import generate_matlab_bundle
from softarm.codegen.matlab import _HELPERS
from softarm.config import (
    ActuationConfig,
    BaseConfig,
    ConstraintConfig,
    IntegrationConfig,
    ModelConfig,
    TendonChannelConfig,
    TendonSpanConfig,
)
from softarm.constraints import derive_constraint
from softarm.derive import derive
from softarm.geometry import pac_transform, pcs_transform

ROOT = Path(__file__).parents[1]


def _small_plant():
    return derive(ModelConfig(
        rod="euler_bernoulli", parameterization="ritz", segments=1,
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5), ritz_y=(0.0, 0.0, 1.5, -0.5),
    ))


def test_ast_round_trip_uses_the_current_dag_schema():
    plant = _small_plant()
    graph = encode_dag([plant.mass[0, 0]])
    assert set(graph) == {"nodes", "roots"}
    assert decode_dag(graph) == [plant.mass[0, 0]]


def test_ast_round_trip_preserves_cosserat_special_functions():
    symbols = sp.symbols("kx ky kz vx vy vz L", real=True)
    expression = pcs_transform(
        sp.Matrix(symbols[:3]), sp.Matrix(symbols[3:6]), symbols[6]
    )[0, 3]
    assert decode_dag(encode_dag([expression])) == [expression]


def test_ast_round_trip_preserves_pac_moment_functions():
    c0, c1, phi, length, xi = sp.symbols("c0 c1 phi length xi", real=True)
    expression = pac_transform(c0, c1, phi, length, xi)[0, 3]
    assert decode_dag(encode_dag([expression])) == [expression]


def test_minimal_manifest_and_fixed_functions(tmp_path):
    plant = _small_plant()
    generate_matlab_bundle(plant, tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert set(manifest) == {"model", "coordinates", "parameters", "actuation", "constraint"}
    assert manifest["coordinates"]["base"] == []
    assert manifest["model"]["rod"] == "euler_bernoulli"
    assert manifest["model"]["parameterization"] == "ritz"
    assert manifest["model"]["base_mode"] == "fixed"
    assert manifest["model"]["mount_xyz"] == [0.0, 0.0, 0.0]
    assert manifest["model"]["mount_rpy"] == [0.0, 0.0, 0.0]
    assert manifest["model"]["centerline"] is True
    assert len(manifest["coordinates"]["arm"]) == 2
    for filename in (
        "softarm_mass.m", "softarm_bias.m", "softarm_kinematics.m",
        "softarm_centerline.m", "softarm_centerline_at.m",
        "softarm_end_jacobian.m", "softarm_applied_force.m",
        "softarm_forward_dynamics.m", "softarm_state_rhs.m", "softarm_model.tex",
    ):
        assert (tmp_path / filename).is_file()
    document = (tmp_path / "softarm_model.tex").read_text(encoding="utf-8")
    assert document.startswith(r"\documentclass[11pt]{article}")
    assert r"\begin{document}" in document
    assert document.endswith("\\end{document}\n")
    assert "Exact Symbolic Appendix" not in document


def test_reference_bundles_include_all_fixed_matlab_helpers():
    manifests = sorted((ROOT / "examples" / "generated").glob("*/manifest.json"))
    assert manifests
    for manifest in manifests:
        for filename, expected in _HELPERS.items():
            path = manifest.parent / filename
            assert path.is_file(), f"{manifest.parent.name} is missing {filename}"
            assert path.read_text(encoding="utf-8") == expected


def test_pac_bundle_contains_moment_helpers_and_public_coordinates(tmp_path):
    plant = derive(ModelConfig(
        rod="extensible_euler_bernoulli", parameterization="pac", segments=1,
        inertia="lumped", integration=IntegrationConfig(),
    ))
    generate_matlab_bundle(plant, tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["coordinates"]["arm"] == ["c0_1", "c1_1", "phi1", "l1"]
    assert manifest["model"]["parameterization"] == "pac"
    for filename in (
        "softarm_affine_cos_moment.m",
        "softarm_affine_sin_moment.m",
        "softarm_affine_moment_parts.m",
    ):
        assert (tmp_path / filename).is_file()
    assert "softarm_affine_" in (tmp_path / "softarm_kinematics.m").read_text()
    document = (tmp_path / "softarm_model.tex").read_text(encoding="utf-8")
    assert "Piecewise-Affine-Curvature" in document
    assert r"R_i(0)=R_z(\phi_i)" in document


def test_bundle_without_material_kinematics_omits_centerline(tmp_path):
    plant = _small_plant()
    generate_matlab_bundle(plant, tmp_path)
    plant._material_coordinate = None
    plant._material_kinematics = None
    generate_matlab_bundle(plant, tmp_path)
    model = json.loads((tmp_path / "manifest.json").read_text())["model"]
    assert model["centerline"] is False
    assert not (tmp_path / "softarm_centerline.m").exists()
    assert not (tmp_path / "softarm_centerline_at.m").exists()


def test_manifest_preserves_nonzero_base_mount(tmp_path):
    plant = _small_plant()
    plant.config = replace(
        plant.config,
        base=BaseConfig(
            mode="fixed",
            mount_xyz=(0.12, -0.23, 0.34),
            mount_rpy=(0.41, -0.32, 0.13),
        ),
    )
    generate_matlab_bundle(plant, tmp_path)
    model = json.loads((tmp_path / "manifest.json").read_text())["model"]
    assert model["mount_xyz"] == [0.12, -0.23, 0.34]
    assert model["mount_rpy"] == [0.41, -0.32, 0.13]


def test_actuated_bundle_keeps_minimal_manifest_and_generic_functions(tmp_path):
    config = load_config(
        ROOT / "examples/config/extensible_euler_bernoulli_pcs_three_tendon_n2.toml"
    )
    plant = derive(config)
    actuation = derive_actuation(plant)
    generate_matlab_bundle(plant, tmp_path, actuation=actuation)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert set(manifest) == {"model", "coordinates", "parameters", "actuation", "constraint"}
    assert [item["name"] for item in manifest["parameters"]][-2:] == [
        "act_t3_s1_radius", "act_t3_s2_radius",
    ]
    for filename in (
        "softarm_actuator_coordinates.m",
        "softarm_actuator_jacobian.m",
            "softarm_actuator_velocity_bias.m",
            "softarm_actuator_force.m",
            "softarm_actuator_acceleration.m",
    ):
        assert (tmp_path / filename).is_file()
    assert not list(tmp_path.glob("*two_tendon*"))
    assert not (tmp_path / "softarm_actuator_info.m").exists()


def test_force_only_actuation_omits_strict_acceleration_function(tmp_path):
    plant = _small_plant()
    channels = tuple(
        TendonChannelConfig(
            f"s{index}", "signed", (TendonSpanConfig(1, 0.02, float(index)),)
        )
        for index in range(3)
    )
    actuation = derive_actuation(plant, ActuationConfig("tendon", "none", channels))
    generate_matlab_bundle(plant, tmp_path, actuation=actuation)
    assert (tmp_path / "softarm_actuator_force.m").is_file()
    assert not (tmp_path / "softarm_actuator_acceleration.m").exists()


def test_constraint_bundle_has_manifest_metadata_and_solver(tmp_path):
    config = ModelConfig(
        rod="euler_bernoulli", parameterization="ritz", segments=1,
        integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5), ritz_y=(0.0, 0.0, 1.5, -0.5),
        constraint=ConstraintConfig("plane_point_contact", {
            "family": "plane_point_contact", "plane_normal": [0.0, 0.0, -1.0],
        }),
    )
    plant = derive(config)
    constraint = derive_constraint(plant)
    generate_matlab_bundle(plant, tmp_path, constraint=constraint)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["constraint"]["family"] == "plane_point_contact"
    assert manifest["constraint"]["channels"] == [
        {"name": "normal_contact", "kind": "unilateral"}
    ]
    for filename in (
        "softarm_constraint_value.m", "softarm_constraint_jacobian.m",
        "softarm_constraint_velocity_bias.m", "softarm_constraint_reaction_map.m",
        "softarm_constraint_acceleration.m",
    ):
        assert (tmp_path / filename).is_file()
