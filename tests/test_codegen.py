import json
from pathlib import Path

from softarm.ast import decode, encode
from softarm.codegen import generate_matlab_bundle
from softarm import derive_actuation, load_config
from softarm.config import (
    ActuationConfig, ConstraintConfig, IntegrationConfig, ModelConfig, TendonChannelConfig,
    TendonSpanConfig,
)
from softarm.derive import derive
from softarm.constraints import derive_constraint


ROOT = Path(__file__).parents[1]


def _small_plant():
    return derive(ModelConfig(
        family="euler", segments=1, integration=IntegrationConfig("analytic"),
        ritz_x=(0.0, 0.0, 1.5, -0.5), ritz_y=(0.0, 0.0, 1.5, -0.5),
    ))


def test_ast_round_trip_has_no_metadata():
    plant = _small_plant()
    node = encode(plant.mass[0, 0])
    assert "version" not in node
    assert decode(node) == plant.mass[0, 0]


def test_minimal_manifest_and_fixed_functions(tmp_path):
    plant = _small_plant()
    generate_matlab_bundle(plant, tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert set(manifest) == {"model", "coordinates", "parameters", "actuation", "constraint"}
    assert manifest["coordinates"]["base"] == []
    assert len(manifest["coordinates"]["arm"]) == 2
    for filename in (
        "softarm_mass.m", "softarm_bias.m", "softarm_kinematics.m",
        "softarm_end_jacobian.m", "softarm_applied_force.m",
        "softarm_forward_dynamics.m", "softarm_state_rhs.m",
    ):
        assert (tmp_path / filename).is_file()


def test_actuated_bundle_keeps_minimal_manifest_and_generic_functions(tmp_path):
    config = load_config(ROOT / "examples/config/pcc_three_tendon_extensible_n2.toml")
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
        family="euler", segments=1, integration=IntegrationConfig("analytic"),
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
