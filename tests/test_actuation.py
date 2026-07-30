from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from softarm import ActuationModel, derive, derive_actuation, load_config, register_actuator
from softarm.config import ActuationConfig, TendonChannelConfig, TendonSpanConfig


ROOT = Path(__file__).parents[1]


def _three_tendon():
    plant = derive(load_config(ROOT / "examples/config/pcc_three_tendon_extensible_n2.toml"))
    actuation = derive_actuation(plant)
    assert actuation is not None
    return plant, actuation


def _values(plant, actuation, q_values):
    return {
        **{symbol: value for symbol, value in zip(plant.q, q_values)},
        **{item.symbol: item.default for item in plant.parameters + actuation.parameters},
    }


def test_three_tendon_axial_mode_jacobian_and_virtual_work():
    plant, actuation = _three_tendon()
    assert actuation.count == 3
    assert actuation.channel_names == ("t1", "t2", "t3")
    assert actuation.channel_kinds == ("unilateral",) * 3
    values = _values(plant, actuation, [0.01, -0.02, 0.46, 0.03, 0.01, 0.51])
    jacobian = np.asarray(actuation.jacobian.subs(values), dtype=float)
    np.testing.assert_allclose(jacobian[:, [2, 5]], np.ones((3, 2)), atol=1e-14)
    np.testing.assert_allclose(jacobian.sum(axis=0)[[0, 1, 3, 4]], 0.0, atol=1e-14)

    step = 1e-7
    numerical = np.empty_like(jacobian)
    for column, symbol in enumerate(plant.q):
        plus = dict(values)
        minus = dict(values)
        plus[symbol] += step
        minus[symbol] -= step
        numerical[:, column] = (
            np.asarray(actuation.coordinates.subs(plus), dtype=float).reshape(3)
            - np.asarray(actuation.coordinates.subs(minus), dtype=float).reshape(3)
        ) / (2 * step)
    np.testing.assert_allclose(jacobian, numerical, rtol=1e-8, atol=1e-9)

    tension = sp.Matrix([2.0, 1.5, 1.0])
    displacement = sp.Matrix([0.003, -0.002, 0.001, 0.004, -0.001, 0.002])
    tau = actuation.generalized_force(tension)
    np.testing.assert_allclose(
        float((tau.T * displacement).subs(values)[0]),
        float((-tension.T * actuation.jacobian * displacement).subs(values)[0]),
        rtol=1e-12,
        atol=1e-14,
    )
    assert actuation.is_feasible(tension)
    assert not actuation.is_feasible(sp.Matrix([-1.0, 0.0, 0.0]))


def test_signed_channels_allow_negative_tension_and_have_no_axial_action():
    plant = derive(load_config(ROOT / "examples/config/pcc_signed_pair_n2.toml"))
    actuation = derive_actuation(plant)
    assert actuation is not None
    assert actuation.channel_kinds == ("signed", "signed")
    assert actuation.jacobian[:, 2] == sp.zeros(2, 1)
    assert actuation.jacobian[:, 5] == sp.zeros(2, 1)
    tension = sp.Matrix([-2.0, 1.5])
    assert actuation.is_feasible(tension)
    assert actuation.generalized_force(tension)[2] == 0
    assert actuation.generalized_force(tension)[5] == 0


def test_euler_uses_ritz_end_slope_without_axial_coordinate():
    plant = derive(load_config(ROOT / "examples/config/euler_ritz_n2.toml"))
    routing = ActuationConfig("tendon", "none", (
        TendonChannelConfig(
            "physical_x", "unilateral", (TendonSpanConfig(1, 0.02, 0.0),)
        ),
        TendonChannelConfig(
            "signed_y", "signed", (TendonSpanConfig(1, 0.02, np.pi / 2),)
        ),
    ))
    actuation = derive_actuation(plant, routing)
    assert actuation is not None
    defaults = {item.symbol: item.default for item in plant.parameters + actuation.parameters}
    reference = {coordinate: 0 for coordinate in plant.q}
    coordinate = np.asarray(actuation.coordinates.subs(defaults).subs(reference), dtype=float).reshape(2)
    np.testing.assert_allclose(coordinate, [0.45, 0.0], atol=1e-14)
    jacobian = np.asarray(actuation.jacobian.subs(defaults), dtype=float)
    np.testing.assert_allclose(jacobian[0, 0], -0.02 * 1.5 / 0.45, rtol=1e-12)
    np.testing.assert_allclose(jacobian[1, 1], -0.02 * 1.5 / 0.45, rtol=1e-12)


def test_strict_acceleration_rejects_dependent_or_excess_channels():
    plant = derive(load_config(ROOT / "examples/config/pcc_lumped_n2.toml"))
    span = (TendonSpanConfig(1, 0.02, 0.0),)
    dependent = ActuationConfig("tendon", "strict", (
        TendonChannelConfig("a", "unilateral", span),
        TendonChannelConfig("b", "unilateral", span),
    ))
    with pytest.raises(ValueError, match="full-row-rank"):
        derive_actuation(plant, dependent)

    channels = tuple(
        TendonChannelConfig(f"t{index}", "signed", (TendonSpanConfig(1, 0.02, index),))
        for index in range(7)
    )
    with pytest.raises(ValueError, match="7 constraints"):
        derive_actuation(plant, ActuationConfig("tendon", "strict", channels))
    assert derive_actuation(plant, ActuationConfig("tendon", "none", channels)).count == 7


def test_custom_actuator_builder_registration():
    plant = derive(load_config(ROOT / "examples/config/pcc_lumped_n2.toml"))

    def builder(symbolic_plant, config):
        coordinate = sp.Matrix([symbolic_plant.q[0]])
        jacobian = coordinate.jacobian(symbolic_plant.q)
        return ActuationModel(
            config.family, config.acceleration, ("custom",), ("signed",), (),
            coordinate, jacobian, sp.zeros(1, 1),
        )

    register_actuator("test_custom_actuator", builder)
    result = derive_actuation(plant, ActuationConfig("test_custom_actuator", "none"))
    assert result is not None
    assert result.channel_names == ("custom",)
