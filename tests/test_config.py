from pathlib import Path

import pytest

from softarm.config import ConfigError, load_config


ROOT = Path(__file__).parents[1]


def test_reference_configs_are_valid():
    for path in (ROOT / "examples" / "config").glob("*.toml"):
        assert load_config(path).segments >= 1


def test_base_and_constraint_validation(tmp_path):
    path = tmp_path / "floating.toml"
    path.write_text(
        '[model]\nfamily="pcc"\nsegments=1\n'
        '[base]\nmode="floating_rpy"\nmount_xyz=[0,0,0]\nmount_rpy=[0,0,0]\n'
        '[constraint]\nfamily="plane_point_contact"\nplane_normal=[0,0,0]\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="plane_normal.*nonzero"):
        load_config(path)


def test_euler_requires_explicit_normalized_ritz(tmp_path):
    path = tmp_path / "bad.toml"
    path.write_text('[model]\nfamily="euler"\nsegments=1\n[ritz]\nx=[0,0,1]\ny=[0,0,2]\n')
    with pytest.raises(ConfigError, match="normalized"):
        load_config(path)


@pytest.mark.parametrize("actuation,match", [
    ('family="tendon"\n', "acceleration"),
    ('family="tendon"\nacceleration="auto"\n', "strict.*none"),
    ('family="tendon"\nacceleration="none"\n', "at least one channel"),
])
def test_tendon_actuation_requires_explicit_valid_structure(tmp_path, actuation, match):
    path = tmp_path / "bad_actuation.toml"
    path.write_text(
        '[model]\nfamily="pcc"\nsegments=1\n[actuation]\n' + actuation,
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match=match):
        load_config(path)


def test_tendon_config_rejects_duplicate_span(tmp_path):
    path = tmp_path / "duplicate.toml"
    path.write_text(
        '[model]\nfamily="pcc"\nsegments=1\n'
        '[actuation]\nfamily="tendon"\nacceleration="none"\n'
        '[[actuation.channels]]\nname="t1"\nkind="unilateral"\n'
        '[[actuation.channels.spans]]\nsection=1\nradius=0.02\nangle=0\n'
        '[[actuation.channels.spans]]\nsection=1\nradius=0.03\nangle=1\n',
        encoding="utf-8",
    )
    with pytest.raises(ConfigError, match="repeats section"):
        load_config(path)
