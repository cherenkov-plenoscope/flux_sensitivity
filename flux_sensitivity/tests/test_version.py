import flux_sensitivity


def test_version_exists():
    assert len(str.split(flux_sensitivity.__version__, ".")) == 3