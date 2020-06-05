import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_array_almost_equal

from barrier3d import Barrier3d, Barrier3dConfiguration, load_inputs
from barrier3d.parameters import _from_xlsx


def test_from_xlsx(datadir):
    params = _from_xlsx(datadir / "barrier3d.xlsx")
    assert isinstance(params, dict)


def test_barrier3d_configuration(datadir, fmt):
    expected = load_inputs(datadir, prefix="barrier3d", fmt="py")
    actual = load_inputs(datadir, prefix="barrier3d", fmt=fmt)

    # expected.update({"StormTimeSeries": 0, "StormSeries": []})

    assert set(actual.keys()) - set(expected.keys()) == set()
    assert set(expected.keys()) - set(actual.keys()) == set()

    actual = dict(actual.items())

    assert_array_almost_equal(actual.pop("InteriorDomain"), expected.pop("InteriorDomain"))
    assert_array_almost_equal(actual.pop("PC"), expected.pop("PC"))
    assert_array_almost_equal(actual.pop("growthparam"), expected.pop("growthparam"))
    assert_array_almost_equal(actual.pop("DuneDomain"), expected.pop("DuneDomain"))
    assert_array_almost_equal(actual.pop("StormSeries"), expected.pop("StormSeries"))

    # These are arrays filled with random numbers
    # assert actual.pop("growthparam").shape == expected.pop("growthparam").shape
    # assert actual.pop("DuneDomain").shape == expected.pop("DuneDomain").shape

    assert actual == expected


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_init(datadir, fmt):
    # params = Barrier3dParameters.from_xlsx(datadir / "barrier3d.xlsx")
    params = load_inputs(datadir, prefix="barrier3d", fmt=fmt)
    barrier3d = Barrier3d(**params)


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_from_path(datadir, fmt):
    barrier3d = Barrier3d.from_path(datadir, fmt=fmt)


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_update(datadir, fmt):
    barrier3d = Barrier3d.from_path(datadir, fmt=fmt)
    barrier3d.update()
