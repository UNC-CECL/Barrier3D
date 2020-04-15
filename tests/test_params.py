import numpy as np
import xlrd
from numpy.testing import assert_array_equal

from barrier3d.parameters import from_xls, load_params
from barrier3d.barrier3d import Barrier3D


def test_from_xls(datadir):
    params = from_xls(datadir / "barrier3d.xlsx")
    assert isinstance(params, dict)


def test_load_params(datadir):
    expected = from_xls(datadir / "barrier3d.xlsx")
    actual = load_params(datadir / "barrier3d.xlsx")

    assert set(actual.keys()) - set(expected.keys()) == set()
    assert set(expected.keys()) - set(actual.keys()) == set()

    assert np.all(actual.pop("InteriorDomain") == expected.pop("InteriorDomain"))
    assert np.all(actual.pop("PC") == expected.pop("PC"))

    # These are arrays filled with random numbers
    assert actual.pop("growthparam").shape == expected.pop("growthparam").shape
    assert actual.pop("DuneDomain").shape == expected.pop("DuneDomain").shape

    assert actual == expected


def test_barrier3d_init(datadir):
    params = load_params(datadir / "barrier3d.xlsx")
    barrier3d = Barrier3D(**params)


def test_barrier3d_from_xls(datadir):
    barrier3d = Barrier3D.from_xls(datadir / "barrier3d.xlsx")
