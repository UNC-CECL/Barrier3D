import importlib

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from barrier3d import Barrier3d, Barrier3dConfiguration, load_inputs
from barrier3d.load_input import as_cwd


def test_parameters_to_expected(datadir):
    with as_cwd(datadir):
        mod = importlib.import_module("barrier3d-parameters-expected")
        expected = dict(
            [(k, v) for k, v in mod.__dict__.items() if not k.startswith("_")]
        )
        actual = load_inputs(datadir, prefix="barrier3d", fmt="yaml")

    for key in actual:
        assert key in expected
        try:
            assert actual[key] == expected[key], key
        except ValueError:
            assert_array_almost_equal(actual[key], expected[key], err_msg=key)
    # assert actual == expected


def test_barrier3d_configuration(datadir):
    expected = load_inputs(datadir, prefix="barrier3d", fmt="py")
    actual = load_inputs(datadir, prefix="barrier3d", fmt="yaml")

    assert set(actual.keys()) - set(expected.keys()) == set()
    assert set(expected.keys()) - set(actual.keys()) == set()

    actual = dict(actual.items())

    assert_array_almost_equal(
        actual.pop("InteriorDomain"), expected.pop("InteriorDomain")
    )
    assert_array_almost_equal(actual.pop("PC"), expected.pop("PC"))
    assert_array_almost_equal(actual.pop("growthparam"), expected.pop("growthparam"))
    assert_array_almost_equal(actual.pop("DuneDomain"), expected.pop("DuneDomain"))
    assert_array_almost_equal(actual.pop("StormSeries"), expected.pop("StormSeries"))

    assert actual == expected


def test_bad_input_format(datadir):
    with pytest.raises(ValueError):
        load_inputs(datadir, prefix="barrier3d", fmt="xlsx")


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_missing_file(datadir, fmt):
    with pytest.raises(FileNotFoundError, match=f"missing-prefix-parameters.{fmt}"):
        load_inputs(datadir, prefix="missing-prefix", fmt=fmt)


def test_missing_folder(datadir):
    with pytest.raises(FileNotFoundError, match="missing-folder"):
        load_inputs("missing-folder", fmt="yaml")


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_init(datadir, fmt):
    params = load_inputs(datadir, prefix="barrier3d", fmt=fmt)
    barrier3d = Barrier3d(**params)


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_from_path(datadir, fmt):
    barrier3d = Barrier3d.from_path(datadir, fmt=fmt)


@pytest.mark.parametrize("fmt", ["yaml", "py"])
def test_barrier3d_update(datadir, fmt):
    barrier3d = Barrier3d.from_path(datadir, fmt=fmt)
    barrier3d.update()
