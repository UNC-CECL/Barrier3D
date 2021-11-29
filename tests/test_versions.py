import numpy as np
from barrier3d import Barrier3dBmi
import pytest


def current_version_check_against_V1():

    # specify data directories with initial conditions
    datadir_V1 = "V1_NoBMI/"
    datadir_current = "tests/test_versions/"

    # create an instance of the new BMI class, which is the model
    barrier3d = Barrier3dBmi()
    input_file = "barrier3d-parameters.yaml"
    barrier3d.initialize(datadir_current + input_file)

    # increase time step
    for time_step in range(1, barrier3d._model._TMAX):
        barrier3d.update()

    # Version 1.0 ------------------------------

    DuneDomain = []
    x_s_TS = []
    ShorelineChangeTS = []
    s_sf_TS = []
    QowTS = []

    # starts running immediately and ends with plots!
    exec(open(datadir_V1 + "Barrier3D.py").read())
    # execfile(datadir_V1 + "Barrier3D.py")

    # ----------------------------- #

    assert np.all(barrier3d._model._DuneDomain == DuneDomain)  # dune height over time
    assert np.all(barrier3d._model._x_s_TS == x_s_TS)  # shoreline change time series
    assert np.all(
        barrier3d._model._ShorelineChangeTS == ShorelineChangeTS
    )  # dune migration time series
    assert np.all(barrier3d._model._s_sf_TS == s_sf_TS)  # shoreface slope time series
    assert np.all(barrier3d._model._QowTS == QowTS)
