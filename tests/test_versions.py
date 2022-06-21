import numpy as np
import pathlib
import os
from barrier3d import Barrier3dBmi

# BMI_DATA_DIR = pathlib.Path(__file__).parent / "test_versions_inputs"
# V1_RUN_FILE_PATH = os.path.abspath(
#     os.path.join(__file__, "../../V1_NoBMI/Barrier3D.py")
# )
# V1_MAIN_DIR = os.path.abspath(os.path.join(__file__, "../../"))

datadir_V1 = "V1_NoBMI/"
bmi_datadir = "tests/test_versions_inputs/"


def test_BMI_against_V1():

    # Version 1.0 (original version written by Ian Reeves) ------------------------------
    DuneDomain = []
    x_s_TS = []
    ShorelineChangeTS = []
    s_sf_TS = []
    QowTS = []

    # os.chdir(V1_MAIN_DIR)
    # exec(open(V1_RUN_FILE_PATH).read())
    exec(open(datadir_V1 + "Barrier3D.py").read())

    # Version 2.0 and beyond (contains a BMI) ------------------------------
    barrier3d = Barrier3dBmi()
    input_file = "barrier3d-parameters.yaml"
    # barrier3d.initialize(str(BMI_DATA_DIR / input_file))
    barrier3d.initialize(bmi_datadir + input_file)

    # increase time step
    for time_step in range(1, barrier3d._model._TMAX):
        barrier3d.update()

    assert np.all(barrier3d._model._DuneDomain == DuneDomain)  # dune height over time
    assert np.all(
        barrier3d._model._x_s_TS == x_s_TS
    )  # shoreline change time series, BMI consistently less
    assert np.all(
        barrier3d._model._ShorelineChangeTS == ShorelineChangeTS  #
    )  # dune migration time series
    assert np.all(barrier3d._model._s_sf_TS == s_sf_TS)  # shoreface slope time series
    assert np.all(
        barrier3d._model._QowTS == QowTS
    )  # because we are producing less overwash
