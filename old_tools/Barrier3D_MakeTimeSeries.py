"""MakeTimeSeries

These functions create, return, and if specified save time series of annual storm parameters, initial dune height, and
dune growth rates for use as inputs in Barrier3D simulations.

References
----------
.. [1] Wahl, T., Plant, N. G., & Long, J. W. (2016). Probabilistic assessment of erosion and flooding risk in the
northern Gulf of Mexico. Journal of Geophysical Research: Oceans, 121(5), 3029-3043.

Notes
---------
The storm time series can be created using one of three functions, all of which require a list of storms generated
using the multivariateSeaStorm.m module, developed following the method of ...[1]:
1) yearly_storms_from_MSSM - specify the mean and standard deviation from a normal distribution, which are used to
   select a random number of storms per year from the MSSM list
2) frequency_storms_from_MSSM - specify the TWL and frequency of a return period storm, find closest match in MSSM list
3) shift_storm_intensity - shifts the TWL distribution created in 1) to the left or right to produce storms of different
   "intensity".

"""
import pathlib

import numpy as np


# Generate dune height start
def gen_dune_height_start(datadir, name, Dstart=0.5, ny=1000):
    datadir = pathlib.Path(datadir)

    # convert to decameters
    Dstart = Dstart / 10

    DuneStart = np.ones([ny]) * (
        Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
    )

    return np.save(datadir / name, DuneStart)


# Generate along-shore varying rmin & rmax
def gen_alongshore_variable_rmin_rmax(datadir, name, rmin=0.35, rmax=0.85, ny=1000):
    datadir = pathlib.Path(datadir)

    growthparam = rmin + (rmax - rmin) * np.random.rand(ny)

    return np.save(datadir / name, growthparam)
