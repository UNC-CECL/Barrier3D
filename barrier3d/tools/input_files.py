"""
These functions create, return, and if specified save time series of
annual storm parameters, initial dune height, and dune growth rates
for use as inputs in Barrier3D simulations.

References
----------
.. [1] Wahl, T., Plant, N. G., & Long, J. W. (2016). Probabilistic assessment of erosion and flooding risk in the
northern Gulf of Mexico. Journal of Geophysical Research: Oceans, 121(5), 3029-3043.

Notes
-----

The storm time series can be created using one of three functions, all of
which require a list of storms generated using the multivariateSeaStorm.m
module, developed following the method of ...[1]:
1. yearly_storms - specify the mean and standard deviation
   from a normal distribution, which are used to select a random number
   of storms per year from the MSSM list
2. frequency_storms - specify the TWL and frequency of a
   return period storm, find closest match in MSSM list
3. shift_storm_intensity - shifts the TWL distribution created in 1) to
   the left or right to produce storms of different "intensity".

"""
import bisect
import pathlib
import random
from distfit import distfit

import matplotlib.pyplot as plt
import numpy as np

from ..load_input import _guess_format


def yearly_storms(
    datadir=".",
    storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
    mean_yearly_storms=8.3,
    SD_yearly_storms=5.9,
    MHW=0.46,  # m NAVD88
    StormStart=2,
    BermEl=1.9,  # m NAVD88, just used for plotting
    model_years=10000,
    bPlot=True,
    bSave=False,
    output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04",
):
    """
    Use a normal distribution (provided the mean and standard deviation)
    to select a random number of storms per year from a list of
    multivariate sea storms, created using the Wahl et al., 2016 method.
    """

    datadir = pathlib.Path(datadir)

    # convert to decameters
    MHW = MHW / 10
    BermEl = BermEl / 10 - MHW  # just for plotting

    # load list of storms (created using multivariateSeaStorm.m)
    fmt = _guess_format(datadir / storm_list_name)
    if fmt == "npy":
        StormList = np.load(datadir / storm_list_name, allow_pickle=True)
    elif fmt == "csv":
        StormList = np.loadtxt(datadir / storm_list_name, delimiter=",", encoding='utf-8-sig')

    # pad with zeros until storms start
    StormSeries = np.zeros([StormStart, 5])

    for t in range(StormStart, model_years):

        # Calculate number of storms in year
        numstorm = round(np.random.normal(mean_yearly_storms, SD_yearly_storms))
        if numstorm < 0:
            numstorm = 0
        stormTS = np.zeros([numstorm, 5])

        # Select storms for year
        for n in range(numstorm):
            storm = random.randint(1, len(StormList) - 1)

            dur = StormList[storm, 1]  # Duration
            Rhigh = StormList[storm, 2]  # TWL
            period = StormList[storm, 4]  # Tp
            Rlow = StormList[storm, 6]  # Rlow

            stormTS[n, 0] = t
            stormTS[n, 1] = Rhigh / 10 - MHW  # make relative to MHW
            stormTS[n, 2] = Rlow / 10 - MHW
            stormTS[n, 3] = period
            stormTS[n, 4] = round(
                dur / 2
            )  # Divided by two assuming TWL only for only half of storm
            # (NOTE from KA: we could probably do better here, like assume that the TWL follows a distribution;
            # or create a time series as in the Wahl 2019 follow up paper for MSSM. Future work!)

        StormSeries = np.vstack([StormSeries, stormTS])

    if bPlot:
        # Plots
        Bin = np.linspace(-1, 4.6, 57)

        plt.figure()
        surgetidecalc = StormSeries[:, 2] * 10
        plt.hist(surgetidecalc, bins=Bin)
        plt.title("StormSeries Rlow")
        plt.xlabel("m MHW")

        plt.figure()
        twlcalc = StormSeries[:, 1] * 10
        plt.hist(twlcalc, bins=Bin)
        plt.title("StormSeries Rhigh")
        plt.xlabel("m MHW")

        plt.figure()
        fig = plt.gcf()
        fig.set_size_inches(16, 4)
        plt.plot(twlcalc)
        plt.plot(
            np.arange(0, len(twlcalc), 1),
            np.ones(len(twlcalc)) * BermEl * 10,
            "r--",
        )
        plt.xlabel("Storm")
        plt.ylabel("TWL (m MHW)")
        plt.legend(["TWL", "BermEl"])

        print("Max TWL (m MHW): ", max(StormSeries[:, 1]) * 10)
        print("Max Rexcess (m above berm): ", (max(StormSeries[:, 1]) - BermEl) * 10)
        print(
            "% Rhigh > BermEl: ", (twlcalc > (BermEl * 10)).sum() / len(twlcalc) * 100
        )
        print(
            "% Rlow  > BermEl: ",
            (surgetidecalc > (BermEl * 10)).sum() / len(surgetidecalc) * 100,
        )

    if bSave:
        np.save(datadir / output_filename, StormSeries)

    return StormSeries


def shift_storm_intensity(
    datadir=".",
    storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
    mean_yearly_storms=8.3,
    SD_yearly_storms=5.9,
    shift=-0.15,  # shift the TWL distribution to change intensity, m NAVD88; [-0.15, 0.15] for Reeves et al., 2021
    MHW=0.46,  # m NAVD88
    StormStart=2,
    BermEl=1.9,  # m NAVD88, just used for plotting
    model_years=1000,
    bPlot=True,
    bSave=False,
    output_filename="StormList_10kyrs_VCR_Berm1pt9m_Slope0pt04-lowIntensity",
):
    """
    Fit a beta distribution to the TWL time series and then shift the beta
    distribution to the left or right to simulate TWLs of higher or lower intensities
    """
    datadir = pathlib.Path(datadir)

    # convert to decameters
    MHW = MHW / 10
    BermEl = BermEl / 10 - MHW  # just for plotting

    # load list of storms (created using multivariateSeaStorm.m)
    fmt = _guess_format(datadir / storm_list_name)
    if fmt == "npy":
        StormList = np.load(datadir / storm_list_name, allow_pickle=True)
    elif fmt == "csv":
        StormList = np.loadtxt(datadir / storm_list_name, delimiter=",")

    # sort the storms based on TWL, from min to max (probably a more elegant way to do this)
    dur = StormList[:, 1]  # Duration
    simTWL = StormList[:, 2]
    period = StormList[:, 4]  # Tp
    Rlow = StormList[:, 6]  # Rlow
    zip_storm_list = list(zip(simTWL, dur, period, Rlow))
    zip_storm_list.sort()
    simTWL_sorted = np.array([simTWL for (simTWL, dur, period, Rlow) in zip_storm_list])
    dur_sorted = np.array([dur for (simTWL, dur, period, Rlow) in zip_storm_list])
    period_sorted = np.array([period for (simTWL, dur, period, Rlow) in zip_storm_list])
    Rlow_sorted = np.array([Rlow for (simTWL, dur, period, Rlow) in zip_storm_list])

    # Fit Distribution - Note: this is typically beta for VCR TWLs, so we assume beta here
    dist = distfit(distr="beta")
    fit = dist.fit_transform(simTWL_sorted)
    dist.plot()
    model_parameters = dist.model  # dictionary of beta distribution parameters
    loc = model_parameters["loc"]
    scale = model_parameters["scale"]
    [a, b] = model_parameters["arg"]

    # Bin Specifications
    BinStart = round(np.min(fit["histdata"][1]), 1)  # round to one decimal
    BinStop = round(np.max(fit["histdata"][1]), 1)
    BinWidth = 0.1  # m
    BinNum = int(((BinStop - BinStart) / BinWidth) + 2)
    Bin = np.linspace(BinStart, BinStop, BinNum)

    # Make storm series
    StormSeries = np.zeros([StormStart, 5])

    for t in range(StormStart, model_years):

        # Calculate number of storms in year
        numstorm = max(0, round(np.random.normal(mean_yearly_storms, SD_yearly_storms)))
        stormTS = np.zeros([numstorm, 5])

        # Select storms for year
        for n in range(numstorm):
            ST = (
                np.random.beta(a, b) * scale + loc + shift
            )  # the TWL from the beta distribution
            if ST < 0:
                ST = 0
            elif ST > simTWL_sorted[-1]:
                ST = simTWL_sorted[-1]
            STmin = Bin[bisect.bisect(Bin, ST) - 1]
            STmax = STmin + 0.1
            indexMin = bisect.bisect(simTWL_sorted, STmin)
            indexMax = bisect.bisect(simTWL_sorted, STmax) - 1
            if indexMin == 0:
                storm = 1
            elif indexMin >= indexMax:
                storm = indexMin
            else:
                storm = random.randint(indexMin, indexMax)

            dur = dur_sorted[storm]
            Rhigh = simTWL_sorted[storm]
            period = period_sorted[storm]
            Rlow = Rlow_sorted[storm]

            stormTS[n, 0] = t
            stormTS[n, 1] = Rhigh / 10 - MHW
            stormTS[n, 2] = Rlow / 10 - MHW
            stormTS[n, 3] = period
            stormTS[n, 4] = round(
                dur / 2
            )  # Divided by two assuming TWL only for only half of storm

        # Save
        StormSeries = np.vstack([StormSeries, stormTS])

    if bPlot:
        # Plots
        Bin = np.linspace(-1, 4.6, 57)

        plt.figure()
        surgetidecalc = StormSeries[:, 2] * 10
        plt.hist(surgetidecalc, bins=Bin)
        plt.title("StormSeries Rlow")
        plt.xlabel("m MHW")

        plt.figure()
        twlcalc = StormSeries[:, 1] * 10
        plt.hist(twlcalc, bins=Bin)
        plt.title("StormSeries Rhigh")
        plt.xlabel("m MHW")

        plt.figure()
        fig = plt.gcf()
        fig.set_size_inches(16, 4)
        plt.plot(twlcalc)
        plt.plot(
            np.arange(0, len(twlcalc), 1),
            np.ones(len(twlcalc)) * BermEl * 10,
            "r--",
        )
        plt.xlabel("Storm")
        plt.ylabel("TWL (m MHW)")
        plt.legend(["TWL", "BermEl"])

        print("Max TWL (m MHW): ", max(StormSeries[:, 1]) * 10)
        print("Max Rexcess (m above berm): ", (max(StormSeries[:, 1]) - BermEl) * 10)
        print(
            "% Rhigh > BermEl: ", (twlcalc > (BermEl * 10)).sum() / len(twlcalc) * 100
        )
        print(
            "% Rlow  > BermEl: ",
            (surgetidecalc > (BermEl * 10)).sum() / len(surgetidecalc) * 100,
        )

    if bSave:
        np.save(datadir / output_filename, StormSeries)

    return StormSeries


def frequency_storms(
    datadir,
    storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
    MHW=0.46,
    return_period=50,  # minimum of 1
    return_period_TWL=2.0,  # in m above MHW (from NOAA Annual Exceedance Probability Curves)
    StormStart=2,
    BermEl=1.9,  # m NAVD88, just used for plotting
    model_years=1000,
    bPlot=True,
    bSave=False,
    output_filename="StormList_50yrRP_2mTWL_1kyrs_VCR_Berm1pt9m_Slope0pt04",
):
    """
    Select a storm from the list of multivariate sea storms -- created using
    the Wahl et al., 2016 method -- that matches the TWL specified for a
    given return period (in m above MHW), and returns the storm time
    series at the specified return period
    """

    datadir = pathlib.Path(datadir)

    # convert to decameters
    MHW = MHW / 10
    BermEl = BermEl / 10 - MHW  # just for plotting
    return_period_TWL = return_period_TWL / 10

    # load list of storms (created using multivariateSeaStorm.m)
    fmt = _guess_format(datadir / storm_list_name)
    if fmt == "npy":
        StormList = np.load(datadir / storm_list_name, allow_pickle=True)
    elif fmt == "csv":
        StormList = np.loadtxt(datadir / storm_list_name, delimiter=",")

    # find storm that has the closest TWL (in m above MHW) to the return period storm
    dur = np.round(
        StormList[:, 1] / 2
    )  # Duration, divided by two assuming TWL only for only half of storm
    Rhigh = StormList[:, 2] / 10 - MHW  # TWL in dam above MHW
    period = StormList[:, 4]  # Tp
    Rlow = StormList[:, 6] / 10 - MHW  # Rlow in dam relative to MHW

    [id, closest_TWL] = min(
        enumerate(Rhigh), key=lambda x: abs(x[1] - return_period_TWL)
    )

    # pad with zeros until storms start
    StormSeries = np.zeros([StormStart, 5])

    for t in range(StormStart, model_years):

        # only allow for one storm per year
        numstorm = 1
        stormTS = np.zeros([numstorm, 5])

        # Select storms for year
        if t % return_period == 0:

            stormTS[0, 0] = t
            stormTS[0, 1] = Rhigh[id]
            stormTS[0, 2] = Rlow[id]
            stormTS[0, 3] = period[id]
            stormTS[0, 4] = dur[id]

        StormSeries = np.vstack([StormSeries, stormTS])

    if bPlot:
        # Plots
        Bin = np.linspace(-1, 4.6, 57)

        plt.figure()
        surgetidecalc = StormSeries[:, 2] * 10
        plt.hist(surgetidecalc, bins=Bin)
        plt.title("StormSeries Rlow")
        plt.xlabel("m MHW")

        plt.figure()
        twlcalc = StormSeries[:, 1] * 10
        plt.hist(twlcalc, bins=Bin)
        plt.title("StormSeries Rhigh")
        plt.xlabel("m MHW")

        plt.figure()
        fig = plt.gcf()
        fig.set_size_inches(16, 4)
        plt.plot(twlcalc)
        plt.plot(
            np.arange(0, len(twlcalc), 1),
            np.ones(len(twlcalc)) * BermEl * 10,
            "r--",
        )
        plt.xlabel("Storm")
        plt.ylabel("TWL (m MHW)")
        plt.legend(["TWL", "BermEl"])

        print("Max TWL (m MHW): ", max(StormSeries[:, 1]) * 10)
        print("Max Rexcess (m above berm): ", (max(StormSeries[:, 1]) - BermEl) * 10)
        print(
            "% Rhigh > BermEl: ", (twlcalc > (BermEl * 10)).sum() / len(twlcalc) * 100
        )
        print(
            "% Rlow  > BermEl: ",
            (surgetidecalc > (BermEl * 10)).sum() / len(surgetidecalc) * 100,
        )

    if bSave:
        np.save(datadir / output_filename, StormSeries)

    return StormSeries


def gen_dune_height_start(datadir, name, Dstart=0.5, ny=1000):
    """Generate dune height start"""
    datadir = pathlib.Path(datadir)

    # convert to decameters
    Dstart = Dstart / 10

    DuneStart = np.ones([ny]) * (
        Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
    )

    return np.save(datadir / name, DuneStart)


def gen_alongshore_variable_rmin_rmax(datadir, name, rmin=0.35, rmax=0.85, ny=1000):
    """Generate along-shore varying rmin & rmax"""
    datadir = pathlib.Path(datadir)

    growthparam = rmin + (rmax - rmin) * np.random.rand(ny)

    return np.save(datadir / name, growthparam)
