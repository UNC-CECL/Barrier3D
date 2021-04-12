"""MakeTimeSeries

These functions are for creating and saving storm time series, initial dune height, and dune growth rates for
use as inputs in Barrier3D simulations.

"""

import numpy as np
import random
import matplotlib.pyplot as plt


def storms_per_year_from_MSSM_output(
    datadir,
    name,
    storm_list_name="VCRStormList.npy",
    mean_storm=8.3,
    SD_storm=5.9,
    MHW=0.46,
    StormStart=2,
    BermEl=1.9,
    model_years=10000,
):
    r"""This function uses a normal distribution (provided the mean and standard deviation) to select a random number
    of storms per year from a list of multivariate sea storms, created using the Wahl et al., 2016 method.

    """

    # convert to decameters
    MHW = MHW / 10
    BermEl = BermEl / 10 - MHW  # just for plotting

    # load list of storms (created using a MSSM model)
    StormList = np.load(datadir + storm_list_name)

    StormSeries = np.zeros([StormStart, 5])  # pad with zeros until storms start

    for t in range(StormStart, model_years):

        # Calculate number of storms in year
        numstorm = round(np.random.normal(mean_storm, SD_storm))

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
            stormTS[n, 1] = Rhigh / 10 - MHW
            stormTS[n, 2] = Rlow / 10 - MHW
            stormTS[n, 3] = period
            stormTS[n, 4] = round(
                dur / 2
            )  # Divided by two assuming TWL only for only half of storm
            # (we could do better here, like assume that the TWL follows a distribution;
            # or create a time series as in the Wahl follow up paper)

        StormSeries = np.vstack([StormSeries, stormTS])

    # Plots
    Bin = np.linspace(-1, 4.6, 57)

    plt.figure()
    surgetidecalc = StormSeries[:, 2] * 10
    plt.hist(surgetidecalc, bins=Bin)
    plt.title("StormSeries Rlow")

    plt.figure()
    twlcalc = StormSeries[:, 1] * 10
    plt.hist(twlcalc, bins=Bin)
    plt.title("StormSeries Rhigh")

    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(16, 4)
    plt.plot(twlcalc)
    plt.xlabel("Storm")
    plt.hlines(BermEl * 10, -5, len(twlcalc) + 20, colors="red", linestyles="dashed")
    plt.show()
    plt.ylabel("TWL (m MHW)")

    print("Max TWL (m):          ", max(StormSeries[:, 1]) * 10)
    print("Max Rexcess (m):      ", (max(StormSeries[:, 1]) - BermEl) * 10)

    print("% Rhigh > BermEl: ", (twlcalc > (BermEl * 10)).sum() / len(twlcalc) * 100)
    print(
        "% Rlow  > BermEl: ",
        (surgetidecalc > (BermEl * 10)).sum() / len(surgetidecalc) * 100,
    )

    return np.save(datadir + name, StormSeries)


# Generate dune height start
def gen_dune_height_start(datadir, name, Dstart=0.5, ny=1000):

    # convert to decameters
    Dstart = Dstart / 10

    DuneStart = np.ones([ny]) * (
        Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
    )

    return np.save(datadir + name, DuneStart)


# Generate along-shore varying rmin & rmax
def gen_alongshore_variable_rmin_rmax(datadir, name, rmin=0.35, rmax=0.85, ny=1000):

    growthparam = rmin + (rmax - rmin) * np.random.rand(ny)

    return np.save(datadir + name, growthparam)
