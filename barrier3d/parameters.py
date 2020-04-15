# Parameter value loading script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 1
# Updated: 22 Jan 2020


# Script loads all input parameters from Barrier_Input.xlsx
# Barrier3D_Input.xlsx must be in same directory as Barrier and Barrier_Functions
# Converts from meters to decameters for simulation


import numpy as np
import xlrd
import yaml


def from_yaml(path_to_yaml):
    with open(path_to_yaml, "r") as fp:
        params = yaml.safe_load(fp)

    return process_params(params)


def load_time_params(param):
    return {
        "TMAX": param.cell_value(4, 3),
        "StormStart": param.cell_value(5, 3),
    }


def load_domain_params(param):
    return {
        "LShoreface": param.cell_value(12, 3),
        "DShoreface": param.cell_value(13, 3),
        "BayDepth": param.cell_value(14, 3),
        # Used as offset to convert given elevations relative to a MHW of 0
        "MHW": param.cell_value(15, 3),
        "BarrierLength": param.cell_value(9, 3),
        "DuneWidth": param.cell_value(11, 3),
    }


def load_dune_params(param):
    return {
        "Dstart": param.cell_value(19, 3),
        "BermEl": param.cell_value(20, 3),
        "rmin": param.cell_value(21, 3),
        "rmax": param.cell_value(22, 3),
        "HdDiffu": param.cell_value(23, 3),
        "Dmaxel": param.cell_value(24, 3),
        "C1": param.cell_value(25, 3),
        "C2": param.cell_value(26, 3),
        "DuneRestart": param.cell_value(27, 3),
    }


def load_alongshore_transport_params(param):
    return {
        "Rat": param.cell_value(31, 3),
        "RSLR": param.cell_value(32, 3),
    }


def load_storm_params(param):
    return {
        "mean_storm": param.cell_value(36, 3),
        "SD_storm": param.cell_value(37, 3),
        "numstorm": param.cell_value(38, 3),
        # Water Level Forcing
        # In meters - converted to decameters after water level calculations
        "surge_tide_m": param.cell_value(39, 3),
        "surge_tide_sd": param.cell_value(40, 3),
        # height_m = param.cell_value(43,3)
        # height_sd = param.cell_value(44,3)
        "height_mu": param.cell_value(41, 3),
        "height_sigma": param.cell_value(42, 3),
        "period_m": param.cell_value(43, 3),
        "period_sd": param.cell_value(44, 3),
        "duration_mu": param.cell_value(45, 3),
        "duration_sigma": param.cell_value(46, 3),
        "beta": param.cell_value(47, 3),
    }


def load_overwash_params(param):
    return {
        "nn": param.cell_value(51, 3),
        "mm": param.cell_value(52, 3),
        "Rin_r": param.cell_value(53, 3),
        "Rin_i": param.cell_value(54, 3),
        "Qs_min": param.cell_value(55, 3),
        "threshold_in": param.cell_value(56, 3),
        # Sediment Transport
        "Kr": param.cell_value(57, 3),
        "Ki": param.cell_value(58, 3),
        "Cbb_r": param.cell_value(59, 3),
        "Cbb_i": param.cell_value(60, 3),
    }


def load_shoreface_dynamics_params(param):
    return {
        "k_sf": param.cell_value(64, 3),
        "s_sf_eq": param.cell_value(65, 3),
    }


def load_shrub_parameters(param):
    return {
        # Dispersal
        "Shrub_ON": param.cell_value(69, 3),
        "Seedmin": param.cell_value(70, 3),
        "Seedmax": param.cell_value(71, 3),
        "disp_mu": param.cell_value(72, 3),
        "disp_sigma": param.cell_value(73, 3),
        # Growth
        "Dshrub": param.cell_value(74, 3),
        "GermRate": param.cell_value(75, 3),
        "TimeFruit": param.cell_value(76, 3),
        "Female": param.cell_value(77, 3),
        "ShrubEl_min": param.cell_value(78, 3),
        "ShrubEl_max": param.cell_value(79, 3),
        # Overwash Interaction
        "BurialLimit": param.cell_value(80, 3),
        "UprootLimit": param.cell_value(81, 3),
        "SalineLimit": param.cell_value(82, 3),
        "Qshrub_max": param.cell_value(83, 3),
    }


def load_params(path_to_xls):
    # Open Excel Files
    wb = xlrd.open_workbook(path_to_xls)
    param = wb.sheet_by_index(0)  # Parameters sheet
    morph = wb.sheet_by_index(1)  # Initial morphology sheet

    params = dict()

    params.update(load_time_params(param))
    params.update(load_domain_params(param))
    params.update(load_initial_elevation(morph))
    params.update(load_dune_params(param))
    params.update(load_alongshore_transport_params(param))
    params.update(load_storm_params(param))
    params.update(load_overwash_params(param))
    params.update(load_shoreface_dynamics_params(param))
    params.update(load_shrub_parameters(param))

    return process_params(params)


def process_params(params):
    params["TMAX"] = int(params["TMAX"]) + 1
    params["StormStart"] = int(params["StormStart"])

    params["LShoreface"] /= 10.0
    params["DShoreface"] /= 10.0
    params["BayDepth"] /= 10.0
    params["MHW"] /= 10.0

    # intElev = load_initial_elevation(morph)

    params["BarrierLength"] = int(params["BarrierLength"] / 10.0)
    params["DuneWidth"] = int(params["DuneWidth"] / 10.0)

    z = params["InteriorDomain"] / 10.0 - params["MHW"]
    z[z <= 0.0] = -params["BayDepth"]
    # above_water = np.any(z > 0, axis=0)
    # z = z[:, above_water]
    u = 1
    while u == 1:  # Remove all rows that have zero subaerial cells
        if all(z[:, -1] <= 0):
            z = z[:, :-1]
        else:
            u = 0
    # while np.all(z[:, -1] <= 0):
    #     z = z[:, :-1]
    params["InteriorDomain"] = np.flipud(np.rot90(z))

    if len(params["InteriorDomain"][0]) > params["BarrierLength"]:
        # Reduce to specified max length
        params["InteriorDomain"] = params["InteriorDomain"][
            :, : params["BarrierLength"]
        ]
    else:
        params["BarrierLength"] = len(params["InteriorDomain"][0])

    # intElev = intElev / 10.0 - params["MHW"]
    # intElev[intElev <= 0.0] = - params["BayDepth"]

    # NOTE: This removes ALL columns that have zero subaerial cells
    # intElev = intElev[:, np.all(intElev <= 0, axis=0)]

    # Flip to correct orientation
    # params["InteriorDomain"] = np.flipud(np.rot90(intElev))

    params["DomainWidth"] = len(params["InteriorDomain"])

    params["Dstart"] /= 10.0
    params["BermEl"] = params["BermEl"] / 10.0 - params["MHW"]

    # Initialize dune crest height domain
    params["DuneDomain"] = np.zeros(
        [params["TMAX"], params["BarrierLength"], params["DuneWidth"]]
    )
    params["DuneDomain"][0, :, 0] = np.ones([1, params["BarrierLength"]]) * (
        params["Dstart"]
        + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1, params["BarrierLength"]))
    )
    params["DuneDomain"][0, :, 1:] = params["DuneDomain"][0, :, 0, None]
    # for col in range(1, DuneWidth):
    #     DuneDomain[0, :, col] = DuneDomain[0, :, 0]

    params["growthparam"] = params["rmin"] + (
        params["rmax"] - params["rmin"]
    ) * np.random.rand(1, params["BarrierLength"])
    params["HdDiffu"] /= 10.0

    params.pop("rmin")
    params.pop("rmax")

    # Maximum dune height
    params["Dmaxel"] /= 10.0

    # Erosion parameters
    params["DuneRestart"] /= 10.0

    # Volume of sediment lost via alongshore transport
    params["Rat"] = params["Rat"] / -10.0
    params["Qat"] = params["Rat"] * params["DShoreface"]
    params.pop("Rat")

    # Relative Sea Level Rise Rate
    params["RSLR"] /= 10.0

    params["surge_tide_m"] -= params["MHW"] * 10.0

    params["Dshrub"] /= 10.0
    params["ShrubEl_min"] = params["ShrubEl_min"] / 10.0 - params["MHW"]
    params["ShrubEl_max"] = params["ShrubEl_max"] / 10.0 - params["MHW"]

    # Overwash Interaction
    params["BurialLimit"] /= 10.0
    params["UprootLimit"] /= 10.0

    # Percent cover change (years 0-9)
    PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
    addend = np.ones(params["TMAX"] + 50)  # (years 10+)
    params["PC"] = np.append(PC, addend)

    return params


def load_initial_elevation(morph):
    # Elevation (decameters)
    x_ind = np.asarray(morph.col_values(0, 1, morph.nrows), dtype=int)
    y_ind = np.asarray(morph.col_values(1, 1, morph.nrows), dtype=int)
    z = np.asarray(morph.col_values(2, 1, morph.nrows), dtype=float)

    shape = (int(x_ind.max() + 1), int(y_ind.max() + 1))
    elevation = np.zeros(shape)
    elevation[(x_ind, y_ind)] = z

    return {"InteriorDomain": elevation}


def from_xls(path_to_xls):
    """Read Barrier3D parameters from an Excel file.

    Parameters
    ----------
    path_to_xls : str
        Path to the Excel paramter file.

    Returns
    -------
    dict
        Dictionary that contains the parameters.
    """
    params = dict()

    # Set-up input files
    # loc = "Barrier3D_Input.xlsx"  # Input File Path

    # Open Excel Files
    wb = xlrd.open_workbook(path_to_xls)
    param = wb.sheet_by_index(0)  # Parameters sheet
    morph = wb.sheet_by_index(1)  # Initial morphology sheet

    # Time Parameters
    TMAX = int(param.cell_value(4, 3)) + 1
    StormStart = int(param.cell_value(5, 3))

    params["TMAX"] = TMAX
    params["StormStart"] = StormStart

    # Computational domain

    # Vertical Dimensions
    LShoreface = param.cell_value(12, 3) / 10
    DShoreface = param.cell_value(13, 3) / 10
    BayDepth = param.cell_value(14, 3) / 10
    MHW = (
        param.cell_value(15, 3) / 10
    )  # Used as offset to convert given elevations relative to a MHW of 0

    params["LShoreface"] = LShoreface
    params["DShoreface"] = DShoreface
    params["BayDepth"] = BayDepth
    params["MHW"] = MHW

    # Elevation (decameters)
    X = []  # Import X coordinates
    for k in range(1, morph.nrows):
        X.append(morph.row_values(k)[0])
    Y = []  # Import Y coordinates
    for k in range(1, morph.nrows):
        Y.append(morph.row_values(k)[1])
    elev = []  # Import elevation values
    for k in range(1, morph.nrows):
        elev.append(morph.row_values(k)[2])

    width = int(max(Y) + 1)
    length = int(max(X) + 1)
    elev_array = np.zeros([length, width])  # Initialize array

    for n in range(len(X)):  # Add elevation values to array
        xx = int(X[n])
        yy = int(Y[n])
        zz = elev[n]
        elev_array[xx, yy] = zz

    intElev = elev_array[:, :] / 10  # Convert to decameters
    intElev = intElev - MHW  # Convert elevations relative to a MHW of 0
    intElev[intElev <= 0] = -BayDepth  # Set all subaerial cells to bay depth
    u = 1
    while u == 1:  # Remove all rows that have zero subaerial cells
        if all(intElev[:, -1] <= 0):
            intElev = intElev[:, :-1]
        else:
            u = 0

    InteriorDomain = np.flipud(np.rot90(intElev))  # Flip to correct orientation

    # Horizontal Dimensions
    BarrierLength = int(param.cell_value(9, 3) / 10)
    if len(InteriorDomain[0]) > BarrierLength:
        InteriorDomain = InteriorDomain[
            :, 0:BarrierLength
        ]  # Reduce to specified max length
    else:
        BarrierLength = len(InteriorDomain[0])
    DomainWidth = len(InteriorDomain)
    DuneWidth = int(param.cell_value(11, 3) / 10)

    params["BarrierLength"] = BarrierLength
    params["InteriorDomain"] = InteriorDomain
    params["DomainWidth"] = DomainWidth
    params["DuneWidth"] = DuneWidth

    # Dunes

    # Dune height refers to heigh of dune above the static berm elevation
    Dstart = param.cell_value(19, 3) / 10
    BermEl = param.cell_value(20, 3) / 10 - MHW

    # Initialize dune crest height domain
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0, :, 0] = np.ones([1, BarrierLength]) * (
        Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1, BarrierLength))
    )
    for w in range(1, DuneWidth):
        DuneDomain[0, :, w] = DuneDomain[0, :, 0]

    # Dune growth parameter
    rmin = param.cell_value(21, 3)
    rmax = param.cell_value(22, 3)
    growthparam = rmin + (rmax - rmin) * np.random.rand(1, BarrierLength)

    # Dune diffusion parameter
    HdDiffu = param.cell_value(23, 3) / 10

    # Maximum dune height
    Dmaxel = param.cell_value(24, 3) / 10

    # Erosion parameters
    C1 = param.cell_value(25, 3)
    C2 = param.cell_value(26, 3)
    DuneRestart = param.cell_value(27, 3) / 10

    params["Dstart"] = Dstart
    params["BermEl"] = BermEl
    params["DuneDomain"] = DuneDomain
    params["growthparam"] = growthparam
    params["HdDiffu"] = HdDiffu
    params["Dmaxel"] = Dmaxel
    params["C1"] = C1
    params["C2"] = C2
    params["DuneRestart"] = DuneRestart

    # Alongshore transport & RSLR

    # Volume of sediment lost via alongshore transport
    Rat = (param.cell_value(31, 3)) * (-1) / 10  # dam
    Qat = Rat * DShoreface  # dam^3/dam

    # Relative Sea Level Rise Rate
    RSLR = param.cell_value(32, 3) / 10

    params["Qat"] = Qat
    params["RSLR"] = RSLR

    # Storm

    # Number Per Year
    mean_storm = param.cell_value(36, 3)
    SD_storm = param.cell_value(37, 3)
    numstorm = param.cell_value(38, 3)

    # Water Level Forcing
    # In meters - converted to decameters after water level calculations
    surge_tide_m = param.cell_value(39, 3) - (MHW * 10)
    surge_tide_sd = param.cell_value(40, 3)
    # height_m = param.cell_value(43,3)
    # height_sd = param.cell_value(44,3)
    height_mu = param.cell_value(41, 3)
    height_sigma = param.cell_value(42, 3)
    period_m = param.cell_value(43, 3)
    period_sd = param.cell_value(44, 3)
    duration_mu = param.cell_value(45, 3)
    duration_sigma = param.cell_value(46, 3)
    beta = param.cell_value(47, 3)

    params["mean_storm"] = mean_storm
    params["SD_storm"] = SD_storm
    params["numstorm"] = numstorm
    params["surge_tide_m"] = surge_tide_m
    params["surge_tide_sd"] = surge_tide_sd
    params["height_mu"] = height_mu
    params["height_sigma"] = height_sigma
    params["period_m"] = period_m
    params["period_sd"] = period_sd
    params["duration_mu"] = duration_mu
    params["duration_sigma"] = duration_sigma
    params["beta"] = beta

    # Overwash

    # Flow Routing
    nn = param.cell_value(51, 3)
    mm = param.cell_value(52, 3)
    Rin_r = param.cell_value(53, 3)
    Rin_i = param.cell_value(54, 3)
    Qs_min = param.cell_value(55, 3)
    threshold_in = param.cell_value(56, 3)

    # Sediment Transport
    Kr = param.cell_value(57, 3)
    Ki = param.cell_value(58, 3)
    Cbb_r = param.cell_value(59, 3)
    Cbb_i = param.cell_value(60, 3)

    params["nn"] = nn
    params["mm"] = mm
    params["Rin_r"] = Rin_r
    params["Rin_i"] = Rin_i
    params["Qs_min"] = Qs_min
    params["threshold_in"] = threshold_in
    params["Kr"] = Kr
    params["Ki"] = Ki
    params["Cbb_r"] = Cbb_r
    params["Cbb_i"] = Cbb_i

    # Shoreface dynamics

    k_sf = param.cell_value(64, 3)
    s_sf_eq = param.cell_value(65, 3)

    params["k_sf"] = k_sf
    params["s_sf_eq"] = s_sf_eq

    # Shrubs

    # Dispersal
    Shrub_ON = param.cell_value(69, 3)
    Seedmin = param.cell_value(70, 3)
    Seedmax = param.cell_value(71, 3)
    disp_mu = param.cell_value(72, 3)
    disp_sigma = param.cell_value(73, 3)

    # Growth
    Dshrub = param.cell_value(74, 3) / 10
    GermRate = param.cell_value(75, 3)
    TimeFruit = param.cell_value(76, 3)
    Female = param.cell_value(77, 3)
    ShrubEl_min = param.cell_value(78, 3) / 10 - MHW
    ShrubEl_max = param.cell_value(79, 3) / 10 - MHW

    # Overwash Interaction
    BurialLimit = param.cell_value(80, 3) / 10
    UprootLimit = param.cell_value(81, 3) / 10
    SalineLimit = param.cell_value(82, 3)
    Qshrub_max = param.cell_value(83, 3)

    # Percent cover change (years 0-9)
    PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
    addend = np.ones(TMAX + 50)  # (years 10+)
    PC = np.append(PC, addend)

    params["Shrub_ON"] = Shrub_ON
    params["Seedmin"] = Seedmin
    params["Seedmax"] = Seedmax
    params["disp_mu"] = disp_mu
    params["disp_sigma"] = disp_sigma
    params["Dshrub"] = Dshrub
    params["GermRate"] = GermRate
    params["TimeFruit"] = TimeFruit
    params["Female"] = Female
    params["ShrubEl_min"] = ShrubEl_min
    params["ShrubEl_max"] = ShrubEl_max
    params["BurialLimit"] = BurialLimit
    params["UprootLimit"] = UprootLimit
    params["SalineLimit"] = SalineLimit
    params["Qshrub_max"] = Qshrub_max
    params["PC"] = PC

    return params
