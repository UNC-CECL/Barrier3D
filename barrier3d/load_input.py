import contextlib
import math
import os
import pathlib

import numpy as np
import pandas

from .configuration import Barrier3dConfiguration


@contextlib.contextmanager
def as_cwd(path):
    prev_cwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(prev_cwd)


def load_inputs(path_to_folder, prefix="barrier3d", fmt="yaml"):
    """Load input parameter and elevation data for Barrier3d.

    Parameters
    ----------
    path_to_folder : str, or Path
        Path to folder that contains barrier3d input files.
    prefix : str, optional
        Filename prefix for input files.
    fmt : {'py', 'yaml'}, optional
        Input file format.

    Returns
    -------
    dict
        Processes input parameters for a barrier3d simulation.
    """
    fmts = ("py", "yaml")
    if fmt not in fmts:
        raise ValueError(
            "unrecognized format for input files ({fmt} not in {fmts})".format(
                fmt=fmt, fmts=", ".join(fmts)
            )
        )
    parameter_file = pathlib.Path(f"{prefix}-parameters.{fmt}")

    with as_cwd(path_to_folder):
        if not parameter_file.is_file():
            raise FileNotFoundError(parameter_file.absolute())

        params = load_parameters(parameter_file, fmt=fmt)

        params["InteriorDomain"] = load_elevation(
            params.pop("elevation_file"), fmt=None
        )  # interior elevations must come from a time series

        params["StormSeries"] = load_storms(
            params.pop("storm_file"), fmt=None
        )  # storms must come from a time series

        if params[
            "DuneParamStart"
        ]:  # user specifies if dune height will come from external file
            params["DuneStart"] = load_dunes(params.pop("dune_file"), fmt=None)
        else:
            params.pop("dune_file")  # get rid of the keyword

        if params[
            "GrowthParamStart"
        ]:  # user specifies if growth parameters will come from external file
            params["GrowthStart"] = load_growth_param(
                params.pop("growth_param_file"), fmt=None
            )
        else:
            params.pop("growth_param_file")

        _process_raw_input(params)

    return params


def load_parameters(path_to_file, fmt="yaml"):
    try:
        loader = getattr(Barrier3dConfiguration, f"from_{fmt}")
    except AttributeError:
        raise ValueError(f"format not understood ({fmt})")
    else:
        return loader(path_to_file).data


def load_elevation(path_to_file, fmt="npy"):
    fmt = fmt or _guess_format(path_to_file)

    # load_elevation(path_to_file, fmt="npy"):
    if fmt == "npy":
        elevation = np.load(path_to_file, allow_pickle=True)
    elif fmt == "csv":
        elevation = np.loadtxt(path_to_file, delimiter=",")
        # data = pandas.read_csv(
        #     path_to_csv,
        #     names=("x", "y", "z"),
        #     dtype={"x": int, "y": int, "z": float},
        #     comment="#",
        #     header=0,
        # )
        # elevation = np.zeros((data.x.max() + 1, data.y.max() + 1), dtype=float)
        # elevation[(data.x, data.y)] = data.z
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(
            f"unrecognized format for storm time series ({fmt} not one of {fmts})"
        )

    return elevation


def _guess_format(path_to_file):
    """Guess the format of a file based on its name.

    Parameters
    ----------
    path_to_file : str, path-like
        Name of the file.

    Returns
    -------
    str
        The format of the file or an empty string if it can't be determined.

    Examples
    --------
    >>> from barrier3d.load_input import _guess_format
    >>> _guess_format("file.csv")
    'csv'
    >>> _guess_format("file.npy")
    'npy'
    >>> _guess_format("file-with-an-extension")
    ''
    """
    return pathlib.Path(path_to_file).suffix[1:]


def load_storms(path_to_file, fmt="npy"):
    fmt = fmt or _guess_format(path_to_file)

    if fmt == "npy":
        data = np.load(path_to_file)
    elif fmt == "csv":
        df = pandas.read_csv(
            path_to_file,
            names=("time", "Rhigh", "Rlow", "period", "duration"),
            dtype={
                "time": int,
                "Rhigh": float,
                "Rlow": float,
                "period": float,
                "duration": int,
            },
            comment="#",
            header=0,
        )
        df.to_numpy()
        data = np.hstack(
            (
                df["time"][:, None],
                df["Rhigh"][:, None],
                df["Rlow"][:, None],
                df["period"][:, None],
                df["duration"][:, None],
            )
        )
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(
            f"unrecognized format for storm time series ({fmt} not one of {fmts})"
        )

    return data


def load_dunes(path_to_file, fmt="npy"):
    fmt = fmt or _guess_format(path_to_file)

    if fmt == "npy":
        data = np.load(path_to_file)
    elif fmt == "csv":
        data = pandas.read_csv(
            path_to_file, names=("DuneStart",), dtype=float, comment="#", header=0
        )["DuneStart"].values
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(f"unrecognized format for dunes ({fmt} not one of {fmts})")

    return data


def load_growth_param(path_to_file, fmt="npy"):
    fmt = fmt or _guess_format(path_to_file)

    if fmt == "npy":
        data = np.load(path_to_file)
    elif fmt == "csv":
        data = pandas.read_csv(
            path_to_file, names=("growth_param",), dtype=float, comment="#", header=0
        )["growth_param"].values
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(
            f"unrecognized format for growth param ({fmt} not one of {fmts})"
        )

    return data


def _process_raw_input(params):
    params["RNG"] = np.random.default_rng(seed=1973)

    params["TMAX"] = int(params["TMAX"])
    params["StormStart"] = int(params["StormStart"])

    params.setdefault("StormSeries", [])

    params["LShoreface"] /= 10.0
    params["DShoreface"] /= 10.0
    params["BayDepth"] /= 10.0
    params["MHW"] /= 10.0

    # intElev = load_initial_elevation(morph)

    params["BarrierLength"] = int(params["BarrierLength"] / 10.0)
    params["DuneWidth"] = int(params["DuneWidth"] / 10.0)

    if params["InteriorDomain"].shape[1] > params["BarrierLength"]:
        params["InteriorDomain"] = params["InteriorDomain"][
            :, : params["BarrierLength"]
        ]
    else:
        params["BarrierLength"] = params["InteriorDomain"].shape[1]

    params["DomainWidth"] = len(params["InteriorDomain"])

    params["Dstart"] /= 10.0
    params["BermEl"] = params["BermEl"] / 10.0 - params["MHW"]

    # Initialize dune crest height domain
    params["DuneDomain"] = np.zeros(
        (params["TMAX"], params["BarrierLength"], params["DuneWidth"])
    )

    if "DuneStart" in params:
        params["DuneDomain"][0, :, 0] = params["DuneStart"][0 : params["BarrierLength"]]
        params.pop("DuneStart")
    else:
        if params["SeededRNG"]:
            params["DuneDomain"][0, :, 0] = np.ones([1, params["BarrierLength"]]) * (
                params["Dstart"]
                + (
                    -0.01
                    + (0.01 - (-0.01))
                    * params["RNG"].random((1, params["BarrierLength"]))
                )
            )
        else:
            params["DuneDomain"][0, :, 0] = np.ones([1, params["BarrierLength"]]) * (
                params["Dstart"]
                + (
                    -0.01
                    + (
                        -0.01
                        + (0.01 - (-0.01)) * np.random.rand(1, params["BarrierLength"])
                    )
                )
            )
    params["DuneDomain"][0, :, 1:] = params["DuneDomain"][0, :, 0, None]

    if "GrowthStart" in params:
        params["growthparam"] = params["GrowthStart"][0 : params["BarrierLength"]]
        params.pop("GrowthStart")
    else:
        if params["SeededRNG"]:
            params["growthparam"] = params["rmin"] + (
                params["rmax"] - params["rmin"]
            ) * params["RNG"].random((1, params["BarrierLength"]))
        else:
            params["growthparam"] = params["rmin"] + (
                params["rmax"] - params["rmin"]
            ) * np.random.rand(1, params["BarrierLength"])

    params["HdDiffu"] /= 10.0

    # params.pop("rmin")
    # params.pop("rmax")

    # Maximum dune height
    # params["Dmaxel"] /= 10.0
    params["Dmaxel"] = params["Dmaxel"] / 10.0 - params["MHW"]

    # Erosion parameters
    params["DuneRestart"] /= 10.0

    # Volume of sediment lost via alongshore transport
    params["Rat"] = params["Rat"] / 10.0
    # params["Rat"] = params["Rat"] / -10.0
    params["Qat"] = params["Rat"] * params["DShoreface"]
    params.pop("Rat")

    # Relative Sea Level Rise Rate
    # params["RSLR"] /= 10.0
    if params["RSLR_Constant"]:
        # Constant RSLR
        params["RSLR_const"] /= 10.0
        params["RSLR"] = [params["RSLR_const"]] * params["TMAX"]
    else:
        # Logistic RSLR rate projection - Rohling et al. (2013)
        params["RSLR"] = []
        alpha = 0.75  # m/yr -- probability maximum = 0.75, 68% upper bound = 2.0
        beta = alpha / 0.003 - 1  # constant
        gamma = 350  # yr -- probability maximum = 350, 68% upper bound = 900
        C = 12  # constant
        for t in range(150, params["TMAX"] + 150):
            delta = (
                alpha / (1 + beta * math.exp(-t / gamma * C)) / 10000 * 10
            )  # Convert from m/cy to dam/yr
            params["RSLR"].append(delta)

    # Shoreface
    params["x_t"] = params["ShorefaceToe"] / 10.0
    params.pop("ShorefaceToe")

    # Shrubs
    params["Dshrub"] /= 10.0
    params["ShrubEl_min"] = params["ShrubEl_min"] / 10.0 - params["MHW"]
    params["ShrubEl_max"] = params["ShrubEl_max"] / 10.0 - params["MHW"]
    # params["TideAmp"] /= 10.0
    params["SprayDist"] /= 10.0
    params["MaxShrubHeight"] /= 10.0

    # # Shrub height (with age as proxy)
    # SH = np.linspace(0, 0.3, 10)
    # addend = np.ones(params["TMAX"] + 50)  # (years 10+)
    # params["SH"] = np.append(SH, addend)

    # Overwash Interaction (convert to decameters)
    params["BurialLimit"] /= 10.0
    params["UprootLimit"] /= 10.0
    params["Qs_min"] /= 1000
    params["Qs_bb_min"] /= 1000

    # Percent cover change (years 0-9)
    PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
    addend = np.ones(params["TMAX"] + 50)  # (years 10+)
    params["PC"] = np.append(PC, addend)

    return params
