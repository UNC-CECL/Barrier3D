import pathlib

import numpy as np
import pandas
import xlrd

from exconfig import Configuration, ArrayField, BooleanField, FloatField, IntegerField
from exconfig.validators import Range


class Barrier3dConfiguration(Configuration):
    TMAX = IntegerField(
        "TMAX",
        default=50,
        units="y",
        description="Duration of simulation",
        validators=[Range(lower=0)],
    )
    StormStart = IntegerField(
        "StormStart",
        default=2,
        units="y",
        description="Year when storm can start occurring",
        validators=[Range(lower=0)],
    )
    BarrierLength = FloatField(
        "BarrierLength",
        default=1000.0,
        units="m",
        description="Static length (alongshore) of island segment",
        validators=[Range(lower=0)],
    )
    DuneWidth = FloatField(
        "DuneWidth",
        default=20.0,
        units="m",
        description="Width (cross-shore) of island dune field; for illustration purposes only",
        validators=[Range(lower=0)],
    )
    LShoreface = FloatField(
        "LShoreface",
        default=1000.0,
        units="m",
        description="Length of shoreface",
        validators=[Range(lower=0)],
    )
    DShoreface = FloatField(
        "DShoreface",
        default=10.0,
        units="m",
        description="Height of shoreface",
        validators=[Range(lower=0)],
    )
    BayDepth = FloatField(
        "BayDepth",
        default=1.5,
        units="m",
        description="Depth of bay benind island segment",
        validators=[Range(lower=0)],
    )
    MHW = FloatField(
        "MHW",
        default=0.46,
        units="m",
        description="Elevation of Mean High Water",
        validators=[Range(lower=0)],
    )
    Dstart = FloatField(
        "Dstart",
        default=0.25,
        units="m",
        description="Initial height of dune domain above berm elevation",
        validators=[Range(lower=0)],
    )
    BermEl = FloatField(
        "BermEl",
        default=1.7,
        units="m",
        description="Static elevation of berm; berm elevation + dune height = dune elevation",
        validators=[Range(lower=0)],
    )
    rmin = FloatField(
        "rmin",
        default=0.35,
        description="Minimum growth rate for logistic dune growth",
    )
    rmax = FloatField(
        "rmax",
        default=0.85,
        description="Maximum growth rate for logistic dune growth",
    )
    HdDiffu = FloatField(
        "HdDiffu",
        default=0.8,
        units="m",
        description="Dune diffusion parameter (i.e. max height offset between adjacent dune cells)",
        validators=[Range(lower=0)],
    )
    Dmaxel = FloatField(
        "Dmaxel",
        default=3.0,
        units="m",
        description="Maximum elevation of dunes",
        validators=[Range(lower=0)],
    )
    C1 = FloatField(
        "C1",
        default=8.8,
        units="m",
        description="Empirical dune erosion parameter",
    )
    C2 = FloatField(
        "C2",
        default=4.6,
        units="m",
        description="Empirical dune erosion parameter",
    )
    DuneRestart = FloatField(
        "DuneRestart",
        default=0.05,
        units="m",
        description="Restart height for dunes lowered to essentially zero",
        validators=[Range(lower=0)],
    )
    Rat = FloatField(
        "Rat",
        default=-1,
        units="m / y",
        description="Rate of shoreline reatreat attributed to alongshore transport; (-) = erosion, (+) = accretion",
    )
    RSLR = FloatField(
        "RSLR",
        default=0.004,
        units="m / y",
        description="Relative sea-level rise rate",
    )
    mean_storm = FloatField(
        "mean_storm",
        default=13.17,
        description="For a random number of storms per year sampled from normal distribution",
        validators=[Range(lower=0)],
    )
    SD_storm = FloatField(
        "SD_storm",
        default=5.16,
        description="For a random number of storms per year sampled from normal distribution",
        validators=[Range(lower=0)],
    )
    numstorm = IntegerField(
        "numstorm",
        default=0,
        description="For a single constant number of storms per year",
        validators=[Range(lower=0)],
    )
    surge_tide_m = FloatField(
        "surge_tide_m",
        default=0.91,
        units="m",
        description="Mean storm surge + tide water levels",
        validators=[Range(lower=0)],
    )
    surge_tide_sd = FloatField(
        "surge_tide_sd",
        default=0.3,
        units="m",
        description="Standard deviation of surge + tide water levels",
        validators=[Range(lower=0)],
    )
    height_mu = FloatField(
        "height_mu",
        default=1.02,
        units="m",
        description="Mu storm wave height (lognormal distribution)",
        validators=[Range(lower=0)],
    )
    height_sigma = FloatField(
        "height_sigma",
        default=0.19,
        units="m",
        description="Sigma storm wave height (lognormal distribution)",
        validators=[Range(lower=0)],
    )
    period_m = FloatField(
        "period_m",
        default=10.74,
        description="Mean storm wave period",
        validators=[Range(lower=0)],
    )
    period_sd = FloatField(
        "period_sd",
        default=2.76,
        description="Standard deviation of storm wave periods",
        validators=[Range(lower=0)],
    )
    duration_mu = FloatField(
        "duration_mu",
        default=2.9249,
        units="hr",
        description="Mu storm duration (lognormal distribution)",
        validators=[Range(lower=0)],
    )
    duration_sigma = FloatField(
        "duration_sigma",
        default=0.638673,
        units="hr",
        description="Sigma storm duration (lognormal distribution)",
        validators=[Range(lower=0)],
    )
    beta = FloatField(
        "beta",
        default=0.04,
        description="Beach slope for runup calculations",
        validators=[Range(lower=0)],
    )
    StormTimeSeries = BooleanField(
        "StormTimeSeries",
        default=0,
        description="Storms will come from a time series",
    )
    StormSeries = ArrayField(
        "StormSeries",
        default=(),
        description="Time series of storms",
    )
    nn = FloatField(
        "nn",
        default=0.5,
        description="Flow routing constant",
    )
    mm = FloatField(
        "mm",
        default=2.0,
        description="Exponent constant for sediment transport",
    )
    Rin_r = FloatField(
        "Rin_r",
        default=6,
        description="Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)",
    )
    Rin_i = FloatField(
        "Rin_i",
        default=1,
        description="Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)",
    )
    Qs_min = FloatField(
        "Qs_min",
        default=0.05,
        units="m^3 / hr",
        description="Minimum discharge needed for sediment transport",
    )
    threshold_in = FloatField(
        "threshold_in",
        default=0.25,
        units="m^3 / hr",
        description="Threshold to determine if in inundation regime",
    )
    Kr = FloatField(
        "Kr",
        default=0.0001,
        description="Sediment flux constant, run-up regime",
    )
    Ki = FloatField(
        "Ki",
        default=0.0004,
        description="Sediment flux constant, inundation regime",
    )
    Cbb_r = FloatField(
        "Cbb_r",
        default=0.65,
        units="",
        description="Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime",
        validators=[Range(lower=0, upper=1)],
    )
    Cbb_i = FloatField(
        "Cbb_i",
        default=0.9,
        units="",
        description="Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime",
        validators=[Range(lower=0, upper=1)],
    )
    Qs_bb_min = FloatField(
        "Qs_bb_min",
        default=0.0001,
        units="",
        description="",
        validators=[Range(lower=0)],
    )
    Cx = FloatField(
        "Cx",
        default=2.0,
        units="",
        description="",
        validators=[Range(lower=0)],
    )
    OWss_i = IntegerField(
        "OWss_i",
        default=2,
        description="Overwash substep",
        validators=[Range(lower=1)],
    )
    OWss_r = IntegerField(
        "OWss_r",
        default=2,
        description="Overwash substep",
        validators=[Range(lower=1)],
    )
    k_sf = FloatField(
        "k_sf",
        default=5000,
        units="m^3 / m / y",
        description="Shoreface flux rate constant",
    )
    s_sf_eq = FloatField(
        "s_sf_eq",
        default=0.01,
        units="",
        description="Equilibrium shoreface slope",
    )
    Shrub_ON = BooleanField(
        "Shrub_ON",
        default=False,
        description="1 = shrubs on in simulation, 0 = shrubs off",
    )
    Seedmin = FloatField(
        "Seedmin",
        default=100,
        units="1 / yr",
        description="Seeds produced per shrub per year (fecundity)",
        validators=[Range(lower=0)],
    )
    Seedmax = FloatField(
        "Seedmax",
        default=1000,
        units="1 / yr",
        description="Seeds produced per shrub per year (fecundity)",
        validators=[Range(lower=0)],
    )
    disp_mu = FloatField(
        "disp_mu",
        default=-0.721891,
        description="For lognormal probability distribution of seed dispersal distance",
    )
    disp_sigma = FloatField(
        "disp_sigma",
        default=1.5,
        description="For lognormal probability distribution of seed dispersal distance",
    )
    Dshrub = FloatField(
        "Dshrub",
        default=2,
        units="m",
        description="Minimum elevation of fronting dune for shrub growth",
        validators=[Range(lower=0)],
    )
    GermRate = FloatField(
        "GermRate",
        default=0.6,
        units="",
        description="Germination rate",
        validators=[Range(lower=0, upper=1)],
    )
    TimeFruit = FloatField(
        "TimeFruit",
        default=5,
        units="yr",
        description="Age shrubs need to be before they start fruiting",
        validators=[Range(lower=0)],
    )
    Female = FloatField(
        "Female",
        default=0.5,
        units="",
        description="Percentage of shrubs that are female",
        validators=[Range(lower=0, upper=1)],
    )
    ShrubEl_min = FloatField(
        "ShrubEl_min",
        default=1.2,
        units="m",
        description="Elevation range for shrub growth, minimum bound",
        validators=[Range(lower=0)],
    )
    ShrubEl_max = FloatField(
        "ShrubEl_max",
        default=2.3,
        units="m",
        description="Elevation range for shrub growth, maximum bound",
        validators=[Range(lower=0)],
    )
    BurialLimit = FloatField(
        "BurialLimit",
        default=0.5,
        units="m",
        description="Shrubs buried beyond this limit killed",
        validators=[Range(lower=0)],
    )
    UprootLimit = FloatField(
        "UprootLimit",
        default=-0.3,
        units="m",
        description="Shrubs eroded beyond this limit killed",
        validators=[Range(upper=0)],
    )
    SalineLimit = FloatField(
        "SalineLimit",
        default=0.05,
        units="m^3 / hr",
        description="Dishcharge limit to determine shrub mortality via saline flooding",
        validators=[Range(lower=0)],
    )
    Qshrub_max = FloatField(
        "Qshrub_max",
        default=0.15,
        units="",
        description="Maximum percentage of overwash reduction through a shrub cell with full percent cover",
        validators=[Range(lower=0, upper=1)],
    )

    @classmethod
    def from_py(cls, path_to_py):
        import importlib

        mod = importlib.import_module(str(pathlib.Path(path_to_py).stem))
        return cls(
            data=dict(
                [(k, v) for k, v in mod.__dict__.items() if not k.startswith("_")]
            )
        )

    @staticmethod
    def from_xlsx(path_to_xlsx):
        from .parameters import from_xlsx

        return from_xlsx(path_to_xlsx)


def load_inputs(path_to_folder, prefix="barrier3d", fmt="yaml"):
    """Load input parameter and elevation data for Barrier3d."""
    base_dir = pathlib.Path(path_to_folder)
    if fmt == "yaml":
        path_to_parameters = f"{prefix}-parameters.yaml"
        path_to_elevations = f"{prefix}-elevations.csv"
        # path_to_parameters = "parameters.yaml"
        # path_to_elevations = "elevations.csv"
    elif fmt == "xlsx":
        path_to_parameters = f"{prefix}.xlsx"
        path_to_elevations = f"{prefix}.xlsx"
    elif fmt == "py":
        path_to_parameters = f"{prefix}-parameters.py"
        path_to_elevations = f"{prefix}-elevations.csv"
    else:
        raise ValueError(
            "unrecognized format for input files ({fmt} not in {fmts})".format(
                fmt=fmt, fmts=", ".join(fmts)
            )

        )

    if base_dir.is_file():
        path_to_parameters = base_dir.name
        base_dir = base_dir.parent

    params = load_configuration(base_dir / path_to_parameters, fmt=fmt)
    params["InteriorDomain"] = load_elevation(
        base_dir / path_to_elevations, fmt="csv" if fmt != "xlsx" else "xlsx"
    )

    # if fmt == "yaml":
    #     fmts = ("yaml", "csv")
    # elif fmt == "xlsx":
    #     fmts = ("xlsx", "xlsx")
    # else:
    #     raise ValueError(
    #         "unrecognized format for input files ({fmt} not in {fmts})".format(
    #             fmt=fmt, fmts=", ".join(fmts)
    #         )
    #     )

    # path_to_config = path_to_folder / (prefix + "." + fmts[0])
    # path_to_elevation = path_to_folder / (prefix + "." + fmts[1])

    # params = load_configuration(path_to_config, fmt=fmts[0])
    # params["InteriorDomain"] = load_elevation(path_to_elevation, fmt=fmts[1])

    return _process(params)


def load_configuration(path_to_file, fmt="yaml"):
    if fmt == "yaml":
        return Barrier3dConfiguration.from_yaml(path_to_file).data
    elif fmt == "py":
        return Barrier3dConfiguration.from_py(path_to_file).data
    else:
        raise ValueError(f"format not understood ({fmt})")
    # return Barrier3dConfiguration.from_path(path_to_file, fmt=fmt).data


def load_storms(path_to_file, fmt="csv"):
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
        )
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


def load_dunes(path_to_file, fmt="csv"):
    if fmt == "npy":
        data = np.load(path_to_file)
    elif fmt == "csv":
        data = pandas.read_csv(
            path_to_file, names=("DuneStart",), dtype=float, comment="#",
        )["DuneStart"].values
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(
            f"unrecognized format for dunes ({fmt} not one of {fmts})"
        )

    return data


def load_growth_param(path_to_file, fmt="csv"):
    if fmt == "npy":
        data = np.load(path_to_file)
    elif fmt == "csv":
        data = pandas.read_csv(
            path_to_file, names=("growth_param",), dtype=float, comment="#",
        )["growth_param"].values
    else:
        fmts = ", ".join(["npy", "csv"])
        raise ValueError(
            f"unrecognized format for growth param ({fmt} not one of {fmts})"
        )

    return data


def load_elevation(path_to_file, fmt="csv"):
    loaders = {
        "csv": _load_initial_elevation_csv,
        "xlsx": _load_initial_elevation_xlsx,
    }
    try:
        loader = loaders[fmt]
    except KeyError:
        fmts = ", ".join(loaders.keys())
        raise ValueError(f"unrecognized format for elevation ({fmt} not one of {fmts})")

    def _process(data):
        shape = data.x.max() + 1, data.y.max() + 1
        elevation = np.zeros(shape, dtype=float)
        elevation[(data.x, data.y)] = data.z

        return elevation
        # return np.flipud(np.rot90(elevation))

    return _process(loader(path_to_file))


def _load_initial_elevation_xlsx(path_to_xlsx):
    wb = xlrd.open_workbook(path_to_xlsx)
    morph = wb.sheet_by_index(1)  # Initial morphology sheet

    return pandas.DataFrame(
        {
            "x": np.asarray(morph.col_values(0, 1, morph.nrows), dtype=int),
            "y": np.asarray(morph.col_values(1, 1, morph.nrows), dtype=int),
            "z": np.asarray(morph.col_values(2, 1, morph.nrows), dtype=float),
        }
    )


def _load_initial_elevation_csv(path_to_csv):
    return pandas.read_csv(
        path_to_csv,
        names=("x", "y", "z"),
        dtype={"x": int, "y": int, "z": float},
        comment="#"
    )


def _process(params):
    params["TMAX"] = int(params["TMAX"]) + 1
    params["StormStart"] = int(params["StormStart"])

    if params["StormTimeSeries"]:
        params["StormSeries"] = load_storms("barrier3d-storms.csv")
    else:
        params["StormSeries"] = []

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
    # NOTE: if we haven't already flipped/rotated
    u = 1
    while u == 1:  # Remove all rows that have zero subaerial cells
        if all(z[:, -1] <= 0):
            z = z[:, :-1]
        else:
            u = 0
    # while np.all(z[:, -1] <= 0):
    #     z = z[:, :-1]
    # NOTE: if we have already flipped/rotated
    # u = 1
    # while u == 1:  # Remove all rows that have zero subaerial cells
    #     if all(z[0, :] <= 0):
    #         z = z[1:, :]
    #     else:
    #         u = 0

    # NOTE: We do this in the process step
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
        (params["TMAX"], params["BarrierLength"], params["DuneWidth"])
    )
    if params["StormTimeSeries"]:
        dune_start = load_dunes("barrier3d-dunes.csv")
        params["DuneDomain"][0, :, 0] = dune_start[0:params["BarrierLength"]]    
    else:
        params["DuneDomain"][0, :, 0] = np.ones([1, params["BarrierLength"]]) * (
            params["Dstart"]
            + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1, params["BarrierLength"]))
        )
    params["DuneDomain"][0, :, 1:] = params["DuneDomain"][0, :, 0, None]

    if params["StormTimeSeries"]:
        growth_param_start = load_growth_param("barrier3d-growthparam.csv")
        params["growthparam"] = growth_param_start[0: params["BarrierLength"]]
    else:
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
