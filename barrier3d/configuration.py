import pathlib

import numpy as np
import pandas
import xlrd
from exconfig import (
    ArrayField,
    BooleanField,
    Configuration,
    Field,
    FloatField,
    IntegerField,
)
from exconfig.validators import Path, Range


class Barrier3dConfiguration(Configuration):
    elevation_file = Field(
        "elevation_file",
        default="barrier3d-elevations.npy",
        description="File that contains initial elevations in [m MHH]",
        validators=[Path(file_okay=True, dir_okay=False)],
    )
    dune_file = Field(
        "dune_file",
        default="barrier3d-dunes.npy",
        description="File that contains initial dune height values [m]",
        validators=[Path(file_okay=True, dir_okay=False)],
    )
    growth_param_file = Field(
        "growth_param_file",
        default="barrier3d-growthparam.npy",
        description="File that contains initial growth parameters",
        validators=[Path(file_okay=True, dir_okay=False)],
    )
    storm_file = Field(
        "storm_file",
        default="barrier3d-storms.npy",
        description="File that contains storm data",
        validators=[Path(file_okay=True, dir_okay=False)],
    )
    TMAX = IntegerField(
        "TMAX",
        default=150,
        units="y",
        description="Duration of simulation",
        validators=[Range(lower=0)],
    )
    StormStart = IntegerField(
        "StormStart",
        default=2,
        units="y",
        description="Year when storm can start occurring (NOTE: if changed, need new storm time series)",
        validators=[Range(lower=0)],
    )
    BarrierLength = FloatField(
        "BarrierLength",
        default=500.0,
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
        default=500.0,
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
        default=3.0,
        units="m",
        description="Depth of bay benind island segment",
        validators=[Range(lower=0)],
    )
    MHW = FloatField(
        "MHW",
        default=0.46,
        units="m NAVD88",
        description="Elevation of Mean High Water (NOTE: if changed, need new storm time series)",
        validators=[Range(lower=0)],
    )
    Dstart = FloatField(
        "Dstart",
        default=0.50,
        units="m",
        description="Initial height of dune domain above berm elevation",
        validators=[Range(lower=0)],
    )
    BermEl = FloatField(
        "BermEl",
        default=1.9,
        units="m NAVD88",
        description="Static elevation of berm; berm elevation + dune height = dune elevation (NOTE: if changed, need new MSSM and storms)",
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
        default=0.75,
        units="m",
        description="Dune diffusion parameter (i.e. max height offset between adjacent dune cells)",
        validators=[Range(lower=0)],
    )
    Dmaxel = FloatField(
        "Dmaxel",
        default=3.4,
        units="m NAVD88",
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
        default=0.075,
        units="m",
        description="Restart height for dunes lowered to essentially zero",
        validators=[Range(lower=0)],
    )
    Rat = FloatField(
        "Rat",
        default=0,
        units="m / y",
        description="Rate of shoreline retreat attributed to alongshore transport; (-) = erosion, (+) = accretion",
    )
    RSLR_Constant = BooleanField(
        "RSLR_Constant",
        default=True,
        description="Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series",
    )
    RSLR_const = FloatField(
        "RSLR_const",
        default=0.004,
        units="m / y",
        description="Relative sea-level rise rate",
    )
    beta = FloatField(
        "beta",
        default=0.04,
        description="Beach slope for runup calculations",
        validators=[Range(lower=0)],
    )
    #    StormTimeSeries = BooleanField(
    #        "StormTimeSeries",
    #        default=True,
    #        description="Storms will come from a time series",
    #    )
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
        default=2,
        description="Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)",
    )
    Rin_i = FloatField(
        "Rin_i",
        default=0.25,
        description="Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)",
    )
    Qs_min = FloatField(
        "Qs_min",
        default=1.0,
        units="m^3 / hr",
        description="Minimum discharge needed for sediment transport",
    )
    MaxUpSlope = FloatField(
        "MaxUpSlope",
        default=0.25,
        units="m / m",
        description="Maximum slope water can flow upward",
    )
    threshold_in = FloatField(
        "threshold_in",
        default=0.25,
        units="m^3 / hr",
        description="Threshold to determine if in inundation regime",
    )
    Kr = FloatField(
        "Kr",
        default=0.000075,
        description="Sediment flux constant, run-up regime",
    )
    Ki = FloatField(
        "Ki",
        default=0.0000075,
        description="Sediment flux constant, inundation regime",
    )
    Cbb_r = FloatField(
        "Cbb_r",
        default=0.5,
        units="",
        description="Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime",
        validators=[Range(lower=0, upper=1)],
    )
    Cbb_i = FloatField(
        "Cbb_i",
        default=0.8,
        units="",
        description="Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime",
        validators=[Range(lower=0, upper=1)],
    )
    Qs_bb_min = FloatField(
        "Qs_bb_min",
        default=1,
        units="m^3 / hr",
        description="Minimum sediment flux in back-barrier bay (below which sediment won't flux)",
        validators=[Range(lower=0)],
    )
    Cx = FloatField(
        "Cx",
        default=10.0,
        units="",
        description="Multiplier with the average slope of the interior for constant C in inundation transport rule",
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
        default=1,
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
        default=0.02,
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
        default=1000,
        units="1 / yr",
        description="Seeds produced per shrub per year (fecundity)",
        validators=[Range(lower=0)],
    )
    Seedmax = FloatField(
        "Seedmax",
        default=10000,
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
        default=2.75,
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
        units="m NAVD88",
        description="Elevation range for shrub growth, minimum bound",
        validators=[Range(lower=0)],
    )
    ShrubEl_max = FloatField(
        "ShrubEl_max",
        default=2.3,
        units="m NAVD88",
        description="Elevation range for shrub growth, maximum bound",
        validators=[Range(lower=0)],
    )
    # TideAmp = FloatField(
    #     "TideAmp",
    #     default=1.2,
    #     units="m",
    #     description="Tidal amplitude",
    #     validators=[Range(lower=0)],
    # )
    SprayDist = FloatField(
        "SprayDist",
        default=170,
        units="m",
        description="Distance from ocean shoreline that shrubs can establish",
        validators=[Range(lower=0)],
    )
    BurialLimit = FloatField(
        "BurialLimit",
        default=0.75,
        units="m",
        description="Maximum percentage of height that a shrub can be buried up to before dying",
        validators=[Range(lower=0)],
    )
    UprootLimit = FloatField(
        "UprootLimit",
        default=-0.2,
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
    DuneParamStart = BooleanField(
        "DuneParamStart",
        default=True,
        description="Dune height will come from external file",
    )
    GrowthParamStart = BooleanField(
        "GrowthParamStart",
        default=True,
        description="Dune growth parameters will come from external file",
    )
    MaxShrubHeight = FloatField(
        "MaxShrubHeight",
        default=5.2,
        units="m",
        description="Maximum shrub height",
        validators=[Range(lower=0)],
    )
    ShorefaceToe = FloatField(
        "ShorefaceToe",
        default=0,
        units="m",
        description="Start location of shoreface toe ",
    )
    SeededRNG = BooleanField(
        "SeededRNG",
        default=True,
        description="Use seeded random number generator for reproducibility",
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
        comment="#",
        header=0,
    )
