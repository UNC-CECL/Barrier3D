import os
import pathlib
import textwrap
from contextlib import suppress

from pydantic import (
    BaseModel,
    Field,
    FilePath,
    NonNegativeFloat,
    NonPositiveFloat,
    PositiveInt,
    confloat,
    conint,
    conlist,
    validator,
)


class Barrier3dConfiguration(BaseModel):
    elevation_file: FilePath = Field(
        "barrier3d-default-elevations.npy",
        description="File that contains initial elevations in [m MHH]",
    )
    dune_file: FilePath = Field(
        "barrier3d-default-dunes.npy",
        description="File that contains initial dune height values [m]",
    )
    growth_param_file: FilePath = Field(
        "barrier3d-default-growthparam.npy",
        description="File that contains initial growth parameters",
    )
    storm_file: FilePath = Field(
        "barrier3d-default-storms.npy", description="File that contains storm data"
    )

    TMAX: PositiveInt = Field(
        150,
        description="Duration of simulation",
        unit="y",
    )
    StormStart: PositiveInt = Field(
        2,
        description=(
            "Year when storm can start occurring (NOTE: if changed, need new storm "
            "time series)"
        ),
        unit="y",
    )
    BarrierLength: NonNegativeFloat = Field(
        500.0,
        description="Static length (alongshore) of island segment",
        unit="m",
    )
    DuneWidth: NonNegativeFloat = Field(
        20.0,
        description=(
            "Width (cross-shore) of island dune field; for illustration purposes only"
        ),
        unit="m",
    )
    LShoreface: NonNegativeFloat = Field(
        500.0,
        description="Length of shoreface",
        unit="m",
    )
    DShoreface: NonNegativeFloat = Field(
        10.0,
        description="Height of shoreface",
        unit="m",
    )
    BayDepth: NonNegativeFloat = Field(
        3.0,
        description="Depth of bay behind island segment",
        unit="m",
    )
    MHW: NonNegativeFloat = Field(
        0.46,
        description=(
            "Elevation of Mean High Water (NOTE: if changed, need new storm time "
            "series)"
        ),
        unit="m NAVD88",
    )
    Dstart: NonNegativeFloat = Field(
        0.5, description="Initial height of dune domain above berm elevation", unit="m"
    )
    BermEl: NonNegativeFloat = Field(
        1.9,
        description=(
            "Static elevation of berm; berm elevation + dune height = dune elevation "
            "(NOTE: if changed, need new MSSM and storms)"
        ),
        unit="m NAVD88",
    )
    rmin: float = Field(
        0.35,
        description="Minimum growth rate for logistic dune growth",
    )
    rmax: float = Field(
        0.85,
        description="Maximum growth rate for logistic dune growth",
    )
    HdDiffu: NonNegativeFloat = Field(
        0.75,
        description=(
            "Dune diffusion parameter (i.e. max height offset between adjacent dune "
            "cells)"
        ),
        unit="m",
    )
    Dmaxel: NonNegativeFloat = Field(
        3.4,
        description="Maximum elevation of dunes",
        unit="m NAVD88",
    )
    C1: float = Field(8.8, description="Empirical dune erosion parameter", unit="m")
    C2: float = Field(4.6, description="Empirical dune erosion parameter", unit="m")
    DuneRestart: NonNegativeFloat = Field(
        0.075,
        description="Restart height for dunes lowered to essentially zero",
        unit="m",
    )
    Rat: float = Field(
        0.0,
        description=(
            "Rate of shoreline retreat attributed to alongshore transport; "
            "(-) = erosion, (+) = accretion"
        ),
        unit="m / y",
    )
    RSLR_Constant: bool = Field(
        True,
        description=(
            "Relative sea-level rise rate will be constant, otherwise logistic growth "
            "function used for time series"
        ),
    )
    RSLR_const: float = Field(
        0.004,
        description="Relative sea-level rise rate",
        unit="m / y",
    )
    beta: NonNegativeFloat = Field(
        0.04,
        description="Beach slope for runup calculations",
    )
    StormSeries: conlist(item_type=NonNegativeFloat) = Field(
        [],
        description="Time series of storms",
    )
    nn: float = Field(
        0.5,
        description="Flow routing constant",
    )
    mm: float = Field(
        2.0,
        description="Exponent constant for sediment transport",
    )
    Rin_r: float = Field(
        2.0,
        description=(
            "Run-up regime infiltration rate (volume of overwash flow lost per m "
            "cross-shore per time step)"
        ),
    )
    Rin_i: float = Field(
        0.25,
        description=(
            "Inundation regime infiltration rate (volume of overwash flow lost per m "
            "cross-shore per time step)"
        ),
    )
    Qs_min: float = Field(
        1.0,
        description="Minimum discharge needed for sediment transport",
        unit="m^3 / hr",
    )
    MaxUpSlope: float = Field(
        0.25,
        description="Maximum slope water can flow upward",
        unit="m / m",
    )
    threshold_in: float = Field(
        0.25,
        description="Threshold to determine if in inundation regime",
        unit="m^3 / hr",
    )
    Kr: float = Field(
        0.000075,
        description="Sediment flux constant, run-up regime",
    )

    Ki: float = Field(
        0.0000075,
        description="Sediment flux constant, inundation regime",
    )
    Cbb_r: confloat(gt=0.0, lt=1.0) = Field(
        0.5,
        description=(
            "Coefficient for exponential decay of sediment load entering back-barrier "
            "bay in run-up regime"
        ),
    )
    Cbb_i: confloat(gt=0.0, lt=1.0) = Field(
        0.8,
        description=(
            "Coefficient for exponential decay of sediment load entering back-barrier "
            "bay in inundation regime"
        ),
    )
    Qs_bb_min: NonNegativeFloat = Field(
        1.0,
        description=(
            "Minimum sediment flux in back-barrier bay (below which sediment won't "
            "flux)"
        ),
        unit="m^3 / hr",
    )
    Cx: NonNegativeFloat = Field(
        10.0,
        description=(
            "Multiplier with the average slope of the interior for constant C in "
            "inundation transport rule"
        ),
    )
    OWss_i: conint(ge=1) = Field(2, description="Overwash substep")
    OWss_r: conint(ge=1) = Field(1, description="Overwash substep")
    k_sf: float = Field(
        5000,
        description="Shoreface flux rate constant",
        unit="m^3 / m / y",
    )
    s_sf_eq: float = Field(
        0.02,
        description="Equilibrium shoreface slope",
    )
    Shrub_ON: bool = Field(
        False, description="1 = shrubs on in simulation, 0 = shrubs off"
    )
    Seedmin: NonNegativeFloat = Field(
        1000.0,
        description="Seeds produced per shrub per year (fecundity)",
        unit="1 / yr",
    )
    Seedmax: NonNegativeFloat = Field(
        10000.0,
        description="Seeds produced per shrub per year (fecundity)",
        unit="1 / yr",
    )
    disp_mu: float = Field(
        -0.721891,
        description="For lognormal probability distribution of seed dispersal distance",
    )
    disp_sigma: float = Field(
        1.5,
        description="For lognormal probability distribution of seed dispersal distance",
    )
    Dshrub: NonNegativeFloat = Field(
        2.75,
        description="Minimum elevation of fronting dune for shrub growth",
        unit="m",
    )
    GermRate: confloat(ge=0.0, le=1.0) = Field(0.6, description="Germination rate")
    TimeFruit: NonNegativeFloat = Field(
        5.0, description="Age shrubs need to be before they start fruiting", unit="y"
    )
    Female: confloat(ge=0.0, le=1.0) = Field(
        0.5, description="Percentage of shrubs that are female"
    )
    ShrubEl_min: NonNegativeFloat = Field(
        1.2,
        description="Elevation range for shrub growth, minimum bound",
        unit="m NAVD88",
    )
    ShrubEl_max: NonNegativeFloat = Field(
        2.3,
        description="Elevation range for shrub growth, maximum bound",
        unit="m NAVD88",
    )
    SprayDist: NonNegativeFloat = Field(
        170.0,
        description="Distance from ocean shoreline that shrubs can establish",
        unit="m",
    )
    BurialLimit: NonNegativeFloat = Field(
        0.75,
        description=(
            "Maximum percentage of height that a shrub can be buried up to before dying"
        ),
        unit="m",
    )
    UprootLimit: NonPositiveFloat = Field(
        -0.2, description="Shrubs eroded beyond this limit killed", unit="m"
    )
    SalineLimit: NonNegativeFloat = Field(
        0.05,
        description="Dishcharge limit to determine shrub mortality via saline flooding",
        unit="m^3 / hr",
    )
    Qshrub_max: confloat(ge=0.0, le=1.0) = Field(
        0.15,
        description=(
            "Maximum percentage of overwash reduction through a shrub cell with full "
            "percent cover"
        ),
    )
    DuneParamStart: bool = Field(
        True, description="Dune crest heights will come from external file; copied for each sequential dune row"
    )
    DuneParamMultipleRows: bool = Field(
        False, description="if DuneParamStart = True, use dune crest heights to populate multiple dune rows with different crest heights"
    )
    GrowthParamStart: bool = Field(
        True, description="Dune growth parameters will come from external file"
    )
    MaxShrubHeight: NonNegativeFloat = Field(
        5.2, description="Maximum shrub height", unit="m"
    )
    ShorefaceToe: float = Field(0.0, description="Start location of shoreface toe [m]")
    SeededRNG: bool = Field(
        True, description="Use seeded random number generator for reproducibility"
    )

    @validator("rmax", always=True)
    def check_rmax(cls, v, values, **kwargs):  # noqa: B902
        min_field = "rmin"

        with suppress(KeyError):
            if v < values[min_field]:
                raise ValueError(
                    f"ensure this value is greater than or equal to {min_field} "
                    f"({values[min_field]})"
                )
        return v

    @validator("Seedmax", always=True)
    def check_seed_max(cls, v, values, **kwargs):  # noqa: B902
        min_field = "Seedmin"

        with suppress(KeyError):
            if v < values[min_field]:
                raise ValueError(
                    f"ensure this value is greater than or equal to {min_field} "
                    f"({values[min_field]})"
                )
        return v

    @validator("ShrubEl_max", always=True)
    def check_shrub_elevation_max(cls, v, values, **kwargs):  # noqa: B902
        min_field = "ShrubEl_min"

        with suppress(KeyError):
            if v < values[min_field]:
                raise ValueError(
                    f"ensure this value is greater than or equal to {min_field} "
                    f"({values[min_field]})"
                )
        return v

    def to_yaml(self):
        lines = []
        for k, meta in sorted(self._meta().items()):
            unit = meta["unit"]
            value = meta["value"]
            desc = meta["description"]

            if desc:
                lines += textwrap.wrap(
                    desc, width=88, initial_indent="# ", subsequent_indent="# "
                )
            lines += [f"{k}: {value!r}" + (f"  # [{unit}]" if unit else "")]

        return os.linesep.join(lines)

    def to_py(self):
        lines = []
        for k, meta in sorted(self._meta().items()):
            unit = meta["unit"]
            value = meta["value"]
            desc = meta["description"]

            if desc:
                lines += textwrap.wrap(
                    desc, width=88, initial_indent="# ", subsequent_indent="# "
                )
            lines += [f"{k} = {value!r}" + (f"  # [{unit}]" if unit else "")]

        return os.linesep.join(lines)

    def _meta(self):
        schema = self.schema()

        meta = {}
        for k, v in sorted(self.dict().items()):
            if isinstance(v, pathlib.Path):
                v = str(v)
            meta[k] = {
                "value": v,
                "unit": schema["properties"][k].get("unit", None),
                "description": schema["properties"][k].get("description", None),
            }

        return meta

    @classmethod
    def from_yaml(cls, filepath):
        import yaml

        with open(filepath) as fp:
            config = yaml.safe_load(fp)

        return cls(**config)

    @classmethod
    def from_py(cls, filepath):
        from importlib.machinery import SourceFileLoader

        valid_keys = set(cls.schema()["properties"])

        mod = SourceFileLoader("config", str(filepath)).load_module()
        config = {k: v for k, v in mod.__dict__.items() if k in valid_keys}

        return cls(**config)
