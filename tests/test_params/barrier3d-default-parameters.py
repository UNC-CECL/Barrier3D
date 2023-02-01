elevation_file = "barrier3d-default-elevations.npy"  # File that contains initial elevations
dune_file = "barrier3d-default-dunes.npy"  # File that contains initial dune values
growth_param_file = (
    "barrier3d-default-growthparam.npy"  # File that contains initial growth parameters
)
storm_file = "barrier3d-default-storms.npy"  # File that contains storm data
TMAX = 51  # [y] Duration of simulation
StormStart = 2  # [y] Year when storm can start occurring (NOTE: if changed, need new storm time series)
BarrierLength = 500.0  # [m] Static length (alongshore) of island segment
DuneWidth = 20.0  # [m] Width (cross-shore) of island dune field
LShoreface = 500.0  # [m] Initial length of shoreface
DShoreface = 10.0  # [m] Height of shoreface
BayDepth = 3.0  # [m] Depth of bay behind island segment
MHW = 0.46  # [m] Elevation of Mean High Water (NOTE: if changed, need new storm time series)
DuneParamStart = True  # Dune height will come from external file
Dstart = 0.5  # [m] Initial height of dune domain above berm elevation
BermEl = 1.9  # [m] Static elevation of berm; berm elevation + dune height = dune elevation (NOTE: if changed, need new MSSM and storms)
GrowthParamStart = True  # Dune growth parameter will come from external file
rmin = 0.35  # Minimum growth rate for logistic dune growth
rmax = 0.85  # Maximum growth rate for logistic dune growth
HdDiffu = 0.75  # [m] Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
Dmaxel = 3.4  # [m] Maximum elevation of dunes
C1 = 8.8  # [m] Empirical dune erosion parameter
C2 = 4.6  # [m] Empirical dune erosion parameter
DuneRestart = 0.075  # [m] Restart height for dunes lowered to essentially zero
Rat = 0.0  # [m / y] Rate of shoreline reatreat attributed to alongshore transport; (-) = erosion, (+) = accretion
RSLR_Constant = True  # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
RSLR_const = 0.004  # [m / y] Relative sea-level rise rate
beta = 0.04  # Beach slope for runup calculations
StormSeries = []  # Time series of storms
nn = 0.5  # Flow routing constant
mm = 2.0  # Exponent constant for sediment transport
Rin_r = 2.0  # Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Rin_i = 0.1  # Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Qs_min = 1.0  # [m^3 / hr] Minimum discharge needed for sediment transport
MaxUpSlope = 0.25  # Maximum slope water can flow upward
threshold_in = 0.25  # [m^3 / hr] Threshold to determine if in inundation regime
Kr = 0.000075  # Sediment flux constant, run-up regime
Ki = 0.0000075  # Sediment flux constant, inundation regime
Cbb_r = 0.7  # Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime
Cbb_i = 0.85  # Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime
Qs_bb_min = 1  # [m^3 / hr] Minimum sediment flux in back-barrier bay (below which sediment won't flux)
Cx = 10.0  # Multiplier with the average slope of the interior for constant "C" in inundation transport rule
OWss_i = 2  # Overwash substep
OWss_r = 1  # Overwash substep
k_sf = 5000.0  # [m^3 / m / y] Shoreface flux rate constant
s_sf_eq = 0.02  # Equilibrium shoreface slope
Shrub_ON = 0  # 1 = shrubs on in simulation, 0 = shrubs off
Seedmin = 100.0  # [1 / yr] Seeds produced per shrub per year (fecundity)
Seedmax = 1000.0  # [1 / yr] Seeds produced per shrub per year (fecundity)
disp_mu = -0.721891  # For lognormal probability distribution of seed dispersal distance
disp_sigma = 1.5  # For lognormal probability distribution of seed dispersal distance
Dshrub = 2.75  # [m] Minimum elevation of fronting dune for shrub growth
GermRate = 0.6  # Germination rate
TimeFruit = 5.0  # [yr] Age shrubs need to be before they start fruiting
Female = 0.5  # Percentage of shrubs that are female
ShrubEl_min = 1.2  # [m] Elevation range for shrub growth, minimum bound
ShrubEl_max = 2.3  # [m] Elevation range for shrub growth, maximum bound
SprayDist = 170  # [m] Distance from ocean shoreline that shrubs can establish
BurialLimit = 0.75  # [m] Maximum percentage of height that a shrub can be buried up to before dying
UprootLimit = -0.2  # [m] Shrubs eroded beyond this limit killed
SalineLimit = (
    5  # [m^3 / hr] Dishcharge limit to determine shrub mortality via saline flooding
)
Qshrub_max = 0.15  # Maximum percentage of overwash reduction through a shrub cell with full percent cover
MaxShrubHeight = 5.3  # [m] Maximum shrub height
ShorefaceToe = 0  # [m] Start location of shoreface toe
SeededRNG = True  # Use seeded random number generator for reproducibility
