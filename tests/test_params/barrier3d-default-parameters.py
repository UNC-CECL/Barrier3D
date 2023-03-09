# Static length (alongshore) of island segment
BarrierLength = 500.0  # [m]
# Depth of bay behind island segment
BayDepth = 3.0  # [m]
# Static elevation of berm; berm elevation + dune height = dune elevation (NOTE: if
# changed, need new MSSM and storms)
BermEl = 1.9  # [m NAVD88]
# Maximum percentage of height that a shrub can be buried up to before dying
BurialLimit = 0.75  # [m]
# Empirical dune erosion parameter
C1 = 8.8  # [m]
# Empirical dune erosion parameter
C2 = 4.6  # [m]
# Coefficient for exponential decay of sediment load entering back-barrier bay in
# inundation regime
Cbb_i = 0.85
# Coefficient for exponential decay of sediment load entering back-barrier bay in run-up
# regime
Cbb_r = 0.7
# Multiplier with the average slope of the interior for constant C in inundation
# transport rule
Cx = 10.0
# Height of shoreface
DShoreface = 10.0  # [m]
# Maximum elevation of dunes
Dmaxel = 3.4  # [m NAVD88]
# Minimum elevation of fronting dune for shrub growth
Dshrub = 2.75  # [m]
# Initial height of dune domain above berm elevation
Dstart = 0.5  # [m]
# Dune crest heights will come from external file; copied for each sequential dune row
DuneParamStart = True
# if DuneParamStart = True, use dune crest heights to populate multiple dune rows with different crest heights
DuneParamMultipleRows: False
# Restart height for dunes lowered to essentially zero
DuneRestart = 0.075  # [m]
# Width (cross-shore) of island dune field; for illustration purposes only
DuneWidth = 20.0  # [m]
# Percentage of shrubs that are female
Female = 0.5
# Germination rate
GermRate = 0.6
# Dune growth parameters will come from external file
GrowthParamStart = True
# Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
HdDiffu = 0.75  # [m]
# Sediment flux constant, inundation regime
Ki = 7.5e-06
# Sediment flux constant, run-up regime
Kr = 7.5e-05
# Length of shoreface
LShoreface = 500.0  # [m]
# Elevation of Mean High Water (NOTE: if changed, need new storm time series)
MHW = 0.46  # [m NAVD88]
# Maximum shrub height
MaxShrubHeight = 5.3  # [m]
# Maximum slope water can flow upward
MaxUpSlope = 0.25  # [m / m]
# Overwash substep
OWss_i = 2
# Overwash substep
OWss_r = 1
# Minimum sediment flux in back-barrier bay (below which sediment won't flux)
Qs_bb_min = 1.0  # [m^3 / hr]
# Minimum discharge needed for sediment transport
Qs_min = 1.0  # [m^3 / hr]
# Maximum percentage of overwash reduction through a shrub cell with full percent cover
Qshrub_max = 0.15
# Relative sea-level rise rate will be constant, otherwise logistic growth function used
# for time series
RSLR_Constant = True
# Relative sea-level rise rate
RSLR_const = 0.004  # [m / y]
# Rate of shoreline retreat attributed to alongshore transport; (-) = erosion, (+) =
# accretion
Rat = 0.0  # [m / y]
# Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore
# per time step)
Rin_i = 0.1
# Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per
# time step)
Rin_r = 2.0
# Dishcharge limit to determine shrub mortality via saline flooding
SalineLimit = 5.0  # [m^3 / hr]
# Use seeded random number generator for reproducibility
SeededRNG = True
# Seeds produced per shrub per year (fecundity)
Seedmax = 1000.0  # [1 / yr]
# Seeds produced per shrub per year (fecundity)
Seedmin = 100.0  # [1 / yr]
# Start location of shoreface toe [m]
ShorefaceToe = 0.0
# Elevation range for shrub growth, maximum bound
ShrubEl_max = 2.3  # [m NAVD88]
# Elevation range for shrub growth, minimum bound
ShrubEl_min = 1.2  # [m NAVD88]
# 1 = shrubs on in simulation, 0 = shrubs off
Shrub_ON = False
# Distance from ocean shoreline that shrubs can establish
SprayDist = 170.0  # [m]
# Time series of storms
StormSeries = []
# Year when storm can start occurring (NOTE: if changed, need new storm time series)
StormStart = 2  # [y]
# Duration of simulation
TMAX = 51  # [y]
# Age shrubs need to be before they start fruiting
TimeFruit = 5.0  # [y]
# Shrubs eroded beyond this limit killed
UprootLimit = -0.2  # [m]
# Beach slope for runup calculations
beta = 0.04
# For lognormal probability distribution of seed dispersal distance
disp_mu = -0.721891
# For lognormal probability distribution of seed dispersal distance
disp_sigma = 1.5
# File that contains initial dune height values [m]
dune_file = "barrier3d-default-dunes.npy"
# File that contains initial elevations in [m MHH]
elevation_file = "barrier3d-default-elevations.npy"
# File that contains initial growth parameters
growth_param_file = "barrier3d-default-growthparam.npy"
# Shoreface flux rate constant
k_sf = 5000.0  # [m^3 / m / y]
# Exponent constant for sediment transport
mm = 2.0
# Flow routing constant
nn = 0.5
# Maximum growth rate for logistic dune growth
rmax = 0.85
# Minimum growth rate for logistic dune growth
rmin = 0.35
# Equilibrium shoreface slope
s_sf_eq = 0.02
# File that contains storm data
storm_file = "barrier3d-default-storms.npy"
# Threshold to determine if in inundation regime
threshold_in = 0.25  # [m^3 / hr]
