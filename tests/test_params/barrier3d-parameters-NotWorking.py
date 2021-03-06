# Parameter value loading script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 4
# Updated: 26 August 2020


# Script sets up all input parameters
# Converts from meters to decameters for simulation


TMAX = 150                       # [y] Duration of simulation
StormStart = 2                   # [y] Year when storm can start occurring
BarrierLength = 500.0            # [m] Static length (alongshore) of island segment
DuneWidth = 20.0                 # [m] Width (cross-shore) of island dune field; for illustration purposes only
LShoreface = 500.0               # [m] Length of shoreface
DShoreface = 10.0                # [m] Height of shoreface
BayDepth = 3.0                   # [m] Depth of bay benind island segment
MHW = 0.46                       # [m] Elevation of Mean High Water
DuneParamStart = True            # Dune height will come from external file
Dstart = 0.50                    # [m] Initial height of dune domain above berm elevation
BermEl = 1.9                     # [m] Static elevation of berm; berm elevation + dune height = dune elevation
GrowthParamStart = True          # Dune growth parameter will come from external file
rmin = 0.35                      # Minimum growth rate for logistic dune growth
rmax = 0.85                      # Maximum growth rate for logistic dune growth
HdDiffu = 0.75                   # [m] Dune diffusion parameter (i.e. max height offset between adjacent dune cells)
Dmaxel = 3.4                     # [m] Maximum elevation of dunes
C1 = 8.8                         # [m] Empirical dune erosion parameter
C2 = 4.6                         # [m] Empirical dune erosion parameter
DuneRestart = 0.075              # [m] Restart height for dunes lowered to essentially zero
Rat = 0.0                        # [m / y] Rate of shoreline reatreat attributed to alongshore transport; (-) = erosion, (+) = accretion
RSLR_Constant = True             # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
RSLR_const = 0.004               # [m / y] Relative sea-level rise rate
# mean_storm = 8.3                 # For a random number of storms per year sampled from normal distribution
# SD_storm = 5.9                   # For a random number of storms per year sampled from normal distribution
# numstorm = 0                     # For a single constant number of storms per year
beta = 0.04                      # Beach slope for runup calculations
# StormTimeSeries = True           # Storms will come from a time series
StormSeries = []                 # Time series of storms
nn = 0.5                         # Flow routing constant
mm = 2.0                         # Exponent constant for sediment transport
Rin_r = 2.0                      # Run-up regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Rin_i = 0.25                     # Inundation regime infiltration rate (volume of overwash flow lost per m cross-shore per time step)
Qs_min = 1.0                     # [m^3 / hr] Minimum discharge needed for sediment transport
MaxUpSlope = 0.25                # Maximum slope water can flow upward
threshold_in = 0.25              # [m^3 / hr] Threshold to determine if in inundation regime
Kr = 0.000075                    # Sediment flux constant, run-up regime
Ki = 0.0000075                   # Sediment flux constant, inundation regime
Cbb_r = 0.5                      # Coefficient for exponential decay of sediment load entering back-barrier bay in run-up regime
Cbb_i = 0.8                      # Coefficient for exponential decay of sediment load entering back-barrier bay in inundation regime
Qs_bb_min = 1                    # [m^3 / hr] Minimum sediment flux in back-barrier bay (below which sediment won't flux)
Cx = 10.0                        # Multiplier with the average slope of the interior for constant "C" in inundation transport rule
OWss_i = 2                       # Overwash substep
OWss_r = 1                       # Overwash substep
k_sf = 5000.0                    # [m^3 / m / y] Shoreface flux rate constant
s_sf_eq = 0.02                   # Equilibrium shoreface slope
Shrub_ON = False                 # 1 = shrubs on in simulation, 0 = shrubs off
Seedmin = 100.0                  # [1 / yr] Seeds produced per shrub per year (fecundity)
Seedmax = 1000.0                 # [1 / yr] Seeds produced per shrub per year (fecundity)
disp_mu = -0.721891              # For lognormal probability distribution of seed dispersal distance
disp_sigma = 1.5                 # For lognormal probability distribution of seed dispersal distance
Dshrub = 2.0                     # [m] Minimum elevation of fronting dune for shrub growth
GermRate = 0.6                   # Germination rate
TimeFruit = 5.0                  # [yr] Age shrubs need to be before they start fruiting
Female = 0.5                     # Percentage of shrubs that are female
ShrubEl_min = 0.6                # [m] Elevation range for shrub growth, minimum bound
ShrubEl_max = 2.3                # [m] Elevation range for shrub growth, maximum bound
TideAmp = 1.2                    # [m] Tidal amplitude
SprayDist = 170                  # [m] Distance from ocean shoreline that shrubs can establish
BurialLimit = 0.5                # [m] Shrubs buried beyond this limit killed
UprootLimit = -0.3               # [m] Shrubs eroded beyond this limit killed
SalineLimit = 0.05               # [m^3 / hr] Dishcharge limit to determine shrub mortality via saline flooding
Qshrub_max = 0.15                # Maximum percentage of overwash reduction through a shrub cell with full percent cover
MaxShrubHeight = 3.5             # [m] Maximum shrub height
