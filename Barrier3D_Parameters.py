# Parameter value loading script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 4
# Updated: 30 April 2021


# Script sets up all input parameters
# Converts from meters to decameters for simulation



import numpy as np
import math



elevfile = 'Parameterization/InitElevHog.npy'
stormfile = 'Parameterization/StormTimeSeries_1000yr.npy'
dunestartfile = 'Parameterization/DuneStart_1000dam.npy'
growthparamfile = 'Parameterization/growthparam_1000dam.npy'

################################
### TIME

TMAX = 50 + 1                 
StormStart = 2




################################
### COMPUTATIONAL DOMAIN 

# Vertical Dimensions
LShoreface = 500 /10
DShoreface = 10 /10
BayDepth = 3 /10
MHW = 0.46 /10 # Used as offset to convert given elevations relative to a MHW of 0

# Elevation (decameters) 
InteriorDomain = np.load(elevfile)

# Horizontal Dimensions
BarrierLength = int(500 /10)
if len(InteriorDomain[0]) > BarrierLength:
    InteriorDomain = InteriorDomain[:,0:BarrierLength] # Reduce to specified max length
else:
    BarrierLength = len(InteriorDomain[0])
DomainWidth = len(InteriorDomain)
DuneWidth = int(20 /10)




################################
### Storm Time Series
StormTimeSeries = True
StormSeries = np.load(stormfile)




################################
### DUNES

# Dune height refers to heigh of dune above the static berm elevation
Dstart = 0.5 /10
BermEl = 1.9 /10 - MHW

# Initialize dune crest height domain
if StormTimeSeries:
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0,:,0] = np.ones([1, BarrierLength]) * (Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1,BarrierLength)))
    for w in range(1,DuneWidth):
        DuneDomain[0,:,w] = DuneDomain[0,:,0]
else:
    DuneStart = np.load(dunestartfile)
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0,:,0] = DuneStart[0:BarrierLength]    
    for w in range(1,DuneWidth):
        DuneDomain[0,:,w] = DuneDomain[0,:,0]  
    
# Dune growth parameter
rmin = 0.35
rmax = 0.85
if StormTimeSeries:
    growthparam = rmin + (rmax-rmin) * np.random.rand(1,BarrierLength)
else:
    growthparamstart = np.load(growthparamfile)
    growthparam = growthparamstart[0:BarrierLength]

# Dune diffusion parameter
HdDiffu = 0.75 /10

# Maximum dune height
Dmaxel = 3.4 /10 - MHW

# Erosion parameters
C1 = 8.8 
C2 = 4.6
DuneRestart = 0.075 /10




################################
### ALONGSHORE TRANSPORT & RSLR

# Volume of sediment lost via alongshore transport
Rat = 0 /10 # dam   Note: Positive value will result in erosion, negative will result in progradation
Qat = Rat * DShoreface # dam^3/dam

# Relative Sea Level Rise Rate
RSLR_Constant = True
if RSLR_Constant:
    # Constant RSLR
    RSLR_const = 0.004 /10 
    RSLR = [RSLR_const] * TMAX 
else:
    # Logistic RSLR rate projection - Rohling et al. (2013)
    RSLR = []
    alpha = 0.75 # m/yr -- probability maximum = 0.75, 68% upper bound = 2.0
    beta = alpha / 0.003 - 1 # constant
    gamma = 350 # yr -- probability maximum = 350, 68% upper bound = 900
    C = 12 # constant
    for t in range(150,TMAX+150):
        delta = alpha / (1 + beta * math.exp(-t / gamma * C)) / 10000 * 10 # Convert from m/cy to dam/yr
        RSLR.append(delta)



################################
### STORM

mean_storm = 8.3 
SD_storm = 5.9 
numstorm = 0
beta = 0.04




################################
### OVERWASH

# Substeps
OWss_r = 1
OWss_i = 2

# Flow Routing
nn = 0.5
mm = 2
Rin_r = 2
Rin_i = 0.1
MaxUpSlope = 0.25
threshold_in = 0.25

# Sediment Transport
Kr = 0.000075
Ki = 0.0000075
Qs_min = 1 /1000 # Convert to dam^3
Qs_bb_min = 1 /1000 # Convert to dam^3
Cx = 10
Cbb_r = 0.7
Cbb_i = 0.85




################################
### SHOREFACE DYNAMICS

k_sf = 5000
s_sf_eq = 0.02




################################
### SHRUBS

# Dispersal
Shrub_ON = 0
Seedmin = 100
Seedmax = 1000
disp_mu = -0.721891
disp_sigma = 1.5

# Growth
Dshrub = 2.75 /10
GermRate = 0.6
TimeFruit = 5
Female = 0.5
ShrubEl_min = 1.2 /10 - MHW
ShrubEl_max = 2.3 /10 - MHW
MaxShrubHeight = 5.3 /10
SprayDist = 170 /10

# Overwash Interaction
BurialLimit = 0.75 # percent
UprootLimit = -0.2 /10 # Convert to dam
SalineLimit = 5 /1000 # Convert to dam^3
Qshrub_max = 0.15 # percent

# Percent cover change (years 0-9)
PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
addend = np.ones(TMAX+50) # (years 10+)
PC = np.append(PC, addend)


SL = 0

    
SimParams = [TMAX,
             RSLR, 
             MHW, 
             BermEl, 
             BarrierLength, 
             s_sf_eq, 
             DShoreface, 
             LShoreface, 
             Shrub_ON,
             Dmaxel,
             beta,
             BayDepth,
             growthparam, 
             DuneWidth, 
             HdDiffu, 
             k_sf, 
             Qat, 
             Dshrub,
             Female, 
             ShrubEl_min, 
             ShrubEl_max,
             BurialLimit, 
             UprootLimit, 
             TimeFruit, 
             Seedmin, 
             Seedmax, 
             GermRate, 
             disp_mu, 
             disp_sigma, 
             SalineLimit, 
             PC,
             Dstart,
             SD_storm, 
             StormStart,
             mean_storm, 
             numstorm, 
             Qshrub_max, 
             C1, 
             C2, 
             DuneRestart, 
             nn, 
             mm,
             threshold_in, 
             Rin_r, 
             Rin_i, 
             Qs_min, 
             Kr, 
             Ki, 
             Cbb_r, 
             Cbb_i, 
             Qs_bb_min, 
             Cx, 
             OWss_i, 
             OWss_r,
             MaxShrubHeight,
             SprayDist,
             SL,
             MaxUpSlope]
          

