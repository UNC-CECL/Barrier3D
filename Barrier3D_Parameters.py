# Parameter value loading script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 3
# Updated: 28 May 2020


# Script sets up all input parameters
# Converts from meters to decameters for simulation



import numpy as np


################################
### TIME

TMAX = 100                  
StormStart = 2




################################
### COMPUTATIONAL DOMAIN 

# Vertical Dimensions
LShoreface = 1000 /10
DShoreface = 10 /10
BayDepth = 2 /10
MHW = 0.46 /10 # Used as offset to convert given elevations relative to a MHW of 0

# Elevation (decameters) 
InteriorDomain = np.load('Parameterization/InitElev.npy')

# Horizontal Dimensions
BarrierLength = int(300 /10)
if len(InteriorDomain[0]) > BarrierLength:
    InteriorDomain = InteriorDomain[:,0:BarrierLength] # Reduce to specified max length
else:
    BarrierLength = len(InteriorDomain[0])
DomainWidth = len(InteriorDomain)
DuneWidth = int(20 /10)




################################
### Storm Time Series
StormTimeSeries = True
StormSeries = np.load('Parameterization/StormTimeSeries_1000yr.npy') # TEMP HARDWIRED




################################
### DUNES

# Dune height refers to heigh of dune above the static berm elevation
Dstart = 0.25 /10
BermEl = 1.7 /10 - MHW

# Initialize dune crest height domain
if StormTimeSeries == 0:
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0,:,0] = np.ones([1, BarrierLength]) * (Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1,BarrierLength)))
    for w in range(1,DuneWidth):
        DuneDomain[0,:,w] = DuneDomain[0,:,0]
else:
    DuneStart = np.load('Parameterization/DuneStart_1000dam.npy') # TEMP HARDWIRED
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0,:,0] = DuneStart[0:BarrierLength]    
    for w in range(1,DuneWidth):
        DuneDomain[0,:,w] = DuneDomain[0,:,0]  
    
# Dune growth parameter
rmin = 0.05
rmax = 0.55
if StormTimeSeries == 0:
    growthparam = rmin + (rmax-rmin) * np.random.rand(1,BarrierLength)
else:
    growthparamstart = np.load('Parameterization/growthparam_1000dam.npy') # TEMP HARDWIRED
    growthparam = growthparamstart[0:BarrierLength]

# Dune diffusion parameter
HdDiffu = 0.45 /10

# Maximum dune height
Dmaxel = 2.30 /10 - MHW

# Erosion parameters
C1 = 8.8 
C2 = 4.6
DuneRestart = 0.05 /10




################################
### ALONGSHORE TRANSPORT & RSLR

# Volume of sediment lost via alongshore transport
Rat = 0 /10 # dam   Note: Positive value will result in erosion, negative will result in progradation
Qat = Rat * DShoreface # dam^3/dam

# Relative Sea Level Rise Rate
RSLR = 0.004 /10 




################################
### STORM

# Number Per Year
mean_storm = 13.17 
SD_storm = 5.16 
numstorm = 0

# Water Level Forcing
# In meters - converted to decameters after water level calculations
surge_tide_m = 0 - (MHW *10)
surge_tide_sd = 0
height_mu = 0
height_sigma = 0
period_m = 0
period_sd = 0
duration_mu = 0
duration_sigma = 0
beta = 0.04




################################
### OVERWASH

# Flow Routing
nn = 0.5
mm = 2
Rin_r = 2
Rin_i = 0.25
Qs_min = 1
threshold_in = 0.25

# Sediment Transport
Kr = 0.0003
Ki = 0.000075
Cbb_r = 0.5
Cbb_i = 0.8
Qs_bb_min = 0.0001
Cx = 2

# Substeps
OWss_i = 2
OWss_r = 2



################################
### SHOREFACE DYNAMICS

k_sf = 5000
s_sf_eq = 0.01




################################
### SHRUBS

# Dispersal
Shrub_ON = 0
Seedmin = 100
Seedmax = 1000
disp_mu = -0.721891
disp_sigma = 1.5

# Growth
Dshrub = 2 /10
GermRate = 0.6
TimeFruit = 5
Female = 0.5
ShrubEl_min = 0.6 /10 - MHW
ShrubEl_max = 2.3 /10 - MHW

# Overwash Interaction
BurialLimit = 0.5 /10
UprootLimit = -0.3 /10
SalineLimit = 0.05
Qshrub_max = 0.15

# Percent cover change (years 0-9)
PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
addend = np.ones(TMAX+50) # (years 50+)
PC = np.append(PC, addend)

