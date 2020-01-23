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


# Set-up input files
loc = ("Barrier3D_Input.xlsx") # Input File Path

# Open Excel Files
wb = xlrd.open_workbook(loc)
param = wb.sheet_by_index(0) # Parameters sheet
morph = wb.sheet_by_index(1) # Initial morphology sheet


################################
### TIME

TMAX = int(param.cell_value(4,3)) + 1
StormStart = int(param.cell_value(5,3))




################################
### COMPUTATIONAL DOMAIN 

# Vertical Dimensions
LShoreface = param.cell_value(12,3) /10
DShoreface = param.cell_value(13,3) /10
BayDepth = param.cell_value(14,3) /10
MHW = param.cell_value(15,3) /10 # Used as offset to convert given elevations relative to a MHW of 0

# Elevation (decameters)
X = [] # Import X coordinates
for k in range(1,morph.nrows):
    X.append(morph.row_values(k)[0])    
Y = [] # Import Y coordinates
for k in range(1,morph.nrows):
    Y.append(morph.row_values(k)[1])   
elev = [] # Import elevation values
for k in range(1,morph.nrows):
    elev.append(morph.row_values(k)[2])
    
width = int(max(Y) + 1)
length = int(max(X) + 1)
elev_array = np.zeros([length,width]) # Initialize array

for n in range(len(X)): # Add elevation values to array
    xx = int(X[n])
    yy = int(Y[n])
    zz = elev[n]
    elev_array[xx,yy] = zz

intElev = elev_array[:,:] /10 # Convert to decameters
intElev = intElev - MHW # Convert elevations relative to a MHW of 0
intElev[intElev <= 0] = -BayDepth # Set all subaerial cells to bay depth

u = 1
while u == 1: # Remove all rows that have zero subaerial cells
    if all(intElev[:,-1] <= 0):
        intElev = intElev[:,:-1]
    else:
        u = 0     

InteriorDomain = np.flipud(np.rot90(intElev)) # Flip to correct orientation

# Horizontal Dimensions
BarrierLength = int(param.cell_value(9,3) /10)
if len(InteriorDomain[0]) > BarrierLength:
    InteriorDomain = InteriorDomain[:,0:BarrierLength] # Reduce to specified max length
else:
    BarrierLength = len(InteriorDomain[0])
DomainWidth = len(InteriorDomain)
DuneWidth = int(param.cell_value(11,3) /10)




################################
### DUNES

# Dune height refers to heigh of dune above the static berm elevation
Dstart = param.cell_value(19,3) /10
BermEl = param.cell_value(20,3) /10 - MHW

# Initialize dune crest height domain
DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
DuneDomain[0,:,0] = np.ones([1, BarrierLength]) * (Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1,BarrierLength)))
for w in range(1,DuneWidth):
    DuneDomain[0,:,w] = DuneDomain[0,:,0]

# Dune growth parameter
rmin = param.cell_value(21,3)
rmax = param.cell_value(22,3)
growthparam = rmin + (rmax-rmin) * np.random.rand(1,BarrierLength)

# Dune diffusion parameter
HdDiffu = param.cell_value(23,3) /10

# Maximum dune height
Dmaxel = param.cell_value(24,3) /10

# Erosion parameters
C1 = param.cell_value(25,3) 
C2 = param.cell_value(26,3)
DuneRestart = param.cell_value(27,3) /10




################################
### ALONGSHORE TRANSPORT & RSLR

# Volume of sediment lost via alongshore transport
Rat = (param.cell_value(31,3)) * (-1) /10 # dam
Qat = Rat * DShoreface # dam^3/dam

# Relative Sea Level Rise Rate
RSLR = param.cell_value(32,3) /10 




################################
### STORM

# Number Per Year
mean_storm = param.cell_value(36,3) 
SD_storm = param.cell_value(37,3) 
numstorm = param.cell_value(38,3)

# Water Level Forcing
# In meters - converted to decameters after water level calculations
surge_tide_m = param.cell_value(39,3) - (MHW *10)
surge_tide_sd = param.cell_value(40,3)
#height_m = param.cell_value(43,3)
#height_sd = param.cell_value(44,3)
height_mu = param.cell_value(41,3)
height_sigma = param.cell_value(42,3)
period_m = param.cell_value(43,3)
period_sd = param.cell_value(44,3)
duration_mu = param.cell_value(45,3)
duration_sigma = param.cell_value(46,3)
beta = param.cell_value(47,3)




################################
### OVERWASH

# Flow Routing
nn = param.cell_value(51,3)
mm = param.cell_value(52,3)
Rin_r = param.cell_value(53,3)
Rin_i = param.cell_value(54,3)
Qs_min = param.cell_value(55,3)
threshold_in = param.cell_value(56,3)

# Sediment Transport
Kr = param.cell_value(57,3)
Ki = param.cell_value(58,3)
Cbb_r = param.cell_value(59,3)
Cbb_i = param.cell_value(60,3)




################################
### SHOREFACE DYNAMICS

k_sf = param.cell_value(64,3)
s_sf_eq = param.cell_value(65,3)




################################
### SHRUBS

# Dispersal
Shrub_ON = param.cell_value(69,3)
Seedmin = param.cell_value(70,3)
Seedmax = param.cell_value(71,3)
disp_mu = param.cell_value(72,3)
disp_sigma = param.cell_value(73,3)

# Growth
Dshrub = param.cell_value(74,3) /10
GermRate = param.cell_value(75,3)
TimeFruit = param.cell_value(76,3)
Female = param.cell_value(77,3)
ShrubEl_min = param.cell_value(78,3) /10 - MHW
ShrubEl_max = param.cell_value(79,3) /10 - MHW

# Overwash Interaction
BurialLimit = param.cell_value(80,3) /10
UprootLimit = param.cell_value(81,3) /10
SalineLimit = param.cell_value(82,3)
Qshrub_max = param.cell_value(83,3)

# Percent cover change (years 0-9)
PC = np.array([0, 0.04, 0.08, 0.10, 0.15, 0.15, 0.20, 0.35, 0.80, 1])
addend = np.ones(TMAX+50) # (years 10+)
PC = np.append(PC, addend)

