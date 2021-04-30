# Script for running batch simulations with

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


# I.R.B. Reeves - 19 Sep 2019

# Runs range of dune growth rates and average storms/year -- with permanent shift in average TWL


# Version Number: 1
# Updated: 30 April 2021


# Imports
import os
import time
import Barrier3D_Batch as batch
from datetime import datetime
import numpy as np
import random
import bisect
#import multiprocessing
from joblib import Parallel, delayed

from Barrier3D_Parameters_Batch import (SD_storm, MHW, StormStart)

Time = time.time()


#==============================================================================================================================================
# Set-up

# Batch Specifications
Num = 25 # Duplicates

# Make new directory for batch
name = datetime.today().strftime('%Y_%m%d_%H%M')
directory = 'Output/BatchSims_' + name
print('BatchSims_' + name)
os.makedirs(directory)

# Batch loops
Rmin = [0.05, 0.20, 0.35, 0.50, 0.65]
Rmax = [0.55, 0.70, 0.85, 1.00, 1.15]
Nstorm = [4, 6, 8, 10, 12]

# Beta Distribution Parameters
loc = -0.23670528422151482
scale = 733.1726404400683
a = 6.168920215856851
b = 2845.158579719937

# Bin Specifications
BinStart = -0.2
BinStop = 5.0
BinWidth = 0.1
BinNum = int(((BinStop - BinStart) / BinWidth) + 2)
Bin = np.linspace(BinStart, BinStop, BinNum)

# Load Storms
StormList = np.load('Parameterization/StormList.npy')
simTWL = StormList[:, 2]

# Shift in Average TWL
shift = -0.15

# Makes sets of storm series
def makeStormSet():    
    StormSet = [None] * len(Nstorm)
    for Ns in range(len(Nstorm)):
        S = Nstorm[Ns]
        # Make Time series
        StormSeries = np.zeros([StormStart, 5])
        for t in range(StormStart, 1000):   
            # Calculate number of storms in year
            numstorm = round(np.random.normal(S, SD_storm))
            if numstorm < 0:
                numstorm = 0
            stormTS = np.zeros([numstorm,5])    
            # Select storms for year
            for n in range(numstorm):       
                ST = np.random.beta(a,b) * scale + loc + shift
                if ST < 0: 
                    ST = 0
                elif ST > simTWL[-1]: 
                    ST = simTWL[-1]      
                STmin = Bin[bisect.bisect(Bin, ST) - 1]
                STmax = STmin + 0.1                    
                indexMin = bisect.bisect(simTWL, STmin)
                indexMax = bisect.bisect(simTWL, STmax) - 1               
                if indexMin == 0:
                    storm = 1
                elif indexMin >= indexMax:
                    storm = indexMin
                else:
                    storm = random.randint(indexMin, indexMax) 
                
                dur = StormList[storm, 1] # Duration          
                Rhigh = StormList[storm, 2] # TWL
                period = StormList[storm, 4] # Tp
                Rlow = StormList[storm, 6] # Rlow
                
                stormTS[n,0] = t
                stormTS[n,1] = Rhigh / 10 - MHW
                stormTS[n,2] = Rlow / 10 - MHW
                stormTS[n,3] = period
                stormTS[n,4] = round(dur / 2)            
            # Save
            StormSeries = np.vstack([StormSeries, stormTS])
            StormSet[Ns] = StormSeries
    return StormSet

# Makes sets of dune growth rates
def makeRSet():
    RSet = [None] * len(Rmin)
    for g in range(len(Rmin)):             
        rmin = Rmin[g] 
        rmax = Rmax[g]
        growthparam = rmin + (rmax-rmin) * np.random.rand(1000)
        RSet[g] = growthparam 
    return RSet

# Batch script
def RunBatch(n):
    
    # Initiate BatchData array
    sims = len(Rmin) * len(Nstorm)
    Data = np.zeros([sims,14])
    
    # Initiate SimNumber - unique ID for each simulation
    SimNumber = 1 + (sims * n)
          
    # Make sets of storm series
    StormSet = makeStormSet()
    
    # Make sets of dune growth rates
    RSet = makeRSet()   
  
    ### Run simulations looping through parameter space
    for Ns in range(len(Nstorm)): 
        
        StormSeries = StormSet[Ns]
        S = Nstorm[Ns]
        
        for g in range(len(Rmin)):             

            growthparam = RSet[g]
            
            rmin = Rmin[g] 
            rmax = Rmax[g]
            growthrate = (rmin + rmax)/2
                     
            # Run simulation
            print('=================================')
            print('Storm Number = ', S, ', GrowthRates = ', growthrate, ', Batch Run ', SimNumber)
            Result, t_Result, InteriorWidth_Avg, ShorelineChange, Hd_avg, IslandArea, AvgInteriorElevation, SCperiod, AvgFastDur, AvgSlowDur, Punc = batch.Barrier3D_Batch(SimNumber, directory, StormSeries, growthparam)
            print(Result, t_Result, InteriorWidth_Avg, ShorelineChange, Hd_avg, IslandArea, AvgInteriorElevation, SCperiod, AvgFastDur, AvgSlowDur)
                       
            # Save data to array
            row = SimNumber -1 - (sims * n)
            Data[row,0] = SimNumber
            Data[row,1] = growthrate
            Data[row,2] = S
            Data[row,3] = Punc
            Data[row,4] = Result
            Data[row,5] = t_Result
            Data[row,6] = InteriorWidth_Avg
            Data[row,7] = ShorelineChange
            Data[row,8] = Hd_avg
            Data[row,9] = AvgInteriorElevation
            Data[row,10] = IslandArea
            Data[row,11] = SCperiod
            Data[row,12] = AvgSlowDur
            Data[row,13] = AvgFastDur
            
            # Increment SimNumber
            SimNumber += 1
    
    savename1 = directory + '/BatchData_' + name + '_PS' + str(n) + '.npy'
    np.save(savename1, Data)
         
    return Data


#==============================================================================================================================================
# Run Simulations
    
#num_cores = multiprocessing.cpu_count()
num_cores = 32
print('Cores: ' + str(num_cores))
print('Running...')

results = Parallel(n_jobs=num_cores)(delayed(RunBatch)(n) for n in range(Num))
BatchData = np.vstack(results)


#==============================================================================================================================================
# End
savename2 = directory + '/BatchData_' + name + '.npy'
np.save(savename2, BatchData)
            
print()
print('Elapsed Time of Batch Run: ', (time.time() - Time) / 60, 'minutes')


