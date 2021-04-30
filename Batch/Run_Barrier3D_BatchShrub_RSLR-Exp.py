# Script for running batch simulations with

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


# I.R.B. Reeves - 19 Sep 2019

# Runs range of [RSLR] with shrubs to explore impact of [RSLR] on shrub expansion


# Version Number: 1
# Updated: 24 November 2020


# Imports
import os
import time
import Barrier3D_BatchShrub as batch
from datetime import datetime
import numpy as np
import random
import math
#import multiprocessing
from joblib import Parallel, delayed

from Barrier3D_Parameters_Batch import (mean_storm, SD_storm, MHW, StormStart, rmin, rmax, TMAX, Dmaxel, BermEl, BarrierLength, DuneWidth)

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
slr = [3,6,9,12,15]

# Make set of dune growth rates
def makeGrowthParam():
    growthparam = rmin + (rmax-rmin) * np.random.rand(1000)  
    return growthparam
    
# Makes sets of storm series
def makeStormSet():    
    Ns = mean_storm
    # Make Time series
    StormList = np.load('Parameterization/StormList.npy')
    StormSeries = np.zeros([StormStart, 5])
    for t in range(StormStart, 1000):   
        # Calculate number of storms in year
        numstorm = round(np.random.normal(Ns, SD_storm))
        if numstorm < 0:
            numstorm = 0
        stormTS = np.zeros([numstorm,5])    
        # Select storms for year
        for n in range(numstorm):       
            storm = random.randint(0,len(StormList)-1) # Pick random row index of each storm in list 
            
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
    return StormSeries

def makeDuneDomain():    
    Dstart = np.random.uniform((0.1 /10), (Dmaxel - BermEl - (0.1/10))) # Random start height   
    DuneDomain = np.zeros([TMAX, BarrierLength, DuneWidth])
    DuneDomain[0,:,0] = np.ones([1, BarrierLength]) * (Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1,BarrierLength)))
    for w in range(1,DuneWidth):
        DuneDomain[0,:,w] = DuneDomain[0,:,0]
    return DuneDomain, Dstart



# Batch script
def RunBatch(n):
    
    Shrub_ON = 1
    
    # Initiate BatchData array
    sims = len(slr)
    Data = np.zeros([sims,14])
    
    # Initiate SimNumber - unique ID for each simulation
    SimNumber = 1 + (sims * n)
          
    # Make new growhparam set
    growthparam = makeGrowthParam()
    
    # Make sets of storm series
    StormSeries = makeStormSet() 
  
    # Make dune domain
    DuneDomain, Dstart = makeDuneDomain()  
  
    ### Run simulations looping through parameter space
        
    for r in range(len(slr)):             
        
        rslr = slr[r]
        RSLR = [(rslr / 10000)] * TMAX # convert to dam/yr 
        
        # Run simulation
        print('=================================')
        print('Shrub = ', Shrub_ON, ', RSLR = ', rslr, ', Batch Run ', SimNumber)
        Result, t_Result, InteriorWidth_Avg, ShorelineChange, Hd_avg, ShrubArea, AvgInteriorElevation, SCperiod, AvgFastDur, AvgSlowDur, Punc = batch.Barrier3D_BatchShrub(SimNumber, directory, StormSeries, growthparam, Shrub_ON, RSLR, DuneDomain, Dstart)
        print(Result, t_Result, InteriorWidth_Avg, ShorelineChange, Hd_avg, ShrubArea, AvgInteriorElevation, SCperiod, AvgFastDur, AvgSlowDur)
                   
        # Save data to array
        row = SimNumber -1 - (sims * n)
        Data[row,0] = SimNumber
        Data[row,1] = rslr
        Data[row,2] = Shrub_ON
        Data[row,3] = Punc
        Data[row,4] = Result
        Data[row,5] = t_Result
        Data[row,6] = InteriorWidth_Avg
        Data[row,7] = ShorelineChange
        Data[row,8] = Hd_avg
        Data[row,9] = AvgInteriorElevation
        Data[row,10] = ShrubArea
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


