# Time series creation for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 3
# Updated: 21 Sep 2020

# Creates and saves storm time series, initial dune height, and dune growth rates

# This version randomly selects a storm from a given storm list

import numpy as np
import random
import matplotlib.pyplot as plt

#####################################################################################################################
### Generate storm time series

from Barrier3D_Parameters import (mean_storm, SD_storm, MHW, StormStart, Dstart, rmin, rmax, BermEl)
print()


# Set-up input files
StormList = np.load('C:/Barrier3D/Parameterization/StormList.npy')

# Time series
StormSeries = np.zeros([StormStart, 5])
for t in range(StormStart, 1000):   
    # Calculate number of storms in year
    numstorm = round(np.random.normal(mean_storm, SD_storm))
   
    if numstorm < 0:
        numstorm = 0
    stormTS = np.zeros([numstorm,5])    
    
    # Select storms for year
    for n in range(numstorm):       
        storm = random.randint(1,len(StormList)-1)
       
        dur = StormList[storm, 1] # Duration          
        Rhigh = StormList[storm, 2] # TWL
        period = StormList[storm, 4] # Tp
        Rlow = StormList[storm, 6] # Rlow
               
        stormTS[n,0] = t
        stormTS[n,1] = Rhigh / 10 - MHW
        stormTS[n,2] = Rlow / 10 - MHW
        stormTS[n,3] = period
        stormTS[n,4] = round(dur / 2) # Divided by two assuming TWL only for only half of storm
    
    # Save
    StormSeries = np.vstack([StormSeries, stormTS])
    
# np.save('C:/Barrier3D/Parameterization/StormTimeSeries_1000yr.npy', StormSeries)

## Plots
Bin = np.linspace(-1, 4.6, 57)

plt.figure()
surgetidecalc = StormSeries[:,2] * 10
plt.hist(surgetidecalc, bins=Bin)
plt.title('StormSeries Rlow')

plt.figure()
twlcalc = StormSeries[:,1] * 10
plt.hist(twlcalc, bins=Bin)
plt.title('StormSeries Rhigh')

plt.figure()
fig = plt.gcf()
fig.set_size_inches(16,4)
plt.plot(twlcalc)
plt.xlabel('Storm')
plt.hlines(BermEl*10, 0, len(twlcalc),colors='red',linestyles='dashed')
plt.ylabel('TWL (m MHW)')

print('Max TWL:          ', max(StormSeries[:,1])*10)
print('Max Rexcess:      ', (max(StormSeries[:,1])-BermEl)*10)

print('% Rhigh > BermEl: ', (twlcalc>(BermEl*10)).sum()/len(twlcalc)*100)
print('% Rlow  > BermEl: ', (surgetidecalc>(BermEl*10)).sum()/len(surgetidecalc)*100)
print()
print('[X] Storm series')


#####################################################################################################################
### Generate dune height start
DuneStart = np.ones([1000]) * (Dstart + (-0.01 + (0.01 - (-0.01)) * np.random.rand(1000)))
np.save('C:/Barrier3D/Parameterization/DuneStart_1000dam.npy', DuneStart)
print('[X] Dune start heights')


#####################################################################################################################
### Generate along-shore varying rmin & rmax
growthparam = rmin + (rmax-rmin) * np.random.rand(1000)
np.save('C:/Barrier3D/Parameterization/growthparam_1000dam.npy', growthparam)
print('[X] Dune growth rates')


print()
print('- Time series generation complete -')





