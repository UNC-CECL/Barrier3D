# Time series creation for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 3
# Updated: 21 Sep 2020

# Creates and saves storm time series, initial dune height, and dune growth rates

# This version randomly selects from actual storms binned by surge/tide level, and can vary the frequency that a storm is chosen from each bin


import numpy as np
import random
import bisect
import matplotlib.pyplot as plt
from distfit import distfit



#####################################################################################################################
### Generate storm time series

from Barrier3D_Parameters import (mean_storm, SD_storm, MHW, StormStart, BermEl)

# Load Storms
StormList = np.load('Parameterization/StormList.npy')
simTWL = StormList[:, 2]

# # Fit Distribution - Note: Beta distribution is best, with parameters listed below
# dist = distfit()
# dist.fit_transform(simTWL)
# dist.plot()

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

# Shift in location of distribution (i.e. SHIFT STORM INTENSITY) 
shift = -0.15

# Make storm series
StormSeries = np.zeros([StormStart, 5])

for t in range(StormStart, 1000):
    
    # Calculate number of storms in year
    numstorm = max(0,round(np.random.normal(mean_storm, SD_storm)))
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
        stormTS[n,4] = round(dur / 2) # Divided by two assuming TWL only for only half of storm
          
    # Save
    StormSeries = np.vstack([StormSeries, stormTS])

np.save('C:/Barrier3D/Parameterization/StormTimeSeries_1000yr.npy', StormSeries)    


### Plots
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

