# Reads WIS data, calculates TWL time series, and counts storms
 
# Data needed: one .onlns WIS file and one .csv tide gauge file

# IRBR 9Jun20


import numpy as np
import csv
import math
import matplotlib.pyplot as plt



################################################################
### Load WIS data 1980-2014

pathpre = ('C:\Barrier3D\Parameterization\Storms\WIS-63183')
name = ('\ST63183_v03.onlns')


filename = pathpre + name

data = np.loadtxt(filename)
DT = data[:,0]
Hs = data[:,9]
Tp = data[:,11]
WavD = data[:,15]

# Date Time
year = []
month = []
day = []
hour = []

for n in range(len(DT)):
    DTstring = str(DT[n])
    year.append(int(DTstring[0:4]))
    month.append(int(DTstring[4:6]))
    day.append(int(DTstring[6:8]))
    hour.append(int(DTstring[8:10]))



################################################################
### Load Water Level from CSV
filename2 = pathpre + '/Tide-8631044.csv'
with open(filename2,newline='') as csvfile:
    sl = list(csv.reader(csvfile))
sl = [float(i) for i in sl[0]]
sl = sl[1:-1] # Trim off first and last readings from tide gauge to match waves
SL = np.asarray(sl)



################################################################
### Running means

N = 24 * 30 # For 30 day running mean

SL_rm = np.convolve(SL, np.ones((N,))/N, mode='same')

Hs_rm = np.convolve(Hs, np.ones((N,))/N, mode='same')

Tp_rm = np.convolve(Tp, np.ones((N,))/N, mode='same')



################################################################
### Storm surge (non-tidal residuals) & astronomical tide

# Get from Matlab t_Tide for:

    # NTR
    # AT



################################################################
### R2%, TWL, Rlow
beta = 0.04 # Beach slope, Hog Island    
L0 = (9.8 * Tp**2) / (2 * math.pi) # Wavelength       
R2 = []
Rlow = []
TWL = []

for n in range(len(L0)):
    # Setup
    setup = 0.35 * beta * math.sqrt(Hs[n] * L0[n]) 
    
    # Incident band swash
    Sin = 0.75 * beta * math.sqrt(Hs[n] * L0[n]) 
    
    # Infragravity band swash
    Sig = 0.06 * math.sqrt(Hs[n] * L0[n])
    
    # Swash
    swash = math.sqrt((Sin**2) + (Sig**2))
    
    # R2%
    r2 = 1.1 * (setup + (swash/2))
    R2.append(r2)
    
    # TWL & Rlow
    twl = SL[n] + r2
    TWL.append(twl)
    Rlow.append(twl - (swash/2))
    
TWL = np.asarray(TWL)
Rlow = np.asarray(Rlow)
R2 = np.asarray(R2)



################################################################
### Plot
plt.figure()
fig = plt.gcf()
fig.set_size_inches(14,18)
plt.rcParams.update({'font.size':14})

x = []                          # TEMP: this is not exact!
for t in range(len(year)):
    yr = 1980 + t/(365*24)
    x.append(yr)

# SL
plt.subplot(6,1,1)
plt.plot(x,SL,color='silver')
plt.plot(x,SL_rm,color='teal')
plt.ylabel('Sea level [mNAVD88]')

# Hs
plt.subplot(6,1,2)
plt.plot(x,Hs,color='silver')
plt.plot(x,Hs_rm,color='teal')
plt.ylabel('Hs [m]')

# Tp
plt.subplot(6,1,3)
plt.plot(x,Tp,color='silver')
plt.plot(x,Tp_rm,color='teal')
plt.ylabel('Tp [m]')

# WavD
plt.subplot(6,1,4)
plt.plot(x,WavD,color='dimgray')
plt.ylabel('Wave Direction [degree]')

# R2
plt.subplot(6,1,5)
plt.plot(x,R2,color='dimgray')
plt.ylabel('R2% [m]')

# TWL
plt.subplot(6,1,6)
plt.plot(x,TWL,color='dimgray')
plt.ylabel('TWL [mNAVD88]')

plt.show() 



################################################################
### Find Storms

# Find Hs threshold
BermEl = 1.7

Hs_over_yearly = []
for y in range(35):
    start = 365 * 24 * y
    stop = 365 * 24 * (y + 1)
    hh = Hs[start:stop]
    tt = TWL[start:stop]
    Hs_over = hh[tt > BermEl]
    Hs_over_yearly.append(np.mean(Hs_over))
   
Hs_min = min(Hs_over_yearly)
Hs_threshold = math.floor(Hs_min / 0.05) * 0.05 # Threshold Hs needed to qualify as storm event, rounded to nearst 0.05 m    
    
# Plot yearly means
plt.figure()
x = np.linspace(1980,2014,35)
plt.plot(x,Hs_over_yearly)
fig = plt.gcf()
fig.set_size_inches(14,5)
plt.xlabel('Year')
plt.ylabel('Hs [m]')
plt.hlines(Hs_threshold, x[0], x[-1],colors='black',linestyles='dashed')
plt.show()

# Find storms
Storms = np.zeros([0,12]) # Storm start, storm stop, datenum start, datenum stop, Hs, dur, TWL, SL, Tp, R2, Rlow, year

t = 0
while t < len(hour):
    if Hs[t] >= Hs_threshold:    # Future improvement: discard storms where simultaneous surge is negative (but need to separate surge from tide first)
        stormStart = t
        height = Hs[t]
        dur = 1
        t += 1
        while len([x for x in Hs[t:t+25] if x > Hs_threshold]) > 0: # If Hs drops below Hs_threshold for only 24 hrs or less, exceedence is assumed part of same weather system (Wahl et al., 2016; Li et al., 2014)
            if Hs[t] > Hs_threshold:
                dur += 1
                t += 1
            else:
                t += 1
        if dur > 8:
            stormStop = t
            datenumStart = DT[stormStart]
            datenumStop =  DT[stormStop]
            storm = [stormStart, stormStop, datenumStart, datenumStop, max(Hs[stormStart:stormStop+1]), dur, max(TWL[stormStart:stormStop+1]), max(SL[stormStart:stormStop+1]), max(Tp[stormStart:stormStop+1]), max(R2[stormStart:stormStop+1]), max(Rlow[stormStart:stormStop+1]), year[stormStart]]
            Storms = np.vstack([Storms, storm])
        t += 1
    else:
        t += 1

# Plot storm TWL histogram
plt.figure()
Bin = np.linspace(0.2,4.7,46)
stormtwl = Storms[:,6]
plt.hist(stormtwl, bins=Bin)
plt.xlabel('Storm TWL [mNAVD88]')
plt.title('1980 - 2014')

print('Max TWL:          ', max(Storms[:,6]))
print('Max Rexcess:      ', (max(Storms[:,6])-BermEl))

print('% Rhigh > BermEl: ', ((Storms[:,6]) > BermEl).sum() / len(Storms) * 100)
print('% Rlow  > BermEl: ', ((Storms[:,10]) > BermEl).sum() / len(Storms) * 100)


################################################################




