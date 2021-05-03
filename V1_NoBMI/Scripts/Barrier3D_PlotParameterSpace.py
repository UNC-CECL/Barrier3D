# Parameter space plotting script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 2
# Updated: 9 Apr 2020

# Creates a parameter space figure from saved numpy BatchData array

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt



#===============================
### Inputs

Group = 'BatchSims_Group30'
BatchDay = '2020_1106'
BatchTime = '2118'

Save = False # Whether or not to save figures

Rep = 25 # Number of duplicates

Rmin = [0.05, 0.20, 0.35, 0.50, 0.65]
Rmax = [0.55, 0.70, 0.85, 1.00, 1.15]
Nstorm = [4, 6, 8, 10, 12]


TMAX = 1000

# Calulations
Rspread = 0.25 #(Rmax[0] - Rmin[0]) / 2
simnum = len(Rmin) * len(Nstorm)



#===============================
### Load Data Into Arrays

loc = ('C:/Barrier3D/Output/' + Group + '/BatchSims_' + BatchDay + '_' + BatchTime + '/BatchData_' + BatchDay + '_' + BatchTime + '.npy') # Input File Path


BatchData = np.load(loc) 

Result = np.zeros([Rep, len(Nstorm), len(Rmin)])
Tresult = np.zeros([Rep, len(Nstorm), len(Rmin)])
Width = np.zeros([Rep, len(Nstorm), len(Rmin)])
SC = np.zeros([Rep, len(Nstorm), len(Rmin)])
Hd = np.zeros([Rep, len(Nstorm), len(Rmin)])
Punctuated = np.zeros([Rep, len(Nstorm), len(Rmin)])
SCperiod = np.zeros([Rep, len(Nstorm), len(Rmin)])
IslArea = np.zeros([Rep, len(Nstorm), len(Rmin)])
IntElev = np.zeros([Rep, len(Nstorm), len(Rmin)])
Pspace1 = np.zeros([Rep, len(Nstorm), len(Rmin)])
Pspace2 = np.zeros([Rep, len(Nstorm), len(Rmin)])
Pspace3 = np.zeros([Rep, len(Nstorm), len(Rmin)])


# Loop through replicate sets
for r in range(Rep):
    add = r*simnum
    # Loop through parameter space
    for sim in range(add, simnum + add):
        
        G = BatchData[sim, 1]
        N = BatchData[sim, 2]
        PP = BatchData[sim, 3]
        R = BatchData[sim, 4]
        T = BatchData[sim, 5]
        W = BatchData[sim, 6]
        S = BatchData[sim, 7]
        H = BatchData[sim, 8]
        P = BatchData[sim, 11]
        A = BatchData[sim, 10]
        E = BatchData[sim, 9]
        V1 = BatchData[sim, 12]
        V2 = BatchData[sim, 13]
       
        if G == 0:
            g = 0
        else:
            g = Rmin.index(round(G - Rspread, 2))
        n = Nstorm.index(N)
                 
        Result[r,n,g] = R       # Correct Version
        Tresult[r,n,g] = T      # Only average values for runs that drowned
        Width[r,n,g] = W
        SC[r,n,g] = S/T         # Rate 
        Hd[r,n,g] = H
        Punctuated[r,n,g] = PP
        SCperiod[r,n,g] = P     # Only average values for runs that are puncuated
        IslArea[r,n,g] = A
        IntElev[r,n,g] = E
        Pspace1[r,n,g] = V1     # Only average if values are greater than 0
        Pspace2[r,n,g] = V2     # Only average if values are greater than 0


# Take average of replicates       
Result2 = np.mean(Result, axis=0)
Tresult = np.nanmean(np.where(Result!=0,Tresult,np.nan),axis=0)
Width = np.mean(Width, axis=0)
SC = np.mean(SC, axis=0)                                    
Hd = np.mean(Hd, axis=0)
Punctuated2 = np.mean(Punctuated, axis=0)
SCperiod = np.true_divide(SCperiod.sum(0),(SCperiod!=0).sum(0))
IslArea = np.mean(IslArea, axis=0)
IntElev = np.mean(IntElev, axis=0)
Pspace1 = np.true_divide(Pspace1.sum(0),(Pspace1!=0).sum(0))
Pspace2 = np.true_divide(Pspace2.sum(0),(Pspace2!=0).sum(0))



#===============================
### Plot

xtic = ['', '0.30', '0.45', '0.60', '0.75', '0.90'] 
ytic =  ['', '4', '6', '8', '10', '12'] 

xlab = 'Mean Dune Growth Rate'
ylab = 'Average Storms Per Year' #'Mean Surge-Tide'


saveloc = 'C:\Barrier3D\Output'

# Result
PSFig = plt.figure(figsize=(14,8))
plt.rcParams.update({'font.size':13})
ax = PSFig.add_subplot(111)
cax = ax.matshow(Result2, origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Probability', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Drowning Probability')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Result') 

# Tresult
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Tresult, origin='lower', cmap='Greys', aspect='auto')
none = np.argwhere(np.isnan(Tresult))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Years', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Time to Drowning')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_TimeToDrown') 

# Width
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Width, origin='lower', cmap='Greys', aspect='auto')#, vmin=0, vmax=450)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Meters', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Island Width')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Width') 

# SC
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(SC, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.25, vmax=1)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Meters $Year^{-1}$', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Shoreline Change')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_ShorelineChange') 

# Hd
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Hd, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.1, vmax=1.5)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Meters', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Mean Dune Height')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_DuneHeight') 

# IntElev
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(IntElev, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.1, vmax=1.7)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Meters', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Mean Interior Elevation')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Elevation') 

# Punctuated
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Punctuated2, origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Probability', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Punctuated Retreat Probability')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Punctuated') 

# Characteristic Timescale
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(SCperiod, origin='lower', cmap='Greys', aspect='auto')#, vmin=80, vmax=420)
none = np.argwhere(np.isnan(SCperiod))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Years', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Characteristic Timescale')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Timescale') 

## IslArea
#PSFig = plt.figure(figsize=(14,8))
#ax = PSFig.add_subplot(111)
#cax = ax.matshow(IslArea, origin='lower', cmap='Greys', aspect='auto')
#ax.xaxis.set_ticks_position('bottom')
#cbar = PSFig.colorbar(cax) 
#cbar.set_label('Cubic Meters', rotation=270, labelpad=20)
#ax.set_xticklabels(xtic)
#ax.set_yticklabels(ytic)
#plt.xlabel(xlab)
#plt.ylabel(ylab)
#plt.title('Island Area')
#if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Area') 

# SlowDur
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Pspace1, origin='lower', cmap='Greys', aspect='auto')#, vmin=10, vmax=90)
none = np.argwhere(np.isnan(Pspace1))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Years', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Average Duration of Slow Periods')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_SlowDur') 

# FastDur
PSFig = plt.figure(figsize=(14,8))
ax = PSFig.add_subplot(111)
cax = ax.matshow(Pspace2, origin='lower', cmap='Greys', aspect='auto')#, vmin=20, vmax=400)
none = np.argwhere(np.isnan(Pspace2))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
cbar = PSFig.colorbar(cax) 
cbar.set_label('Years', rotation=270, labelpad=20)
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title('Average Duration of Fast Periods')
if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_FastDur') 


        
        
        



