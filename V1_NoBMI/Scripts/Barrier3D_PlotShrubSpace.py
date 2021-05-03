# Parameter space plotting script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 2
# Updated: 9 Apr 2020

# Creates a parameter space figure from saved numpy BatchData array for shrub/no shrub runs

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt



#===============================
### Inputs

Save = False # Whether or not to save figures


paths = ['C:/Barrier3D/Output/BatchSims_2021_0220_1452/BatchData_2021_0220_1452.npy',
          'C:/Barrier3D/Output/BatchSims_2021_0220_2005/BatchData_2021_0220_2005.npy',
          'C:/Barrier3D/Output/BatchSims_2021_0225_1615/BatchData_2021_0225_1615.npy',
          'C:/Barrier3D/Output/BatchSims_2021_0225_1617/BatchData_2021_0225_1617.npy']



Rep = [25, 25, 25, 25]


sealevel = False
slr = [3, 6, 9, 12, 15] # mm/yr
ravg = [0.3, 0.45, 0.6, 0.75, 0.9]
shrub = [0,1]


TMAX = 1000


Reps = np.sum(Rep)
simnum = len(slr) * len(shrub) 



#===============================
### Load Data Into Arrays

Result = np.zeros([Reps, len(slr), len(shrub)])
Tresult = np.zeros([Reps, len(slr), len(shrub)])
Width = np.zeros([Reps, len(slr), len(shrub)])
SC = np.zeros([Reps, len(slr), len(shrub)])
Hd = np.zeros([Reps, len(slr), len(shrub)])
Punctuated = np.zeros([Reps, len(slr), len(shrub)])
SCperiod = np.zeros([Reps, len(slr), len(shrub)])
ShrubArea = np.zeros([Reps, len(slr), len(shrub)])
IntElev = np.zeros([Reps, len(slr), len(shrub)])
Pspace1 = np.zeros([Reps, len(slr), len(shrub)])
Pspace2 = np.zeros([Reps, len(slr), len(shrub)])
Pspace3 = np.zeros([Reps, len(slr), len(shrub)])


# Loop through replicate sets
for n in range(len(Rep)):
    loc = paths[n]
    BatchData = np.load(loc, allow_pickle = True) 
    
    for r in range(Rep[n]):
        add = r*simnum
        # Loop through parameter space
        for sim in range(add, simnum + add):
            
            Z_1 = BatchData[sim, 1]
            SH = BatchData[sim, 2]
            PP = BatchData[sim, 3]
            Z_2 = BatchData[sim, 4]
            T = BatchData[sim, 5]
            W = BatchData[sim, 6]
            S = BatchData[sim, 7]
            H = BatchData[sim, 8]
            E = BatchData[sim, 9]
            A = BatchData[sim, 10]
            P = BatchData[sim, 11]        
            SD = BatchData[sim, 12]
            FD = BatchData[sim, 13]
    
            if sealevel:
                z = slr.index(Z_1)
            else:
                z = ravg.index(round(Z_2,2))
            
            s = int(SH)
            
            radd = int(r +  np.sum(Rep[:n]))      
            if T >= 1000:
                Result[radd,z,s] = 0
            else:
                Result[radd,z,s] = 1
            Tresult[radd,z,s] = T      # Only average values for runs that drowned
            Width[radd,z,s] = W
            SC[radd,z,s] = S/T         # Rate 
            Hd[radd,z,s] = H
            Punctuated[radd,z,s] = PP
            SCperiod[radd,z,s] = P     # Only average values for runs that are puncuated
            ShrubArea[radd,z,s] = A
            IntElev[radd,z,s] = E
            Pspace1[radd,z,s] = SD     # Only average if values are greater than 0
            Pspace2[radd,z,s] = FD     # Only average if values are greater than 0


# Take average of replicates       
Result2 = np.mean(Result, axis=0)
Tresult = np.nanmean(np.where(Result!=0,Tresult,np.nan),axis=0)
Width = np.mean(Width, axis=0)
SC = np.mean(SC, axis=0)                                    
Hd = np.mean(Hd, axis=0)
Punctuated2 = np.mean(Punctuated, axis=0)
SCperiod = np.true_divide(SCperiod.sum(0),(SCperiod!=0).sum(0))
ShrubArea = np.mean(ShrubArea, axis=0)
IntElev = np.mean(IntElev, axis=0)

Pspace1copy = Pspace1
Pspace2copy = Pspace2
Pspace1 = np.true_divide(Pspace1.sum(0),(Pspace1!=0).sum(0))
Pspace2 = np.true_divide(Pspace2.sum(0),(Pspace2!=0).sum(0))

# Standard Deviations
Pspace1copy[Pspace1copy == 0] = np.nan
Pspace2copy[Pspace2copy == 0] = np.nan
Pspace1_StD = np.nanstd(Pspace1copy, axis=0)
Pspace2_StD = np.nanstd(Pspace2copy, axis=0)




#===============================
### Plot
xtic = ['', 'No Shrub', 'Shrub'] 
xlab = ''
if sealevel:
    ytic =  ['', '3', '6', '9', '12', '15']   
    ylab = 'RSLR (mm/yr)' 
else:
    ytic = ['', '0.3', '0.45', '0.6', '0.75', '0.9']
    ylab = 'Mean Dune Growth Rate' 

saveloc = 'C:\Barrier3D\Output'


# #===============================
# # Phase Spaces

# # Result
# PSFig = plt.figure(figsize=(14,8))
# plt.rcParams.update({'font.size':13})
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Result2, origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Probability', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Drowning Probability')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Result') 

# # Tresult
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Tresult, origin='lower', cmap='Greys', aspect='auto')
# none = np.argwhere(np.isnan(Tresult))
# for j in none:
#     y, x = j
#     ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Years', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Time to Drowning')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_TimeToDrown') 

# # Width
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Width, origin='lower', cmap='Greys', aspect='auto')#, vmin=0, vmax=450)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Meters', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Island Width')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Width') 

# # SC
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(SC, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.25, vmax=1)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Meters $Year^{-1}$', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Shoreline Change')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_ShorelineChange') 

# # Hd
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Hd, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.1, vmax=1.5)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Meters', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Mean Dune Height')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_DuneHeight') 

# # IntElev
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(IntElev, origin='lower', cmap='Greys', aspect='auto')#, vmin=0.1, vmax=1.7)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Meters', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Mean Interior Elevation')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Elevation') 

# # ShrubArea
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(ShrubArea, origin='lower', cmap='Greys', aspect='auto')
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Cubic Meters', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Shrub Area')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Area') 

# # Punctuated
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Punctuated2, origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Probability', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Punctuated Retreat Probability')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Punctuated') 

# # Characteristic Timescale
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(SCperiod, origin='lower', cmap='Greys', aspect='auto')#, vmin=80, vmax=420)
# none = np.argwhere(np.isnan(SCperiod))
# for j in none:
#     y, x = j
#     ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Years', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Characteristic Timescale')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_Timescale') 

# # SlowDur
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Pspace1, origin='lower', cmap='Greys', aspect='auto')#, vmin=10, vmax=90)
# none = np.argwhere(np.isnan(Pspace1))
# for j in none:
#     y, x = j
#     ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Years', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Average Duration of Slow Periods')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_SlowDur') 

# # FastDur
# PSFig = plt.figure(figsize=(14,8))
# ax = PSFig.add_subplot(111)
# cax = ax.matshow(Pspace2, origin='lower', cmap='Greys', aspect='auto')#, vmin=20, vmax=400)
# none = np.argwhere(np.isnan(Pspace2))
# for j in none:
#     y, x = j
#     ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
# ax.xaxis.set_ticks_position('bottom')
# cbar = PSFig.colorbar(cax) 
# cbar.set_label('Years', rotation=270, labelpad=20)
# ax.set_xticklabels(xtic)
# ax.set_yticklabels(ytic)
# plt.xlabel(xlab)
# plt.ylabel(ylab)
# plt.title('Average Duration of Fast Periods')
# if Save == True: PSFig.savefig(saveloc + 'PhaseSpace_FastDur') 


        

#===============================
# Line Plots

if sealevel:
    xtic =  ['', '3', '6', '9', '12', '15']   
    xlab = 'RSLR (mm/yr)' 
else:
    xtic = ['', '0.3', '0.45', '0.6', '0.75', '0.9']
    xlab = 'Mean Dune Growth Rate' 

        
fig = plt.figure(figsize=(20,20))
plt.rcParams.update({'font.size':18})

# 1
# ax = fig.add_subplot(4,4,1)
# ax.scatter(SCrate_rm, aHd_rm, marker='o', s=10, c='black', edgecolors='none')
# ax.set_ylim(0, (Dmaxel - BermEl) * 10 + 0.1)
# plt.ylabel('DuneHeight (m)')


# Interior Width
ax = fig.add_subplot(4,2,1)
ax.plot(Width[:,0], color='r')
ax.plot(Width[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Island Width (m)')


# Shoreline Change
ax = fig.add_subplot(4,2,2)
plt.plot(SC[:,0], color='r')
plt.plot(SC[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Shoreline Change (m/yr)')


# Dune Height
ax = fig.add_subplot(4,2,3)
plt.plot(Hd[:,0], color='r')
plt.plot(Hd[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Dune Height (m)')


# Interior Elevation
ax = fig.add_subplot(4,2,4)
plt.plot(IntElev[:,0], color='r')
plt.plot(IntElev[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Interior Elevation (m)')


# Shrub Area
ax = fig.add_subplot(4,2,5)
plt.plot(ShrubArea[:,0], color='r')
plt.plot(ShrubArea[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Shrub Area (m^2)')


# Punctuated
ax = fig.add_subplot(4,2,6)
plt.plot(Punctuated2[:,0], color='r')
plt.plot(Punctuated2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Punctuated Probability')


# Slow Dur
ax = fig.add_subplot(4,2,7)
plt.plot(Pspace1[:,0], color='r')
plt.plot(Pspace1[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Slow Duration (yrs)')


# Fast Dur
ax = fig.add_subplot(4,2,8)
plt.plot(Pspace2[:,0], color='r')
plt.plot(Pspace2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Fast Duration (yrs)')

plt.show()

if Save == True: fig.savefig(saveloc + 'ShrubLine') 



fig = plt.figure(figsize=(15,10))
plt.rcParams.update({'font.size':18})

### Punctuated Retreat
# Timescale
ax = fig.add_subplot(2,2,2)
plt.plot(SCperiod[:,0], color='b')
plt.plot(SCperiod[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Characteristic Timescale (yrs)')

# Punctuated
ax = fig.add_subplot(2,2,1)
plt.plot(Punctuated2[:,0], color='b')
plt.plot(Punctuated2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Punctuated Probability')

# Slow Dur
ax = fig.add_subplot(2,2,3)
plt.plot(Pspace1[:,0], color='b')
plt.plot(Pspace1[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Slow Duration (yrs)')

# Fast Dur
ax = fig.add_subplot(2,2,4)
plt.plot(Pspace2[:,0], color='b')
plt.plot(Pspace2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Fast Duration (yrs)')

plt.show()



# Select
fig = plt.figure(figsize=(20,4))
plt.rcParams.update({'font.size':15})

# Punctuated
ax = fig.add_subplot(1,3,1)
plt.plot(Punctuated2[:,0], color='b')
plt.plot(Punctuated2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Punctuated Probability')

# Slow Dur
ax = fig.add_subplot(1,3,2)
plt.plot(Pspace1[:,0], color='b')
plt.plot(Pspace1[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Slow Duration (yrs)')

# Fast Dur
ax = fig.add_subplot(1,3,3)
plt.plot(Pspace2[:,0], color='b')
plt.plot(Pspace2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Fast Duration (yrs)')
plt.tight_layout()

plt.show()




# Select Error
if sealevel:
    x = [3,6,9,12,15]
else:
    x = [0.3,0.45,0.6,0.75,0.9]


fig = plt.figure(figsize=(20,4))
plt.rcParams.update({'font.size':15})

# Punctuated
ax = fig.add_subplot(1,3,1)
plt.plot(Punctuated2[:,0], color='b')
plt.plot(Punctuated2[:,1], color='g')
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Punctuated Probability')

# Slow Dur
ax = fig.add_subplot(1,3,2)
plt.errorbar(x, Pspace1[:,0], yerr=Pspace1_StD[:,0], color='b', capsize=5, capthick=2)
plt.errorbar(x, Pspace1[:,1], yerr=Pspace1_StD[:,0], color='g', capsize=5, capthick=2)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Slow Duration (yrs)')

# Fast Dur
ax = fig.add_subplot(1,3,3)
plt.errorbar(x, Pspace2[:,0], yerr=Pspace2_StD[:,0], color='b', capsize=5, capthick=2)
plt.errorbar(x, Pspace2[:,1], yerr=Pspace2_StD[:,0], color='g', capsize=5, capthick=2)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xticklabels(xtic)
plt.xlabel(xlab)
plt.ylabel('Fast Duration (yrs)')
plt.tight_layout()

plt.show()


