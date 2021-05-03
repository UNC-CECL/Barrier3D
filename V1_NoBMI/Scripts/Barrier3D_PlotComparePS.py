# Parameter space plotting script for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 1
# Updated: 17 Aug 2020

# Creates a figures of three different parameter spaces from batch run spreadsheet data

import numpy as np
import xlrd
import matplotlib as mpl
import matplotlib.pyplot as plt



#%% Inputs

# Save Figures
Save = False 
saveloc = 'C:/Barrier3D/Output'

# Three parameter spaces - Ksf
Group = [35,30,36] # Give in order of plotting, left to right
Rep = [100,100,100] # Number of duplicates
title1 = '$k_{sf}$ = 1000'
title2 = '$k_{sf}$ = 5000'
title3 = '$k_{sf}$ = 10000'

# # Three parameter spaces - RSLR
# Group = [32,30,31] # Give in order of plotting, left to right
# Rep = [100,100,100] # Number of duplicates
# title1 = 'RSLR = 2'
# title2 = 'RSLR = 4'
# title3 = 'RSLR = 8'

# # Three parameter spaces - Average TWL
# Group = [34,30,33] # Give in order of plotting, left to right
# Rep = [100,100,100] # Number of duplicates
# title1 = 'Average TWL - 0.15 m'
# title2 = 'Average TWL'
# title3 = 'Average TWL + 0.15 m'

# # Three parameter spaces - Punctuated min duration
# Group = [32,32,32] # Give in order of plotting, left to right
# Rep = [100,100,100] # Number of duplicates
# title1 = 'MinDur = 20 yrs'
# title2 = 'MinDur = 30 yrs'
# title3 = 'MinDur = 40 yrs'

# # Three parameter spaces - Punctuated rate threshold
# Group = [32,32,32] # Give in order of plotting, left to right
# Rep = [100,100,100] # Number of duplicates
# title1 = 'Threhold = 0.25 m/yr'
# title2 = 'Threhold = 0.50 m/yr'
# title3 = 'Threhold = 0.75 m/yr'

# Range
Rmin = [0.05, 0.20, 0.35, 0.50, 0.65]
Rmax = [0.55, 0.70, 0.85, 1.00, 1.15]
Nstorm = [4, 6, 8, 10, 12]

# Plotting Labels
xtic = ['', '0.30', '0.45', '0.60', '0.75', '0.90'] 
ytic =  ['', '4', '6', '8', '10', '12'] 

xlab = 'Mean Dune Growth Rate'
ylab = 'Average Storms Per Year'

# Calulations
Rspread = 0.25 #(Rmax[0] - Rmin[0]) / 2
simnum = len(Rmin) * len(Nstorm)

# Initialize
Result = []
Tresult = []
Width = []
SC = [] 
Hd = []
Punctuated = []
SCperiod = []
IslArea = []
IntElev = []
Pspace1 = []
Pspace2 = []



#%% Load Data Into Arrays

for x in range(3):
    
    if x == 1:
        loc = ('C:\Barrier3D\Output\BatchSims_Group' + str(Group[x]) + '\BatchData_Combined_Group' + str(Group[x]) + '\Data_BatchSims_Combined_Group' + str(Group[x]) + '.xlsx')
    else:
        loc = ('E:\Reeves\Backups\Barrier3D_Data\Runs_122-156\BatchData_Combined_Group' + str(Group[x]) + '\Data_BatchSims_Combined_Group' + str(Group[x]) + '.xlsx')
    
    
    
    wb = xlrd.open_workbook(loc)
    s1 = wb.sheet_by_index(0) 
    
    allResult = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allTresult = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allWidth = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allSC = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allHd = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allPunctuated = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allSCperiod = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allIslArea = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allIntElev = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allPspace1 = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allPspace2 = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    allPspace3 = np.zeros([Rep[x], len(Rmin), len(Nstorm)])
    
    # Loop through replicate sets
    for r in range(Rep[x]):
        add = r*simnum
        # Loop through parameter space
        for sim in range(1 + add, simnum + 1 + add):
            
            G = s1.cell_value(sim, 1) 
            N = s1.cell_value(sim, 2)
            PP = s1.cell_value(sim, 3) #3 / 15 / 19 / 23 / 27 / 31 / 35
            R = s1.cell_value(sim, 4) 
            T = s1.cell_value(sim, 5)
            W = s1.cell_value(sim, 6)
            S = s1.cell_value(sim, 7)
            H = s1.cell_value(sim, 8)
            P = s1.cell_value(sim, 11) #11 / 16 / 20 / 24 / 28 / 32 / 36
            A = s1.cell_value(sim, 10)
            E = s1.cell_value(sim, 9)
            V1 = s1.cell_value(sim, 12) #12 / 17 / 21 / 25 / 29 / 33 / 37
            V2 = s1.cell_value(sim, 13) #13 / 18 / 22 / 26 / 30 / 34 / 38
            
            
            if G == 0:
                g = 0
            else:
                g = Rmin.index(round(G - Rspread, 2))
            n = Nstorm.index(N)
         
            allResult[r,n,g] = R
            allTresult[r,n,g] = T      # Only average values for runs that drowned
            allWidth[r,n,g] = W
            allSC[r,n,g] = S/T         # Rate
            allHd[r,n,g] = H
            allPunctuated[r,n,g] = PP
            allSCperiod[r,n,g] = P     # Only average values for runs that are puncuated
            allIslArea[r,n,g] = A
            allIntElev[r,n,g] = E
            allPspace1[r,n,g] = V1     # Only average if values are greater than 0
            allPspace2[r,n,g] = V2     # Only average if values are greater than 0
    
    # Take average of replicates       
    Result.append(np.mean(allResult, axis=0))
    Tresult.append(np.nanmean(np.where(allResult!=0,allTresult,np.nan),axis=0))
    Width.append(np.mean(allWidth, axis=0))
    SC.append(np.mean(allSC, axis=0) )
    Hd.append(np.mean(allHd, axis=0))
    Punctuated.append(np.mean(allPunctuated, axis=0))
    SCperiod.append(np.true_divide(allSCperiod.sum(0),(allSCperiod!=0).sum(0)))
    IslArea.append(np.mean(allIslArea, axis=0))
    IntElev.append(np.mean(allIntElev, axis=0))
    Pspace1.append(np.true_divide(allPspace1.sum(0),(allPspace1!=0).sum(0)))
    Pspace2.append(np.true_divide(allPspace2.sum(0),(allPspace2!=0).sum(0)))



#%% Plot
    

#=====================================
# Result
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Result[0], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Result[1], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Result[2], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Probability', rotation=270, labelpad=20)
PSFig.suptitle('Drowning Probability', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_Result', bbox_inches='tight') 


#=====================================
# Tresult
VMin = np.nanmin(Tresult)
VMax = np.nanmax(Tresult)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Tresult[0], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult[0]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Tresult[1], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult[1]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Tresult[2], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult[2]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Years', rotation=270, labelpad=20)
PSFig.suptitle('Time to Drowning', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_TimeToDrown', bbox_inches='tight')


#=====================================
# Width
VMin = np.nanmin(Width)
VMax = np.nanmax(Width)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Width[0], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Width[1], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Width[2], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Meters', rotation=270, labelpad=20)
PSFig.suptitle('Island Width', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_Width', bbox_inches='tight')


#=====================================
# SC
VMin = np.nanmin(SC)
VMax = np.nanmax(SC)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(SC[0], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(SC[1], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(SC[2], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Meters $Year^{-1}$', rotation=270, labelpad=20)
PSFig.suptitle('Shoreline Change', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_ShorelineChange', bbox_inches='tight')


#=====================================
# Hd
VMin = np.nanmin(Hd)
VMax = np.nanmax(Hd)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})
cmap0 = 'Reds'

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Hd[0], origin='lower', cmap=cmap0, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Hd[1], origin='lower', cmap=cmap0, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Hd[2], origin='lower', cmap=cmap0, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Meters', rotation=270, labelpad=20)
PSFig.suptitle('Mean Dune Height', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_DuneHeight', bbox_inches='tight')

PSFig.savefig(saveloc + 'PhaseSpace_DuneHeight.eps', format='eps', bbox_inches='tight')


#=====================================
# IntElev
VMin = np.nanmin(IntElev)
VMax = np.nanmax(IntElev)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(IntElev[0], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(IntElev[1], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(IntElev[2], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Tresult))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Meters', rotation=270, labelpad=20)
PSFig.suptitle('Mean Interior Elevation', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_Elevation', bbox_inches='tight')


#=====================================
# Punctuated
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Punctuated[0], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Punctuated[1], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Punctuated[2], origin='lower', cmap='Greys', aspect='auto', vmin=0, vmax=1)
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Probability', rotation=270, labelpad=20)
PSFig.suptitle('Punctuated Retreat Probability', y=1.07, size=22)

plt.show()
if Save: 
    PSFig.savefig(saveloc + 'PhaseSpace_Punctuated', bbox_inches='tight') 

PSFig.savefig(saveloc + 'PhaseSpace_Punctuated.eps', format='eps', bbox_inches='tight')

#=====================================
# SCperiod
VMin = np.nanmin(SCperiod)
VMax = np.nanmax(SCperiod)
PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(SCperiod[0], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(SCperiod[0]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(SCperiod[1], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(SCperiod[1]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(SCperiod[2], origin='lower', cmap='Greys', aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(SCperiod[2]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Years', rotation=270, labelpad=20)
PSFig.suptitle('Characteristic Timescale', y=1.07, size=22)

plt.show()
if Save: PSFig.savefig(saveloc + 'PhaseSpace_Timescale', bbox_inches='tight')


#=====================================
# SloWDur
VMin = np.nanmin(Pspace1)
VMax = np.nanmax(Pspace1)

# VMin = 36
# VMax = 367
print('Vmin-S: ', VMin)
print('Vmax-S: ', VMax)
cm2 = 'Purples'

PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Pspace1[0], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace1[0]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Pspace1[1], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace1[1]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Pspace1[2], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace1[2]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Years', rotation=270, labelpad=20)
PSFig.suptitle('Average Duration of Slow Periods', y=1.07, size=22)

plt.show()
if Save: 
    PSFig.savefig(saveloc + 'PhaseSpace_SlowDur', bbox_inches='tight')

    PSFig.savefig(saveloc + 'PhaseSpace_SlowDur.eps', format='eps', bbox_inches='tight')

#=====================================
# FastDur
VMin = np.nanmin(Pspace2)
VMax = np.nanmax(Pspace2)

# VMin = 36
# VMax = 367
print('Vmin-F: ', VMin)
print('Vmax-F: ', VMax)

PSFig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size':16})

ax = PSFig.add_subplot(1,3,1)
Cax = ax.matshow(Pspace2[0], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace1[0]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.ylabel(ylab, labelpad=10)
plt.title(title1)

ax = PSFig.add_subplot(1,3,2)
Cax = ax.matshow(Pspace2[1], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace2[1]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.xlabel(xlab, labelpad=15)
plt.title(title2)

ax = PSFig.add_subplot(1,3,3)
Cax = ax.matshow(Pspace2[2], origin='lower', cmap=cm2, aspect='auto', vmin=VMin, vmax=VMax)
none = np.argwhere(np.isnan(Pspace2[2]))
for j in none:
    y, x = j
    ax.add_patch(mpl.patches.Rectangle((x-0.5, y-0.5), 1, 1, hatch='x', fill=False, snap=False, linewidth=0, alpha=0.25))
ax.xaxis.set_ticks_position('bottom')
ax.set_xticklabels(xtic)
ax.set_yticklabels(ytic)
plt.title(title3)

PSFig.subplots_adjust(right=0.82)
cbar_ax = PSFig.add_axes([0.85, 0.12, 0.02, 0.76])
cbar = PSFig.colorbar(Cax, cax=cbar_ax)
cbar.set_label('Years', rotation=270, labelpad=20)
PSFig.suptitle('Average Duration of Fast Periods', y=1.07, size=22)

plt.show()
if Save: 
    PSFig.savefig(saveloc + 'PhaseSpace_FastDur', bbox_inches='tight')

    PSFig.savefig(saveloc + 'PhaseSpace_FastDur.eps', format='eps', bbox_inches='tight')


