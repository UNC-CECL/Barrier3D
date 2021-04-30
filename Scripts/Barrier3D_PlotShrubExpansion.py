# Script for comparing shrub vs no shrub simulations

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 1
# Updated: 13 November 2020


import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt


#%% Initialize

Plot = True
Save = False

#=============================
paths = ['C:\Barrier3D\Output\BatchSims_2021_0312_0016',
            'C:\Barrier3D\Output\BatchSims_2021_0312_0018',
            'C:\Barrier3D\Output\BatchSims_2021_0319_1101',
            'C:\Barrier3D\Output\BatchSims_2021_0319_1106']

#=============================
PSNums = [25,25,25,25]

#=============================
Var = 5 # Number of variables
SimDur = 150

#=============================
SimNum = np.sum(PSNums) * Var

Front = np.zeros([np.sum(PSNums),5,SimDur+1])
Area_Norm = np.zeros([np.sum(PSNums),5,SimDur+1])
Area = np.zeros([np.sum(PSNums),5,SimDur+1])
Dune = np.zeros([np.sum(PSNums),5,SimDur+1])
Width = np.zeros([np.sum(PSNums),5,SimDur+1])
Arrested = np.zeros([np.sum(PSNums),5])
Arrested_tstart = np.zeros([np.sum(PSNums),5])

T_Full = np.zeros([np.sum(PSNums),5])
AvgTimeFull = np.zeros([np.sum(PSNums),5])
AvgExpRate = np.zeros([np.sum(PSNums),5])
Qow_cumul_Full = np.zeros([np.sum(PSNums),5])
avgHd_Full = np.zeros([np.sum(PSNums),5])
avgShrubArea_Full = np.zeros([np.sum(PSNums),5])


#%% LOAD DATA 

for q in range(len(PSNums)):
    
    print(paths[q])
    
    for n in range(PSNums[q]):
    
        for Sim in range(Var):
            
            add = Var * n
           
            # Create File Path
            filename = paths[q] + '/SimData_' + str(Sim+1+add) + '.npz'
            
            # Load File
            SimData = np.load(filename, allow_pickle = True)
            
            # Load  Data
            x_s_TS = SimData['x_s_TS']
            x_b_TS = SimData['x_b_TS']
            InteriorWidth_AvgTS = SimData['InteriorWidth_AvgTS']
            AvgInteriorElevationTS = SimData['AvgInteriorElevationTS']
            DomainTS = SimData['DomainTS']
            QowTS = SimData['QowTS']
            Hd_AverageTS = SimData['Hd_AverageTS']
            SimParams = SimData['SimParams']
            RSLR = SimParams[1]
            BermEl = SimParams[3] 
            DShoreface = SimParams[6] 
            LShoreface = SimParams[7] 
            Shrub_ON = SimParams[8]
            TMAX = len(x_s_TS)
            PercentCoverTS = SimData['PercentCoverTS']
            ShrubArea = SimData['ShrubArea']
            ShrubDomainAll = SimData['ShrubDomainAll']
            DeadPercentCoverTS = SimData['DeadPercentCoverTS']
            ShrubDeadTS = SimData['ShrubDeadTS']

            
                
    #%% CALCULATE & STORE DATA    
            
            if q > 0:
                addy = np.sum(PSNums[:q])
                row = addy + n
            else:
                row = n
    
            # Shrub Dead Area
            ShrubArea_Dead = []
            for t in range(TMAX):
                DeadArea = np.count_nonzero(ShrubDeadTS[t])
                ShrubArea[t] += DeadArea # Add dead shrubs area to alive shrub area!
                ShrubArea_Dead.append(DeadArea)
    
            # Shrub Front: Defined as the left 95% of all shrub cells where age > 1        
            percentile = 0.95    
            # Initialize
            ShrubFront = []
            FrontLoc = 0 
            ArrestedYears = 0
            ArrestedYears_tstart = 0
            # Loop through each year
            for t in range(TMAX):        
                shrub = PercentCoverTS[t]
                shrub[shrub < 0.05] = 0 # Only counts shrubs older than 1 year
                ShrubTotal = np.count_nonzero(shrub) # Total amount of shrub plants in given year 
                domainW = len(shrub)
                domainL = len(shrub[0])
                endshrub = int(ShrubTotal * percentile)
                # Loop through each shrub cell by cross-shore column first to identify which column the the 95th% shrub is in
                count = 0 # Initialize
                for l in range(domainL):
                    for w in range(domainW):
                        if shrub[w,l] > 0:
                            count +=1
                            if count == endshrub:
                                FrontLoc = l # Column where 95th% shrub is in
                

                if t > 0 and FrontLoc <= ShrubFront[-1] and FrontLoc <= 450:
                    ArrestedYears += 1
                    if FrontLoc >= 10:
                        ArrestedYears_tstart += 1
                ShrubFront.append(FrontLoc) # Store in array           
            ShrubFront = [i * 10 for i in ShrubFront]      
            # Store
            Front[row,Sim,:] = ShrubFront
            Arrested[row,Sim] = ArrestedYears
            Arrested_tstart[row,Sim] = ArrestedYears_tstart
            
            # Shrub area normalized by barrier area
            ShrubArea_Norm = []
            AreaShrub = []
            DuneAvg = []
            Wid = []
            for t in range(TMAX):
                ShrubAr = ShrubArea[t] 
                Barr = DomainTS[t]
                Barr[Barr < 0] = 0
                BarrierAr = np.count_nonzero(Barr)
                ShrubArea_Norm.append(ShrubAr / BarrierAr)
                AreaShrub.append(ShrubAr)
                DuneAvg.append(Hd_AverageTS[t])
                Wid.append(InteriorWidth_AvgTS[t])
            # Store
            Area_Norm[row,Sim,:] = ShrubArea_Norm   
            Area[row,Sim,:] = AreaShrub
            Dune[row,Sim,:] = DuneAvg
            Width[row,Sim,:] = Wid
           
            #Cumulative overwash flux
            Qow_cumul = np.cumsum(QowTS)
            
            # Average Dune Height
            aHd = [a * 10 for a in Hd_AverageTS] # Convert to m
            
            # Average Expansion Rate Until Full (or Sim End)
            t_start = np.argmax(np.asarray(ShrubFront) > 100) # Could also calculate by area
            t_full = np.argmax(np.asarray(ShrubFront) > 4500)
            if t_full == 0:
                t_full = SimDur
            
            T_Full[row,Sim] = t_full
            AvgTimeFull[row,Sim] = ShrubFront[t_full]/t_full
            AvgExpRate[row,Sim] = ShrubFront[t_full]/(t_full-t_start)
            Qow_cumul_Full[row,Sim] = Qow_cumul[t_full]
            avgHd_Full[row,Sim] = np.mean(aHd[:t_full])
            avgShrubArea_Full[row,Sim] = np.mean(ShrubArea[:t_full])
            
           
#%% PLOT 

lab = ['0.3', '0.45', '0.6', '0.75', '0.9']
labtype = 'Characteristic Dune Growth Rate'
# lab = ['3', '6', '9', '12', '15'] 
# labtype = 'RSLR (mm/yr)'

# Box plots
ExpansionBox = [AvgExpRate[:,0], AvgExpRate[:,1], AvgExpRate[:,2], AvgExpRate[:,3], AvgExpRate[:,4]]
TimeRateBox = [AvgTimeFull[:,0], AvgTimeFull[:,1], AvgTimeFull[:,2], AvgTimeFull[:,3], AvgTimeFull[:,4]]
QowBox = [Qow_cumul_Full[:,0], Qow_cumul_Full[:,1], Qow_cumul_Full[:,2], Qow_cumul_Full[:,3], Qow_cumul_Full[:,4]]
HdBox = [avgHd_Full[:,0], avgHd_Full[:,1], avgHd_Full[:,2], avgHd_Full[:,3], avgHd_Full[:,4]]
ShrubCoverBox = [avgShrubArea_Full[:,0], avgShrubArea_Full[:,1], avgShrubArea_Full[:,2], avgShrubArea_Full[:,3], avgShrubArea_Full[:,4]]
ArrestedBox = [Arrested[:,0], Arrested[:,1], Arrested[:,2], Arrested[:,3], Arrested[:,4]]
ArrestedBox_tstart = [Arrested_tstart[:,0], Arrested_tstart[:,1], Arrested_tstart[:,2], Arrested_tstart[:,3], Arrested_tstart[:,4]]



#%%
# For Publication



ExpRateAll = np.zeros([np.sum(PSNums)*Var])
OwAll = np.zeros([np.sum(PSNums)*Var])
HdAll = np.zeros([np.sum(PSNums)*Var])
for r in range(Var):
    minny = np.sum(PSNums)*r
    maxxy = minny + np.sum(PSNums)
    ExpRateAll[minny:maxxy] = AvgExpRate[:,r]
    OwAll[minny:maxxy] = Qow_cumul_Full[:,r]
    HdAll[minny:maxxy] = avgHd_Full[:,r]


# All
plt.figure()   
fig = plt.gcf()
fig.set_size_inches(16,15)
plt.rcParams.update({'font.size':15})
color = ['darkviolet', 'blue', 'green', 'gold', 'red']

ax = fig.add_subplot(3,2,1)
for s in range(Var): 
    for n in range(np.sum(PSNums)):
        plt.scatter(AvgExpRate[n,s], Qow_cumul_Full[n,s], c=color[s])  
    plt.xlabel('Average Expansion Rate (m/yr)')
    plt.ylabel('Cumulative Overwash (m^3/m)')
# # Linear Regression
# slope, intercept, r_value, p_value, std_err = st.linregress(ExpRateAll, OwAll)
# xl = [min(ExpRateAll), max(ExpRateAll)]
# yl = [slope*xx + intercept for xx in xl]
# plt.plot(xl,yl, '-k')
# print('R-squared', r_value**2)    
    
ax = fig.add_subplot(3,2,2)
for s in range(Var): 
    for n in range(np.sum(PSNums)):
        plt.scatter(AvgExpRate[n,s], avgHd_Full[n,s], c=color[s])  
    plt.xlabel('Average Expansion Rate (m/yr)')
    plt.ylabel('Average Dune Height (m)')
# # Linear Regression
# slope, intercept, r_value, p_value, std_err = st.linregress(ExpRateAll, HdAll)
# xl = [min(ExpRateAll), max(ExpRateAll)]
# yl = [slope*xx + intercept for xx in xl]
# plt.plot(xl,yl, '-k')
# print('R-squared', r_value**2)  


ax = fig.add_subplot(3,2,3)
plt.boxplot(ExpansionBox, labels=lab)
plt.ylabel('Average Expansion Rate (m/yr)')
plt.xlabel(labtype)

ax = fig.add_subplot(3,2,6)
plt.boxplot(QowBox, labels=lab)
plt.ylabel('Cumulative Overwash Flux (m^3/m)')
plt.xlabel(labtype)

ax = fig.add_subplot(3,2,5)
plt.boxplot(HdBox, labels=lab)
plt.ylabel('Average Dune Height (m)')
plt.xlabel(labtype)

ax = fig.add_subplot(3,2,4)
plt.boxplot(ArrestedBox_tstart, labels=lab)
plt.ylabel('Years of No Expansion')
plt.xlabel(labtype)

plt.tight_layout(pad=1.5)
plt.show()


