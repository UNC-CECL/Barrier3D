# Script for comparing time series of shrub vs no shrub simulations

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions

# Ian R.B. Reeves

# Version Number: 1
# Updated: 13 November 2020


import numpy as np
import matplotlib.pyplot as plt


#%% Initialize

Plot = True
Save = False
sealevel = True


#=============================
paths = ['C:/Barrier3D/Output/BatchSims_2021_0220_1452',
            'C:/Barrier3D/Output/BatchSims_2021_0220_2005',
            'C:/Barrier3D/Output/BatchSims_2021_0225_1615',
            'C:/Barrier3D/Output/BatchSims_2021_0225_1617']


#============================
PSNums = [25,25,25,25]


SimNum = np.sum(PSNums) * 10
Stats = np.zeros([SimNum, 4]) 


#============================
# Initialize

Drown_NS = 0
Drown_S = 0

# RSLR0
HdA0 = []
SCRA0 = []
bbSCRA0 = []
WidthA0 = []
SAA0 = []
xsA0 = []
ElevA0 = []
QowA0 = []
QowcA0 = []
VolA0 = []
xbA0 = []

HdB0 = []
SCRB0 = []
bbSCRB0 = []
WidthB0 = []
SAB0 = []
xsB0 = []
ElevB0 = []
QowB0 = []
QowcB0 = []
VolB0 = []
xbB0 = []

# RSLR1
HdA1 = []
SCRA1 = []
bbSCRA1 = []
WidthA1 = []
SAA1 = []
xsA1 = []
ElevA1 = []
QowA1 = []
QowcA1 = []
VolA1 = []
xbA1 = []

HdB1 = []
SCRB1 = []
bbSCRB1 = []
WidthB1 = []
SAB1 = []
xsB1 = []
ElevB1 = []
QowB1 = []
QowcB1 = []
VolB1 = []
xbB1 = []

# RSLR2
HdA2 = []
SCRA2 = []
bbSCRA2 = []
WidthA2 = []
SAA2 = []
xsA2 = []
ElevA2 = []
QowA2 = []
QowcA2 = []
VolA2 = []
xbA2 = []

HdB2 = []
SCRB2 = []
bbSCRB2 = []
WidthB2 = []
SAB2 = []
xsB2 = []
ElevB2 = []
QowB2 = []
QowcB2 = []
VolB2 = []
xbB2 = []

# RSLR3
HdA3 = []
SCRA3 = []
bbSCRA3 = []
WidthA3 = []
SAA3 = []
xsA3 = []
ElevA3 = []
QowA3 = []
QowcA3 = []
VolA3 = []
xbA3 = []

HdB3 = []
SCRB3 = []
bbSCRB3 = []
WidthB3 = []
SAB3 = []
xsB3 = []
ElevB3 = []
QowB3 = []
QowcB3 = []
VolB3 = []
xbB3 = []

# RSLR4
HdA4 = []
SCRA4 = []
bbSCRA4 = []
WidthA4 = []
SAA4 = []
xsA4 = []
ElevA4 = []
QowA4 = []
QowcA4 = []
VolA4 = []
xbA4 = []

HdB4 = []
SCRB4 = []
bbSCRB4 = []
WidthB4 = []
SAB4 = []
xsB4 = []
ElevB4 = []
QowB4 = []
QowcB4 = []
VolB4 = []
xbB4 = []


#%% LOAD DATA 

for q in range(len(paths)):
    
    print(paths[q])
    
    SimN = PSNums[q] * 10 + 1
    
    for Sim in range(1,SimN):
           
        # Create File Path
        filename = paths[q] + '/SimData_' + str(Sim) + '.npz'
    
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
        BarrierLength = SimParams[4] 
        DShoreface = SimParams[6] 
        LShoreface = SimParams[7] 
        Shrub_ON = SimParams[8]
        TMAX = len(x_s_TS)
    
        
        if Shrub_ON == 1:
            PercentCoverTS = SimData['PercentCoverTS']
            ShrubArea = SimData['ShrubArea']
            DeadPercentCoverTS = SimData['DeadPercentCoverTS']
            ShrubDeadTS = SimData['ShrubDeadTS']
        else:
            PercentCoverTS = 0
            ShrubArea = 0
            DeadPercentCoverTS = 0
            
                
        #%% STORE DATA    
            
        N = Sim % 10 - 1
        if N < 0: N = 9

        # Average Dune Height
        aHd = [a * 10 for a in Hd_AverageTS] # Convert to m
    
        # Average Island Width
        aW = [a * 10 for a in InteriorWidth_AvgTS] # Convert to m
        
        # Average Island Elevation
        aE = [a * 10 for a in AvgInteriorElevationTS] # Convert to m
    
        # Shorline Change & Change Rate
        scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS] # Convert to m
        bbscts = [(x - x_b_TS[0]) * 10 for x in x_b_TS] # Convert to m
        SCrate = [0]
        bbSCrate = [0]
        for k in range(1,len(scts)):
            SCrate.append(scts[k]- scts[k-1])   
        for k in range(1,len(bbscts)):
            bbSCrate.append(bbscts[k]- bbscts[k-1])
            
        # Island Volume
        volumeTS = []
        for t in range(TMAX):
            domain = DomainTS[t]
            vol = (np.sum(domain[domain > 0])*1000) / (BarrierLength*10) # Convert to m^3, Normalize alongshore
            volumeTS.append(vol)
        # volumeTS = [a * 1000 for a in volumeTS] # Convert to m^3
        
        # Cumulative overwash
        Qow_cumul = np.cumsum(QowTS)
        
        # Shrub Dead Area
        if Shrub_ON == 1:
            ShrubArea_Dead = []
            for t in range(TMAX):
                DeadArea = np.count_nonzero(ShrubDeadTS[t])
                ShrubArea[t] += DeadArea
                ShrubArea_Dead.append(DeadArea)
        
        # Count drowning
        if TMAX < 1000:
            if Shrub_ON == 1:
                Drown_S += 1
            else:
                Drown_NS += 1            
        
            AddZeros = [0] * (1001 - TMAX)
            volumeTS.extend(AddZeros)
            aE.extend(AddZeros)
            aW.extend(AddZeros)
            aHd.extend(AddZeros)
            SCrate.extend(AddZeros)
            bbSCrate.extend(AddZeros)
            if Shrub_ON == 1:
                AreaTemp = np.ndarray.tolist(ShrubArea)
                AreaTemp.extend(AddZeros)
                ShrubArea = np.asarray(AreaTemp)
            QowTemp = np.ndarray.tolist(QowTS)
            QowTemp.extend(AddZeros)
            QowTS = np.asarray(QowTemp) 
            
            # Cumulative stats
            Add = [scts[-1]] * (1001 - TMAX)
            scts.extend(Add)
            
            Add = [bbscts[-1]] * (1001 - TMAX)
            bbscts.extend(Add)
            
            Add = [Qow_cumul[-1]] * (1001 - TMAX)
            QowCumulTemp = np.ndarray.tolist(Qow_cumul)
            QowCumulTemp.extend(Add)
            Qow_cumul = np.asarray(QowCumulTemp)
            

        # RSLR0
        if N == 0:
            HdA0.append(aHd)
            SCRA0.append(SCrate)
            bbSCRA0.append(bbSCrate)
            WidthA0.append(aW)
            SAA0.append(0)
            xsA0.append(scts)
            ElevA0.append(aE)
            QowA0.append(QowTS)
            QowcA0.append(Qow_cumul)
            VolA0.append(volumeTS)
            xbA0.append(bbscts)
        
        elif N == 5: 
            HdB0.append(aHd)
            SCRB0.append(SCrate)
            bbSCRB0.append(bbSCrate)
            WidthB0.append(aW)
            SAB0.append(ShrubArea)
            xsB0.append(scts)
            ElevB0.append(aE)
            QowB0.append(QowTS)
            QowcB0.append(Qow_cumul)
            VolB0.append(volumeTS)
            xbB0.append(bbscts)
            
        # RSLR1
        elif N == 1:
            HdA1.append(aHd)
            SCRA1.append(SCrate)
            bbSCRA1.append(bbSCrate)
            WidthA1.append(aW)
            SAA1.append(0)
            xsA1.append(scts)
            ElevA1.append(aE)
            QowA1.append(QowTS)
            QowcA1.append(Qow_cumul)
            VolA1.append(volumeTS)
            xbA1.append(bbscts)
        
        elif N == 6: 
            HdB1.append(aHd)
            SCRB1.append(SCrate)
            bbSCRB1.append(bbSCrate)
            WidthB1.append(aW)
            SAB1.append(ShrubArea)
            xsB1.append(scts)
            ElevB1.append(aE)
            QowB1.append(QowTS)
            QowcB1.append(Qow_cumul)
            VolB1.append(volumeTS)
            xbB1.append(bbscts)
        
        # RSLR2
        elif N == 2:
            HdA2.append(aHd)
            SCRA2.append(SCrate)
            bbSCRA2.append(bbSCrate)
            WidthA2.append(aW)
            SAA2.append(0)
            xsA2.append(scts)
            ElevA2.append(aE)
            QowA2.append(QowTS)
            QowcA2.append(Qow_cumul)
            VolA2.append(volumeTS)
            xbA2.append(bbscts)
        
        elif N == 7: 
            HdB2.append(aHd)
            SCRB2.append(SCrate)
            bbSCRB2.append(bbSCrate)
            WidthB2.append(aW)
            SAB2.append(ShrubArea)
            xsB2.append(scts)
            ElevB2.append(aE)
            QowB2.append(QowTS)
            QowcB2.append(Qow_cumul)
            VolB2.append(volumeTS)
            xbB2.append(bbscts)
        
        # RSLR3
        elif N == 3:
            HdA3.append(aHd)
            SCRA3.append(SCrate)
            bbSCRA3.append(bbSCrate)
            WidthA3.append(aW)
            SAA3.append(0)
            xsA3.append(scts)
            ElevA3.append(aE)
            QowA3.append(QowTS)
            QowcA3.append(Qow_cumul)
            VolA3.append(volumeTS)
            xbA3.append(bbscts)
        
        elif N == 8: 
            HdB3.append(aHd)
            SCRB3.append(SCrate)
            bbSCRB3.append(bbSCrate)
            WidthB3.append(aW)
            SAB3.append(ShrubArea)
            xsB3.append(scts)
            ElevB3.append(aE)
            QowB3.append(QowTS)
            QowcB3.append(Qow_cumul)
            VolB3.append(volumeTS)
            xbB3.append(bbscts)
        
        # RSLR 4
        elif N == 4:
            HdA4.append(aHd)
            SCRA4.append(SCrate)
            bbSCRA4.append(bbSCrate)
            WidthA4.append(aW)
            SAA4.append(0)
            xsA4.append(scts)
            ElevA4.append(aE)
            QowA4.append(QowTS)
            QowcA4.append(Qow_cumul)
            VolA4.append(volumeTS)
            xbA4.append(bbscts)
        
        elif N == 9: 
            HdB4.append(aHd)
            SCRB4.append(SCrate)
            bbSCRB4.append(bbSCrate)
            WidthB4.append(aW)
            SAB4.append(ShrubArea)
            xsB4.append(scts)
            ElevB4.append(aE)
            QowB4.append(QowTS)
            QowcB4.append(Qow_cumul)
            VolB4.append(volumeTS)
            xbB4.append(bbscts)


#%% CALCULATE MEANS

# RSLR0
HdA0_Avg = np.average(HdA0, axis=0)
SCRA0_Avg = np.average(SCRA0, axis=0)
bbSCRA0_Avg = np.average(bbSCRA0, axis=0)
WidthA0_Avg = np.average(WidthA0, axis=0)
xsA0_Avg = np.average(xsA0, axis=0)
ElevA0_Avg = np.average(ElevA0, axis=0)
QowA0_Avg = np.average(QowA0, axis=0)
QowcA0_Avg = np.average(QowcA0, axis=0)
VolA0_Avg = np.average(VolA0, axis=0)
xbA0_Avg = np.average(xbA0, axis=0)

HdB0_Avg = np.average(HdB0, axis=0)
SCRB0_Avg = np.average(SCRB0, axis=0)
bbSCRB0_Avg = np.average(bbSCRB0, axis=0)
WidthB0_Avg = np.average(WidthB0, axis=0)
xsB0_Avg = np.average(xsB0, axis=0)
SAB0_Avg = np.average(SAB0, axis=0)
ElevB0_Avg = np.average(ElevB0, axis=0)
QowB0_Avg = np.average(QowB0, axis=0)
QowcB0_Avg = np.average(QowcB0, axis=0)
VolB0_Avg = np.average(VolB0, axis=0)
xbB0_Avg = np.average(xbB0, axis=0)
    
# RSLR1
HdA1_Avg = np.average(HdA1, axis=0)
SCRA1_Avg = np.average(SCRA1, axis=0)
bbSCRA1_Avg = np.average(bbSCRA1, axis=0)
WidthA1_Avg = np.average(WidthA1, axis=0)
xsA1_Avg = np.average(xsA1, axis=0)
ElevA1_Avg = np.average(ElevA1, axis=0)
QowA1_Avg = np.average(QowA1, axis=0)
QowcA1_Avg = np.average(QowcA1, axis=0)
VolA1_Avg = np.average(VolA1, axis=0)
xbA1_Avg = np.average(xbA1, axis=0)

HdB1_Avg = np.average(HdB1, axis=0)
SCRB1_Avg = np.average(SCRB1, axis=0)
bbSCRB1_Avg = np.average(bbSCRB1, axis=0)
WidthB1_Avg = np.average(WidthB1, axis=0)
xsB1_Avg = np.average(xsB1, axis=0)
SAB1_Avg = np.average(SAB1, axis=0)
ElevB1_Avg = np.average(ElevB1, axis=0)
QowB1_Avg = np.average(QowB1, axis=0)
QowcB1_Avg = np.average(QowcB1, axis=0)
VolB1_Avg = np.average(VolB1, axis=0)
xbB1_Avg = np.average(xbB1, axis=0)

# RSLR2
HdA2_Avg = np.average(HdA2, axis=0)
SCRA2_Avg = np.average(SCRA2, axis=0)
bbSCRA2_Avg = np.average(bbSCRA2, axis=0)
WidthA2_Avg = np.average(WidthA2, axis=0)
xsA2_Avg = np.average(xsA2, axis=0)
ElevA2_Avg = np.average(ElevA2, axis=0)
QowA2_Avg = np.average(QowA2, axis=0)
QowcA2_Avg = np.average(QowcA2, axis=0)
VolA2_Avg = np.average(VolA2, axis=0)
xbA2_Avg = np.average(xbA2, axis=0)

HdB2_Avg = np.average(HdB2, axis=0)
SCRB2_Avg = np.average(SCRB2, axis=0)
bbSCRB2_Avg = np.average(bbSCRB2, axis=0)
WidthB2_Avg = np.average(WidthB2, axis=0)
xsB2_Avg = np.average(xsB2, axis=0)
SAB2_Avg = np.average(SAB2, axis=0)
ElevB2_Avg = np.average(ElevB2, axis=0)
QowB2_Avg = np.average(QowB2, axis=0)
QowcB2_Avg = np.average(QowcB2, axis=0)
VolB2_Avg = np.average(VolB2, axis=0)
xbB2_Avg = np.average(xbB2, axis=0)

# RSLR3
HdA3_Avg = np.average(HdA3, axis=0)
SCRA3_Avg = np.average(SCRA3, axis=0)
bbSCRA3_Avg = np.average(bbSCRA3, axis=0)
WidthA3_Avg = np.average(WidthA3, axis=0)
xsA3_Avg = np.average(xsA3, axis=0)
ElevA3_Avg = np.average(ElevA3, axis=0)
QowA3_Avg = np.average(QowA3, axis=0)
QowcA3_Avg = np.average(QowcA3, axis=0)
VolA3_Avg = np.average(VolA3, axis=0)
xbA3_Avg = np.average(xbA3, axis=0)

HdB3_Avg = np.average(HdB3, axis=0)
SCRB3_Avg = np.average(SCRB3, axis=0)
bbSCRB3_Avg = np.average(bbSCRB3, axis=0)
WidthB3_Avg = np.average(WidthB3, axis=0)
xsB3_Avg = np.average(xsB3, axis=0)
SAB3_Avg = np.average(SAB3, axis=0)
ElevB3_Avg = np.average(ElevB3, axis=0)
QowB3_Avg = np.average(QowB3, axis=0)
QowcB3_Avg = np.average(QowcB3, axis=0)
VolB3_Avg = np.average(VolB3, axis=0)
xbB3_Avg = np.average(xbB3, axis=0)

# RSLR4
HdA4_Avg = np.average(HdA4, axis=0)
SCRA4_Avg = np.average(SCRA4, axis=0)
bbSCRA4_Avg = np.average(bbSCRA4, axis=0)
WidthA4_Avg = np.average(WidthA4, axis=0)
xsA4_Avg = np.average(xsA4, axis=0)
ElevA4_Avg = np.average(ElevA4, axis=0)
QowA4_Avg = np.average(QowA4, axis=0)
QowcA4_Avg = np.average(QowcA4, axis=0)
VolA4_Avg = np.average(VolA4, axis=0)
xbA4_Avg = np.average(xbA4, axis=0)

HdB4_Avg = np.average(HdB4, axis=0)
SCRB4_Avg = np.average(SCRB4, axis=0)
bbSCRB4_Avg = np.average(bbSCRB4, axis=0)
WidthB4_Avg = np.average(WidthB4, axis=0)
xsB4_Avg = np.average(xsB4, axis=0)
SAB4_Avg = np.average(SAB4, axis=0)
ElevB4_Avg = np.average(ElevB4, axis=0)
QowB4_Avg = np.average(QowB4, axis=0)
QowcB4_Avg = np.average(QowcB4, axis=0)
VolB4_Avg = np.average(VolB4, axis=0)
xbB4_Avg = np.average(xbB4, axis=0)



#%% PLOT       
    
if Plot:
    
    c1 = 'darkviolet'
    c2 = 'blue'
    c3 = 'green'
    c4 = 'gold'
    c5 = 'red'
    
    
    Fig = plt.figure(figsize=(20,20))
    plt.rcParams.update({'font.size':10})

    ax = Fig.add_subplot(4,2,1)
    plt.plot(xsA0_Avg, c1, ls='-')
    plt.plot(xsB0_Avg, c1, ls='--')
    plt.plot(xsA1_Avg, c2, ls='-')
    plt.plot(xsB1_Avg, c2, ls='--')
    plt.plot(xsA2_Avg, c3, ls='-')
    plt.plot(xsB2_Avg, c3, ls='--')
    plt.plot(xsA3_Avg, c4, ls='-')
    plt.plot(xsB3_Avg, c4, ls='--')
    plt.plot(xsA4_Avg, c5, ls='-')
    plt.plot(xsB4_Avg, c5, ls='--')
    plt.ylabel('Shoreline Position (m)')
    plt.xlabel('Year')
    if sealevel:
        plt.legend(['RSLR = 3', '', 'RSLR = 6', '', 'RSLR = 9', '', 'RSLR = 12', '', 'RSLR = 15', ''])
    else:
        plt.legend(['r = 0.3', '', 'r = 0.45', '', 'r = 0.6', '', 'r = 0.75', '', 'r = 0.9', ''])
    
    ax = Fig.add_subplot(4,2,4)
    plt.plot(WidthA0_Avg, c1, ls='-')
    plt.plot(WidthB0_Avg, c1, ls='--')
    plt.plot(WidthA1_Avg, c2, ls='-')
    plt.plot(WidthB1_Avg, c2, ls='--')
    plt.plot(WidthA2_Avg, c3, ls='-')
    plt.plot(WidthB2_Avg, c3, ls='--')
    plt.plot(WidthA3_Avg, c4, ls='-')
    plt.plot(WidthB3_Avg, c4, ls='--')
    plt.plot(WidthA4_Avg, c5, ls='-')
    plt.plot(WidthB4_Avg, c5, ls='--')
    plt.ylabel('Island Width (m)')
    plt.xlabel('Year') 
    
    ax = Fig.add_subplot(4,2,3)
    plt.plot(QowcA0_Avg, c1, ls='-')
    plt.plot(QowcB0_Avg, c1, ls='--')
    plt.plot(QowcA1_Avg, c2, ls='-')
    plt.plot(QowcB1_Avg, c2, ls='--')
    plt.plot(QowcA2_Avg, c3, ls='-')
    plt.plot(QowcB2_Avg, c3, ls='--')
    plt.plot(QowcA3_Avg, c4, ls='-')
    plt.plot(QowcB3_Avg, c4, ls='--')
    plt.plot(QowcA4_Avg, c5, ls='-')
    plt.plot(QowcB4_Avg, c5, ls='--')
    plt.ylabel('Cumulative Overwash Flux (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(4,2,5)
    plt.plot(ElevA0_Avg, c1, ls='-')
    plt.plot(ElevB0_Avg, c1, ls='--')
    plt.plot(ElevA1_Avg, c2, ls='-')
    plt.plot(ElevB1_Avg, c2, ls='--')
    plt.plot(ElevA2_Avg, c3, ls='-')
    plt.plot(ElevB2_Avg, c3, ls='--')
    plt.plot(ElevA3_Avg, c4, ls='-')
    plt.plot(ElevB3_Avg, c4, ls='--')
    plt.plot(ElevA4_Avg, c5, ls='-')
    plt.plot(ElevB4_Avg, c5, ls='--')
    plt.ylabel('Interior Elevation (m)')
    plt.xlabel('Year')  
    
    ax = Fig.add_subplot(4,2,6)
    plt.plot(VolA0_Avg, c1, ls='-')
    plt.plot(VolB0_Avg, c1, ls='--')
    plt.plot(VolA1_Avg, c2, ls='-')
    plt.plot(VolB1_Avg, c2, ls='--')
    plt.plot(VolA2_Avg, c3, ls='-')
    plt.plot(VolB2_Avg, c3, ls='--')
    plt.plot(VolA3_Avg, c4, ls='-')
    plt.plot(VolB3_Avg, c4, ls='--')
    plt.plot(VolA4_Avg, c5, ls='-')
    plt.plot(VolB4_Avg, c5, ls='--')
    plt.ylabel('Island Volume (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(4,2,2)
    plt.plot(xbA0_Avg, c1, ls='-')
    plt.plot(xbB0_Avg, c1, ls='--')
    plt.plot(xbA1_Avg, c2, ls='-')
    plt.plot(xbB1_Avg, c2, ls='--')
    plt.plot(xbA2_Avg, c3, ls='-')
    plt.plot(xbB2_Avg, c3, ls='--')
    plt.plot(xbA3_Avg, c4, ls='-')
    plt.plot(xbB3_Avg, c4, ls='--')
    plt.plot(xbA4_Avg, c5, ls='-')
    plt.plot(xbB4_Avg, c5, ls='--')
    plt.ylabel('Back-Barrier Shoreline Position (m)')
    plt.xlabel('Year') 
    
    ax = Fig.add_subplot(4,2,7)
    plt.plot(HdA0_Avg, c1, ls='-')
    plt.plot(HdB0_Avg, c1, ls='--')
    plt.plot(HdA1_Avg, c2, ls='-')
    plt.plot(HdB1_Avg, c2, ls='--')
    plt.plot(HdA2_Avg, c3, ls='-')
    plt.plot(HdB2_Avg, c3, ls='--')
    plt.plot(HdA3_Avg, c4, ls='-')
    plt.plot(HdB3_Avg, c4, ls='--')
    plt.plot(HdA4_Avg, c5, ls='-')
    plt.plot(HdB4_Avg, c5, ls='--')
    plt.ylabel('Average Dune Height (m)')
    plt.xlabel('Year') 
    
    ax = Fig.add_subplot(4,2,8)
    plt.plot(SAB0_Avg, c1)
    plt.plot(SAB1_Avg, c2)
    plt.plot(SAB2_Avg, c3)
    plt.plot(SAB3_Avg, c4)
    plt.plot(SAB4_Avg, c5)
    plt.ylabel('ShrubArea (dam^2)')
    plt.xlabel('Year')
    
    plt.show()
    

    # For publication
    
    Fig = plt.figure(figsize=(16,15))
    plt.rcParams.update({'font.size':14})

    ax = Fig.add_subplot(3,2,3)
    plt.plot(xsA0_Avg, c1, ls='-')
    plt.plot(xsB0_Avg, c1, ls='--')
    plt.plot(xsA1_Avg, c2, ls='-')
    plt.plot(xsB1_Avg, c2, ls='--')
    plt.plot(xsA2_Avg, c3, ls='-')
    plt.plot(xsB2_Avg, c3, ls='--')
    plt.plot(xsA3_Avg, c4, ls='-')
    plt.plot(xsB3_Avg, c4, ls='--')
    plt.plot(xsA4_Avg, c5, ls='-')
    plt.plot(xsB4_Avg, c5, ls='--')
    plt.ylabel('Shoreline Position (m)')
    plt.xlabel('Year')
    # if sealevel:
    #     plt.legend(['RSLR = 3', '', 'RSLR = 6', '', 'RSLR = 9', '', 'RSLR = 12', '', 'RSLR = 15', ''])
    # else:
    #     plt.legend(['r = 0.3', '', 'r = 0.45', '', 'r = 0.6', '', 'r = 0.75', '', 'r = 0.9', ''])
    
    ax = Fig.add_subplot(3,2,1)
    plt.plot(WidthA0_Avg, c1, ls='-')
    plt.plot(WidthB0_Avg, c1, ls='--')
    plt.plot(WidthA1_Avg, c2, ls='-')
    plt.plot(WidthB1_Avg, c2, ls='--')
    plt.plot(WidthA2_Avg, c3, ls='-')
    plt.plot(WidthB2_Avg, c3, ls='--')
    plt.plot(WidthA3_Avg, c4, ls='-')
    plt.plot(WidthB3_Avg, c4, ls='--')
    plt.plot(WidthA4_Avg, c5, ls='-')
    plt.plot(WidthB4_Avg, c5, ls='--')
    plt.ylabel('Barrier Width (m)')
    plt.xlabel('Year') 
    
    ax = Fig.add_subplot(3,2,4)
    plt.plot(QowcA0_Avg, c1, ls='-')
    plt.plot(QowcB0_Avg, c1, ls='--')
    plt.plot(QowcA1_Avg, c2, ls='-')
    plt.plot(QowcB1_Avg, c2, ls='--')
    plt.plot(QowcA2_Avg, c3, ls='-')
    plt.plot(QowcB2_Avg, c3, ls='--')
    plt.plot(QowcA3_Avg, c4, ls='-')
    plt.plot(QowcB3_Avg, c4, ls='--')
    plt.plot(QowcA4_Avg, c5, ls='-')
    plt.plot(QowcB4_Avg, c5, ls='--')
    plt.ylabel('Cumulative Overwash Flux (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(3,2,2)
    plt.plot(VolA0_Avg, c1, ls='-')
    plt.plot(VolB0_Avg, c1, ls='--')
    plt.plot(VolA1_Avg, c2, ls='-')
    plt.plot(VolB1_Avg, c2, ls='--')
    plt.plot(VolA2_Avg, c3, ls='-')
    plt.plot(VolB2_Avg, c3, ls='--')
    plt.plot(VolA3_Avg, c4, ls='-')
    plt.plot(VolB3_Avg, c4, ls='--')
    plt.plot(VolA4_Avg, c5, ls='-')
    plt.plot(VolB4_Avg, c5, ls='--')
    plt.ylabel('Barrier Volume (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(3,2,5)
    plt.plot(SAB0_Avg, c1)
    plt.plot(SAB1_Avg, c2)
    plt.plot(SAB2_Avg, c3)
    plt.plot(SAB3_Avg, c4)
    plt.plot(SAB4_Avg, c5)
    plt.ylabel('Shrub Cover (Count))')
    plt.xlabel('Year')
    
    plt.tight_layout(pad=1.5)
    plt.show()
    
    
    
    
    # For publication - Without No Shrub lines
    
    Fig = plt.figure(figsize=(16,15))
    plt.rcParams.update({'font.size':14})

    ax = Fig.add_subplot(3,2,3)
    plt.plot(xsB0_Avg, c1, ls='--')
    plt.plot(xsB1_Avg, c2, ls='--')
    plt.plot(xsB2_Avg, c3, ls='--')
    plt.plot(xsB3_Avg, c4, ls='--')
    plt.plot(xsB4_Avg, c5, ls='--')
    plt.ylabel('Shoreline Position (m)')
    plt.xlabel('Year')
    # if sealevel:
    #     plt.legend(['RSLR = 3', '', 'RSLR = 6', '', 'RSLR = 9', '', 'RSLR = 12', '', 'RSLR = 15', ''])
    # else:
    #     plt.legend(['r = 0.3', '', 'r = 0.45', '', 'r = 0.6', '', 'r = 0.75', '', 'r = 0.9', ''])
    
    ax = Fig.add_subplot(3,2,1)
    plt.plot(WidthB0_Avg, c1, ls='--')
    plt.plot(WidthB1_Avg, c2, ls='--')
    plt.plot(WidthB2_Avg, c3, ls='--')
    plt.plot(WidthB3_Avg, c4, ls='--')
    plt.plot(WidthB4_Avg, c5, ls='--')
    plt.ylabel('Barrier Width (m)')
    plt.xlabel('Year') 
    
    ax = Fig.add_subplot(3,2,4)
    plt.plot(QowcB0_Avg, c1, ls='--')
    plt.plot(QowcB1_Avg, c2, ls='--')
    plt.plot(QowcB2_Avg, c3, ls='--')
    plt.plot(QowcB3_Avg, c4, ls='--')
    plt.plot(QowcB4_Avg, c5, ls='--')
    plt.ylabel('Cumulative Overwash Flux (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(3,2,2)
    plt.plot(VolB0_Avg, c1, ls='--')
    plt.plot(VolB1_Avg, c2, ls='--')
    plt.plot(VolB2_Avg, c3, ls='--')
    plt.plot(VolB3_Avg, c4, ls='--')
    plt.plot(VolB4_Avg, c5, ls='--')
    plt.ylabel('Barrier Volume (m^3/m)')
    plt.xlabel('Year')   
    
    ax = Fig.add_subplot(3,2,5)
    plt.plot(SAB0_Avg, c1)
    plt.plot(SAB1_Avg, c2)
    plt.plot(SAB2_Avg, c3)
    plt.plot(SAB3_Avg, c4)
    plt.plot(SAB4_Avg, c5)
    plt.ylabel('Shrub Cover (Count))')
    plt.xlabel('Year')
    
    plt.tight_layout(pad=1.5)
    plt.show()
    
    
    
    
    