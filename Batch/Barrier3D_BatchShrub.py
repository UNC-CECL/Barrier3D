# Main runfile for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- Copyright (C) 2020 Ian R.B. Reeves (current developer)                                                                                                                                                                                          -
- Current developer can be contacted by email (reevesi@live.unc.edu) and paper mail (104 South Road, Mitchell Hall CB #3315, Chapel Hill, NC 27599, USA                                                                                           -
- This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. -
- This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.        -
- You should have received a copy of the GNU General Public License along with this program; if not, see <https://www.gnu.org/licenses/>.                                                                                                         -
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""

# Version Number: 1
# Updated: 30 April 2021


# Barrier3D_Parameters.py and Barrier3D_Functions.py must be located in same directory as this main runfile

# All units are in DECAMETERS for simulation; convert to meters after simulation end for presentation

# Processes Included: 
#   RSLR
#   Overwash flow routing with different sediment transport rules depending on regime
#   Ocean and back-barrier shoreline change
#   Extension of back-barrier shoreline where overwash occurs
#   LTA14 shoreface dynamics control shoreline change
#   Active dune field changes location depending on shoreline change
#   Exponential decay of sediment drop-out in back-barrier bay
#   Uses initial morphology from lidar
#   Dune field has width of multiple cells
#   Has option for using a prefab storm time series
#   Includes shurb expansion/mortality module  
   

def Barrier3D_BatchShrub(SimNumber, directory, StormSeries, growthparam, Shrub_ON, RSLR, DuneDomain, Dstart):

    #==============================================================================================================================================
    # IMPORTS       
    #==============================================================================================================================================
    
    import Barrier3D_Functions_Batch as func
    import numpy as np
    import math 
    import time
    import warnings
    warnings.filterwarnings("ignore")
    
    
    Time = time.time()
     
    #==============================================================================================================================================
    # SET PARAMETERS       
    #==============================================================================================================================================
    
    ### Load Input Parameters
    from Barrier3D_Parameters_Batch import ( 
                                             BarrierLength, 
                                             BayDepth, 
                                             BermEl, 
                                             DomainWidth,
                                             DuneWidth, 
                                             DShoreface, 
                                             InteriorDomain, 
                                             LShoreface, 
                                             MHW, 
                                             SD_storm, 
                                             StormStart, 
                                             StormTimeSeries,
                                             TMAX, 
                                             mean_storm, 
                                             numstorm, 
                                             Qshrub_max, 
                                             C1, 
                                             C2, 
                                             DuneRestart, 
                                             nn, 
                                             mm,
                                             MaxUpSlope,
                                             threshold_in, 
                                             Rin_r, 
                                             Rin_i, 
                                             Qs_min, 
                                             Kr, 
                                             Ki, 
                                             Cbb_r, 
                                             Cbb_i, 
                                             Qs_bb_min, 
                                             Cx, 
                                             OWss_i, 
                                             OWss_r,
                                             SimParams
                                             )
    
    
    ### Set variables
    SL = 0 # Does not change (Lagrangian frame of reference)
    MaxAvgSlope = BermEl / 10
    fluxLimit = 100 # Initialize
    
    SimParams[1] = RSLR
    SimParams[8] = Shrub_ON
    
    
    ### Initialize Shrubs
    if Shrub_ON ==1: print('Shrubs ON')
    PercentCoverTS = [None] * TMAX
    PercentCoverTS[0] = np.zeros([DomainWidth, BarrierLength])
    DeadPercentCoverTS = [None] * TMAX
    DeadPercentCoverTS[0] = np.zeros([DomainWidth, BarrierLength])
    ShrubFemaleTS = [None] * TMAX
    ShrubFemaleTS[0] = np.zeros([DomainWidth, BarrierLength])
    ShrubMaleTS = [None] * TMAX
    ShrubMaleTS[0] = np.zeros([DomainWidth, BarrierLength])
    ShrubDeadTS = [None] * TMAX
    ShrubDeadTS[0] = np.zeros([DomainWidth, BarrierLength])
    ShrubDomainFemale = np.zeros([DomainWidth, BarrierLength])
    ShrubDomainMale = np.zeros([DomainWidth, BarrierLength])
    ShrubDomainDead = np.zeros([DomainWidth, BarrierLength])
    ShrubPercentCover = PercentCoverTS[0]
    DeadPercentCover = DeadPercentCoverTS[0]
    BurialDomain = np.zeros([DomainWidth, BarrierLength])
    ShrubArea = [0]
    
    
    ### Initialize variables 
    DomainTS = [None] * TMAX # Stores the elevation domain for each timestep
    DomainTS[0] = InteriorDomain
    SCR = 0 # (dam/yr) Change rate of ocean-fronting shoreline: (+) = progradation, (-) = erosion
    SCRagg = 0 # Counter variable for shoreline change that is less than 1 cell size per year
    ShorelineChangeTS = [0] 
    ShorelineChange = 0 # (dam) Total amount of shoreline change 
    StormCount = [0]
    InundationCount = 0
    RunUpCount = 0
    InteriorWidth_AvgTS = [DomainWidth]
    QowTS = [0] # (m^3/m)
    x_t = 0  # (dam) Start location of shoreface toe
    x_s = x_t + LShoreface # (dam) Start location of shoreline
    x_t_TS = [x_t] # (dam) Shoreface toe locations for each time step
    x_s_TS = [x_s] # (dam) Shoreline locations for each time step
    x_b_TS = [(x_s + DomainWidth)] # (dam) Bay shoreline locations for each time step
    s_sf_TS = [(DShoreface/LShoreface)]
    Hd_AverageTS = [Dstart]
    QsfTS = [0] # (m^3/m)
    Hd_Loss_TS = np.zeros([TMAX,BarrierLength]) # Dune height loss exclusively from vertical storm erosion
    
    # Batch version additions
    IslandAreaTS = []
    AvgInteriorElevationTS = []
    Hd_Loss_TS = np.zeros([TMAX,BarrierLength]) # Dune height loss exclusively from vertical storm erosion
    Result = 0
    
    
    
    #==============================================================================================================================================
    # RUN MODEL      
    #==============================================================================================================================================
    
    for t in range(1, TMAX): # Yearly time steps - actual time = t + 1
        
        ### Print time step to screen
        print("\r", 'Time Step: ', t, end = "")
        
        
        ### RSLR
        InteriorDomain, DuneDomain = func.SeaLevelBatch(InteriorDomain, DuneDomain, t, RSLR)
        
        
        ### Find DomainWidth and InteriorWidth
        DomainWidth, InteriorWidth, InteriorWidth_Avg = func.FindWidths(InteriorDomain, SL)
        
        
        ### Check for drowning
        # Barrier drowns if interior width thins to 10 m or fewer or all interior cells are below SL
        if max(InteriorWidth) <= 1:
            Result = 1
            TMAX = t-1
            break
        if all(j <= SL for j in InteriorDomain[0,:]):
            Result = 2
            TMAX = t-1
            break
    
    
        ### Grow Dunes
        DuneDomain, Dmax, Qdg = func.DuneGrowthBatch(DuneDomain, t, growthparam)
    
        DuneDomainCrest = DuneDomain[t,:,:].max(axis=1) # Maximum height of each row in DuneDomain    
        DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
            
        Hd_AverageTS.append(np.mean(DuneDomainCrest)) # Store average pre-storms dune-heigh for time step
        
        
        
        ###########################################      
        ### Shrubs 
        
        if Shrub_ON == 1:
    
            ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, BurialDomain, ShrubDomainDead = func.Shrubs(InteriorDomain, DuneDomainCrest, t, ShrubDomainFemale, ShrubDomainMale, ShrubPercentCover, 
                                                                                                        DeadPercentCover, BurialDomain, InteriorWidth_Avg, DomainWidth, ShrubDomainDead, SL)
        
        
        
        ###########################################      
        ### Storms
    
        OWloss = 0
        DuneLoss = 0    
        
        if t >= StormStart:
            # Select number of storms for this time step from normal distribution
            if StormTimeSeries == True:
                TSloc = np.argwhere(StormSeries[:,0] == t)
                numstorm = int(len(TSloc)) # analysis:ignore
            else:        
                numstorm = int(numstorm)
                numstorm = round(np.random.normal(mean_storm, SD_storm)) # analysis:ignore # Comment out this line if using pre-specified static number of storms per year        
            
            if numstorm > 0:
                if StormTimeSeries == True:
                    TSloc = np.argwhere(StormSeries[:,0] == t)
                    start = TSloc[0,0]
                    stop = TSloc[-1,0] + 1
                    Rhigh = StormSeries[start:stop,1] 
                    Rlow = StormSeries[start:stop,2] 
                    period = StormSeries[start:stop,3]
                    dur = np.array(StormSeries[start:stop,4], dtype='int')
                    
                else:                    
                    # Draw statistics for each storm           
                    Rhigh, Rlow, period, dur = func.WaterLevels(numstorm, MHW)
                
                ### Individual Storm Impacts       
                for n in range(numstorm): # Loop through each individual storm
                                    
                    ###########################################
                    ### Dune Erosion      
                    
                    DuneChange = np.zeros([BarrierLength, DuneWidth]) # Vector storing dune height change for this storm
                    
                    # Find overwashed dunes and gaps           
                    Dow = [index for index, value in enumerate((DuneDomainCrest + BermEl)) if value < Rhigh[n]]
                    gaps = func.DuneGaps(DuneDomainCrest, Dow, BermEl, Rhigh[n]) # Finds location and Rexcess of continuous gaps in dune ridge                
                    
                    for d in range(len(Dow)): # Loop through each overwashed dune cell
                        for w in range(DuneWidth):
    
                            # Calculate dune elevation loss
                            Rnorm = Rhigh[n] / (DuneDomain[t, Dow[d], w] + BermEl) # Rhigh relative to pre-storm dune elevation
                            Dloss = Rnorm / (C1 + (Rnorm * (Rnorm - C2))) # Amount of dune crest elevation change normalized by pre-storm dune elevation (i.e. a percent change), from Goldstein and Moore (2016)
                            
                            # Set new dune height
                            InitDElev = DuneDomain[t, Dow[d], w] + BermEl
                            NewDElev = InitDElev * (1-Dloss) # Calculate new dune elevation from storm lowering
                            if NewDElev < BermEl:
                                NewDElev = BermEl
                            DuneDomain[t, Dow[d], w] = NewDElev - BermEl # Convert elevation to height above berm
                           
                            DuneChange[Dow[d], w] = InitDElev - NewDElev 
        
                            # If dune is lowered to essentially zero, allow for chance of regrowth by raising dune height to 5 cm                                    
                            if DuneDomain[t, Dow[d], w] < DuneRestart:
                                if DuneRestart < Dmax:
                                    DuneDomain[t, Dow[d], w] = DuneRestart
                                else:
                                    DuneDomain[t, Dow[d], w] = Dmax # Restart height can't be greater than Dmax
    
                    # Dune Height Diffusion
                    DuneDomain = func.DiffuseDunes(DuneDomain, t)
                    DuneDomain[DuneDomain < SL] = DuneRestart
                    
                    DuneLoss = np.sum(DuneChange)/BarrierLength
                    Hd_TSloss = DuneChange.max(axis=1) / dur[n] # Average height of dune loss for each substep during storm
    
                    Hd_Loss_TS[t,:] = Hd_Loss_TS[t,:] + DuneChange.max(axis=1)
    
                    ###########################################                
                    ### Overwash
           
                    Iow = 0 # Count of dune gaps in inundation regime
                    Dunes_prestorm = DuneDomainCrest               
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2] # (m)                    
                        gapwidth = stop - start + 1
                        meandune = (sum(Dunes_prestorm[start:stop+1]) / gapwidth) + BermEl # Average elevation of dune gap
                        
                        # Determine number of gaps in inundation regime
                        if Rlow[n] > meandune:
                            Iow += 1
                            
                    # Dertermine Sediment And Water Routing Rules
                            
                    if len(gaps) > 0 and Iow / len(gaps) >= threshold_in: # If greater than threshold % of dune gaps are in unundation regime, use inundation regime routing
                        inundation = 1                                           
                        substep = OWss_i                                                                 
                        InundationCount += 1
                    else:
                        inundation = 0                                                                 
                        substep = OWss_r                                                                        
                        RunUpCount += 1
                                    
                    # Set Domain
                    add = 10
                    duration = int(dur[n] * substep) 
                    width = np.shape(InteriorDomain)[0] + 1 + add # (dam) Add one for Dunes and 25 for bay
                    Elevation = np.zeros([duration, width, BarrierLength])                
                    Dunes = Dunes_prestorm + BermEl
                    Bay = np.ones([add, BarrierLength]) * -BayDepth
                    Elevation[0,:,:] = np.vstack([Dunes, InteriorDomain, Bay])         
    
                    # Initialize Memory Storage Arrays
                    Discharge = np.zeros([duration, width, BarrierLength])
                    SedFluxIn = np.zeros([duration, width, BarrierLength])
                    SedFluxOut = np.zeros([duration, width, BarrierLength])
    
                    Rin = 0 # (dam^3/t) Infiltration Rate, volume of overwash flow lost per m cross-shore per time 
                                    
                    # Set Water at Dune Crest
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2] # (m)
                        gapwidth = stop - start + 1
                        meandune = (sum(Dunes_prestorm[start:stop+1]) / gapwidth) + BermEl # Average elevation of dune gap
                        
                        # Calculate discharge through each dune cell
                        Vdune = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10 # (dam/s)                      
                        Qdune =  Vdune * Rexcess * 3600 # (dam^3/hr)
    
                        # Set discharge at dune gap
                        Discharge[:,0,start:stop] = Qdune
                        
                        if inundation == 1: # Inundation regime       
                            Rin = Rin_i
                            
                            # Find average slope of interior
                            AvgSlope = BermEl / 20 #InteriorWidth_Avg #<-------------------TEMP!!! Hold C constant
                            
                            # # Enforce max average interior slope
                            # AvgSlope = min(MaxAvgSlope, AvgSlope)
                            
                            C = Cx * AvgSlope # Momentum constant
    
                        else: # Run-up regime
                            Rin = Rin_r 
                    
                    ### Run Flow Routing Algorithm                                      
                    for TS in range(duration):
                        
                        ShrubDomainWidth = np.shape(ShrubDomainFemale)[0]
                        DeadDomainWidth = np.shape(ShrubDomainDead)[0]
                                 
                        if TS > 0:                      
                            Elevation[TS,1:,:] = Elevation[TS-1,1:,:] # Begin timestep with elevation from end of last                                    
                            Elevation[TS,0,:] = Dunes - ((Hd_TSloss/substep) * TS) # Reduce dune in height linearly over course of storm
                            
                        for d in range(width-1):                                              
                            # Reduce discharge across row via infiltration  
                            if d > 0:  
                                Discharge[TS,d,:][Discharge[TS,d,:] > 0] -= Rin                                         
                            Discharge[TS,d,:][Discharge[TS,d,:] < 0] = 0
                                                    
                            for i in range(BarrierLength):                                       
                                if Discharge[TS,d,i] > 0:
                    
                                    Q0 = Discharge[TS,d,i]
                                                 
                                    ### Calculate Slopes                                 
                                    if i > 0:
                                        S1 = (Elevation[TS,d,i] - Elevation[TS,d+1,i-1]) / (math.sqrt(2))
                                        S1 = np.nan_to_num(S1)
                                    else:
                                        S1 = 0
                                        
                                    S2 = Elevation[TS,d,i] - Elevation[TS,d+1,i]
                                    S2 = np.nan_to_num(S2)
                                    
                                    if i < (BarrierLength-1):
                                        S3 = (Elevation[TS,d,i] - Elevation[TS,d+1,i+1]) / (math.sqrt(2))
                                        S3 = np.nan_to_num(S3)
                                    else:
                                        S3 = 0                               
                                    
                                    ### Calculate Discharge To Downflow Neighbors
                                    # One or more slopes positive
                                    if S1 > 0 or S2 > 0 or S3 > 0:
                                        
                                        if S1 < 0:
                                            S1 = 0
                                        if S2 < 0:
                                            S2 = 0    
                                        if S3 < 0:
                                            S3 = 0 
                                        
                                        Q1 = Q0 * S1**nn / (S1**nn + S2**nn + S3**nn)
                                        Q2 = Q0 * S2**nn / (S1**nn + S2**nn + S3**nn)
                                        Q3 = Q0 * S3**nn / (S1**nn + S2**nn + S3**nn)
                                        
                                        Q1 = np.nan_to_num(Q1)
                                        Q2 = np.nan_to_num(Q2)
                                        Q3 = np.nan_to_num(Q3)                                                                        
                                    
                                    # No slopes positive, one or more equal to zero
                                    elif S1 == 0 or S2 == 0 or S3 == 0:
                                        
                                        pos = 0
                                        if S1 == 0:
                                            pos += 1
                                        if S2 == 0:
                                            pos += 1    
                                        if S3 == 0:
                                            pos += 1 
                                    
                                        Qx = Q0 / pos
                                        Qx = np.nan_to_num(Qx)  
                                                                            
                                        if S1 == 0 and i > 0:
                                            Q1 = Qx
                                        else:
                                            Q1 = 0
                                        if S2 == 0:
                                            Q2 = Qx
                                        else:
                                            Q2 = 0
                                        if S3 == 0 and i < (BarrierLength-1):
                                            Q3 = Qx
                                        else:
                                            Q3 = 0
                                                                                  
                                    # All slopes negative
                                    else:                            
                                            
                                        Q1 = Q0 * abs(S1)**(-nn) / (abs(S1)**(-nn) + abs(S2)**(-nn) + abs(S3)**(-nn))
                                        Q2 = Q0 * abs(S2)**(-nn) / (abs(S1)**(-nn) + abs(S2)**(-nn) + abs(S3)**(-nn))
                                        Q3 = Q0 * abs(S3)**(-nn) / (abs(S1)**(-nn) + abs(S2)**(-nn) + abs(S3)**(-nn))                                  
                                                                         
                                        Q1 = np.nan_to_num(Q1)
                                        Q2 = np.nan_to_num(Q2)
                                        Q3 = np.nan_to_num(Q3)  
                                    
                                        #MaxUpSlope = 0.25 # dam
    
                                        if S1 > MaxUpSlope:
                                            Q1 = 0
                                        else:
                                            Q1 = Q1 * (1 - (abs(S1) / MaxUpSlope))
                                    
                                        if S2 > MaxUpSlope:
                                            Q2 = 0
                                        else:
                                            Q2 = Q2 * (1 - (abs(S2) / MaxUpSlope))                                
                                    
                                        if S3 > MaxUpSlope:
                                            Q3 = 0
                                        else:
                                            Q3 = Q3 * (1 - (abs(S3) / MaxUpSlope))                                
                                    
                                    
                                    ### Reduce Overwash Through Shrub Cells and Save Discharge                                
                                    if Shrub_ON == 1:             
                                        # Cell 1
                                        if i > 0:
                                            if d < ShrubDomainWidth and ShrubPercentCover[d,i-1] > 0:
                                                Q1 =  Q1 * ((1 - Qshrub_max) * ShrubPercentCover[d,i-1])   
                                            elif d < (DeadDomainWidth) and DeadPercentCover[d,i-1] > 0: 
                                                Q1 =  Q1 * ((1 - Qshrub_max * 0.66) * DeadPercentCover[d,i-1]) # Dead shrubs block 2/3 of what living shrubs block
                                            Discharge[TS,d+1,i-1] = Discharge[TS,d+1,i-1] + Q1
                                            
                                        # Cell 2                                    
                                        if d < ShrubDomainWidth and ShrubPercentCover[d,i] > 0:
                                            Q2 =  Q2 * ((1 - Qshrub_max) * ShrubPercentCover[d,i])
                                        elif d < (DeadDomainWidth) and DeadPercentCover[d,i] > 0:
                                            Q2 =  Q2 * ((1 - Qshrub_max * 0.66) * DeadPercentCover[d,i])
                                        Discharge[TS,d+1,i] = Discharge[TS,d+1,i] + Q2         
                                        
                                        # Cell 3
                                        if i < (BarrierLength-1):
                                            if d < ShrubDomainWidth and ShrubPercentCover[d,i+1] > 0:
                                                Q3 =  Q3 * ((1 - Qshrub_max) * ShrubPercentCover[d,i+1])  
                                            elif d < (DeadDomainWidth) and DeadPercentCover[d,i+1] > 0:
                                                Q3 =  Q3 * ((1 - Qshrub_max * 0.66) * DeadPercentCover[d,i+1])
                                            Discharge[TS,d+1,i+1] = Discharge[TS,d+1,i+1] + Q3                                   
                                    else:
                                        # Cell 1
                                        if i > 0:
                                            Discharge[TS,d+1,i-1] = Discharge[TS,d+1,i-1] + Q1
                                            
                                        # Cell 2    
                                        Discharge[TS,d+1,i] = Discharge[TS,d+1,i] + Q2 
                                        
                                        # Cell 3
                                        if i < (BarrierLength-1):
                                            Discharge[TS,d+1,i+1] = Discharge[TS,d+1,i+1] + Q3  
                                                                   
                                    
                                    ### Calculate Sed Movement  
                                    fluxLimit = Dmax                                                                                      
                                    
                                    # Run-up Regime  
                                    if inundation == 0:
                                        if Q1 > Qs_min and S1 >= 0:
                                            Qs1 = Kr * Q1
                                            if Qs1 > fluxLimit:
                                                Qs1 = fluxLimit
                                        else:
                                            Qs1 = 0
                                            
                                        if Q2 > Qs_min and S2 >= 0:
                                            Qs2 = Kr * Q2
                                            if Qs2 > fluxLimit:
                                                Qs2 = fluxLimit
                                        else:
                                            Qs2 = 0
                                            
                                        if Q3 > Qs_min and S3 >= 0:
                                            Qs3 = Kr * Q3
                                            if Qs3 > fluxLimit:
                                                Qs3 = fluxLimit
                                        else:
                                            Qs3 = 0                           
                                    
                                    # Inundation Regime - Murray and Paola (1994, 1997) Rule 3 with flux limiter
                                    else:   
                                        if Q1 > Qs_min:
                                            Qs1 = Ki * (Q1*(S1 + C))**mm
                                            if Qs1 < 0: 
                                                Qs1 = 0
                                            elif Qs1 > fluxLimit:
                                                Qs1 = fluxLimit
                                        else:
                                            Qs1 = 0
                                            
                                        if Q2 > Qs_min:
                                            Qs2 = Ki * (Q2*(S2 + C))**mm
                                            if Qs2 < 0: 
                                                Qs2 = 0
                                            elif Qs2 > fluxLimit:
                                                Qs2 = fluxLimit
                                        else:
                                            Qs2 = 0
                                            
                                        if Q3 > Qs_min:
                                            Qs3 = Ki * (Q3*(S3 + C))**mm
                                            if Qs3 < 0: 
                                                Qs3 = 0
                                            elif Qs3 > fluxLimit:
                                                Qs3 = fluxLimit
                                        else:
                                            Qs3 = 0                                
                                                                        
                                    Qs1 = np.nan_to_num(Qs1)
                                    Qs2 = np.nan_to_num(Qs2)
                                    Qs3 = np.nan_to_num(Qs3)          
                                            
                                    ### Calculate Net Erosion/Accretion
                                    # if Elevation[TS,d,i] > SL: # If cell is subaerial, elevation change is determined by difference between flux in vs. flux out
                                    if Elevation[TS,d,i] > SL or any(z > SL for z in Elevation[TS,d+1:d+6,i]):
                                        if i > 0:
                                            SedFluxIn[TS,d+1,i-1] += Qs1
                                            
                                        SedFluxIn[TS,d+1,i] += Qs2
                                        
                                        if i < (BarrierLength-1):
                                            SedFluxIn[TS,d+1,i+1] += Qs3
                                        
                                        Qs_out = Qs1 + Qs2 + Qs3
                                        SedFluxOut[TS,d,i] = Qs_out
                                                                      
                                    else: # If cell is subaqeous, exponentially decay deposition of remaining sediment across bay
    
                                        if inundation == 0: 
                                            Cbb = Cbb_r
                                        else: 
                                            Cbb = Cbb_i                            
                                        
                                        Qs0 = SedFluxIn[TS,d,i] * Cbb
                                        
                                        Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                        Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                        Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)
    
                                        Qs1 = np.nan_to_num(Qs1)
                                        Qs2 = np.nan_to_num(Qs2)
                                        Qs3 = np.nan_to_num(Qs3) 
    
                                        if Qs1 < Qs_bb_min:
                                            Qs1 = 0
                                        elif Qs1 > fluxLimit:
                                            Qs1 = fluxLimit
                                        if Qs2 < Qs_bb_min:
                                            Qs2 = 0
                                        elif Qs2 > fluxLimit:
                                            Qs2 = fluxLimit
                                        if Qs3 < Qs_bb_min:
                                            Qs3 = 0
                                        elif Qs3 > fluxLimit:
                                            Qs3 = fluxLimit
    
                                        if i > 0:
                                            SedFluxIn[TS,d+1,i-1] += Qs1
                                        
                                        SedFluxIn[TS,d+1,i] += Qs2
                                        
                                        if i < (BarrierLength-1):
                                            SedFluxIn[TS,d+1,i+1] += Qs3                           
                                        
                                        Qs_out = Qs1 + Qs2 + Qs3
                                        SedFluxOut[TS,d,i] = Qs_out
                                        
                                            
                                    ### Saline Flooding 
                                    ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead = func.SalineFlooding(ShrubDomainWidth, ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead, d, i, Q0)  
                                       
                                    
                        ### Update Elevation After Every Storm Hour                      
                        if inundation == 1:
                            ElevationChange = (SedFluxIn[TS,:,:] - SedFluxOut[TS,:,:]) / substep
                        else:
                            ElevationChange = (SedFluxIn[TS,:,:] - SedFluxOut[TS,:,:]) / substep
                        Elevation[TS,:,:] = Elevation[TS,:,:] + ElevationChange                    
                                            
                        # Calculate and save volume of sediment deposited on/behind the island for every hour (all four methods below should equal the same!)
                        #OWloss = OWloss + np.sum(ElevationChange[1:,:])
                        #OWloss = OWloss + (np.sum(SedFluxIn[TS,1:,:]) - np.sum(SedFluxOut[TS,1:,:])) / substep  
                        #OWloss = OWloss + np.sum(SedFluxIn[TS,1,:]) / substep
                        OWloss = OWloss + np.sum(SedFluxOut[TS,0,:]) / substep
                        
                        # Update amount of burial/erosion for each shrub
                        BurialDomain = func.UpdateBurial(BurialDomain, ElevationChange, ShrubDomainWidth, ShrubDomainAll)
    
                    ### Update Interior Domain After Every Storm
                    InteriorUpdate = Elevation[-1,1:,:]
                    
                    # Remove all rows of bay without any deposition from the domain
                    check = 1
                    while check == 1:
                        if all(x <= -BayDepth for x in InteriorUpdate[-1,:]):
                            InteriorUpdate = np.delete(InteriorUpdate,(-1), axis=0)
                        else:
                            check = 0        
                  
                    # Update interior domain
                    InteriorDomain = InteriorUpdate                
                    
                    # Update Domain widths
                    DomainWidth = np.shape(InteriorDomain)[0]
                    ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain, ShrubDomainDead, DeadPercentCover = func.UpdateShrubDomains(DomainWidth, ShrubDomainWidth, ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, 
                                                                                                                                                                     ShrubPercentCover, BurialDomain, ShrubDomainDead, DeadPercentCover)
                                      
        # Record storm data
        StormCount.append(numstorm)
    
    
        # Add noise to flat parts of interior domain
        if t > StormStart:
            InteriorDomain = func.InteriorNoise(InteriorDomain, InteriorWidth)
    
    
        ###########################################      
        ### Ocean Shoreline Change
        
        ### Calculate shoreline/shoreface position change following Lorenzo-Trueba and Ashton (2014)  
        SCR, x_s, x_t, x_s_TS, x_t_TS, x_b_TS, s_sf_TS, QowTS, QsfTS = func.LTA_SC(InteriorDomain, OWloss, Qdg, DuneLoss, x_s, x_s_TS, x_t, x_t_TS, x_b_TS, s_sf_TS, InteriorWidth_Avg, SL, QowTS, QsfTS, t)
    
        
        ### Erode/Prograde Ocean Shoreline
        drown_break = 0
             
        SCRagg = SCRagg + SCR # Account for any residual shoreline change (less than cell size) from previous time step        
        
        if abs(SCRagg) >= 1:
            sc = math.floor(abs(SCRagg))
            if SCRagg > 0: # Positive = prograde, add row(s) to front of interior domain
                for d in range(sc):
                    InteriorDomain = np.vstack([DuneDomain[t,:,-1] + BermEl, InteriorDomain]) # New interior row added with elevation of previous dune field (i.e. previous dune row now part of interior)
                    if StormTimeSeries == 0:
                        newDuneHeight = np.ones([BarrierLength]) * (0.005 + (-0.005 + (0.005 - (-0.005)) * np.random.rand(BarrierLength)))
                    else:
                        newDuneHeight = np.ones([BarrierLength]) * 0.01
                    DuneDomain[t,:,:] = np.roll(DuneDomain[t,:,:], 1, axis=1)
                    DuneDomain[t,:,0] = newDuneHeight
                    
                    # Update shrub domains too
                    ShrubDomainFemale = np.concatenate((np.zeros([1,BarrierLength]), ShrubDomainFemale))
                    ShrubDomainMale = np.concatenate((np.zeros([1,BarrierLength]), ShrubDomainMale))
                    ShrubDomainDead = np.concatenate((np.zeros([1,BarrierLength]), ShrubDomainDead))
                    ShrubPercentCover = np.concatenate((np.zeros([1,BarrierLength]), ShrubPercentCover))
                    DeadPercentCover = np.concatenate((np.zeros([1,BarrierLength]), DeadPercentCover))
                    BurialDomain = np.concatenate((np.zeros([1,BarrierLength]), BurialDomain))
                DomainWidth = DomainWidth + sc # Update width of interior domain
                ShorelineChange = ShorelineChange + sc
                ShorelineChangeTS.append(+sc)
                SCRagg = SCRagg - sc # Reset, leaving residual
            elif SCRagg < 0: # Negative = erode, remove front row(s) of interior domain
                for d in range(sc):                 
                    newDuneElev = InteriorDomain[0,:] # First row of interior will now become new part of active dune field
                    newDuneHeight = newDuneElev - BermEl
                    newDuneHeight[newDuneHeight < DuneRestart] = DuneRestart
                    InteriorDomain = np.delete(InteriorDomain,(0), axis=0)
                    if np.shape(InteriorDomain)[0] <= 0:
                        drown_break = 1
                        break                     
                    DuneDomain[t,:,:] = np.roll(DuneDomain[t,:,:], -1, axis=1)
                    DuneDomain[t,:,-1] = newDuneHeight
    
                    # Update shrub domains too
                    ShrubDomainFemale = np.delete(ShrubDomainFemale,(0), axis=0)
                    ShrubDomainMale = np.delete(ShrubDomainMale,(0), axis=0)
                    ShrubDomainDead = np.delete(ShrubDomainDead,(0), axis=0)
                    ShrubPercentCover = np.delete(ShrubPercentCover,(0), axis=0)
                    DeadPercentCover = np.delete(DeadPercentCover,(0), axis=0)
                    BurialDomain = np.delete(BurialDomain,(0), axis=0)
                DomainWidth = DomainWidth - sc # Update width of interior domain
                ShorelineChange = ShorelineChange - sc
                ShorelineChangeTS.append(-sc)
                SCRagg = SCRagg + sc # Reset, leaving residual
        else:
            ShorelineChangeTS.append(0)          
        
        ### Check for drowning    
        if drown_break == 1:
            Result = 1
            TMAX = t-1
            break
        elif all(j <= SL for j in InteriorDomain[0,:]):
            Result = 2
            TMAX = t-1
            break
        
        ### Recalculate and save DomainWidth and InteriorWidth
        DomainWidth, InteriorWidth, InteriorWidth_Avg = func.FindWidths(InteriorDomain, SL)
        InteriorWidth_AvgTS.append(InteriorWidth_Avg)
        
    
        ########################################### 
        ### Save domains of this timestep
        DomainTS[t] = InteriorDomain 
        zero = np.zeros([DomainWidth, BarrierLength])    
        ShrubFemaleTS[t] = ShrubDomainFemale + zero 
        ShrubMaleTS[t] = ShrubDomainMale + zero
        ShrubDeadTS[t] = ShrubDomainDead + zero
        
        # Calculate Percent Cover and Area
        ShrubDomainAll = ShrubDomainMale + ShrubDomainFemale
        ShrubPercentCover, PercentCoverTS, ShrubArea, DeadPercentCover, DeadPercentCoverTS = func.CalcPC(ShrubDomainAll, PercentCoverTS, ShrubDomainDead, DeadPercentCoverTS, ShrubArea, t)
    
        # Find Subaerial Area of Island and Average Interior Elevation
        subaerial = (InteriorDomain >= SL).sum()
        IslandAreaTS.append(subaerial)
        AvgIntTemp = InteriorDomain[InteriorDomain < 1]
        AvgInteriorElevationTS.append(np.average(AvgIntTemp[AvgIntTemp >= SL]))
    
    ###########################################     
    ### End of model run
    
    SimDuration = time.time() - Time
    
    print()
    print('Elapsed Time: ', SimDuration, 'sec') # Print elapsed time of simulation
    
    
    
    #==============================================================================================================================================
    # CALC STATS      
    #==============================================================================================================================================


    # 23: Calculate shoreline change periodicity
    SCperiod, AvgFastDur, AvgSlowDur, Punc = func.calc_ShorelinePeriodicity(TMAX, x_s_TS)
    
    InteriorWidth_Avg = InteriorWidth_Avg * 10
    
    ShorelineChange = (x_s_TS[-1] - LShoreface) * 10
    
    Hd_avg = np.mean(Hd_AverageTS) * 10
    
    AvgInteriorElevation = AvgInteriorElevationTS[-1] * 10
    
    IslandArea = IslandAreaTS[-1] * 100

    t_Result = TMAX


    #==============================================================================================================================================
    # RETURN & SAVE RESULTS      
    #==============================================================================================================================================


    # Save sim data
    filename = directory + '/SimData_' + str(SimNumber) + '.npz'
    
    if Shrub_ON == 1:
        np.savez(filename, DuneDomain = DuneDomain, DomainTS = DomainTS, x_s_TS = x_s_TS, x_b_TS = x_b_TS, x_t_TS = x_t_TS, s_sf_TS = s_sf_TS, InteriorWidth_AvgTS = InteriorWidth_AvgTS, QowTS = QowTS, QsfTS = QsfTS, Hd_AverageTS = Hd_AverageTS, 
             PercentCoverTS = PercentCoverTS, DeadPercentCoverTS = DeadPercentCoverTS, ShrubArea = ShrubArea, ShrubDomainAll = ShrubDomainAll, ShrubDeadTS = ShrubDeadTS, StormCount = StormCount, t = t, RunUpCount = RunUpCount, 
             InundationCount = InundationCount, ShorelineChange = ShorelineChange, Hd_Loss_TS = Hd_Loss_TS, SimDuration = SimDuration, StormSeries = StormSeries, Dmax = Dmax, SL = SL, InteriorWidth_Avg = InteriorWidth_Avg, Hd_avg = Hd_avg, 
              IslandArea = IslandArea, t_Result = t_Result, Result = Result, SCperiod = SCperiod, AvgFastDur = AvgFastDur, AvgSlowDur = AvgSlowDur, Punc = Punc, AvgInteriorElevationTS = AvgInteriorElevationTS, 
              AvgInteriorElevation = AvgInteriorElevation, IslandAreaTS = IslandAreaTS, MaxAvgSlope = MaxAvgSlope, fluxLimit = fluxLimit, Dstart = Dstart, RSLR = RSLR, growthparam = growthparam, SimParams = SimParams)
    else: 
        np.savez(filename, DuneDomain = DuneDomain, DomainTS = DomainTS, x_s_TS = x_s_TS, x_b_TS = x_b_TS, x_t_TS = x_t_TS, s_sf_TS = s_sf_TS, InteriorWidth_AvgTS = InteriorWidth_AvgTS, QowTS = QowTS, QsfTS = QsfTS, Hd_AverageTS = Hd_AverageTS, 
             StormCount = StormCount, t = t, RunUpCount = RunUpCount, InundationCount = InundationCount, ShorelineChange = ShorelineChange, Hd_Loss_TS = Hd_Loss_TS, SimDuration = SimDuration, StormSeries = StormSeries, Dmax = Dmax, SL = SL,
             InteriorWidth_Avg = InteriorWidth_Avg, Hd_avg = Hd_avg, IslandArea = IslandArea, t_Result = t_Result, Result = Result, SCperiod = SCperiod, AvgFastDur = AvgFastDur, AvgSlowDur = AvgSlowDur, Punc = Punc,
             AvgInteriorElevationTS = AvgInteriorElevationTS, AvgInteriorElevation = AvgInteriorElevation, IslandAreaTS = IslandAreaTS, MaxAvgSlope = MaxAvgSlope, fluxLimit = fluxLimit, Dstart = Dstart, RSLR = RSLR, growthparam = growthparam, SimParams = SimParams)
    
    

    return Result, t_Result, InteriorWidth_Avg, ShorelineChange, Hd_avg, ShrubArea[-1], AvgInteriorElevation, SCperiod, AvgFastDur, AvgSlowDur, Punc
