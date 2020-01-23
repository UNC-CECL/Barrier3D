# Main runfile for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- Copyright (C) 2020 Ian R.B. Reeves (current developer)                                                                                                                                                                                          -
- Current developer can be contacted by email (reevesi@live.unc.edu) and paper mail (104 South Road, Mitchell Hall CB #3315, Chapel Hill, NC 27599, USA                                                                                           -
- This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. -
- This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.        -
- You should have received a copy of the GNU General Public License along with this program; if not, see <https://www.gnu.org/licenses/>.                                     -
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""

# Version Number: 1
# Updated: 22 Jan 2020

# All units are in decameters for simulation; convert to meters after simulation end for presentation

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
   



#==============================================================================================================================================
# IMPORTS       
#==============================================================================================================================================

import Barrier3D_Functions as func
import numpy as np
import math 
from mpl_toolkits.mplot3d import Axes3D # analysis:ignore
import time
import sys # analysis:ignore
from sys import platform
import os
import imageio # analysis:ignore
import warnings
warnings.filterwarnings("ignore")


Time = time.time()

print('TEMP!: Salt spray shrub impact turned OFF')

#==============================================================================================================================================
# SET PARAMETERS       
#==============================================================================================================================================

### Load Input Parameters
from Barrier3D_Parameters import (BarrierLength, BayDepth, BermEl, DomainWidth, DuneDomain, Dstart, DuneWidth, DShoreface, InteriorDomain, LShoreface, MHW, SD_storm, StormStart,
                                      TMAX, mean_storm, numstorm, Shrub_ON, Qshrub_max, C1, C2, DuneRestart, nn, mm, threshold_in, Rin_r, Rin_i, Qs_min, Kr, Ki, Cbb_r, Cbb_i)

### Sea Level
SL = 0 # Does not change (Lagrangian frame of reference)


### Initialize Shrubs
if Shrub_ON ==1: print('Shrubs ON')
PercentCoverTS = [None] * TMAX
PercentCoverTS[0] = np.zeros([DomainWidth, BarrierLength])
ShrubFemaleTS = [None] * TMAX
ShrubFemaleTS[0] = np.zeros([DomainWidth, BarrierLength])
ShrubMaleTS = [None] * TMAX
ShrubMaleTS[0] = np.zeros([DomainWidth, BarrierLength])
ShrubDomainFemale = ShrubFemaleTS[0]
ShrubDomainMale = ShrubMaleTS[0]
ShrubPercentCover = PercentCoverTS[0]
BurialDomain = np.zeros([DomainWidth, BarrierLength])
ShrubArea = [0]


### Initialize variables 
DomainTS = [None] * TMAX # Stores the elevation domain for each timestep
DomainTS[0] = InteriorDomain
SCR = 0 # (dam/yr) Change rate of ocean-fronting shoreline: (+) = progradation, (-) = erosion
SCRagg = 0 # Counter variable for shoreline change that is less than 1 cell size per year
ShorelineChangeTS = [0] 
ShorelineChange = 0 # (dam) Total amount of shoreline change 
BBShorelineChangeTS = [0]
BBShorelineChange = 0 # (dam) Total amount of back-barrier shoreline extension 
StormCount = []
InundationCount = 0
RunUpCount = 0
InteriorWidth_AvgTS = [DomainWidth]
QowTS = [] # (m^3/m)
x_t = 0  # (dam) Start location of shoreface toe
x_s = x_t + LShoreface # (dam) Start location of shoreline
x_t_TS = [x_t] # (dam) Shoreface toe locations for each time step
x_s_TS = [x_s] # (dam) Shoreline locations for each time step
x_b_TS = [(x_s + DomainWidth)] # (dam) Bay shoreline locations for each time step
h_b_TS = [BermEl] # (dam) Berm Elevation over time
s_sf_TS = [(DShoreface/LShoreface)]
Hd_AverageTS = [Dstart]
QsfTS = [0] # (m^3/m)




#==============================================================================================================================================
# RUN MODEL      
#==============================================================================================================================================

for t in range(1, TMAX): # Yearly time steps - actual time = t + 1
    
    ### Print time step to screen
    print("\r", 'Time Step: ', t, end = "")
    
    
    ### RSLR
    InteriorDomain, DuneDomain = func.SeaLevel(InteriorDomain, DuneDomain, t)
    
    
    ### Find DomainWidth and InteriorWidth
    DomainWidth, InteriorWidth, InteriorWidth_Avg = func.FindWidths(InteriorDomain, SL)
    
    
    ### Check for drowning
    # Barrier drowns if interior width thins to 10 m or fewer or all interior cells are below SL
    if InteriorWidth_Avg < 1: # if max(InteriorWidth) <= 1:
        print('Barrier has WIDTH DROWNED at t = ' + str(t) + ' years!')
        TMAX = t-1
        break
    if all(j <= SL for j in InteriorDomain[0,:]):
        print('Barrier has HEIGHT DROWNED at t = ' + str(t) + ' years!')
        TMAX = t-1
        break


    ### Grow Dunes
    DuneDomain, Dmax, Qdg = func.DuneGrowth(DuneDomain, t)

    DuneDomainCrest = DuneDomain[t,:,:].max(axis=1) # Maximum height of each row in DuneDomain
    Hd_AverageTS.append(np.mean(DuneDomainCrest)) # Store average pre-storms dune-heigh for time step
    
    
    ###########################################      
    ### Shrubs 
    
    if Shrub_ON == 1:

        ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, BurialDomain = func.Shrubs(InteriorDomain, DuneDomainCrest, t, ShrubDomainFemale, ShrubDomainMale, BurialDomain, InteriorWidth_Avg, DomainWidth)
    
    
    
    ###########################################      
    ### Storms

    OWloss = 0
    
    if t >= StormStart:
        # Select number of storms for this time step from normal distribution
        numstorm = int(numstorm)
        numstorm = round(np.random.normal(mean_storm, SD_storm)) # analysis:ignore # Comment out this line if using pre-specified static number of storms per year
        
        if numstorm > 0:
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
                
                Hd_TSloss = DuneChange.max(axis=1) / dur[n] # Average height of dune loss for each substep during storm


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
                if len(gaps) > 0 and Iow / len(gaps) >= threshold_in: # If greater than 25% of dune gaps are in unundation regime, use inundation regime routing
                    inundation = 1
#                    #### Calc substep
#                    Smax = 0
#                    Smin = 0
#                    for l in range(BarrierLength):
#                        for w in range(InteriorWidth[l]):
#                            slope = InteriorDomain[w,l] - InteriorDomain[w+1,l]
#                            if slope > Smax:
#                                Smax = slope
#                    calc maxQ
#                    accrete = f(maxQ, smax)
#                    ceil(accrete / (Smax *2)) = substep
#                    ####
                    substep = 15 # pre 14Jan20: 5, ~20 needed?
                    InundationCount += 1
                else:
                    inundation = 0
                    surges = (1 / period[n]) * 60 * 60 # Number of surges per hour
                    substep = 1 
                    RunUpCount += 1

                # Set Domain
                add = 25
                duration = dur[n] * substep
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
                inundation = 0 
                
                # Set Water and Sediment at Dune Crest
                for q in range(len(gaps)):
                    start = gaps[q][0]
                    stop = gaps[q][1]
                    Rexcess = gaps[q][2] # (m)
                    gapwidth = stop - start + 1
                    meandune = (sum(Dunes_prestorm[start:stop+1]) / gapwidth) + BermEl # Average elevation of dune gap
                    
                    # Calculate discharge through each dune cell
                    Vdune = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10 # (dam/surge)                      
                    Qdune =  Vdune * Rexcess * 3600 # (dam^3/hr)

                    # Set discharge at dune gap
                    Discharge[:,0,start:stop] = Qdune
                    
                    if inundation == 1: # Inundation regime       
                        Rin = Rin_i
                        slopes = []
                        # Find average slope of interior
                        for u in range(1,width-1):
                            for v in range(BarrierLength):                            
                                if Elevation[0,u,v] > SL and Elevation[0,u+1,v] > SL:
                                    SS = Elevation[0,u,v] - Elevation[0,u+1,v]
                                    slopes.append(SS)
                        AvgSlope = sum(slopes) / len(slopes)
                        C = 2 * AvgSlope                             
                        # Set sediment at dune gap
                        SedFluxIn[:,0,start:stop] = Ki * (Qdune*(AvgSlope + C))**mm        
                    else: # Run-up regime
                        Rin = Rin_r
                        # Set sediment at dune gap
                        SedFluxIn[:,0,start:stop] = Kr * Qdune           

  
                ### Run Flow Routing Algorithm                                      
                for TS in range(duration):
                    
                    ShrubDomainWidth = np.shape(ShrubDomainFemale)[0]
                             
                    if TS > 0:                      
                        Elevation[TS,1:,:] = Elevation[TS-1,1:,:] # Begin timestep with elevation from end of last                                    
                        Elevation[TS,0,:] = Dunes - Hd_TSloss * TS # Reduce dune in height linearly over course of storm

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
                                                                                             
                                
                                ### Reduce Overwash Through Shrub Cells and Save Discharge                                
                                if Shrub_ON == 1:             
                                    # Cell 1
                                    if i > 0:
                                        if d < (ShrubDomainWidth):
                                            if ShrubPercentCover[d,i-1] > 0: 
                                                Q1 =  Q1 * (Qshrub_max * ShrubPercentCover[d,i-1])                                    
                                        Discharge[TS,d+1,i-1] = Discharge[TS,d+1,i-1] + Q1
                                        
                                    # Cell 2                                    
                                    if d < (ShrubDomainWidth - 1):
                                        if ShrubPercentCover[d,i] > 0:
                                            Q2 =  Q2 * (Qshrub_max * ShrubPercentCover[d,i])
                                    Discharge[TS,d+1,i] = Discharge[TS,d+1,i] + Q2         
                                    
                                    # Cell 3
                                    if i < (BarrierLength-1):
                                        if d < (ShrubDomainWidth):
                                            if ShrubPercentCover[d,i+1] > 0:
                                                Q3 =  Q3 * (Qshrub_max * ShrubPercentCover[d,i+1])                                         
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
                                # Run-up regime
                                if inundation == 0:
                                    if Q1 > Qs_min and S1 >= 0:
                                        Qs1 = Kr * Q1
                                    else:
                                        Qs1 = 0
                                        
                                    if Q2 > Qs_min and S2 >= 0:
                                        Qs2 = Kr * Q2
                                    else:
                                        Qs2 = 0
                                        
                                    if Q3 > Qs_min and S3 >= 0:
                                        Qs3 = Kr * Q3
                                    else:
                                        Qs3 = 0
                                
                                # Inundation Regime (Murray and Paola (1994, 1997) Rule 3)
                                else:                                    
                                    if Q1 > Qs_min:
                                        Qs1 = Ki * (Q1*(S1 + C))**mm
                                    else:
                                        Qs1 = 0
                                        
                                    if Q2 > Qs_min:
                                        Qs2 = Ki * (Q2*(S2 + C))**mm
                                    else:
                                        Qs2 = 0
                                        
                                    if Q3 > Qs_min:
                                        Qs3 = Ki * (Q3*(S3 + C))**mm
                                    else:
                                        Qs3 = 0
                                
#                                # Inundation Regime (Murray and Paola (1994, 1997) Rule 5)
#                                else:                                    
#                                    Qt = np.mean(Discharge[0,0,:]) # Typical discharge
#                                    St = AvgSlope # Typical slope
#                                    ThD = 2 # Threshold denomenator, threshold therefore around half of the typical stream power
#                                    Th = Qt * St / ThD 
#                                    if Q1 > Qs_min:
#                                        Qs1 = Ki * (Q1 * (S1 + C) - Th)**mm
#                                    else:
#                                        Qs1 = 0
#                                        
#                                    if Q2 > Qs_min:
#                                        Qs2 = Ki * (Q2 * (S2 + C) - Th)**mm
#                                    else:
#                                        Qs2 = 0
#                                        
#                                    if Q3 > Qs_min:
#                                        Qs3 = Ki * (Q3 * (S3 + C) - Th)**mm
#                                    else:
#                                        Qs3 = 0
                                        
                                        
                                ### Calculate Net Erosion/Accretion
                                if Elevation[TS,d,i] > SL: # If cell is subaerial, elevation change is determined by difference between flux in vs. flux out
                                    if i > 0:
                                        SedFluxIn[TS,d+1,i-1] += Qs1
                                        
                                    SedFluxIn[TS,d+1,i] += Qs2
                                    
                                    if i < (BarrierLength-1):
                                        SedFluxIn[TS,d+1,i+1] += Qs3
                                    
                                    Qs_out = Qs1 + Qs2 + Qs3
                                    SedFluxOut[TS,d,i] = Qs_out
                                    
                                else: # If cell is subaqeous, exponentially decay deposition of remaining sediment across bay
                                    Qs_bb_min = 0.0001 # Assumes sediment flux less than this threshold can be ignored
                                    if inundation == 1: # Inundation regime

                                        SedFluxOut[TS,d,i] = SedFluxIn[TS,d,i] * Cbb_i
                                        Qs0 = SedFluxOut[TS,d,i]
                                        
                                        Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                        Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                        Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                        if i > 0 and Qs1 > Qs_bb_min:
                                            SedFluxIn[TS,d+1,i-1] += Qs1
                                        
                                        if Qs1 > Qs_bb_min:    
                                            SedFluxIn[TS,d+1,i] += Qs2
                                        
                                        if i < (BarrierLength-1) and Qs1 > Qs_bb_min:
                                            SedFluxIn[TS,d+1,i+1] += Qs3
                                    
                                    else: # Run-up regime

                                        SedFluxOut[TS,d,i] = SedFluxIn[TS,d,i] * Cbb_r
                                        Qs0 = SedFluxOut[TS,d,i]
                                        
                                        Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                        Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                        Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                        if i > 0 and Qs1 > Qs_bb_min:
                                            SedFluxIn[TS,d+1,i-1] += Qs1
                                        
                                        if Qs1 > Qs_bb_min:
                                            SedFluxIn[TS,d+1,i] += Qs2
                                        
                                        if i < (BarrierLength-1) and Qs1 > Qs_bb_min:
                                            SedFluxIn[TS,d+1,i+1] += Qs3                                             

                                                                      
                                ### Saline Flooding 
                                ShrubDomainFemale, ShrubDomainMale = func.SalineFlooding(ShrubDomainWidth, ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, d, i, Q0)  

                                           
                    ### Update Elevation After Every Storm Hour                      
                    if inundation == 1:
                        ElevationChange = (SedFluxIn[TS,:,:] - SedFluxOut[TS,:,:]) / substep
                    else:
                        ElevationChange = (SedFluxIn[TS,:,:] - SedFluxOut[TS,:,:]) / substep #/ surges #?
                    Elevation[TS,:,:] = Elevation[TS,:,:] + ElevationChange                    
                    
                    # Calculate and save volume of sedimen deposited on/behind the island for every hour
                    OWloss = OWloss + np.sum(ElevationChange[1:,:])
                    
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
                ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain = func.UpdateShrubDomains(DomainWidth, ShrubDomainWidth, ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain)
                           
                    
    # Record storm data
    StormCount.append(numstorm)


    ###########################################      
    ### Ocean Shoreline Change
    
    ### Calculate shoreline/shoreface position change following Lorenzo-Trueba and Ashton (2014)  
    SCR, x_s, x_t, x_s_TS, x_t_TS, x_b_TS, s_sf_TS, QowTS, QsfTS = func.LTA_SC(InteriorDomain, OWloss, x_s, x_s_TS, x_t, x_t_TS, x_b_TS, s_sf_TS, InteriorWidth_Avg, SL, QowTS, QsfTS)

    
    ### Erode/Prograde Ocean Shoreline
    drown_break = 0
         
    SCRagg = SCRagg + SCR # Account for any residual shoreline change (less than cell size) from previous time step        
    
    if abs(SCRagg) >= 1:
        sc = math.floor(abs(SCRagg))
        if SCRagg > 0: # Positive = prograde, add row(s) to front of interior domain
            for d in range(sc):
                InteriorDomain = np.vstack([DuneDomain[t,:,-1] + BermEl, InteriorDomain]) # New interior row added with elevation of previous dune field (i.e. previous dune row now part of interior)
                
                newDuneHeight = np.ones([BarrierLength]) * (0.01 + (-0.005 + (0.005 - (-0.005)) * np.random.rand(BarrierLength)))
                DuneDomain[t,:,:] = np.roll(DuneDomain[t,:,:], 1, axis=1)
                DuneDomain[t,:,0] = newDuneHeight
                
                # Update shrub domains too
                ShrubDomainFemale = np.concatenate((np.zeros([1,BarrierLength]), ShrubDomainFemale))
                ShrubDomainMale = np.concatenate((np.zeros([1,BarrierLength]), ShrubDomainMale))
                ShrubPercentCover = np.concatenate((np.zeros([1,BarrierLength]), ShrubPercentCover))
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
                ShrubPercentCover = np.delete(ShrubPercentCover,(0), axis=0)
                BurialDomain = np.delete(BurialDomain,(0), axis=0)
            DomainWidth = DomainWidth - sc # Update width of interior domain
            ShorelineChange = ShorelineChange - sc
            ShorelineChangeTS.append(-sc)
            SCRagg = SCRagg + sc # Reset, leaving residual
    else:
        ShorelineChangeTS.append(0)        
    
    ### Check for drowning    
    if drown_break == 1:
        print('Barrier has WIDTH DROWNED at t = ' + str(t) + ' years')
        TMAX = t-1
        break
    elif all(j <= SL for j in InteriorDomain[0,:]):
        print('Barrier has HEIGHT DROWNED at t = ' + str(t) + ' years')
        TMAX = t-1
        break
    
    ### Recalculate and save DomainWidth and InteriorWidth
    DomainWidth, InteriorWidth, InteriorWidth_Avg = func.FindWidths(InteriorDomain, SL)
    InteriorWidth_AvgTS.append(InteriorWidth_Avg)
    

    ########################################### 
    ### Save domains of this timestep
    DomainTS[t] = InteriorDomain 
    
    ShrubFemaleTS[t] = ShrubDomainFemale 
    ShrubMaleTS[t] = ShrubDomainMale
    ShrubDomainAll = ShrubDomainMale + ShrubDomainFemale
    
    # Calculate Percent Cover and Area
    ShrubPercentCover, PercentCoverTS, ShrubArea = func.CalcPC(ShrubDomainAll, PercentCoverTS, ShrubArea, t)


###########################################     
### End of model run
print()
print('Elapsed Time: ', time.time() - Time, 'sec') # Print elapsed time of simulation
if platform == "win32":
    import winsound
    winsound.Beep(800,200)
    winsound.Beep(800,200)
elif platform == "darwin":
    print('\a')
    os.system('say "barrier shrub simulation complete"')
elif platform == "linux" or platform == "linux2":
    print('\a')
    os.system('spd-say "barrier shrub simulation complete"')

 



#==============================================================================================================================================
# PLOT RESULTS      
#==============================================================================================================================================


# 1: Dune Height Over Time
func.plot_DuneHeight(DuneDomain)


# 2: Elevation Domain For Last Time Step
func.plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS)


## 3: Elevation Domain Frames
#func.plot_ElevFrames(TMAX, DomainTS)


## 4: Animation Frames of Barrier and Dune Elevation
func.plot_ElevAnimation(InteriorWidth_AvgTS, ShorelineChange, DomainTS, DuneDomain, SL, x_s_TS, Shrub_ON, PercentCoverTS, TMAX)


# 5: Cross-Shore Transect Every 100 m Alongshore For Last Time Step
func.plot_XShoreTransects(InteriorDomain, DuneDomain, SL, TMAX)


## 6: Shoreline Change Per Year Over Time <-------------- REDO!
#func.plot_ShorelineChange(ShorelineChangeTS, BBShorelineChangeTS)
      

## 7: Shoreline Change Rate Over Time <-------------- REDO!
#func.plot_ShorelineChangeRate(ShorelineChangeTS)
    

## 8: OW vs IN count
#func.plot_RuInCount(RunUpCount, InundationCount)


## 9: Shoreface LTA14 transects over time
#func.plot_LTATransects(SL, TMAX, x_b_TS, x_t_TS, x_s_TS)


## 10: Average Island Elevation Over Time <-------------- BROKEN!
#func.plot_AvgIslandElev(h_b_TS)


## 11: Shoreface Slope Over Time
#func.plot_ShorefaceSlope(s_sf_TS)


## 12: Average Interior Width Over Time
#func.plot_AvgInteriorWidth(InteriorWidth_AvgTS)


## 13: Shoreface Overwash Flux Over Time
#func.plot_OverwashFlux(QowTS)


# 14: Stats Summary: Width, Berm Elevation, SF Slope, Shoreline Change, and Overwash Flux Over Time (all in one)
func.plot_StatsSummary(s_sf_TS, x_s_TS, TMAX, InteriorWidth_AvgTS, QowTS, QsfTS, Hd_AverageTS) 


## 15: 3D Plot of Island Domain For Last Time Step
#func.plot_3DElevTMAX(TMAX, t, SL, DuneDomain, DomainTS)


## 16: 3D Animation Frames of Island Elevation (no translation) <-------------- BROKEN!
#func.plot_3DElevFrames(DomainTS, SL, TMAX, DuneDomain)
       

## 17: 3D Animation Frames of Island Evolution (with translation)
#func.plot_3DElevAnimation(DomainTS, SL, TMAX, DuneDomain, DomainWidth, x_s_TS, ShorelineChange)


## 18: Shrub Age Domain at Simulation End
#func.plot_ShrubAgeTMAX(ShrubDomainAll)


## 19: Percent Cover Domain at Simulation End
#func.plot_ShrubPercentCoverTMAX(PercentCoverTS, TMAX)


# 20: Shrub Area Over Time
func.plot_ShrubArea(ShrubArea)


## 21: Storm count over time
#func.plot_StormCount(StormCount)



