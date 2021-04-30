# Simulation and plotting functions for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 5
# Updated: 30 April 2021


# Simulation Functions Included: 
#   DuneGaps            Finds location and widths of OW throats in dune line
#   SeaLevel            Accounts for relative sea level rise, returns updated elevation domains
#   DuneGrowth          Grows dune cells
#   DiffuseDunes        Diffuses dune heights in alongshore direction
#   FindWidths          Finds DomainWidth and InteriorWidth
#   LTA_SCR             Calculates amounted of shoreline change for year following Lorenzo-Trueba and Ashton (2014)
#   Shrubs              Main function for shrub expansion and mortality 
#   UpdateShrubDomains  Updates size of shrub domains
#   SalineFlooding      Kills young shrubs flooded beyond threshold discharge (Tolliver et al., 1997)
#   UpdateBurial        Updates amount of burial/erosion for each shrub
#   CalcPC              Calculates shrub percent cover and shrub area (dam^2)



import numpy as np
import math
import random
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import imageio 
from scipy import signal

from Barrier3D_Parameters import (
                                  RSLR,
                                  BayDepth,
                                  Dmaxel, 
                                  BermEl, 
                                  growthparam, 
                                  DuneWidth, 
                                  BarrierLength,  
                                  HdDiffu,  
                                  DShoreface,
                                  LShoreface,
                                  k_sf, 
                                  s_sf_eq, 
                                  Qat, 
                                  Dshrub, 
                                  Female, 
                                  ShrubEl_min, 
                                  ShrubEl_max,
                                  MaxShrubHeight,
                                  SprayDist,
                                  BurialLimit, 
                                  UprootLimit, 
                                  TimeFruit, 
                                  Seedmin, 
                                  Seedmax, 
                                  GermRate, 
                                  disp_mu, 
                                  disp_sigma, 
                                  SalineLimit, 
                                  PC,
                                  Shrub_ON
                                  )

#==============================================================================================================================================
# SIMULATION FUNCTIONS      
#==============================================================================================================================================




#===================================================
# DUNEGAPS

# Returns tuple of [gap start index, stop index, avg Rexcess of each gap, alpha: ratio of TWL to dune height]

def DuneGaps(DuneDomain, Dow, bermel, Rhigh):
    
    gaps = []
    start = 0
    stop = 0
    i = start
    while i < (len(Dow)-1):
        adjacent = Dow[i+1] - Dow[i]
        if adjacent == 1:
            i = i+1
        else:
            stop = i
                
            x = DuneDomain[Dow[start]:(Dow[stop]+1)]
            Hmean = sum(x) / float(len(x))
            Rexcess = Rhigh - (Hmean + bermel)
            alpha = Rhigh/(Hmean + bermel)
            gaps.append([Dow[start], Dow[stop], Rexcess, alpha]) 
            
            start = stop+1
            i = start
    if i > 0:
        stop = i-1
        
        x = DuneDomain[Dow[start]:(Dow[stop]+1)]
        if len(x) > 0:
            Hmean = sum(x) / float(len(x))
            Rexcess = Rhigh - (Hmean + bermel)
            alpha = Rhigh/(Hmean + bermel)
            gaps.append([Dow[start], Dow[stop], Rexcess, alpha])
    return gaps





    
#===================================================
# SeaLevel
    
# Accounts for relative sea level rise, returns updated elevation domains

def SeaLevel(InteriorDomain, DuneDomain, t):
    
    # Decrease all elevation this year by RSLR increment
    InteriorDomain = InteriorDomain - RSLR[t]
    DuneDomain[t-1] = DuneDomain[t-1] - RSLR[t]    
    InteriorDomain[InteriorDomain < -BayDepth] = -BayDepth # Bay can't be deeper than BayDepth (roughly equivalent to constant back-barrier slope)
        
    return InteriorDomain, DuneDomain





#===================================================
# DuneGrowth
    
# Grows dune cells

def DuneGrowth(DuneDomain, t):
        
    # Set max dune height 
    Dmax = Dmaxel - BermEl # (dam) Max dune height 
    if Dmax < 0:
        Dmax = 0  
    
    Cf = 3 # Decay coefficient     
    Qdg = 0
    # Grow dune
    for q in range(DuneWidth):
        reduc = 1/(Cf**q)
        G = growthparam * DuneDomain[t-1,:,q] * (1 - DuneDomain[t-1,:,q] / Dmax) * reduc
        DuneDomain[t,:,q] = G + DuneDomain[t-1,:,q]
        Qdg = Qdg + (np.sum(G) / BarrierLength) # Volume of sediment lost from beach/shoreface from dune growth
    
    return DuneDomain, Dmax, Qdg





#===================================================
# DuneGrowthBatch
    
# Grows dune cells (Version for batch runs takes growthparam as direct input)

def DuneGrowthBatch(DuneDomain, t, growthparam):

    growthparam = growthparam[0:BarrierLength]
    
    # Set max dune height
    Dmax = Dmaxel - BermEl # (dam) Max dune height 
    if Dmax < 0:
        Dmax = 0  
    
    Cf = 3 # Decay coefficient
    Qdg = 0
    # Grow dune
    for q in range(DuneWidth):
        reduc = 1/(Cf**q)
        G = growthparam * DuneDomain[t-1,:,q] * (1 - DuneDomain[t-1,:,q] / Dmax) * reduc
        DuneDomain[t,:,q] = G + DuneDomain[t-1,:,q]
        Qdg = Qdg + (np.sum(G) / BarrierLength) # Volume of sediment lost from beach/shoreface from dune growth
    
    return DuneDomain, Dmax, Qdg





#===================================================
# DiffuseDunes
    
# Dune height diffusion in alongshore direction: smoothes dune height by redistributing sand from high dune to neighboring low dune(s)
    
def DiffuseDunes(DuneDomain, t):
    
    for w in range(DuneWidth):
        # Alongshore direction 
        for d in range(2,BarrierLength): # Loop L to R
            Ldiff = DuneDomain[t,d,w] - DuneDomain[t,d-1,w] # Calculate height difference
            # Subtract sand from taller dune cell and add to adjacent shorter one (assumes sand within dunes is conserved)
            if Ldiff > HdDiffu:
                sub = (Ldiff - HdDiffu)/2 # Height difference in excess of maximum, divided by two
                DuneDomain[t,d,w] = DuneDomain[t,d,w] - sub # Lower dune that's too tall
                DuneDomain[t,d-1,w] = DuneDomain[t,d-1,w] + sub # Raise dune that's too short
        for d in range((BarrierLength-2),0,-1): # Loop R to L
            Rdiff = DuneDomain[t,d,w] - DuneDomain[t,d+1,w]
            if Rdiff > HdDiffu:
                sub = (Rdiff - HdDiffu)/2
                DuneDomain[t,d,w] = DuneDomain[t,d,w] - sub
                DuneDomain[t,d+1,w] = DuneDomain[t,d+1,w] + sub           

    return DuneDomain



#===================================================
# FindWidths
    
# Finds DomainWidth and InteriorWidth
# DomainWidth is wide enough to fit widest part of island; InteriorWidth stores the width of the island at each row of the domain where elevation is greater than 0      
    
def FindWidths(InteriorDomain, SL):
    
    DomainWidth = np.shape(InteriorDomain)[0] # (dam) analysis:ignore
    InteriorWidth = [0] * BarrierLength
    for l in range(BarrierLength):
        width = next((index for index, value in enumerate(InteriorDomain[:,l]) if value <= SL), DomainWidth)
        width = width - 1
        if width < 0: width = 0
        InteriorWidth[l] = width
    # Average width of island
    InteriorWidth_Avg = np.nanmean(InteriorWidth)

    return DomainWidth, InteriorWidth, InteriorWidth_Avg




#===================================================
# InteriorNoise
    
# Adds elevation noise to flat parts of the interior to represent aeolian re-working during interstorm periods     
    
def InteriorNoise(InteriorDomain, InteriorWidth):
    
    zthresh = 0.005 # dam
    ythresh = 2 # dam
    gap = 1 # dam
    
    # Loop through each column and 
    for x in range(BarrierLength):
        col = InteriorDomain[0:InteriorWidth[x],x]
    
        # find continuous sections that are too flat
        der = np.diff(col)
        
        der_under = np.where(der < zthresh)[0]
        
        if len(der_under) > 0:
            gaps = np.diff(der_under) > gap
            flat_start = np.insert(der_under[1:][gaps], 0, der_under[0])
            flat_stop = np.append(der_under[:-1][gaps], der_under[-1])
        
            # Add noise to flat sections
            for n in range(1,len(flat_stop)):
                flat_length = flat_stop[n] - flat_start[n] 
                if flat_length >= ythresh:
                    # noise = np.random.rand(flat_length) * zthresh              
                    noise = np.random.uniform(-1, 1, flat_length) * zthresh*2
                    InteriorDomain[flat_start[n]:flat_stop[n],x] += noise

    return InteriorDomain





#===================================================
# LTA_SC

# Finds shoreline change for modeled year following Lorenzo-Trueba and Ashton (2014)

def LTA_SC(InteriorDomain, OWloss, Qdg, DuneLoss, x_s, x_s_TS, x_t, x_t_TS, x_b_TS, s_sf_TS, InteriorWidth_Avg, SL, QowTS, QsfTS, t):
    
    # Find volume of shoreface/beach/dune sand deposited in island interior and back-barrier
    Qow = (OWloss) / (BarrierLength) # (dam^3/dam) Volume of sediment lost from shoreface/beach by overwash
    if Qow < 0:
        Qow = 0
    QowTS.append(Qow * 100) # Save in m^3/m   
    if DuneLoss < Qow:
        Qow = Qow - DuneLoss # Account for dune contribution to overwash volume                     
    else:
        Qow = 0 # Excess DuneLoss assumed lost offshore
    
    if Qdg < 0:
        Qdg = 0
    
    # DefineParams
    d_sf = DShoreface 
    h_b = np.average(InteriorDomain[InteriorDomain >= SL]) 
    x_b = x_s + InteriorWidth_Avg
        
    # Shoreface Flux
    s_sf = d_sf / (x_s - x_t)
    Qsf = (k_sf / 100) * (s_sf_eq - s_sf) # Convert Ksf from m^3/m/yr to dam^3/dam/yr
    QsfTS.append(Qsf * 100) # Save in m^3/m    

    # Toe, Shoreline, and island base elevation changes     
    x_t_dt = (4 * Qsf * (h_b + d_sf) / (d_sf * (2 * h_b + d_sf))) + (2 * RSLR[t] / s_sf)   
    x_s_dt = 2 * (Qow + Qdg + Qat) / ((2 * h_b) + d_sf) - (4 * Qsf * (h_b + d_sf) / (((2 * h_b) + d_sf)**2)) # Dune growth and alongshore transport added to LTA14 formulation
    
    # Record changes
    x_t = x_t + x_t_dt
    x_s = x_s + x_s_dt    
    x_t_TS.append(x_t)
    x_s_TS.append(x_s)
    x_b_TS.append(x_b)
    s_sf_TS.append(s_sf)
    
    SCR = x_s_dt * -1 # (dam) Positive increase in x-location means shoreline erosion


    return SCR, x_s, x_t, x_s_TS, x_t_TS, x_b_TS, s_sf_TS, QowTS, QsfTS

    



#===================================================
# Shrubs

# Main shrub expansion and mortality function
    
def Shrubs(InteriorDomain, DuneDomainCrest, t, ShrubDomainFemale, ShrubDomainMale, ShrubPercentCover, DeadPercentCover, BurialDomain, InteriorWidth_Avg, DomainWidth, ShrubDomainDead, SL):
    
    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
    
    BeachWidth = 50 /10 # Calculated based on berm (dune toe) elevation of 1.9, and beach slope of 0.04 (i.e., 1.9/0.04 = 47.5 ~= 50 m)
     
    
    ### Burial / Uprooting   
    ShrubHeight = ShrubPercentCover * MaxShrubHeight
    DeadHeight = DeadPercentCover * MaxShrubHeight

    for l in range(len(ShrubPercentCover[0])):
        for w in range(len(ShrubPercentCover)):
            if ShrubDomainAll[w,l] > 0:           
                if BurialDomain[w,l] > ShrubHeight[w,l] or BurialDomain[w,l] < UprootLimit: 
                    ShrubDomainFemale[w,l] = 0
                    ShrubDomainMale[w,l] = 0
                                        
                elif BurialDomain[w,l] > (BurialLimit * ShrubHeight[w,l]):
                    ShrubDomainDead[w,l] = ShrubDomainAll[w,l]
                    ShrubDomainFemale[w,l] = 0
                    ShrubDomainMale[w,l] = 0 
                                        
            elif ShrubDomainDead[w,l] > 0:
                if BurialDomain[w,l] > DeadHeight[w,l] or BurialDomain[w,l] < UprootLimit:
                    ShrubDomainDead[w,l] = 0
    
    tempAll = ShrubDomainFemale + ShrubDomainMale + ShrubDomainDead
    BurialDomain[tempAll == 0] = 0 # Reset burial domain

    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale # Recalculate


    ### Inundation
    # Kill shrubs that have fallen below minimum elevation, remove shrubs that have fallen below Mean High Water (i.e. elevation of 0) (passive loss from rising back-barrier water elevations)
    for w in range(DomainWidth):
        for l in range(BarrierLength):
            if InteriorDomain[w,l] < ShrubEl_min:
                if InteriorDomain[w,l] > SL:
                    if ShrubDomainFemale[w,l] > 0 or ShrubDomainMale[w,l] > 0:
                        ShrubDomainDead[w,l] = ShrubDomainFemale[w,l] + ShrubDomainMale[w,l] # Shrub remains as dead if above MHW
                        ShrubDomainFemale[w,l] = 0
                        ShrubDomainMale[w,l] = 0       
                else:
                    ShrubDomainDead[w,l] = 0  # Remove dead shrubs below MHW
                    ShrubDomainFemale[w,l] = 0
                    ShrubDomainMale[w,l] = 0    
                    
    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale # Recalculate                          
    
            
    ### Age the shrubs in years
    ShrubDomainFemale[ShrubDomainFemale > 0] += 1            
    ShrubDomainMale[ShrubDomainMale > 0] += 1
    
    
    #================================================
    # Drop Seed
    # Randomly drop a seed onto the island each time step
    randX = np.random.randint(0,BarrierLength)
    randY = np.random.randint(0,max(1,InteriorWidth_Avg))
    if ShrubDomainFemale[randY,randX] == 0 and ShrubDomainMale[randY,randX] == 0 and ShrubDomainDead[randY,randX] == 0 and DuneDomainCrest[randX] + BermEl >= Dshrub and InteriorDomain[randY,randX] >= ShrubEl_min \
        and InteriorDomain[randY,randX] <= ShrubEl_max:
        if random.random() > Female:
            ShrubDomainFemale[randY, randX] = 1
        else:
            ShrubDomainMale[randY, randX] = 1
     
    # Randomly drop a seed onto the LEFT SIDE of ther island each time step   
    # randX = np.random.randint(0,BarrierLength*0.02) # Can drop anywhere with the first 2% of interior columns    
    # randY = np.random.randint(0,max(1,InteriorWidth_Avg))
    # if ShrubDomainFemale[randY,randX] == 0 and ShrubDomainMale[randY,randX] == 0 and ShrubDomainDead[randY,randX] == 0 and DuneDomainCrest[randX] + BermEl >= Dshrub and InteriorDomain[randY,randX] >= ShrubEl_min \
    #     and InteriorDomain[randY,randX] <= ShrubEl_max:
    #     if random.random() > Female:
    #         ShrubDomainFemale[randY, randX] = 1
    #     else:
    #         ShrubDomainMale[randY, randX] = 1
    #================================================
    
    
    ### Disperse seeds
    for k in range(BarrierLength): # Loop through each row of island width (i.e. from ocean to mainland side of island)
        if 0 in ShrubDomainAll:                
          
            # For all cells with a shrub
            FemaleShrubs = ShrubDomainFemale[:,k]
            I = [index for index, value in enumerate(FemaleShrubs) if value >= TimeFruit]
            numShrubCells = len(I)
            # Determine how many seeds in each cell
            Seedspercell = np.random.randint(Seedmin,high = Seedmax, size = numShrubCells)
            # For each shrub cell, determine survival rate for the cell in this year
            SeedSurvYear = GermRate * np.random.rand(numShrubCells)
    
            for i in range(numShrubCells):
                # For each shrub cell producing seeds, generate a random # of random numbers each representing a single seed
                randseeds = np.random.rand(Seedspercell[i])
                # Find how many seeds produced in each cell actually survive
                Survivors = len(randseeds[randseeds < SeedSurvYear[i]])
    
                # Determine distance, rounding to nearest integer
                if Survivors > 0:
                    DispDist = np.round(np.random.lognormal(disp_mu,disp_sigma,Survivors))
                else:
                    DispDist = []
    
                # If there are seeds to disperse
                if len(DispDist) > 0:
                    for j in range(len(DispDist+1)): # Loop through each individual seed to disperse
                        # Create a meshgrid to find coordinates of all points that are dropdistance from origin 
                        if DispDist[j] > 0:
                            gridsize = int(DispDist[j] * 2 + 2)
                            X, Y = np.meshgrid(range(1, gridsize), range(1, gridsize))
                            originX = I[i] # Sets coordinates of plant where seed is coming from
                            originY = k
                            matOriginX = math.ceil(len(X)/2) 
                            matOriginY = matOriginX
                            distMat = np.round(np.sqrt((X-matOriginX)**2 + (Y-matOriginY)**2)) # Find the distance from origin to every other point on island
                            # Find coordinates of each point on island that is dropdistance away from origin
                            coords = np.where(distMat == DispDist[j]) 
                            row = coords[0]
                            col = coords[1] 
                            
                            # Randomly selct one of those points as the target point - this means equal probability of dispersal in every direction (valid assumption for avian dispersal)
                            if len(col) == 0:
                                targetY = originX
                                targetX = originY
                            else:
                                xx = random.randint(0, (len(col))-1)
                                matTargetX = col[xx]
                                matTargetY = row[xx]
                                
                                targetY = originX + (matTargetX-matOriginX)
                                targetX = originY + (matTargetY-matOriginY)
                                                     
                            # Put a plant in the ground if 
                            #   -the dropdistance>0, 
                            #   -the target drop location is within the island domain, 
                            #   -there is no shrub at the receiving cell (male, female, or dead),
                            #   -the receiving cell has a tall enough fronting dune
                            #   -the receiving cell is within elevation range
                            if  targetY >= 0 and targetY < DomainWidth and targetX >= 0 and targetX < BarrierLength and ShrubDomainFemale[targetY,targetX] == 0 and ShrubDomainMale[targetY,targetX] == 0 \
                                and InteriorDomain[targetY,targetX] >= ShrubEl_min and InteriorDomain[targetY,targetX] <= ShrubEl_max and ShrubDomainDead[targetY,targetX] < 1:
                                if DuneDomainCrest[targetX] + BermEl >= Dshrub: # Shrubs can establish if height of fronting dune greater than threshold
                                    # Decide if the tree wll be a female or male
                                    if random.random() > Female:
                                        ShrubDomainFemale[targetY,targetX] = 1
                                    else:
                                        ShrubDomainMale[targetY,targetX] = 1
                                elif BeachWidth + DuneWidth + targetY > SprayDist: # Shrubs can establish without dune if locaed a certain distance from shoreline...
                                    if random.random() > 0.5: # ... but have about 50% chance of survival 
                                        # Decide if the tree wll be a female or male
                                        if random.random() > Female:
                                            ShrubDomainFemale[targetY,targetX] = 1
                                        else:
                                            ShrubDomainMale[targetY,targetX] = 1
    
    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
  
    return ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, BurialDomain, ShrubDomainDead

    



#===================================================
# UpdateShrubDomains

# Updates size of shrub domains
    
def UpdateShrubDomains(DomainWidth, ShrubDomainWidth, ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain, ShrubDomainDead, DeadPercentCover):

    if DomainWidth > ShrubDomainWidth:
        AddRows = np.zeros([DomainWidth-ShrubDomainWidth,BarrierLength])
        ShrubDomainFemale = np.vstack([ShrubDomainFemale, AddRows])
        ShrubDomainMale = np.vstack([ShrubDomainMale, AddRows])
        ShrubDomainDead = np.vstack([ShrubDomainDead, AddRows])
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
        ShrubPercentCover = np.vstack([ShrubPercentCover, AddRows])
        DeadPercentCover = np.vstack([DeadPercentCover, AddRows])
        BurialDomain = np.vstack([BurialDomain, AddRows])
    elif DomainWidth < ShrubDomainWidth:
        RemoveRows = ShrubDomainWidth - DomainWidth
        ShrubDomainFemale = ShrubDomainFemale[0:-RemoveRows,:]
        ShrubDomainMale = ShrubDomainMale[0:-RemoveRows,:]
        ShrubDomainDead = ShrubDomainDead[0:-RemoveRows,:]
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
        ShrubPercentCover = ShrubPercentCover[0:-RemoveRows,:]
        DeadPercentCover = DeadPercentCover[0:-RemoveRows,:]
        BurialDomain = BurialDomain[0:-RemoveRows,:] 
        
    return ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain, ShrubDomainDead, DeadPercentCover





#===================================================
# SalineFlooding

# Kill all immature (< 1 yr-old) shrubs that have been flooded beyond a threshold discharge (Tolliver et al., 1997)
    
def SalineFlooding(ShrubDomainWidth, ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead, d, i, Q0):
    
    if d < (ShrubDomainWidth-1) and i < (BarrierLength - 1):
        if Q0 >= SalineLimit and ShrubDomainAll[d,i] == 1:
            ShrubDomainDead[d,i] = ShrubDomainMale[d,i] + ShrubDomainFemale[d,i] # Transfer to dead shrub domain
            ShrubDomainFemale[d,i] = 0
            ShrubDomainMale[d,i] = 0
               
    return ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead             

                   
                   

                   
#===================================================
# UpdateBurial

# Updates amount of burial/erosion for each shrub
    
def UpdateBurial(BurialDomain, ElevationChange, ShrubDomainWidth, ShrubDomainAll):
                                       
    BurialDomain = BurialDomain + ElevationChange[1:ShrubDomainWidth+1,:]                                
    BurialDomain[ShrubDomainAll == 0] = 0            
                   
    return BurialDomain     

    



#===================================================
# CalcPC

# Calculates percent cover of shrub domain and shrub coverage area (dam^2)
    
def CalcPC(ShrubDomainAll, PercentCoverTS, ShrubDomainDead, DeadPercentCoverTS, ShrubArea, t):
    
    Allshrub_t = ShrubDomainAll.astype('int64')
    Deadshrub_t = ShrubDomainDead.astype('int64')
    ShrubPercentCover = PC.take(Allshrub_t)
    DeadPercentCover = PC.take(Deadshrub_t)

#    ShrubPercentCover[ShrubPercentCover == 0] = DeadPercentCover[ShrubPercentCover == 0]
    PercentCoverTS[t] = ShrubPercentCover
    DeadPercentCoverTS[t] = DeadPercentCover
    ShrubArea.append(np.count_nonzero(ShrubDomainAll))

    return ShrubPercentCover, PercentCoverTS, ShrubArea, DeadPercentCover, DeadPercentCoverTS





#==============================================================================================================================================
# PLOTTING & CALCULATION FUNCTIONS      
#==============================================================================================================================================
    

#===================================================
# 1: Dune Height Over Time

def plot_DuneHeight(DuneDomain, Dmax):
        
    DuneCrest = DuneDomain.max(axis=2)
    duneFig = plt.figure(figsize=(14,8))
    plt.rcParams.update({'font.size':13})
    ax = duneFig.add_subplot(111)
    ax.matshow((DuneCrest)*10, origin='lower', cmap='bwr', aspect='auto', vmin=0, vmax=Dmax*10)
    cax = ax.xaxis.set_ticks_position('bottom') # analysis:ignore
    #cbar = duneFig.colorbar(cax)
    #cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Year')
    plt.title('Dune Height (m)')
    plt.show()
    name = 'Output/Dunes'
    duneFig.savefig(name)


    

#===================================================
# 2: Elevation Domain For Last Time Step

def plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS, Shrub_ON, PercentCoverTS, DeadPercentCoverTS):

    if TMAX > t:
        TMAX = t
    Dunes = (DuneDomain[TMAX,:,:] + BermEl) * 10
    Dunes = np.rot90(Dunes)
    Dunes = np.flipud(Dunes)
    DuneWidth = int(len(Dunes))
    Domain = DomainTS[TMAX] * 10
    Domain = np.vstack([Dunes, Domain])  
    if Shrub_ON == 1:
        Shrubs = PercentCoverTS[t]
        Dead = DeadPercentCoverTS[t]
        Sy, Sx = np.argwhere(Shrubs > 0).T
        Sz = Shrubs[Sy,Sx]*80 #30
        Dy, Dx = np.argwhere(Dead > 0).T
        Dz = Dead[Dy,Dx]*80 #22 
    elevFig1 = plt.figure(figsize=(14,8))
    ax = elevFig1.add_subplot(111)
    cax = ax.matshow(Domain, origin='lower', cmap='terrain', vmin=-1.1, vmax=4.0)#, interpolation='gaussian') # analysis:ignore
    if Shrub_ON == 1:
        ax.scatter(Sx, Sy + DuneWidth, marker='$*$', s=Sz, c='black', alpha=0.7, edgecolors='none')
        ax.scatter(Dx, Dy + DuneWidth, marker='$*$', s=Dz, c='red', alpha=0.7, edgecolors='none')  
    ax.xaxis.set_ticks_position('bottom')
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Cross-Shore Diatance (dam)')
    plt.title('Interior Elevation (m)')
    timestr = 'Time = ' + str(TMAX) + ' yrs'
    plt.text(1, 1, timestr)
    plt.tight_layout()
    plt.show()
    name = 'Output/FinalElevation'
    elevFig1.savefig(name)




#===================================================
# 3: Elevation Domain Frames

def plot_ElevFrames(TMAX, DomainTS):

    for t in range(TMAX):
        elevFig1 = plt.figure(figsize=(14,5))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(DomainTS[t], origin='lower', cmap='terrain', vmin=-1.1, vmax=4.0)
        ax.xaxis.set_ticks_position('bottom')
        cbar = elevFig1.colorbar(cax) # analysis:ignore
        #cbar.set_label('Elevation (m)', rotation=270)
        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Cross-Shore Diatance (dam)')
        plt.title('Interior Elevation')
        timestr = 'Time = ' + str(t) + ' yrs'
        plt.text(1, 1, timestr)
        plt.show()
        name = 'Output/SimFrames/elev_' + str(t)
        elevFig1.savefig(name)




#===================================================
# 4: Animation Frames of Barrier and Dune Elevation

def plot_ElevAnimation(InteriorWidth_AvgTS, ShorelineChange, DomainTS, DuneDomain, SL, x_s_TS, Shrub_ON, PercentCoverTS, TMAX, DeadPercentCoverTS):
        
    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange/10) + OriginY + 25)
    
    for t in range(TMAX):
        # Build beach elevation domain
        BeachDomain = np.zeros([BeachWidth, BarrierLength])    
        berm = math.ceil(BeachWidth*0.65)
        BeachDomain[berm:BeachWidth+1,:] = BermEl
        add = (BermEl-SL) / berm
        for i in range(berm):
            BeachDomain[i,:] = SL + add * i 
    
        # Make animation frame domain
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t,:,:] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength]) *-1
        widthTS = len(Domain)
        scts = [(x - x_s_TS[0]) for x in x_s_TS]
        if scts[t] >= 0:
            OriginTstart = OriginY + math.floor(scts[t])
        else:
            OriginTstart = OriginY + math.ceil(scts[t])        
        OriginTstop = OriginTstart + widthTS
        AnimateDomain[OriginTstart:OriginTstop, 0:BarrierLength] = Domain    
        if Shrub_ON == 1:
            Shrubs = PercentCoverTS[t]
            Dead = DeadPercentCoverTS[t]
            wid = np.zeros([BeachWidth + DuneWidth + OriginTstart, BarrierLength])
            Shrubs = np.vstack([wid, Shrubs])
            Dead = np.vstack([wid, Dead])
            Sy, Sx = np.argwhere(Shrubs > 0).T
            Sz = Shrubs[Sy,Sx]*80 #30
            Dy, Dx = np.argwhere(Dead > 0).T
            Dz = Dead[Dy,Dx]*80 #22
            
        # Plot and save
        elevFig1 = plt.figure(figsize=(15,13))
        ax = elevFig1.add_subplot(1,1,1)
        cax = ax.matshow(AnimateDomain, origin='lower', cmap='terrain', vmin=-1.1, vmax=4.0)#, interpolation='gaussian') # analysis:ignore
        if Shrub_ON == 1:
            ax.scatter(Sx, Sy, marker='$*$', s=Sz, c='black', alpha=0.7, edgecolors='none')
            ax.scatter(Dx, Dy, marker='$*$', s=Dz, c='red', alpha=0.7, edgecolors='none')
        ax.xaxis.set_ticks_position('bottom')
        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Cross-Shore Diatance (dam)')
        plt.title('Interior Elevation')
        plt.tight_layout()
        timestr = 'Time = ' + str(t) + ' yrs'
        newpath = 'Output/SimFrames/'
        if not os.path.exists(newpath):
            os.makedirs(newpath)            
        plt.text(1, 1, timestr)
        elevFig1.subplots_adjust(right=.85)
        # cbar_ax = elevFig1.add_axes([0.78, 0.15, 0.04, 0.70])
        # cbar = elevFig1.colorbar(cax, cax=cbar_ax)
        # cbar.set_label('Elevation (m)', rotation=270, labelpad=20)       
        name = 'Output/SimFrames/elev_' + str(t)
        elevFig1.savefig(name) # dpi=200
        plt.close(elevFig1)
        
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/elev_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/elev.gif', frames, 'GIF-FI')
    print()
    print('[ * GIF successfully generated * ]')
        
      
    
    
#===================================================        
# 5: Cross-Shore Transect Every 100 m Alongshore For Last Time Step
  
def plot_XShoreTransects(InteriorDomain, DuneDomain, SL, TMAX):

    # Build beach elevation
    BW = 6 # beach width (dam) for illustration purposes
    BeachX = np.zeros(BW)
    berm = math.ceil(BW*0.5)
    BeachX[berm:BW+1] = BermEl
    add = (BermEl-SL) / berm
    for i in range(berm):
        BeachX[i] = SL + add * i 
    # Plot full tranects    
    plt.figure()
    for v in range(0, BarrierLength, 10):
        CrossElev = InteriorDomain[:, v]
        Dunes = DuneDomain[TMAX-1, v, :] + BermEl
        CrossElev1 = np.insert(CrossElev,0,Dunes)
        CrossElev2 = np.insert(CrossElev1,0,BeachX)
        CrossElev = CrossElev2 * 10 # Convert to meters
        plt.plot(CrossElev)
    fig = plt.gcf()
    fig.set_size_inches(14,6)
    plt.hlines(SL,-1,len(CrossElev+1),colors='dodgerblue')
    plt.xlabel('Cross-Shore Distance (dam)')
    plt.ylabel('Elevation (m)')
    plt.title('Cross-shore Topo Transects')
    plt.show()
    name = 'Output/Profiles'
    fig.savefig(name)       
    
   
    
    
#===================================================    
# 6: Shoreline Positions Over Time   

def plot_ShorelinePositions(x_s_TS, x_b_TS):
    
    scts = [(x - x_s_TS[0]) * -10 for x in x_s_TS]
    bscts = [(x - x_s_TS[0]) * -10 for x in x_b_TS]
    shorelinefig = plt.figure()
    plt.plot(scts, 'b')
    plt.plot(bscts, 'g')
    fig = plt.gcf()
    fig.set_size_inches(14,8)
    plt.ylabel('Shoreline Position (m)')
    plt.xlabel('Year')
    plt.show()   
    name = 'Output/Shorelines'
    shorelinefig.savefig(name)   
   
    
    
#===================================================    
# 7: Shoreline Change Rate Over Time

def plot_ShorelineChangeRate(x_s_TS):
    
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
    rate = [0]
    for k in range(1,len(scts)):
        rate.append(scts[k]- scts[k-1])
    ratefig = plt.figure()
    plt.plot(rate)
    fig = plt.gcf()
    fig.set_size_inches(14 ,5)
    plt.xlabel('Year')
    plt.ylabel('Shoreline Erosion Rate(m/yr)')
    plt.show() 
    name = 'Output/ShorelineRate'
    ratefig.savefig(name)  
    
    
#===================================================
# 8: Run-up vs Inundation count

def plot_RuInCount(RunUpCount, InundationCount):

    objects = ('Run Up', 'Inundation')
    y_pos = np.arange(len(objects))
    plt.figure()
    plt.bar(y_pos, [RunUpCount, InundationCount])
    plt.xticks(y_pos, objects)
    plt.ylabel('Count')




#===================================================
# 9: Shoreface LTA14 transects over time

def plot_LTATransects(SL, TMAX, x_b_TS, x_t_TS, x_s_TS):
    

    xmax = x_b_TS[TMAX-1] + 20
    
    SFfig = plt.figure(figsize=(20,5))
    colors = plt.cm.jet(np.linspace(0,1,TMAX))
    
    for t in range(0,TMAX,25): # Plots one transect every 25 years
        # Create data points
        Tx = x_t_TS[t]
        Ty = ((SL + (t * RSLR)) - DShoreface) *10
        Sx = x_s_TS[t]
        Sy = (SL + (t * RSLR)) *10
        Bx = x_b_TS[t]
        By = ((SL + (t * RSLR)) - BayDepth) *10
        Hx1 = Sx
        Hy1 = ((t * RSLR) + BermEl) *10
        Hx2 = Bx
        Hy2 = Hy1
        Mx = xmax
        My = By
        
        x = [Tx, Sx, Hx1, Hx2, Bx, Mx]
        y = [Ty, Sy, Hy1, Hy2, By, My]
        
        # Plot
        plt.plot(x,y,color = colors[t])
        
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Elevation (m)')
    plt.title('Shoreface Evolution')
    plt.show()
    
    # Save   
    name = 'Output/Shoreface'
    SFfig.savefig(name)
    
    
    
    
#===================================================
# 10: Average Island Elevation Over Time

def plot_AvgIslandElev(AvgInteriorElevationTS):
    
    aE = [a * 10 for a in AvgInteriorElevationTS] # Convert to m
    plt.figure()
    plt.plot(aE)
    fig = plt.gcf()
    fig.set_size_inches(14 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Average Island Elevation (m)')
    plt.show()



 
#===================================================
# 11: Shoreface Slope Over Time

def plot_ShorefaceSlope(s_sf_TS):

    ssfTS = s_sf_TS
    plt.figure()
    plt.plot(ssfTS)
    fig = plt.gcf()
    fig.set_size_inches(14 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Shoreface Slope')
    plt.title('Shoreface Slope Over Time')
    plt.show()




#===================================================
# 12: Average Interior Width Over Time

def plot_AvgInteriorWidth(InteriorWidth_AvgTS):

    aiw = InteriorWidth_AvgTS
    plt.figure()
    plt.plot(aiw)
    fig = plt.gcf()
    fig.set_size_inches(14,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Average Interior Width (dam)')
    plt.title('Average Interior Width Over Time')
    plt.show()





#===================================================
# 13: Shoreface Overwash Flux Over Time

def plot_OverwashFlux(QowTS):

    plt.figure()
    plt.plot(QowTS)
    fig = plt.gcf()
    fig.set_size_inches(14,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Qow (m^3)')
    plt.title('Overwash Flux')
    #plt.ylim(0,210)
    plt.show()





#===================================================
# 14: Width, Berm Elevation, SF Slope, Shoreline Change, and Overwash Flux Over Time (all in one)

def plot_StatsSummary(s_sf_TS, x_s_TS, TMAX, InteriorWidth_AvgTS, QowTS, QsfTS, Hd_AverageTS):
    
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(14,18)
    plt.rcParams.update({'font.size':17})
    
    
    # Shoreface Slope
    plt.subplot(6,1,1)
    ssfTS = s_sf_TS
    plt.plot(ssfTS)
    plt.hlines(s_sf_eq,0,TMAX-1,colors='black',linestyles='dashed')
    
    plt.ylabel('Shoreface Slope')
    
    # Interior Width
    plt.subplot(6,1,2)
    aiw = [a * 10 for a in InteriorWidth_AvgTS]
    plt.plot(aiw)
    plt.ylabel('Avg. Width (m)') # Avergae Interior Width
    
    # Shoreline Change
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
    plt.subplot(6,1,3)
    plt.plot(scts)
    plt.ylabel('Shoreline Position (m)')
    
    # Overwash Flux
    plt.subplot(6,1,4)
    plt.plot(QowTS)
    plt.ylabel('Qow (m^3/m)')
    
    # Shoreface Flux
    plt.subplot(6,1,5)
    plt.plot(QsfTS)
    plt.ylabel('Qsf (m^3/m)')
    
    # Dune Height
    aHd = [a * 10 for a in Hd_AverageTS]
    plt.subplot(6,1,6)
    plt.plot(aHd)
    plt.xlabel('Year')
    plt.ylabel('Avg. Dune Height (m)') # Average Dune Height
    
    plt.show() 
    name = 'Output/Stats'
    fig.savefig(name)  




#===================================================
# 15: 3D Plot of Island Domain For Last Time Step

def plot_3DElevTMAX(TMAX, t, SL, DuneDomain, DomainTS):
    
    if TMAX > t:
       TMAX = t   
    # Build beach elevation domain
    BW = 6
    Beach = np.zeros([BW, BarrierLength])    
    berm = math.ceil(BW*0.65)
    Beach[berm:BW+1,:] = BermEl
    add = (BermEl-SL) / berm
    for i in range(berm):
        Beach[i,:] = SL + add * i 
    
    # Construct frame     
    Dunes = [(DuneDomain[TMAX] + BermEl) * 10] * DuneWidth
    Water = np.zeros([3,BarrierLength])
    Domain = DomainTS[TMAX] * 10
    Domain[Domain < 0] = 0
    Domain = np.vstack([Water, Beach, Dunes, Domain, Water])
    Dlen = np.shape(Domain)[1]
    Dwid = np.shape(Domain)[0]
    fig = plt.figure(figsize=(12,9))
    #fig.set_size_inches(12,7)
    ax = fig.add_subplot(111, projection='3d')
    scale_x = 1
    scale_y = Dwid / Dlen
    scale_z = 4 / Dlen * 3
    ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
    X = np.arange(Dlen)
    Y = np.arange(Dwid)
    X, Y = np.meshgrid(X, Y)
    Z = Domain
    ax.plot_surface(X, Y, Z, cmap='terrain', alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True)
    ax.set_zlim(0, 4)
    
    # Plot shrubs - Broken?
    # Shrubs = PercentCoverTS[TMAX]
    # Shrubs[Shrubs>0] = 1   
    # Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
    # Shrubs = Shrubs * Domain
    # Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1 
    # Shrubs[Shrubs<1] = None     
    # ax.scatter(X, Y+1, Shrubs, s=30, c='black')
    
    ax.view_init(10,155)
    plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) #mostly centered
    plt.show()
    name = 'Output/Domain3D'
    fig.savefig(name, dpi=200)




#===================================================
# 16: 3D Animation Frames of Island Elevation and Shrubs (no translation)

def plot_3DElevFrames(DomainTS, SL, TMAX, DuneDomain):
    
    for t in range(0,len(DomainTS)):
        # Build beach elevation domain
        BW = 6
        Beach = np.zeros([BW, BarrierLength])    
        berm = math.ceil(BW*0.65)
        Beach[berm:BW+1,:] = BermEl
        add = (BermEl-SL) / berm
        for i in range(berm):
            Beach[i,:] = SL + add * i 
        
        # Construct frame     
        Dunes = [(DuneDomain[TMAX] + BermEl) * 10] * DuneWidth
        Water = np.zeros([3,BarrierLength])
        Domain = DomainTS[TMAX] * 10
        Domain = np.vstack([Water, Beach, Dunes, Domain, Water])
        Dlen = np.shape(Domain)[1]
        Dwid = np.shape(Domain)[0]
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111, projection='3d')
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 4
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = Domain
        ax.plot_surface(X, Y, Z, cmap='terrain', alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True)
        ax.set_zlim(0, 4)
       
        # Plot shrubs - Broken?
        # Shrubs = PercentCoverTS[t]
        # Shrubs[Shrubs>0] = 1   
        # Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
        # Shrubs = Shrubs * Domain
        # Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1 
        # Shrubs[Shrubs<1] = None     
        # ax.scatter(X, Y+1, Shrubs, s=30, c='black')
       
        timestr = 'Time = ' + str(t) + ' yrs'
        ax.set_ylabel(timestr)
        ax.view_init(20,155)
        plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) #mostly centered
        plt.show()
        name = 'Output/SimFrames/3D_' + str(t)
        fig.savefig(name, dpi=150)
    
    


#===================================================
# 17: 3D Animation Frames of Island Evolution (with barrier translation)

def plot_3DElevAnimation(DomainTS, SL, TMAX, DuneDomain, DomainWidth, x_s_TS, ShorelineChange):
    
    BW = 6   
    AniDomainWidth = DomainWidth + round(BW) + 12 + abs(ShorelineChange)
    OriginY = 5
    for t in range(0,TMAX):
        
        # Build beach elevation domain
        BeachDomain = np.zeros([BW, BarrierLength])    
        berm = math.ceil(BW*0.65) 
        BeachDomain[berm:BW+1,:] = BermEl
        add = (BermEl-SL) / berm
        for i in range(berm):
            BeachDomain[i,:] = SL + add * i 
        
        # Make animation frame
        Domain = DomainTS[t] * 10
        Domain[Domain < 0] = 0
        Dunes = (DuneDomain[t,:,:] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10    
        Domain = np.vstack([Beach, Dunes, Domain]) 
        
        AnimateDomain = np.zeros([AniDomainWidth + 1, BarrierLength])
        widthTS = len(Domain)
        scts = [(x - x_s_TS[0]) for x in x_s_TS]
        if scts[t] >= 0:
            OriginTstart = OriginY + math.floor(scts[t])
        else:
            OriginTstart = OriginY + math.ceil(scts[t])        
        OriginTstop = OriginTstart + widthTS
        AnimateDomain[OriginTstart:OriginTstop, 0:BarrierLength] = Domain
        
        
        Dlen = np.shape(AnimateDomain)[1]
        Dwid = np.shape(AnimateDomain)[0]
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111, projection='3d')
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 1 #4 / Dlen * 4
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = AnimateDomain
        ax.plot_surface(X, Y, Z, cmap='terrain', alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True)
        ax.set_zlim(0, 4)
    
        timestr = 'Time = ' + str(t) + ' yrs'
        ax.set_ylabel(timestr, labelpad=50)
        ax.view_init(20,150) #ax.view_init(20,155)
        # plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        # plt.subplots_adjust(left=-0.7, right=1.3, top=2.2, bottom=-0.3) # mostly centered previous
        plt.subplots_adjust(left=-0.25, right=1.05, top=2.3, bottom=-0.2)
        # plt.show()
        name = 'Output/SimFrames/3DAni_' + str(t)
        fig.savefig(name, dpi=150)
        plt.close(fig)
        
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/3DAni_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/3DAni.gif', frames, 'GIF-FI')    
    
    
    
    
#===================================================
# 18: Shrub Age Domain at Simulation End

def plot_ShrubAgeTMAX(ShrubDomainAll, ShrubDeadTS):

    if Shrub_ON == 1:
        ShrubAll = ShrubDomainAll
        shrubFig1 = plt.figure(figsize=(14,5))
        ax = shrubFig1.add_subplot(111)
        cax = ax.matshow(ShrubAll, origin='lower', cmap='afmhot_r', vmin=0, vmax=10) # analysis:ignore
        ax.xaxis.set_ticks_position('bottom')
        cbar = shrubFig1.colorbar(cax)
        cbar.set_label('Shrub Age', rotation=270,labelpad=20)
        
        Dead = ShrubDeadTS[-1]
        Dy, Dx = np.where(Dead > 0)
        ax.scatter(Dx, Dy, marker='x', s=20, color='black', alpha=0.25)
        
        plt.xlabel('Alongshore Distance (dm)')
        plt.ylabel('Cross-Shore Diatance (dm)')
        plt.title('Final Shrub Age')
        plt.show()    
    
    
    
    
#===================================================
# 19: Percent Cover Domain at Simulation End

def plot_ShrubPercentCoverTMAX(PercentCoverTS, TMAX, DeadPercentCoverTS):
    
    if Shrub_ON == 1:
        ShrubPC = PercentCoverTS[TMAX-1]
        shrubFig2 = plt.figure(figsize=(14,5))
        ax = shrubFig2.add_subplot(111)
        cax = ax.matshow(ShrubPC, origin='lower', cmap='YlGn', vmin=0, vmax=1) # analysis:ignore
        ax.xaxis.set_ticks_position('bottom')
        cbar = shrubFig2.colorbar(cax)
        cbar.set_label('Shrub Percent Cover', rotation=270, labelpad=20)
        
        Dead = DeadPercentCoverTS[TMAX-1]
        Dy, Dx = np.argwhere(Dead > 0).T
        Dz = Dead[Dy,Dx]*20       
        ax.scatter(Dx, Dy, marker='x', s=Dz, color='k', alpha=0.25)
          
        plt.xlabel('Alongshore Distance (dm)')
        plt.ylabel('Cross-Shore Diatance (dm)')
        plt.title('Final Shrub Percent Cover')
        plt.show() 
        name = 'Output/PercentCover'
        shrubFig2.savefig(name)    
    
    
    
    
#===================================================
# 20: Shrub Area Over Time
    
def plot_ShrubArea(ShrubArea):
    
    if Shrub_ON == 1:
        area = ShrubArea
        plt.figure()
        plt.plot(area)
        fig = plt.gcf()
        fig.set_size_inches(14,3)
        plt.xlabel('Year')
        plt.ylabel('Shrub Area (dam^2)')
        plt.show()
        name = 'Output/ShrubArea'
        fig.savefig(name)    
    
    
    
    
#===================================================
# 21: Storm Count Over Time

def plot_StormCount(StormCount):

    plt.figure()
    plt.plot(StormCount)
    fig = plt.gcf()
    fig.set_size_inches(14 ,5)
    plt.xlabel('Year')
    plt.ylabel('Number of Storms')
    plt.title('Storm Count')
    plt.show()          
        
    

    
#===================================================
# 22: Alongshore Dune Height Over Time

def plot_AlongshoreDuneHeight(DuneDomain):
    
    Dunes = DuneDomain.max(axis=2)

    # Plot full tranects    
    plt.figure()
    for x in range(0, BarrierLength, 10):        
        Hd_TS = Dunes[:,x]
        Hd_TS = Hd_TS * 10 # Convert to meters
        plt.plot(Hd_TS)    
    fig = plt.gcf()
    fig.set_size_inches(14,6)
    plt.xlabel('Year')
    plt.ylabel('Dune Height (m)')
    plt.title('Dune Height Alongshore')
    plt.show()
    name = 'Output/Dunes_Alongshore'
    fig.savefig(name)        
     



#===================================================
# 23: Calculate Discontinuous Retreat

def calc_ShorelinePeriodicity(TMAX, x_s_TS):

    # Shoreline Change & Change Rate Over Time    
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
    
    # Filter
    win = 31
    poly = 3
    der1 = (signal.savgol_filter(scts, win, poly, deriv=1))
    
    HitDown = [] # Slow-downs
    HitUp = [] # Speed-ups
    
    window1 = 3 # Maximum allowed length for gaps in slow periods
    window2 = 30 # Minimum length required for slow periods, including gaps
    buffer = 3
    thresh1 = 0.5 # Max slope for slow periods
    thresh2 = 1
        
    # Find slow periods
    der_under = np.where(der1 < thresh1)[0] 
     
    if len(der_under) > 0:
                                        
        gaps = np.diff(der_under) > window1
        peak_start = np.insert(der_under[1:][gaps], 0, der_under[0])
        peak_stop = np.append(der_under[:-1][gaps], der_under[-1])
    
        # for n in range(len(peak_stop)):
        for n in range(1,len(peak_stop)):
            if peak_stop[n] - peak_start[n] > window2:
                if len(HitDown) == 0:
                    if peak_start[n] > buffer:
                        HitDown.append(peak_start[n])
                    if peak_stop[n] < len(scts) - buffer:
                        HitUp.append(peak_stop[n])
                else:
                    a = HitUp[-1]
                    b = peak_start[n]   
                    gap_length = b - a
                    gap_slope = np.average(der1[a:b])
                    # gap_slope = np.max(der1[a:b])
                    if gap_length < window2 and gap_slope < thresh2: 
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp[-1] = peak_stop[n]
                        else:
                            del HitUp[-1]
                    else:
                        if peak_start[n] > buffer:
                            HitDown.append(peak_start[n])
                        if peak_stop[n] < len(scts) - buffer:
                            HitUp.append(peak_stop[n])
                 
    ### CALCULATE STATS      
            
    Jumps = len(HitDown)
    Slows = len(HitUp)
    SlowDur = []
    FastDur = []
    
    if Jumps > 0 and Slows > 0:
        DownFirst = HitDown[0] < HitUp[0]
    elif Jumps == 0 and Slows > 1:
        DownFirst = True
    else:
        DownFirst = False
    
    if Jumps >= 2 or Slows >= 2:
        if Jumps >= 2 and Slows >= 2:
            Periodicity = (np.mean(np.diff(HitDown)) + np.mean(np.diff(HitUp))) / 2
        elif Jumps >= Slows:
            Periodicity = np.mean(np.diff(HitDown))  
        else:
            Periodicity = np.mean(np.diff(HitUp))
        if DownFirst:
            for n in range(Slows):            
                SlowDur.append(HitUp[n] - HitDown[n])
            for n in range(Jumps - 1):
                FastDur.append(HitDown[n+1] - HitUp[n])        
        else:
            for n in range(Slows-1):
                SlowDur.append(HitUp[n+1] - HitDown[n])
            for n in range(Jumps):            
                FastDur.append(HitDown[n] - HitUp[n])   
    else:
        Periodicity = 0
        if Jumps == 1 and Slows == 1:
            if DownFirst:
                SlowDur.append(HitUp[0] - HitDown[0])
            else:
                FastDur.append(HitDown[0] - HitUp[0])  
    
    AvgFastDur = np.mean(FastDur)
    if np.isnan(AvgFastDur):
        AvgFastDur = 0
    AvgSlowDur = np.mean(SlowDur)
    if np.isnan(AvgSlowDur):
        AvgSlowDur = 0
            
    if len(SlowDur) >= 2 and len(FastDur) >= 2:
        Punc = 1
    else:
        Punc = 0
        
    
    return Periodicity, AvgFastDur, AvgSlowDur, Punc
    



#===================================================
# 24: Average Dune Height and Storm TWL

def plot_DuneStorm(Hd_AverageTS, StormSeries, TMAX):
    
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(14,6)
    
    # Dune Height
    aHd = [a * 10 for a in Hd_AverageTS] # Use this for dune height
    # aHd = [(a + BermEl) * 10 for a in Hd_AverageTS] # Use this for dune elevation
    plt.plot(aHd, color='teal')
    plt.xlabel('Year')
    plt.ylabel('Avg. Dune Height (m)') # Use this for dune height
    # plt.ylabel('Avg. Dune Elevation (m)') # Use this for dune elevation
    
    # Storms
    if TMAX >= 1000:
        stop = len(StormSeries)
    else:
        stop = np.where(StormSeries[:,0] >= TMAX)[0][0]
    
    stormX = StormSeries[0:stop,0] # Yr
    stormY = StormSeries[0:stop,1] # Rhigh
    
    stormX = stormX[stormY > BermEl] # Remove storms with TWL lower than berm
    stormY = stormY[stormY > BermEl]
    stormY = [(a - BermEl) * 10 for a in stormY] # Height relative to berm # Use this for dune height
    # stormY = [(a) * 10 for a in stormY] # Use this for dune elevation
    
    plt.scatter(stormX, stormY, c='r', marker='*')
    
    plt.show() 
    name = 'Output/DuneStorm'
    fig.savefig(name)  

    


#===================================================
# 25: Seabed Profile

def plot_SeabedProfile(SL, TMAX, x_t_TS):

    Tx = []
    Ty = []
    for t in range(TMAX):          
        Tx.append(x_t_TS[t] * 10)
        Ty.append(((SL + np.sum(RSLR[0:t])) - DShoreface) *10)

    Sby = np.linspace(Ty[0], Ty[-1], num=TMAX) 
    Sbx = np.linspace(0, Tx[-1], num=TMAX)
    
    # Plot
    plt.figure(figsize=(14,5))
    plt.plot(Sbx, Sby,'slategray',linestyle='--')
    plt.plot(Tx,Ty, 'chocolate')
    plt.xlabel('Distance (m)')
    plt.ylabel('Elevation (m)')
    plt.title('Seabed Profile')
    plt.show()
    



#===================================================
# 26: Animation Frames of Shrub Percent Cover and Barrier Migration

def plot_ShrubAnimation(InteriorWidth_AvgTS, ShorelineChange, DomainTS, DuneDomain, SL, x_s_TS, Shrub_ON, PercentCoverTS, TMAX, DeadPercentCoverTS):
        
    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange) + OriginY + 35) # was +15
    
    for t in range(TMAX):
        # Build beach elevation domain
        BeachDomain = np.zeros([BeachWidth, BarrierLength])    
        berm = math.ceil(BeachWidth*0.65)
        BeachDomain[berm:BeachWidth+1,:] = BermEl
        add = (BermEl-SL) / berm
        for i in range(berm):
            BeachDomain[i,:] = SL + add * i 
    
        # Make animation frame domain
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t,:,:] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1
        Domain[Domain > 0] = 1
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength]) *-1
        widthTS = len(Domain)
        scts = [(x - x_s_TS[0]) for x in x_s_TS]
        if scts[t] >= 0:
            OriginTstart = OriginY + math.floor(scts[t])
        else:
            OriginTstart = OriginY + math.ceil(scts[t])        
        OriginTstop = OriginTstart + widthTS
        AnimateDomain[OriginTstart:OriginTstop, 0:BarrierLength] = Domain    
        if Shrub_ON == 1:
            Shrubs = PercentCoverTS[t]
            Dead = DeadPercentCoverTS[t]
            wid = np.zeros([BeachWidth + DuneWidth + OriginTstart, BarrierLength])
            Shrubs = np.vstack([wid, Shrubs])
            Dead = np.vstack([wid, Dead])
            Sy, Sx = np.argwhere(Shrubs > 0).T
            Sz = Shrubs[Sy,Sx] 
            Dy, Dx = np.argwhere(Dead > 0).T
            
        # Plot and save
        shrubfig = plt.figure(figsize=(10,12))
        ax = shrubfig.add_subplot(111)
        cax = ax.matshow(AnimateDomain, origin='lower', cmap='Blues_r')
        im = ax.scatter(Sx, Sy, marker='o', s=32, c=Sz, cmap='YlGn', vmin=0, vmax=1, alpha=1, edgecolors='none')  
        ax.scatter(Dx, Dy, marker='o', s=22, facecolors='none', edgecolors='maroon', alpha=0.4)
        ax.xaxis.set_ticks_position('bottom')
        cbar = shrubfig.colorbar(im, ax=ax)
        cbar.set_label('Shrub Percent Cover', rotation=270, labelpad=20)      
        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Cross-Shore Diatance (dam)')
        plt.title('Interior Elevation')
        plt.tight_layout()
        timestr = 'Time = ' + str(t) + ' yrs'
        newpath = 'Output/SimFrames/'
        if not os.path.exists(newpath):
            os.makedirs(newpath)            
        plt.text(1, 1, timestr)
        name = 'Output/SimFrames/shrubani_' + str(t)
        shrubfig.savefig(name) # dpi=200
        plt.close(shrubfig)
        
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/shrubani_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/shrubani.gif', frames, 'GIF-FI')
    print()
    print('[ * GIF successfully generated * ]')




#===================================================
# 27: 3D Animation of Island Evolution With Shoreface and Bay 

# WARNING 30Apr21: BROKEN, DOES NOT WORK

def plot_3DElevAnimation_Super(DomainTS, SL, TMAX, DuneDomain, DomainWidth, x_s_TS, x_t_TS, ShorelineChange):
    
    BW = 6   
    AniDomainWidth = int(DomainWidth + round(BW) + 12 + abs(ShorelineChange) + LShoreface)
    OriginY = 5
    for t in range(TMAX):
        
        # Build beach elevation domain
        BeachDomain = np.zeros([BW, BarrierLength])    
        berm = math.ceil(BW*0.65) 
        BeachDomain[berm:BW+1,:] = BermEl
        add = (BermEl-SL) / berm
        for i in range(berm):
            BeachDomain[i,:] = SL + add * i 
        
        # Build shoreface elevation domain
        l_sf = int(x_s_TS[t]-x_t_TS[t])
        SFDomain = np.zeros([l_sf, BarrierLength])    
        add = (DShoreface / l_sf)
        for i in range(l_sf):
            SFDomain[i,:] = -DShoreface + add * i
        
        # Make animation frame
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t,:,:] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        SF = SFDomain * 10
        Domain = np.vstack([SF, Beach, Dunes, Domain]) 
        
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength]) * BayDepth * 10
        widthTS = len(Domain)
        scts = [(x - x_s_TS[0]) for x in x_s_TS]
        # if scts[t] >= 0:
        #     OriginTstart = OriginY + math.floor(scts[t])
        # else:
        #     OriginTstart = OriginY + math.ceil(scts[t])  
        xtoe = [(x - x_t_TS[0]) for x in x_t_TS]
        if xtoe[t] >= 0:
            OriginTstart = OriginY + math.floor(xtoe[t])
        else:
            OriginTstart = OriginY + math.ceil(xtoe[t])        
        OriginTstop = OriginTstart + widthTS
        AnimateDomain[0,:] = DShoreface*10
        AnimateDomain[0:OriginTstart, :] = DShoreface * 10
        AnimateDomain[OriginTstop:-1, :] = BayDepth * 10
        AnimateDomain[OriginTstart:OriginTstop, :] = Domain

        
        
        Dlen = np.shape(AnimateDomain)[1]
        Dwid = np.shape(AnimateDomain)[0]
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111, projection='3d')
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 1 #4 / Dlen * 4
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = AnimateDomain
        ax.plot_surface(X, Y, Z, cmap='terrain', alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True)
        ax.set_zlim(0, 4)
    
        timestr = 'Time = ' + str(t) + ' yrs'
        ax.set_ylabel(timestr, labelpad=50)
        ax.view_init(20,150) #ax.view_init(20,155)
        # plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        # plt.subplots_adjust(left=-0.7, right=1.3, top=2.2, bottom=-0.3) # mostly centered previous
        # plt.subplots_adjust(left=-0.25, right=1.05, top=2.3, bottom=-0.2)
        # plt.show()
        name = 'Output/SimFrames/3DAni_' + str(t)
        fig.savefig(name, dpi=150)
        plt.close(fig)
        
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/3DAni_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/3DAni.gif', frames, 'GIF-FI')    
    




#===================================================
# 28: Shrub front location and/or expansion rate

def plot_ShrubFront(Shrub_ON, TMAX, PercentCoverTS):
    
# Defined as the left 95% of all shrub cells where age > 1

    if Shrub_ON:

        # Initialize
        percentile = 0.95
        ShrubFront = []
        FrontLoc = 0 
        
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
             
            ShrubFront.append(FrontLoc) # Store in array 
        
        # Plot front location
        ShrubFront = [i * 10 for i in ShrubFront]
        plt.figure()
        plt.plot(ShrubFront)
        fig = plt.gcf()
        fig.set_size_inches(14,7)
        plt.xlabel('Year')
        plt.ylabel('Shrub Front (m alongshore)')
        plt.show()
     
        # Plot Rate
        # der1 = (signal.savgol_filter(ShrubFront, 7, 3, deriv=1)) #7,3
        der1 = np.diff(ShrubFront)
        der1_rm = np.convolve(der1, np.ones((5,))/5, mode='same')   
        plt.figure()
        plt.plot(der1)
        plt.plot(der1_rm, c='r')
        fig = plt.gcf()
        fig.set_size_inches(14,7)
        plt.xlabel('Year')
        plt.ylabel('Expansion Rate (m/yr)')
        plt.show()   
            
            