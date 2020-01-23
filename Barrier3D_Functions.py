# Simulation and plotting functions for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 1
# Updated: 22 Jan 2020


# Simulation Functions Included: 
#   DuneGaps            Finds location and widths of OW throats in dune line
#   WaterLevels         Draws probabilistically determined Rhigh and Rlow from storm statistics
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
import matplotlib.pyplot as plt
from matplotlib import cm # analysis:ignore
from mpl_toolkits.mplot3d import Axes3D # analysis:ignore
import imageio # analysis:ignore




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
# WATER LEVELS 
    
# Returns statistics for a given storm
    
def WaterLevels(numstorm, MHW):   
    
    from Barrier3D_Parameters import (surge_tide_m, surge_tide_sd, height_mu, height_sigma, period_m, period_sd, beta, duration_mu, duration_sigma)
    
    ### Generate storm water levels
    Rhigh = []
    Rlow = []
    period = []
    duration = []
    
    for n in range(numstorm):
        ### Draw statistics for this storm from normal distributions
        
        surge_tide = np.random.normal(surge_tide_m, surge_tide_sd)
        if surge_tide < 0:
            surge_tide = 0.01
        
        height = np.random.lognormal(height_mu, height_sigma)
        if height < 0:
            height = 0.01
    
        T = np.random.normal(period_m, period_sd)        
        if T < 0:
            T = 0.01
        L0 = (9.8 * T**2) / (2 * math.pi) # Wavelength       
        if L0 < 0:
            L0 = 0.01

        dur = round(np.random.lognormal(duration_mu, duration_sigma))
        if dur < 8:
            dur = 8

        ### Calculate swash runup

        # Setup
        setup = 0.35 * beta * math.sqrt(height * L0) 
        
        # Incident band swash
        Sin = 0.75 * beta * math.sqrt(height * L0) 
        
        # Infragravity band swash
        Sig = 0.06 * math.sqrt(height * L0)
        
        # Swash
        swash = math.sqrt((Sin**2) + (Sig**2))
        
        # Rhigh (m NAVD88)
        R2 = 1.1 * (setup + (swash/2))
        Rh = surge_tide + R2
        
        # Rlow (m NAVD88)
#        Rl = Rh - (swash/2) # Which calculation is best?
#        Rl = Rh - swash
        Rl = surge_tide

        ### Convert to decameters relative to MHW of 0
        Rh = (Rh / 10) - MHW
        Rl = (Rl / 10) - MHW
        
        # Save to lists
        Rhigh.append(Rh)
        Rlow.append(Rl)
        period.append(T)
        duration.append(dur)
                
    return Rhigh, Rlow, period, duration




    
#===================================================
# SeaLevel
    
# Accounts for relative sea level rise, returns updated elevation domains

def SeaLevel(InteriorDomain, DuneDomain, t):
    
    from Barrier3D_Parameters import (RSLR, BayDepth)
    
    # Decrease all elevation this year by RSLR increment
    InteriorDomain = InteriorDomain - RSLR
    DuneDomain[t-1] = DuneDomain[t-1] - RSLR    
    InteriorDomain[InteriorDomain < -BayDepth] = -BayDepth # Bay can't be deeper than BayDepth (assumes equilibrium depth maintained)
        
    return InteriorDomain, DuneDomain





#===================================================
# DuneGrowth
    
# Grows dune cells

def DuneGrowth(DuneDomain, t):
    
    from Barrier3D_Parameters import (Dmaxel, BermEl, growthparam, BarrierLength, DuneWidth)
    
    # Set max dune height - depends on beach width according to relationship gathered from VCR data
    Dmax = Dmaxel - BermEl # (dam) Max dune height 
    if Dmax < 0:
        Dmax = 0  
    
    Qdg = 0    
    # Grow dune
    for q in range(DuneWidth):
        G = growthparam * DuneDomain[t-1,:,q] * (1 - DuneDomain[t-1,:,q] / Dmax)
        DuneDomain[t,:,q] = G + DuneDomain[t-1,:,q]
        Qdg = Qdg + sum(G[0,:]) / BarrierLength # Volume of sediment lost from beach/shoreface from dune growth
   
    return DuneDomain, Dmax, Qdg






#===================================================
# DiffuseDunes
    
# Dune height diffusion in alongshore direction: smoothes dune height by redistributing sand from high dune to neighboring low dune(s)
    
def DiffuseDunes(DuneDomain, t):
    
    from Barrier3D_Parameters import (HdDiffu, BarrierLength, DuneWidth)
    
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
    
    from Barrier3D_Parameters import (BarrierLength)
    
    DomainWidth = np.shape(InteriorDomain)[0] # (dam) analysis:ignore
    InteriorWidth = [0] * BarrierLength
    for l in range(BarrierLength):
        width = next((index for index, value in enumerate(InteriorDomain[:,l]) if value <= SL), DomainWidth)
        if width < 0: width = 0
        InteriorWidth[l] = width -1
    # Average width of island
    InteriorWidth_Avg = np.nanmean(InteriorWidth)

#    for l in range(BarrierLength):
#        if InteriorDomain[0,l] < SL:
#            front = next((index for index, value in enumerate(InteriorDomain[:,l]) if value >= SL), DomainWidth)
#            back = next((index for index, value in enumerate(InteriorDomain[front:,l]) if value < SL), DomainWidth)
#            width = back - front
#            if width < 0: width = 0
#            InteriorWidth[l] = width
#        else:
#            width = next((index for index, value in enumerate(InteriorDomain[:,l]) if value <= SL), DomainWidth)          
#            if width < 0: width = 0
#            InteriorWidth[l] = width

    return DomainWidth, InteriorWidth, InteriorWidth_Avg





#===================================================
# LTA_SCR

# Finds shoreline change for modeled year following Lorenzo-Trueba and Ashton (2014)

def LTA_SC(InteriorDomain, OWloss, x_s, x_s_TS, x_t, x_t_TS, x_b_TS, s_sf_TS, InteriorWidth_Avg, SL, QowTS, QsfTS):

    from Barrier3D_Parameters import (BarrierLength, DShoreface, k_sf, s_sf_eq, RSLR, Qat)
    
    # Find volume of shoreface/beach sand deposited on island
    Qow = (OWloss) / (BarrierLength) # (dam^3/dam) Volume of sediment lost from shoreface/beach by overwash
    QowTS.append(Qow * 100) # Save in m^3/m    
    
    # DefineParams
    d_sf = DShoreface 
    h_b = np.average(InteriorDomain[InteriorDomain >= SL]) 
    x_b = x_s + InteriorWidth_Avg
        
    # Shoreface Flux
    s_sf = d_sf / (x_s - x_t)
    Qsf = (k_sf / 100) * (s_sf_eq - s_sf) # Convert Ksf from m^3/m/yr to dam^3/dam/yr
    QsfTS.append(Qsf * 100) # Save in m^3/m

    # Toe, Shoreline, and island base elevation changes     
    x_t_dt = (4 * Qsf * (h_b + d_sf) / (d_sf * (2 * h_b + d_sf))) + (2 * RSLR / s_sf)   
    x_s_dt = 2 * (Qow + Qat) / ((2 * h_b) + d_sf) - (4 * Qsf * (h_b + d_sf) / (((2 * h_b) + d_sf)**2)) # Dune growth and alongshore transport added to LTA14 formulation
       
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
    
def Shrubs(InteriorDomain, DuneDomainCrest, t, ShrubDomainFemale, ShrubDomainMale, BurialDomain, InteriorWidth_Avg, DomainWidth):

    from Barrier3D_Parameters import (BarrierLength, Dshrub, BermEl, Female, ShrubEl_min, ShrubEl_max, BurialLimit, UprootLimit, TimeFruit, Seedmin, Seedmax, GermRate, disp_mu, disp_sigma)
    
    ### Age the shrubs in years
    ShrubDomainFemale[ShrubDomainFemale > 0] += 1            
    ShrubDomainMale[ShrubDomainMale > 0] += 1

    ### Randomly disperse a seed onto the island each time step
    randX = np.random.randint(0,BarrierLength)
    randY = np.random.randint(0,InteriorWidth_Avg)
    if  ShrubDomainFemale[randY,randX] == 0 and ShrubDomainMale[randY,randX] == 0 and DuneDomainCrest[randX] + BermEl >= Dshrub and InteriorDomain[randY,randX] >= ShrubEl_min and InteriorDomain[randY,randX] <= ShrubEl_max:
        if random.random() > Female:
            ShrubDomainFemale[randY, randX] = 1
        else:
            ShrubDomainMale[randY, randX] = 1
    
    ### Inundation
    # Remove shrubs that have fallen below SL (passively inundated by rising back-barrier water elevations)
    ShrubDomainFemale[InteriorDomain < 0] = 0
    ShrubDomainMale[InteriorDomain < 0] = 0
    
    ### Burial
    # Kill all shrubs buried by over depth limit or eroded 
    ShrubDomainFemale[BurialDomain > BurialLimit] = 0
    ShrubDomainFemale[BurialDomain < UprootLimit] = 0
    ShrubDomainMale[BurialDomain > BurialLimit] = 0
    ShrubDomainMale[BurialDomain < UprootLimit] = 0
    BurialDomain[BurialDomain > BurialLimit] = 0 # Reset
    BurialDomain[BurialDomain < UprootLimit] = 0 # Reset
    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
    
#        ### Salt spray
#        # Kill first row of shrubs exposed to ocean (i.e. dune below threhold height) from last time step
#        for j in range(BarrierLength):
#            shrubs = np.nonzero(ShrubDomainAll[:,j])[0] # Gather indices of all shrubs in column
#            if len(shrubs) != 0: # If column not empty
#                first = shrubs[0] # Get index of fronting shrub cell
#                if DuneDomainCrest[j] + BermEl < Dshrub and first < np.shape(ShrubDomainAll)[1] - 1: # Crashing here (out of index)
#                    ShrubDomainFemale[first,j] = 0
#                    ShrubDomainMale[first,j] = 0
#        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale

    ### Disperse seeds
    for k in range(BarrierLength): # Loop through each row of island width (i.e. from ocean to mainland side of island)
        if 0 in ShrubDomainAll:                  #<--------------- Working??            
          
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
                            #   -there is no shrub at the receiving cell (male or female),
                            #   -the receiving cell has a tall enough fronting dune
                            #   -the receiving cell is within elevation range
                            if  targetY >= 0 and targetY < DomainWidth and targetX >= 0 and targetX < BarrierLength and ShrubDomainFemale[targetY,targetX] == 0 and ShrubDomainMale[targetY,targetX] == 0 \
                                and DuneDomainCrest[targetX] + BermEl >= Dshrub and InteriorDomain[targetY,targetX] >= ShrubEl_min and InteriorDomain[targetY,targetX] <= ShrubEl_max:
                                # Decide if the tree wll be a female or male
                                if random.random() > Female:
                                    ShrubDomainFemale[targetY,targetX] = 1
                                else:
                                    ShrubDomainMale[targetY,targetX] = 1                                    
    
    ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale

    
    return ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, BurialDomain

    



#===================================================
# UpdateShrubDomains

# Updates size of shrub domains
    
def UpdateShrubDomains(DomainWidth, ShrubDomainWidth, ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain):
    
    from Barrier3D_Parameters import (BarrierLength)


    if DomainWidth > ShrubDomainWidth:
        AddRows = np.zeros([DomainWidth-ShrubDomainWidth,BarrierLength])
        ShrubDomainFemale = np.vstack([ShrubDomainFemale, AddRows])
        ShrubDomainMale = np.vstack([ShrubDomainMale, AddRows])
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
        ShrubPercentCover = np.vstack([ShrubPercentCover, AddRows])   
        BurialDomain = np.vstack([BurialDomain, AddRows])
    elif DomainWidth < ShrubDomainWidth:
        RemoveRows = ShrubDomainWidth - DomainWidth
        ShrubDomainFemale = ShrubDomainFemale[0:-RemoveRows,:]
        ShrubDomainMale = ShrubDomainMale[0:-RemoveRows,:]
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
        ShrubPercentCover = ShrubPercentCover[0:-RemoveRows,:]   
        BurialDomain = BurialDomain[0:-RemoveRows,:] 
        
    return ShrubDomainFemale, ShrubDomainMale, ShrubDomainAll, ShrubPercentCover, BurialDomain





#===================================================
# SalineFlooding

# Kill all immature (< 1 yr-old) shrubs that have been flooded beyond a threshold discharge (Tolliver et al., 1997)
    
def SalineFlooding(ShrubDomainWidth, ShrubDomainAll, ShrubDomainFemale, ShrubDomainMale, d, i, Q0):
                   
    from Barrier3D_Parameters import (BarrierLength, SalineLimit)
    
    if d < (ShrubDomainWidth-1) and i < (BarrierLength - 1):
        if Q0 >= SalineLimit and ShrubDomainAll[d,i] == 1:
            ShrubDomainFemale[d,i] = 0
            ShrubDomainMale[d,i] = 0
            #ShrubDomainMale[d,i] = 0    #???               
                   
    return ShrubDomainFemale, ShrubDomainMale                

                   
                   

                   
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
    
def CalcPC(ShrubDomainAll, PercentCoverTS, ShrubArea, t):

    from Barrier3D_Parameters import PC
    
    Allshrub_t = ShrubDomainAll.astype('int64')
    ShrubPercentCover = PC.take(Allshrub_t) 
    PercentCoverTS[t] = ShrubPercentCover
    ShrubArea.append(np.count_nonzero(ShrubDomainAll))

    return ShrubPercentCover, PercentCoverTS, ShrubArea

    


#==============================================================================================================================================
# PLOTTING FUNCTIONS      
#==============================================================================================================================================
    

#===================================================
# 1: Dune Height Over Time

def plot_DuneHeight(DuneDomain):
    DuneCrest = DuneDomain.max(axis=2)
    duneFig = plt.figure(figsize=(15,8))
    plt.rcParams.update({'font.size':13})
    ax = duneFig.add_subplot(111)
    ax.matshow((DuneCrest)*10, origin='lower', cmap='bwr', aspect='auto')#, vmin=0, vmax=Dmax)
    cax = ax.xaxis.set_ticks_position('bottom') # analysis:ignore
    #cbar = duneFig.colorbar(cax)
    #cbar.set_label('Dune Height Above Berm Elevation (m)', rotation=270)
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Time (yr)')
    plt.title('Dune Height (m)')
    name = 'Output/Dunes'
    duneFig.savefig(name)


    

#===================================================
# 2: Elevation Domain For Last Time Step

def plot_ElevTMAX(TMAX, t, DuneDomain, DomainTS):
    
    from Barrier3D_Parameters import (BermEl)
    
    if TMAX > t:
        TMAX = t
    Dunes = (DuneDomain[TMAX,:,:] + BermEl) * 10
    Dunes = np.rot90(Dunes)
    Dunes = np.flipud(Dunes)
    Domain = DomainTS[TMAX] * 10
    Domain = np.vstack([Dunes, Domain])
    elevFig1 = plt.figure(figsize=(15,5))
    ax = elevFig1.add_subplot(111)
    cax = ax.matshow(Domain, origin='lower', cmap='terrain', vmin=-1.1, vmax=4.0)#, interpolation='gaussian') # analysis:ignore
    ax.xaxis.set_ticks_position('bottom')
    #cbar = elevFig1.colorbar(cax)
    #cbar.set_label('Elevation (m)', rotation=270)
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Cross-Shore Diatance (dam)')
    plt.title('Interior Elevation (m)')
    timestr = 'Time = ' + str(TMAX) + ' yrs'
    plt.text(1, 1, timestr)
    name = 'Output/FinalElevation'
    elevFig1.savefig(name)




#===================================================
# 3: Elevation Domain Frames

def plot_ElevFrames(TMAX, DomainTS):

    for t in range(TMAX):
        elevFig1 = plt.figure(figsize=(15,5))
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
        #plt.show()
        name = 'Output/SimFrames/elev_' + str(t)
        elevFig1.savefig(name)




#===================================================
# 4: Animation Frames of Barrier and Dune Elevation

def plot_ElevAnimation(InteriorWidth_AvgTS, ShorelineChange, DomainTS, DuneDomain, SL, x_s_TS, Shrub_ON, PercentCoverTS, TMAX):

    from Barrier3D_Parameters import (BarrierLength, BermEl, DuneWidth)
    
    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange) + OriginY + 15)
    
#    for t in range(0,len(DomainTS)):
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
            wid = np.zeros([BeachWidth + DuneWidth + OriginTstart, BarrierLength])
            Shrubs = np.vstack([wid, Shrubs])       
            Sy, Sx = np.argwhere(Shrubs > 0).T
        
        # Plot and save
        elevFig1 = plt.figure(figsize=(15,12))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(AnimateDomain, origin='lower', cmap='terrain', vmin=-1.1, vmax=4.0)#, interpolation='gaussian') # analysis:ignore
        if Shrub_ON == 1:
            ax.scatter(Sx, Sy, marker='o', s=12, c='black', alpha=0.35, edgecolors='none')
        ax.xaxis.set_ticks_position('bottom')
        #cbar = elevFig1.colorbar(cax)
        plt.xlabel('Alongshore Distance (dam)')
        plt.ylabel('Cross-Shore Diatance (dam)')
        plt.title('Interior Elevation')
        timestr = 'Time = ' + str(t) + ' yrs'
        plt.text(1, 1, timestr)
        name = 'Output/SimFrames/elev_' + str(t)
        elevFig1.savefig(name) # dpi=200
    
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/elev_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/elev.gif', frames, 'GIF-FI')
        
      
    
    
#===================================================        
# 5: Cross-Shore Transect Every 100 m Alongshore For Last Time Step
  
def plot_XShoreTransects(InteriorDomain, DuneDomain, SL, TMAX):
    
    from Barrier3D_Parameters import (BarrierLength, BermEl)

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
    for v in range(0, BarrierLength, 20):
        CrossElev = InteriorDomain[:, v]
        Dunes = DuneDomain[TMAX-1, v, :] + BermEl
        CrossElev1 = np.insert(CrossElev,0,Dunes)
        CrossElev2 = np.insert(CrossElev1,0,BeachX)
        CrossElev = CrossElev2 * 10 # Convert to meters
        plt.plot(CrossElev)
    fig = plt.gcf()
    fig.set_size_inches(15,6)
    plt.hlines(SL,-1,len(CrossElev+1),colors='dodgerblue')
    plt.xlabel('Cross-Shore Distance (dam)')
    plt.ylabel('Elevation (m)')
    plt.title('Cross-shore Topo Transects')
    plt.show()
    name = 'Output/Profiles'
    fig.savefig(name)       
    
   
    
    
#===================================================    
# 6: Shoreline Change Per Year Over Time   

def plot_ShorelineChange(ShorelineChangeTS, BBShorelineChangeTS):
    
    shorelinesum = 0
    BBsum = 0
    scts = []
    bbscts = []
    for i in ShorelineChangeTS:
        shorelinesum += i
        scts.append(shorelinesum)
    for q in BBShorelineChangeTS:
        BBsum += q
        bbscts.append(BBsum)
    plt.figure()
    plt.plot(bbscts)
    plt.plot(scts)
    fig = plt.gcf()
    fig.set_size_inches(25 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Shoreline Change (dam)')
    plt.title('Shoreline Change')
    plt.show()   
    
   
    
    
#===================================================    
# 7: Shoreline Change Rate Over Time

def plot_ShorelineChangeRate(ShorelineChangeTS):

    scRts = ShorelineChangeTS
    plt.figure()
    plt.plot(scRts)
    fig = plt.gcf()
    fig.set_size_inches(25 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Shoreline Change Rate(dam/yr)')
    plt.title('Shoreline Change Rate Over Time')
    plt.show() 
    
    
    
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

    from Barrier3D_Parameters import (RSLR, BermEl, DShoreface, BayDepth)
    
    ymin = -11
    ymax = SL + RSLR * TMAX + 4
    xmin = -5
    xmax = x_b_TS[TMAX-1] + 20
    
    SFfig = plt.figure(figsize=(20,5))
    axes = plt.gca()
    colors = plt.cm.jet(np.linspace(0,1,TMAX))
    
    for t in range(0,TMAX,5): # Plots one transect every 5 years
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
        
    axes.set_ylim(ymin,ymax) # Doesn't work - why not?
    axes.set_xlim(xmin,xmax)    
    plt.axis('equal')
    plt.xlabel('Alongshore Distance (dam)')
    plt.ylabel('Elevation (m)')
    plt.title('Shoreface Evolution')
    plt.show()
    
    # Save   
    name = 'Output/Shoreface'
    SFfig.savefig(name)
    
    
    
    
#===================================================
# Average Island Elevation Over Time

def plot_AvgIslandElev(h_b_TS):
    
    beTS = h_b_TS
    plt.figure()
    plt.plot(beTS)
    fig = plt.gcf()
    fig.set_size_inches(20 ,5)
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
    fig.set_size_inches(20 ,5)
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
    fig.set_size_inches(20 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Average Interior Width (dam)')
    plt.title('Average Interior Width Over Time')
    plt.show()





#===================================================
# 13: Shoreface Overwash Flux Over Time

def plot_OverwashFlux(QowTS):

    qow = [i * 1000 for i in QowTS]
    plt.figure()
    plt.plot(qow)
    fig = plt.gcf()
    fig.set_size_inches(20 ,5)
    plt.xlabel('Time (yrs)')
    plt.ylabel('Qow (m^3)')
    plt.title('Shoreface Overwash Flux')
    plt.show()





#===================================================
# 14: Width, Berm Elevation, SF Slope, Shoreline Change, and Overwash Flux Over Time (all in one)

def plot_StatsSummary(s_sf_TS, x_s_TS, TMAX, InteriorWidth_AvgTS, QowTS, QsfTS, Hd_AverageTS):

    from Barrier3D_Parameters import (s_sf_eq)
    
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(14,20)
    plt.rcParams.update({'font.size':11})
    
    
    # Shoreface Slope
    plt.subplot(6,1,1)
    ssfTS = s_sf_TS
    plt.plot(ssfTS)
    plt.hlines(s_sf_eq,0,TMAX,colors='black',linestyles='dashed')
    
    plt.ylabel('Shoreface Slope')
    
    # Interior Width
    plt.subplot(6,1,2)
    aiw = [a * 10 for a in InteriorWidth_AvgTS]
    plt.plot(aiw)
    plt.ylabel('Average Interior Width (m)')
    
    # Shoreline Change
    scts = [(x - x_s_TS[0]) * -10 for x in x_s_TS]
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
    plt.xlabel('Time (yrs)')
    plt.ylabel('Average Dune Height (m)')
    
    plt.show() 
    name = 'Output/Stats'
    fig.savefig(name)  




#===================================================
# 15: 3D Plot of Island Domain For Last Time Step

def plot_3DElevTMAX(TMAX, t, SL, DuneDomain, DomainTS):

    from Barrier3D_Parameters import (BarrierLength, BermEl, DuneWidth)
    
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
    
    # Plot shrubs
    #Shrubs = PercentCoverTS[TMAX]
    #Shrubs[Shrubs>0] = 1   
    #Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
    #Shrubs = Shrubs * Domain
    #Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1 
    #Shrubs[Shrubs<1] = None     
    #ax.scatter(X, Y+1, Shrubs, s=30, c='black')
    
    ax.view_init(10,155)
    plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) #mostly centered
    plt.show()
    name = 'Output/Domain3D'
    fig.savefig(name, dpi=200)




#===================================================
# 16: 3D Animation Frames of Island Elevation and Shrubs (no translation)

def plot_3DElevFrames(DomainTS, SL, TMAX, DuneDomain):

    from Barrier3D_Parameters import (BarrierLength, BermEl, DuneWidth)
    
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
       
        # Plot shrubs
    #    Shrubs = PercentCoverTS[t]
    #    Shrubs[Shrubs>0] = 1   
    #    Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
    #    Shrubs = Shrubs * Domain
    #    Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1 
    #    Shrubs[Shrubs<1] = None     
    #    ax.scatter(X, Y+1, Shrubs, s=30, c='black')
       
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

    from Barrier3D_Parameters import (BarrierLength, BermEl, DuneWidth)
    
    BW = 6   
    AniDomainWidth = DomainWidth + round(BW) + 12 + abs(ShorelineChange)
    OriginY = 5
    for t in range(0,len(DomainTS)):
        
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
        Dunes = [(DuneDomain[t] + BermEl) * 10] * DuneWidth
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
        scale_z = 4 / Dlen * 4
        ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = AnimateDomain
        ax.plot_surface(X, Y, Z, cmap='terrain', alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True)
        ax.set_zlim(0, 4)
    
        timestr = 'Time = ' + str(t) + ' yrs'
        ax.set_ylabel(timestr)
        ax.view_init(20,155)
    #    plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        plt.subplots_adjust(left=-0.7, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        plt.show()
        name = 'Output/SimFrames/3DAni_' + str(t)
        fig.savefig(name, dpi=150)
        
        
    frames = []
    for filenum in range(TMAX):
        filename = 'Output/SimFrames/3DAni_' + str(filenum) + '.png'
        frames.append(imageio.imread(filename))
    imageio.mimsave('Output/SimFrames/3DAni.gif', frames, 'GIF-FI')    
    
    
    
#===================================================
# 18: Shrub Age Domain at Simulation End

def plot_ShrubAgeTMAX(ShrubDomainAll):
    
    from Barrier3D_Parameters import Shrub_ON

    if Shrub_ON == 1:
        ShrubAll = ShrubDomainAll
        shrubFig1 = plt.figure(figsize=(40,5))
        ax = shrubFig1.add_subplot(111)
        cax = ax.matshow(ShrubAll, origin='lower', cmap='RdYlGn', vmin=0, vmax=10) # analysis:ignore
        ax.xaxis.set_ticks_position('bottom')
        #cbar = shrubFig1.colorbar(cax)
        #cbar.set_label('Shrub Age', rotation=270)
        plt.xlabel('Alongshore Distance (dm)')
        plt.ylabel('Cross-Shore Diatance (dm)')
        plt.title('Final Shrub Age')
        plt.show()    
    
    
    
    
#===================================================
# 19: Percent Cover Domain at Simulation End

def plot_ShrubPercentCoverTMAX(PercentCoverTS, TMAX):

    from Barrier3D_Parameters import Shrub_ON
    
    if Shrub_ON == 1:
        ShrubPC = PercentCoverTS[TMAX-1]
        shrubFig2 = plt.figure(figsize=(15,5))
        ax = shrubFig2.add_subplot(111)
        cax = ax.matshow(ShrubPC, origin='lower', cmap='YlGn', vmin=0, vmax=1) # analysis:ignore
        ax.xaxis.set_ticks_position('bottom')
        #cbar = shrubFig2.colorbar(cax)
        #cbar.set_label('Shrub Percent Cover', rotation=270)
        plt.xlabel('Alongshore Distance (dm)')
        plt.ylabel('Cross-Shore Diatance (dm)')
        plt.title('Final Shrub Percent Cover')
        plt.show() 
        name = 'Output/PercentCover'
        shrubFig2.savefig(name)    
    
    
    
    
#===================================================
# 20: Shrub Area Over Time
    
def plot_ShrubArea(ShrubArea):    
    
    from Barrier3D_Parameters import Shrub_ON
    
    if Shrub_ON == 1:
        area = ShrubArea
        plt.figure()
        plt.plot(area)
        fig = plt.gcf()
        fig.set_size_inches(15,3)
        plt.xlabel('Time (yrs)')
        plt.ylabel('Shrub Area (dam^2)')
        plt.title('Shrub Coverage Over Time')
        plt.show()
        name = 'Output/ShrubArea'
        fig.savefig(name)    
    
    
    
#===================================================
# 21: Storm count over time

def plot_StormCount(StormCount):

    plt.figure()
    plt.plot(StormCount)
    fig = plt.gcf()
    fig.set_size_inches(12 ,5)
    plt.xlabel('Year')
    plt.ylabel('Number of Storms')
    plt.title('Storm Count')
    plt.show()          
        
        
        
        