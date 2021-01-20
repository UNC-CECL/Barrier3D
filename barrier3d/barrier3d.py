import math
# import random   # KA: mcflugen added a seeded number generator so we can reproduce runs for testing

import numpy as np

# from scipy import signal

# from .parameters import Barrier3dParameters
from .load_input import load_inputs


# KA: called within Barrier3d - does this force an exit from the program? replaces former "warnings" package?
class Barrier3dError(Exception):
    pass


class Barrier3d:
    def SeaLevel(self, InteriorDomain, DuneDomain, t):
        """Accounts for relative sea level rise, returns updated elevation domains."""

        # Decrease all elevation this year by RSLR increment
        InteriorDomain = InteriorDomain - self._RSLR[t]
        DuneDomain[t - 1] = DuneDomain[t - 1] - self._RSLR[t]

        # Bay can't be deeper than BayDepth (roughly equivalent to constant back-barrier slope)
        InteriorDomain[InteriorDomain < -self._BayDepth] = -self._BayDepth

        return InteriorDomain, DuneDomain

    def FindWidths(self, InteriorDomain, SL):
        """Finds DomainWidth and InteriorWidth

        DomainWidth is wide enough to fit widest part of island; InteriorWidth
        stores the width of the island at each row of the domain where elevation
        is greater than 0
        """

        DomainWidth = np.shape(InteriorDomain)[0]  # (dam) analysis:ignore
        InteriorWidth = [0] * self._BarrierLength
        for bl in range(self._BarrierLength):
            width = next(
                (
                    index
                    for index, value in enumerate(InteriorDomain[:, bl])
                    if value <= SL
                ),
                DomainWidth,
            )
            width = width - 1
            if width < 0:
                width = 0
            InteriorWidth[bl] = width
        # Average width of island
        InteriorWidth_Avg = np.nanmean(InteriorWidth)

        return DomainWidth, InteriorWidth, InteriorWidth_Avg

    def DuneGrowth(self, DuneDomain, t):
        """Grows dune cells"""

        # Set max dune height
        Dmax = self._Dmaxel - self._BermEl  # (dam) Max dune height
        if Dmax < 0:
            Dmax = 0

        Cf = 3  # Decay coefficient
        Qdg = 0
        # Grow dune
        for q in range(self._DuneWidth):
            reduc = 1 / (Cf ** q)
            G = (
                    self._growthparam
                    * DuneDomain[t - 1, :, q]
                    * (1 - DuneDomain[t - 1, :, q] / Dmax)
                    * reduc
            )
            DuneDomain[t, :, q] = G + DuneDomain[t - 1, :, q]
            Qdg = Qdg + (
                    np.sum(G) / self._BarrierLength
            )  # Volume of sediment lost from beach/shoreface from dune growth

        return DuneDomain, Dmax, Qdg

    def Shrubs(
            self,
            InteriorDomain,
            DuneDomainCrest,
            ShrubDomainFemale,
            ShrubDomainMale,
            BurialDomain,
            InteriorWidth_Avg,
            DomainWidth,
            ShrubDomainDead,
    ):
        """Main shrub expansion and mortality function"""

        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale

        BeachWidth = 60 / 10  # Temp

        # ### Burial
        # Find heights of shrub plants
        ShrubHeight = ShrubDomainAll * self._MaxShrubHeight

        # Kill all shrubs buried by over depth limit or eroded
        ShrubDomainDead[BurialDomain > (self._BurialLimit * ShrubHeight)] = ShrubDomainAll[
            BurialDomain > (self._BurialLimit * ShrubHeight)]  # Transfer to dead domain
        ShrubDomainDead[
            BurialDomain >= 0.3] = 0  # If buried by more than 3 m, considered shrub to be gone(dead or alive)
        ShrubDomainDead[BurialDomain < self._UprootLimit] = ShrubDomainAll[
            BurialDomain < self._UprootLimit]  # Transfer to dead domain

        ShrubDomainFemale[BurialDomain > (self._BurialLimit * ShrubHeight)] = 0
        ShrubDomainFemale[BurialDomain < self._UprootLimit] = 0
        ShrubDomainMale[BurialDomain > (self._BurialLimit * ShrubHeight)] = 0
        ShrubDomainMale[BurialDomain < self._UprootLimit] = 0

        BurialDomain[BurialDomain > (self._BurialLimit * ShrubHeight)] = 0  # Reset burial domain
        BurialDomain[BurialDomain < self._UprootLimit] = 0  # Reset burial domain
        BurialDomain[BurialDomain >= 0.3] = 0  # Reset burial domain
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale  # Recalculate

        # ### Inundation
        # Remove shrubs that have fallen below Mean High Water (i.e. elevation of 0)
        # (passively inundated by rising back-barrier water elevations)

        for w in range(DomainWidth):
            for bl in range(self._BarrierLength):
                if InteriorDomain[w, bl] < 0 and ShrubDomainAll[w, bl] > 0:
                    if InteriorDomain[w, bl] > -self._TideAmp:
                        ShrubDomainDead[w, bl] = ShrubDomainAll[
                            w, bl
                        ]  # Shrub remains as dead if within tidal range (i.e. MHW - tidal range)
                    ShrubDomainFemale[w, bl] = 0
                    ShrubDomainMale[w, bl] = 0
                elif InteriorDomain[w, bl] < 0 and ShrubDomainDead[w, bl] > 0 \
                        and InteriorDomain[w, bl] < -self._TideAmp:  # Remove dead shrubs below tidal range
                    ShrubDomainDead[w, bl] = 0

        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale  # Recalculate

        # ### Age the shrubs in years
        ShrubDomainFemale[ShrubDomainFemale > 0] += 1
        ShrubDomainMale[ShrubDomainMale > 0] += 1

        # ### Randomly drop a seed onto the island each time step
        # randX = np.random.randint(0, self._BarrierLength)
        # randY = np.random.randint(0, InteriorWidth_Avg)
        randX = self._RNG.integers(0, self._BarrierLength)
        randY = self._RNG.integers(0, max(1, InteriorWidth_Avg))
        if (
                ShrubDomainFemale[randY, randX] == 0
                and ShrubDomainMale[randY, randX] == 0
                and DuneDomainCrest[randX] + self._BermEl >= self._Dshrub
                and self._ShrubEl_min <= InteriorDomain[randY, randX] <= self._ShrubEl_max
        ):
            # if random.random() > self._Female:
            if self._RNG.uniform() > self._Female:
                ShrubDomainFemale[randY, randX] = 1
            else:
                ShrubDomainMale[randY, randX] = 1

        # ### Disperse seeds
        for k in range(
                self._BarrierLength
        ):  # Loop through each row of island width (i.e. from ocean to mainland side of island)
            if 0 in ShrubDomainAll:

                # For all cells with a shrub
                FemaleShrubs = ShrubDomainFemale[:, k]
                I = [
                    index
                    for index, value in enumerate(FemaleShrubs)
                    if value >= self._TimeFruit
                ]
                numShrubCells = len(I)
                # Determine how many seeds in each cell
                # Seedspercell = np.random.randint(
                Seedspercell = self._RNG.integers(
                    self._Seedmin, high=self._Seedmax, size=numShrubCells
                )
                # For each shrub cell, determine survival rate for the cell in this year
                SeedSurvYear = self._GermRate * self._RNG.random(numShrubCells)
                # SeedSurvYear = self._GermRate * np.random.rand(numShrubCells)

                for i in range(numShrubCells):
                    # For each shrub cell producing seeds, generate a random # of random numbers each rep. a single seed
                    randseeds = self._RNG.random(Seedspercell[i])
                    # randseeds = np.random.rand(Seedspercell[i])
                    # Find how many seeds produced in each cell actually survive
                    Survivors = len(randseeds[randseeds < SeedSurvYear[i]])

                    # Determine distance, rounding to nearest integer
                    if Survivors > 0:
                        DispDist = np.round(
                            # np.random.lognormal(
                            self._RNG.lognormal(
                                self._disp_mu, self._disp_sigma, Survivors
                            )
                        )
                    else:
                        DispDist = []

                    # If there are seeds to disperse
                    if len(DispDist) > 0:
                        for j in range(
                                len(DispDist + 1)
                        ):  # Loop through each individual seed to disperse
                            # Create a meshgrid to find coordinates of all points that are dropdistance from origin
                            if DispDist[j] > 0:
                                gridsize = int(DispDist[j] * 2 + 2)
                                X, Y = np.meshgrid(
                                    range(1, gridsize), range(1, gridsize)
                                )
                                originX = I[
                                    i
                                ]  # Sets coordinates of plant where seed is coming from
                                originY = k
                                matOriginX = math.ceil(len(X) / 2)
                                matOriginY = matOriginX
                                distMat = np.round(
                                    np.sqrt(
                                        (X - matOriginX) ** 2 + (Y - matOriginY) ** 2
                                    )
                                )  # Find the distance from origin to every other point on island
                                # Find coordinates of each point on island that is dropdistance away from origin
                                coords = np.where(distMat == DispDist[j])
                                row = coords[0]
                                col = coords[1]

                                # Randomly selct one of those points as the target point - this means equal probability
                                # of dispersal in every direction (valid assumption for avian dispersal)
                                if len(col) == 0:
                                    targetY = originX
                                    targetX = originY
                                else:
                                    # xx = random.randint(0, (len(col)) - 1)
                                    xx = self._RNG.integers(0, (len(col)) - 1)
                                    matTargetX = col[xx]
                                    matTargetY = row[xx]

                                    targetY = originX + (matTargetX - matOriginX)
                                    targetX = originY + (matTargetY - matOriginY)

                                # Put a plant in the ground if
                                #   -the dropdistance>0,
                                #   -the target drop location is within the island domain,
                                #   -there is no shrub at the receiving cell (male, female, or dead),
                                #   -the receiving cell has a tall enough fronting dune
                                #   -the receiving cell is within elevation range
                                if (
                                        not (not (0 <= targetY < DomainWidth) or not (
                                                0 <= targetX < self._BarrierLength) or not (
                                                ShrubDomainFemale[targetY, targetX] == 0))
                                        and ShrubDomainMale[targetY, targetX] == 0
                                        and self._ShrubEl_min <= InteriorDomain[targetY, targetX] <= self._ShrubEl_max
                                        and ShrubDomainDead[targetY, targetX] < 1
                                ):
                                    # Shrubs can establish if height of fronting dune greater than threshold
                                    if DuneDomainCrest[targetX] + self._BermEl >= self._Dshrub:
                                        # Decide if the tree wll be a female or male
                                        if self._RNG.uniform() > self._Female:
                                            ShrubDomainFemale[targetY, targetX] = 1
                                        else:
                                            ShrubDomainMale[targetY, targetX] = 1

                                    # Shrubs can establish without dune if locaed a certain distance from shoreline...
                                    elif BeachWidth + self._DuneWidth + targetY > self._SprayDist:
                                        if self._RNG.uniform() > 0.5:  # ... but have about 50% chance of survival
                                            # Decide if the tree wll be a female or male
                                            if self._RNG.uniform() > self._Female:
                                                ShrubDomainFemale[targetY, targetX] = 1
                                            else:
                                                ShrubDomainMale[targetY, targetX] = 1

        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale

        return (
            ShrubDomainAll,
            ShrubDomainFemale,
            ShrubDomainMale,
            BurialDomain,
            ShrubDomainDead,
        )

    def DuneGaps(self, DuneDomain, Dow, bermel, Rhigh):
        """Returns tuple of [gap start index, stop index, avg Rexcess of each gap,
        alpha: ratio of TWL / dune height]"""
        gaps = []
        start = 0
        #        stop = 0
        i = start
        while i < (len(Dow) - 1):
            adjacent = Dow[i + 1] - Dow[i]
            if adjacent == 1:
                i = i + 1
            else:
                stop = i

                x = DuneDomain[Dow[start]: (Dow[stop] + 1)]
                Hmean = sum(x) / float(len(x))
                Rexcess = Rhigh - (Hmean + bermel)
                alpha = Rhigh / (Hmean + bermel)
                gaps.append([Dow[start], Dow[stop], Rexcess, alpha])

                start = stop + 1
                i = start
        if i > 0:
            stop = i - 1

            x = DuneDomain[Dow[start]: (Dow[stop] + 1)]
            if len(x) > 0:
                Hmean = sum(x) / float(len(x))
                Rexcess = Rhigh - (Hmean + bermel)
                alpha = Rhigh / (Hmean + bermel)
                gaps.append([Dow[start], Dow[stop], Rexcess, alpha])
        return gaps

    def DiffuseDunes(self, DuneDomain, t):
        """Dune height diffusion in alongshore direction: smoothes dune height by redistributing
        sand from high dune to neighboring low dune(s)"""

        # from Barrier3D_Parameters import (HdDiffu, BarrierLength, DuneWidth)

        for w in range(self._DuneWidth):
            # Alongshore direction
            for d in range(2, self._BarrierLength):  # Loop L to R
                Ldiff = (
                        DuneDomain[t, d, w] - DuneDomain[t, d - 1, w]
                )  # Calculate height difference
                # Subtract sand from taller dune cell and add to adjacent shorter one (assumes sand in dunes is consrvd)
                if Ldiff > self._HdDiffu:
                    sub = (
                                  Ldiff - self._HdDiffu
                          ) / 2  # Height difference in excess of maximum, divided by two
                    DuneDomain[t, d, w] = (
                            DuneDomain[t, d, w] - sub
                    )  # Lower dune that's too tall
                    DuneDomain[t, d - 1, w] = (
                            DuneDomain[t, d - 1, w] + sub
                    )  # Raise dune that's too short
            for d in range((self._BarrierLength - 2), 0, -1):  # Loop R to L
                Rdiff = DuneDomain[t, d, w] - DuneDomain[t, d + 1, w]
                if Rdiff > self._HdDiffu:
                    sub = (Rdiff - self._HdDiffu) / 2
                    DuneDomain[t, d, w] = DuneDomain[t, d, w] - sub
                    DuneDomain[t, d + 1, w] = DuneDomain[t, d + 1, w] + sub

        return DuneDomain

    def SalineFlooding(
            self,
            ShrubDomainWidth,
            ShrubDomainAll,
            ShrubDomainFemale,
            ShrubDomainMale,
            ShrubDomainDead,
            d,
            i,
            Q0,
    ):
        """Kill all immature (< 1 yr-old) shrubs that have been flooded
        beyond a threshold discharge (Tolliver et al., 1997)"""

        if d < (ShrubDomainWidth - 1) and i < (self._BarrierLength - 1):
            if Q0 >= self._SalineLimit and ShrubDomainAll[d, i] == 1:
                ShrubDomainDead[d, i] = (
                        ShrubDomainMale[d, i] + ShrubDomainFemale[d, i]
                )  # Transfer to dead shrub domain
                ShrubDomainFemale[d, i] = 0
                ShrubDomainMale[d, i] = 0

        return ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead

    def UpdateBurial(self, BurialDomain, ElevationChange, ShrubDomainWidth, ShrubDomainAll):
        """Updates amount of burial/erosion for each shrub"""

        BurialDomain = BurialDomain + ElevationChange[1: ShrubDomainWidth + 1, :]
        BurialDomain[ShrubDomainAll == 0] = 0

        return BurialDomain

    def UpdateShrubDomains(
            self,
            DomainWidth,
            ShrubDomainWidth,
            ShrubDomainFemale,
            ShrubDomainMale,
            ShrubDomainAll,
            ShrubPercentCover,
            BurialDomain,
            ShrubDomainDead,
            DeadPercentCover,
    ):
        """Updates size of shrub domains"""

        if DomainWidth > ShrubDomainWidth:
            AddRows = np.zeros([DomainWidth - ShrubDomainWidth, self._BarrierLength])
            ShrubDomainFemale = np.vstack([ShrubDomainFemale, AddRows])
            ShrubDomainMale = np.vstack([ShrubDomainMale, AddRows])
            ShrubDomainDead = np.vstack([ShrubDomainDead, AddRows])
            ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
            ShrubPercentCover = np.vstack([ShrubPercentCover, AddRows])
            DeadPercentCover = np.vstack([DeadPercentCover, AddRows])
            BurialDomain = np.vstack([BurialDomain, AddRows])
        elif DomainWidth < ShrubDomainWidth:
            RemoveRows = ShrubDomainWidth - DomainWidth
            ShrubDomainFemale = ShrubDomainFemale[0:-RemoveRows, :]
            ShrubDomainMale = ShrubDomainMale[0:-RemoveRows, :]
            ShrubDomainDead = ShrubDomainDead[0:-RemoveRows, :]
            ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
            ShrubPercentCover = ShrubPercentCover[0:-RemoveRows, :]
            DeadPercentCover = DeadPercentCover[0:-RemoveRows, :]
            BurialDomain = BurialDomain[0:-RemoveRows, :]

        return (
            ShrubDomainFemale,
            ShrubDomainMale,
            ShrubDomainAll,
            ShrubPercentCover,
            BurialDomain,
            ShrubDomainDead,
            DeadPercentCover,
        )

    def LTA_SC(
            self,
            InteriorDomain,
            OWloss,
            Qdg,
            DuneLoss,
            x_s,
            x_s_TS,
            x_t,
            x_t_TS,
            x_b_TS,
            s_sf_TS,
            InteriorWidth_Avg,
            SL,
            QowTS,
            QsfTS,
            t,
            h_b_TS,
    ):
        """Finds shoreline change for modeled year following Lorenzo-Trueba and Ashton (2014)"""

        # Find volume of shoreface/beach/dune sand deposited in island interior and back-barrier
        Qow = OWloss / (
            self._BarrierLength
        )  # (dam^3/dam) Volume of sediment lost from shoreface/beach by overwash
        if Qow < 0:
            Qow = 0
        QowTS.append(Qow * 100)  # Save in m^3/m
        if DuneLoss < Qow:
            Qow = Qow - DuneLoss  # Account for dune contribution to overwash volume
        else:
            Qow = 0  # Excess DuneLoss assumed lost offshore

        if Qdg < 0:
            Qdg = 0

        # DefineParams
        d_sf = self._DShoreface
        h_b = np.average(InteriorDomain[InteriorDomain >= SL])
        x_b = x_s + InteriorWidth_Avg

        # Shoreface Flux
        s_sf = d_sf / (x_s - x_t)
        Qsf = (self._k_sf / 100) * (
                self._s_sf_eq - s_sf
        )  # Convert Ksf from m^3/m/yr to dam^3/dam/yr
        QsfTS.append(Qsf * 100)  # Save in m^3/m

        # Toe, Shoreline, and island base elevation changes
        x_t_dt = (4 * Qsf * (h_b + d_sf) / (d_sf * (2 * h_b + d_sf))) + (
                2 * self._RSLR[t] / s_sf
        )
        x_s_dt = 2 * (Qow + Qdg + self._Qat) / ((2 * h_b) + d_sf) - (
                4 * Qsf * (h_b + d_sf) / (((2 * h_b) + d_sf) ** 2)
        )  # Dune growth and alongshore transport added to LTA14 formulation

        # Record changes
        x_t = x_t + x_t_dt
        x_s = x_s + x_s_dt
        x_t_TS.append(x_t)
        x_s_TS.append(x_s)
        x_b_TS.append(x_b)
        s_sf_TS.append(s_sf)
        h_b_TS.append(h_b)

        return x_s, x_t, x_s_TS, x_t_TS, x_b_TS, s_sf_TS, QowTS, QsfTS, h_b_TS

    def CalcPC(self, ShrubDomainAll, PercentCoverTS, ShrubDomainDead, DeadPercentCoverTS, ShrubArea, t):
        """Calculates percent cover of shrub domain and shrub coverage area (dam^2)"""

        Allshrub_t = ShrubDomainAll.astype("int64")
        Deadshrub_t = ShrubDomainDead.astype("int64")
        ShrubPercentCover = self._PC.take(Allshrub_t)
        DeadPercentCover = self._PC.take(Deadshrub_t)

        #    ShrubPercentCover[ShrubPercentCover == 0] = DeadPercentCover[ShrubPercentCover == 0]
        PercentCoverTS[t] = ShrubPercentCover
        DeadPercentCoverTS[t] = DeadPercentCover
        ShrubArea.append(np.count_nonzero(ShrubDomainAll))

        return ShrubPercentCover, PercentCoverTS, ShrubArea, DeadPercentCover, DeadPercentCoverTS

    def __init__(self, **kwds):
        self._RNG = kwds.pop("RNG")

        self._SL = 0  # Does not change (Lagrangian frame of reference)

        self._BarrierLength = kwds.pop("BarrierLength")
        self._BayDepth = kwds.pop("BayDepth")
        self._BermEl = kwds.pop("BermEl")
        self._beta = kwds.pop("beta")
        self._BurialLimit = kwds.pop("BurialLimit")
        self._C1 = kwds.pop("C1")
        self._C2 = kwds.pop("C2")
        self._Cbb_i = kwds.pop("Cbb_i")
        self._Cbb_r = kwds.pop("Cbb_r")
        self._Cx = kwds.pop("Cx")
        self._disp_mu = kwds.pop("disp_mu")
        self._disp_sigma = kwds.pop("disp_sigma")
        self._Dmaxel = kwds.pop("Dmaxel")
        self._DomainWidth = kwds.pop("DomainWidth")
        self._DShoreface = kwds.pop("DShoreface")
        self._Dshrub = kwds.pop("Dshrub")
        Dstart = kwds.pop("Dstart")
        self._DuneDomain = kwds.pop("DuneDomain")
        self._DuneRestart = kwds.pop("DuneRestart")
        self._DuneWidth = kwds.pop("DuneWidth")
        self._Female = kwds.pop("Female")
        self._GermRate = kwds.pop("GermRate")
        self._growthparam = kwds.pop("growthparam")
        self._HdDiffu = kwds.pop("HdDiffu")
        self._InteriorDomain = kwds.pop("InteriorDomain")
        self._k_sf = kwds.pop("k_sf")
        self._Ki = kwds.pop("Ki")
        self._Kr = kwds.pop("Kr")
        LShoreface = kwds.pop("LShoreface")
        self._mm = kwds.pop("mm")
        self._MHW = kwds.pop("MHW")
        self._nn = kwds.pop("nn")
        self._OWss_i = kwds.pop("OWss_i")
        self._OWss_r = kwds.pop("OWss_r")
        self._PC = kwds.pop("PC")
        self._Qat = kwds.pop("Qat")
        self._Qs_bb_min = kwds.pop("Qs_bb_min")
        self._Qs_min = kwds.pop("Qs_min")
        self._Qshrub_max = kwds.pop("Qshrub_max")
        self._Rin_i = kwds.pop("Rin_i")
        self._Rin_r = kwds.pop("Rin_r")
        self._RSLR = kwds.pop("RSLR")
        self._s_sf_eq = kwds.pop("s_sf_eq")
        self._SalineLimit = kwds.pop("SalineLimit")
        self._Seedmax = kwds.pop("Seedmax")
        self._Seedmin = kwds.pop("Seedmin")
        self._Shrub_ON = kwds.pop("Shrub_ON")
        self._ShrubEl_max = kwds.pop("ShrubEl_max")
        self._ShrubEl_min = kwds.pop("ShrubEl_min")
        self._StormSeries = kwds.pop("StormSeries")
        self._StormStart = kwds.pop("StormStart")
        self._threshold_in = kwds.pop("threshold_in")
        self._TimeFruit = kwds.pop("TimeFruit")
        self._TMAX = kwds.pop("TMAX")
        self._UprootLimit = kwds.pop("UprootLimit")
        self._TideAmp = kwds.pop("TideAmp")
        self._SprayDist = kwds.pop("SprayDist")
        self._MaxUpSlope = kwds.pop("MaxUpSlope")
        self._DuneParamStart = kwds.pop("DuneParamStart")
        self._GrowthParamStart = kwds.pop("GrowthParamStart")
        self._MaxShrubHeight = kwds.pop("MaxShrubHeight")
        self._RSLR_Constant = kwds.pop("RSLR_Constant")
        self._RSLR_const = kwds.pop("RSLR_const")
        self._SH = kwds.pop("SH")
        self._x_t = kwds.pop("x_t")

        if len(kwds) > 0:
            raise ValueError(
                "unrecognized keywords ({0})".format(", ".join(kwds.keys()))
            )

        # Note, the following variables are not used
        self._MaxAvgSlope = self._BermEl / 10
        # fluxLimit = 100  # Initialize

        # ### Initialize Shrubs
        # if Shrub_ON ==1: print('Shrubs ON')

        self._PercentCoverTS = [None] * self._TMAX
        self._PercentCoverTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._DeadPercentCoverTS = [None] * self._TMAX
        self._DeadPercentCoverTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubFemaleTS = [None] * self._TMAX
        self._ShrubFemaleTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubMaleTS = [None] * self._TMAX
        self._ShrubMaleTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        # self._ShrubDomainFemale = self._ShrubFemaleTS[0]
        # self._ShrubDomainMale = self._ShrubMaleTS[0]

        self._ShrubDeadTS = [None] * self._TMAX
        self._ShrubDeadTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubDomainFemale = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubDomainMale = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubDomainDead = np.zeros([self._DomainWidth, self._BarrierLength])

        self._ShrubPercentCover = self._PercentCoverTS[0]
        self._DeadPercentCover = self._DeadPercentCoverTS[0]
        self._BurialDomain = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubArea = [0]

        # ### Initialize variables
        self._DomainTS = [
                             None
                         ] * self._TMAX  # Stores the elevation domain for each timestep
        self._DomainTS[0] = self._InteriorDomain
        self._SCRagg = 0  # Counter variable for shoreline change that is less than 1 cell size per year
        self._ShorelineChangeTS = [0]
        self._ShorelineChange = 0  # (dam) Total amount of shoreline change
        # NOTE: BBShorelineChangeTS and BBShorelineChange are unused
        # BBShorelineChangeTS = [0]
        # BBShorelineChange = 0  # (dam) Total amount of back-barrier shoreline extension
        # self._StormCount = []
        self._StormCount = [0]
        self._InundationCount = 0
        self._RunUpCount = 0

        # Added by KA: use FindWidths to calculate the average interior width for setting initial back barrier shoreline
        _, _, InteriorWidth_Avg = self.FindWidths(
            self._InteriorDomain, self._SL
        )

        # self._InteriorWidth_AvgTS = [self._DomainWidth]
        self._InteriorWidth_AvgTS = [InteriorWidth_Avg]
        self._QowTS = [0]  # (m^3/m)
        # self._x_t = 0  # (dam) Start location of shoreface toe
        self._x_s = self._x_t + LShoreface  # (dam) Start location of shoreline
        self._x_t_TS = [self._x_t]  # (dam) Shoreface toe locations for each time step
        self._x_s_TS = [self._x_s]  # (dam) Shoreline locations for each time step
        # self._x_b_TS = [
        #     (self._x_s + self._DomainWidth)
        # ]  # (dam) Bay shoreline locations for each time step
        self._x_b_TS = [
            (self._x_s + InteriorWidth_Avg)
        ]  # (dam) Bay shoreline locations for each time step
        self._h_b_TS = [
            (np.average(self._InteriorDomain[self._InteriorDomain >= self._SL]))
        ]  # (dam) average height of barrier for each time step
        self._s_sf_TS = [(self._DShoreface / LShoreface)]
        self._Hd_AverageTS = [Dstart]
        self._QsfTS = [0]  # (m^3/m)

        # Dune height loss exclusively from vertical storm erosion
        self._Hd_Loss_TS = np.zeros([self._TMAX, self._BarrierLength])

        self._time_index = 1

    @classmethod
    def from_path(cls, path_to_folder, fmt="yaml"):
        return cls(**load_inputs(path_to_folder, prefix="barrier3d", fmt=fmt))
        # return cls(**Barrier3dParameters.from_xlsx(path_to_xlsx))

    @classmethod
    def from_xlsx(cls, path_to_xlsx):
        return cls(**load_inputs(path_to_xlsx, prefix="barrier3d", fmt="xlsx"))
        # return cls(**Barrier3dParameters.from_xlsx(path_to_xlsx))

    @classmethod
    def from_yaml(cls, path_to_yaml):
        return cls(**load_inputs(path_to_yaml, prefix="barrier3d", fmt="yaml"))
        # return cls(**Barrier3dParameters.from_yaml(path_to_yaml))

    def update(self):

        if self._time_index > 1:
            # update dune domain (erode/prograde) based on shoreline change from last time step
            drown_break = 0

            SCR = (
                    (self._x_s_TS[-1] - self._x_s_TS[-2]) * -1
            )  # (dam) Change rate of ocean-fronting shoreline: (+) = progradation, (-) = erosion

            self._SCRagg = (
                   self._SCRagg + SCR
            )  # Account for any residual shoreline change (less than cell size) from previous time step

            if abs(self._SCRagg) >= 1:
                sc = math.floor(abs(self._SCRagg))

                if (
                       self._SCRagg > 0
                ):  # Positive = prograde, add row(s) to front of interior domain
                    for d in range(sc):
                        # New interior row added with elev. of previous dune field
                        # (i.e. previous dune row now part of int.)
                        self._InteriorDomain = np.vstack(
                            [
                                self._DuneDomain[self._time_index-1, :, -1] + self._BermEl,
                                self._InteriorDomain,
                            ]
                        )
                        if self._DuneParamStart == 0:
                            newDuneHeight = np.ones([self._BarrierLength]) * (
                                    0.005
                                    + (
                                            -0.005
                                            + (0.005 - (-0.005))
                                            * self._RNG.random(self._BarrierLength)
                                    )
                            )
                        else:
                            newDuneHeight = np.ones([self._BarrierLength]) * 0.01

                        self._DuneDomain[self._time_index-1, :, :] = np.roll(
                            self._DuneDomain[self._time_index-1, :, :], 1, axis=1
                        )
                        self._DuneDomain[self._time_index-1, :, 0] = newDuneHeight

                        # Update shrub domains too
                        self._ShrubDomainFemale = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._ShrubDomainFemale)
                        )
                        self._ShrubDomainMale = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._ShrubDomainMale)
                        )
                        self._ShrubDomainDead = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._ShrubDomainDead)
                        )
                        self._ShrubPercentCover = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._ShrubPercentCover)
                        )
                        self._DeadPercentCover = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._DeadPercentCover)
                        )
                        self._BurialDomain = np.concatenate(
                            (np.zeros([1, self._BarrierLength]), self._BurialDomain)
                        )

                    # DomainWidth = DomainWidth + sc  # Update width of interior domain
                    self._ShorelineChange = self._ShorelineChange + sc
                    self._ShorelineChangeTS.append(+sc)
                    self._SCRagg = self._SCRagg - sc  # Reset, leaving residual

                elif (
                        self._SCRagg < 0
                ):  # Negative = erode, remove front row(s) of interior domain
                    for d in range(sc):
                        newDuneElev = self._InteriorDomain[
                                      0, :
                                      ]  # First row of interior will now become new part of active dune field
                        newDuneHeight = newDuneElev - self._BermEl
                        newDuneHeight[newDuneHeight < self._DuneRestart] = self._DuneRestart
                        self._InteriorDomain = np.delete(self._InteriorDomain, 0, axis=0)
                        if np.shape(self._InteriorDomain)[0] <= 0:
                            drown_break = 1
                            break
                        self._DuneDomain[self._time_index-1, :, :] = np.roll(
                            self._DuneDomain[self._time_index-1, :, :], -1, axis=1
                        )
                        self._DuneDomain[self._time_index-1, :, -1] = newDuneHeight

                        # Update shrub domains too
                        self._ShrubDomainFemale = np.delete(
                            self._ShrubDomainFemale, 0, axis=0
                        )
                        self._ShrubDomainMale = np.delete(
                            self._ShrubDomainMale, 0, axis=0
                        )
                        self._ShrubDomainDead = np.delete(
                            self._ShrubDomainDead, 0, axis=0
                        )
                        self._ShrubPercentCover = np.delete(
                            self._ShrubPercentCover, 0, axis=0
                        )
                        self._DeadPercentCover = np.delete(
                            self._DeadPercentCover, 0, axis=0
                        )
                        self._BurialDomain = np.delete(self._BurialDomain, 0, axis=0)

                    # DomainWidth = DomainWidth - sc  # Update width of interior domain
                    self._ShorelineChange = self._ShorelineChange - sc
                    self._ShorelineChangeTS.append(-sc)
                    self._SCRagg = self._SCRagg + sc  # Reset, leaving residual
            else:
                self._ShorelineChangeTS.append(0)

            # ### Check for drowning
            if drown_break == 1:
                print(
                    "Barrier has WIDTH DROWNED at t = " + str(self._time_index-1) + " years"
                )
                self._TMAX = self._time_index - 2
                raise Barrier3dError(
                    "Barrier has WIDTH DROWNED at t = {time} years".format(
                        time=self._time_index - 1
                    )
                )
            elif all(j <= self._SL for j in self._InteriorDomain[0, :]):
                print(
                    "Barrier has HEIGHT DROWNED at t = " + str(self._time_index-1) + " years"
                )
                self._TMAX = self._time_index - 2
                raise Barrier3dError(
                    "Barrier has HEIGHT DROWNED at t = {time} years".format(
                        time=self._time_index - 1
                    )
                )

            # ### Recalculate and save DomainWidth and InteriorWidth from previous time step
            DomainWidth, InteriorWidth, InteriorWidth_Avg = self.FindWidths(
                self._InteriorDomain, self._SL
            )
            self._InteriorWidth_AvgTS.append(InteriorWidth_Avg)

            # Added by KA: replace value of x_b_TS from the last time step with new InteriorWidth_Avg
            # (also from the last time step): note that this is just for completeness
            self._x_b_TS[-1] = self._x_s_TS[-1] + self._InteriorWidth_AvgTS[-1]

            # ###########################################
            # ### Save domains of the previous timestep
            # ###########################################

            self._DomainTS[self._time_index-1] = self._InteriorDomain

            zero = np.zeros([DomainWidth, self._BarrierLength])
            self._ShrubFemaleTS[self._time_index-1] = (
                    self._ShrubDomainFemale + zero
            )  # Correctly saves domain to list only if array of zeros are added to domain first, no clue why - 4Mar20
            self._ShrubMaleTS[self._time_index-1] = self._ShrubDomainMale + zero
            self._ShrubDeadTS[self._time_index-1] = self._ShrubDomainDead + zero

            # Calculate Percent Cover and Area
            self._ShrubDomainAll = self._ShrubDomainMale + self._ShrubDomainFemale
            (
                self._ShrubPercentCover,
                self._PercentCoverTS,
                self._ShrubArea,
                self._DeadPercentCover,
                self._DeadPercentCoverTS
            ) = self.CalcPC(
                self._ShrubDomainAll,
                self._PercentCoverTS,
                self._ShrubDomainDead,
                self._DeadPercentCoverTS,
                self._ShrubArea,
                self._time_index-1,
            )

        # ###########################################
        # ### increase sea level and start new time step
        # ###########################################

        # ### RSLR
        self._InteriorDomain, self._DuneDomain = self.SeaLevel(
            self._InteriorDomain, self._DuneDomain, self._time_index
        )

        # ### Find DomainWidth and InteriorWidth
        DomainWidth, InteriorWidth, InteriorWidth_Avg = self.FindWidths(
            self._InteriorDomain, self._SL
        )

        # # ### Check for drowning
        # # Barrier drowns if interior width thins to 10 m or fewer or all interior cells are below SL
        # if max(InteriorWidth) <= 1:
        #     print(
        #         "Barrier has WIDTH DROWNED at t = " + str(self._time_index) + " years!"
        #     )
        #     self._TMAX = self._time_index - 1
        #     raise Barrier3dError(
        #         "Barrier has WIDTH DROWNED at t = {time} years".format(
        #             time=self._time_index
        #         )
        #     )
        # if all(j <= self._SL for j in self._InteriorDomain[0, :]):
        #     print(
        #         "Barrier has HEIGHT DROWNED at t = " + str(self._time_index) + " years!"
        #     )
        #     self._TMAX = self._time_index - 1
        #     raise Barrier3dError(
        #         "Barrier has HEIGHT DROWNED at t = {time} years".format(
        #             time=self._time_index
        #         )
        #     )

        # ### Grow Dunes
        self._DuneDomain, self._Dmax, Qdg = self.DuneGrowth(
            self._DuneDomain, self._time_index
        )

        DuneDomainCrest = self._DuneDomain[self._time_index, :, :].max(
            axis=1
        )  # Maximum height of each row in DuneDomain
        DuneDomainCrest[DuneDomainCrest < self._DuneRestart] = self._DuneRestart

        self._Hd_AverageTS.append(
            np.mean(DuneDomainCrest)
        )  # Store average pre-storms dune-heigh for time step

        # ###########################################
        # ### Shrubs

        # if self._Shrub_ON == 1:
        if self._Shrub_ON:
            (
                self._ShrubDomainAll,
                self._ShrubDomainFemale,
                self._ShrubDomainMale,
                self._BurialDomain,
                self._ShrubDomainDead,
            ) = self.Shrubs(
                self._InteriorDomain,
                DuneDomainCrest,
                self._ShrubDomainFemale,
                self._ShrubDomainMale,
                self._BurialDomain,
                InteriorWidth_Avg,
                DomainWidth,
                self._ShrubDomainDead,
            )

        # ###########################################
        # ### Storms

        OWloss = 0
        DuneLoss = 0
        numstorm = 0

        if self._time_index >= self._StormStart:
            # Select number of storms for this time step from normal distribution
            TSloc = np.argwhere(self._StormSeries[:, 0] == self._time_index)
            numstorm = int(len(TSloc))  # analysis:ignore

            if numstorm > 0:
                # TSloc = np.argwhere(self._StormSeries[:, 0] == self._time_index)
                start = TSloc[0, 0]
                stop = TSloc[-1, 0] + 1
                Rhigh = self._StormSeries[start:stop, 1]
                Rlow = self._StormSeries[start:stop, 2]
                # Note, period not used
                # period = self._StormSeries[start:stop, 3]
                dur = np.array(self._StormSeries[start:stop, 4], dtype="int")

                # ### Individual Storm Impacts
                for n in range(numstorm):  # Loop through each individual storm

                    # ###########################################
                    # ### Dune Erosion

                    DuneChange = np.zeros(
                        [self._BarrierLength, self._DuneWidth]
                    )  # Vector storing dune height change for this storm

                    # Find overwashed dunes and gaps
                    Dow = [
                        index
                        for index, value in enumerate((DuneDomainCrest + self._BermEl))
                        if value < Rhigh[n]
                    ]
                    gaps = self.DuneGaps(
                        DuneDomainCrest, Dow, self._BermEl, Rhigh[n]
                    )  # Finds location and Rexcess of continuous gaps in dune ridge

                    for d in range(len(Dow)):  # Loop through each overwashed dune cell
                        for w in range(self._DuneWidth):

                            # Calculate dune elevation loss
                            Rnorm = Rhigh[n] / (
                                    self._DuneDomain[self._time_index, Dow[d], w]
                                    + self._BermEl
                            )  # Rhigh relative to pre-storm dune elevation
                            Dloss = Rnorm / (
                                    self._C1 + (Rnorm * (Rnorm - self._C2))
                            )  # Amount of dune crest elevation change normalized by pre-storm dune elevation
                            # (i.e. a percent change), from Goldstein and Moore (2016)

                            # Set new dune height
                            InitDElev = (
                                    self._DuneDomain[self._time_index, Dow[d], w]
                                    + self._BermEl
                            )
                            NewDElev = InitDElev * (
                                    1 - Dloss
                            )  # Calculate new dune elevation from storm lowering
                            if NewDElev < self._BermEl:
                                NewDElev = self._BermEl
                            self._DuneDomain[self._time_index, Dow[d], w] = (
                                    NewDElev - self._BermEl
                            )  # Convert elevation to height above berm

                            DuneChange[Dow[d], w] = InitDElev - NewDElev

                            # If dune is lowered to ~ zero, allow for chance of regrowth by raising dune height to 5 cm
                            if (
                                    self._DuneDomain[self._time_index, Dow[d], w]
                                    < self._DuneRestart
                            ):
                                if self._DuneRestart < self._Dmax:
                                    self._DuneDomain[
                                        self._time_index, Dow[d], w
                                    ] = self._DuneRestart
                                else:
                                    self._DuneDomain[
                                        self._time_index, Dow[d], w
                                    ] = self._Dmax  # Restart height can't be greater than Dmax

                    # Dune Height Diffusion
                    self._DuneDomain = self.DiffuseDunes(
                        self._DuneDomain, self._time_index
                    )
                    self._DuneDomain[self._DuneDomain < self._SL] = self._DuneRestart

                    DuneLoss = np.sum(DuneChange) / self._BarrierLength
                    Hd_TSloss = (
                            DuneChange.max(axis=1) / dur[n]
                    )  # Average height of dune loss for each substep during storm

                    self._Hd_Loss_TS[self._time_index, :] = self._Hd_Loss_TS[
                                                            self._time_index, :
                                                            ] + DuneChange.max(axis=1)

                    # ###########################################
                    # ### Overwash

                    Iow = 0  # Count of dune gaps in inundation regime
                    Dunes_prestorm = DuneDomainCrest
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        # Note, Rexcess not used
                        # Rexcess = gaps[q][2]  # (m)
                        gapwidth = stop - start + 1
                        meandune = (
                                           sum(Dunes_prestorm[start: stop + 1]) / gapwidth
                                   ) + self._BermEl  # Average elevation of dune gap

                        # Determine number of gaps in inundation regime
                        if Rlow[n] > meandune:
                            Iow += 1

                    # Determine Sediment And Water Routing Rules

                    if (
                            len(gaps) > 0 and Iow / len(gaps) >= self._threshold_in
                    ):  # If greater than threshold % of dune gaps are inunundation regime, use inun. regime routing
                        inundation = 1
                        substep = self._OWss_i
                        self._InundationCount += 1
                    else:
                        inundation = 0
                        substep = self._OWss_r
                        self._RunUpCount += 1

                    # Set Domain
                    add = 10
                    duration = dur[n] * substep
                    width = (
                            np.shape(self._InteriorDomain)[0] + 1 + add
                    )  # (dam) Add one for Dunes and 25 for bay
                    Elevation = np.zeros([duration, width, self._BarrierLength])
                    Dunes = Dunes_prestorm + self._BermEl
                    Bay = np.ones([add, self._BarrierLength]) * -self._BayDepth
                    Elevation[0, :, :] = np.vstack([Dunes, self._InteriorDomain, Bay])

                    # Initialize Memory Storage Arrays
                    Discharge = np.zeros([duration, width, self._BarrierLength])
                    SedFluxIn = np.zeros([duration, width, self._BarrierLength])
                    SedFluxOut = np.zeros([duration, width, self._BarrierLength])

                    Rin = 0  # (dam^3/t) Infiltration Rate, volume of overwash flow lost per m cross-shore per time

                    # Set Water at Dune Crest
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2]  # (m)
                        #                        gapwidth = stop - start + 1
                        #                        meandune = (
                        #                            sum(Dunes_prestorm[start : stop + 1]) / gapwidth
                        #                        ) + self._BermEl  # Average elevation of dune gap

                        # Calculate discharge through each dune cell
                        Vdune = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/s)
                        Qdune = Vdune * Rexcess * 3600  # (dam^3/hr)

                        # Set discharge at dune gap
                        Discharge[:, 0, start:stop] = Qdune

                        if inundation == 1:  # Inundation regime
                            Rin = self._Rin_i

                            # Find average slope of interior
                            AvgSlope = self._BermEl / InteriorWidth_Avg

                            # Enforce max average interior slope
                            AvgSlope = min(self._MaxAvgSlope, AvgSlope)

                            C = self._Cx * AvgSlope  # Momentum constant

                        else:  # Run-up regime
                            Rin = self._Rin_r

                    # ### Run Flow Routing Algorithm
                    for TS in range(duration):

                        ShrubDomainWidth = np.shape(self._ShrubDomainFemale)[0]
                        DeadDomainWidth = np.shape(self._ShrubDomainDead)[0]

                        if TS > 0:
                            Elevation[TS, 1:, :] = Elevation[
                                                   TS - 1, 1:, :
                                                   ]  # Begin timestep with elevation from end of last
                            Elevation[TS, 0, :] = Dunes - (
                                    Hd_TSloss / substep * TS
                            )  # Reduce dune in height linearly over course of storm

                        for d in range(width - 1):
                            # Reduce discharge across row via infiltration
                            if d > 0:
                                Discharge[TS, d, :][Discharge[TS, d, :] > 0] -= Rin
                            Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                            for i in range(self._BarrierLength):
                                if Discharge[TS, d, i] > 0:

                                    Q0 = Discharge[TS, d, i]

                                    # ### Calculate Slopes
                                    if i > 0:
                                        S1 = (
                                                     Elevation[TS, d, i]
                                                     - Elevation[TS, d + 1, i - 1]
                                             ) / (math.sqrt(2))
                                        S1 = np.nan_to_num(S1)
                                    else:
                                        S1 = 0

                                    S2 = Elevation[TS, d, i] - Elevation[TS, d + 1, i]
                                    S2 = np.nan_to_num(S2)

                                    if i < (self._BarrierLength - 1):
                                        S3 = (
                                                     Elevation[TS, d, i]
                                                     - Elevation[TS, d + 1, i + 1]
                                             ) / (math.sqrt(2))
                                        S3 = np.nan_to_num(S3)
                                    else:
                                        S3 = 0

                                    # ### Calculate Discharge To Downflow Neighbors
                                    # One or more slopes positive
                                    if S1 > 0 or S2 > 0 or S3 > 0:

                                        if S1 < 0:
                                            S1 = 0
                                        if S2 < 0:
                                            S2 = 0
                                        if S3 < 0:
                                            S3 = 0

                                        Q1 = (
                                                Q0
                                                * S1 ** self._nn
                                                / (
                                                        S1 ** self._nn
                                                        + S2 ** self._nn
                                                        + S3 ** self._nn
                                                )
                                        )
                                        Q2 = (
                                                Q0
                                                * S2 ** self._nn
                                                / (
                                                        S1 ** self._nn
                                                        + S2 ** self._nn
                                                        + S3 ** self._nn
                                                )
                                        )
                                        Q3 = (
                                                Q0
                                                * S3 ** self._nn
                                                / (
                                                        S1 ** self._nn
                                                        + S2 ** self._nn
                                                        + S3 ** self._nn
                                                )
                                        )

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
                                        if S3 == 0 and i < (self._BarrierLength - 1):
                                            Q3 = Qx
                                        else:
                                            Q3 = 0

                                    # All slopes negative
                                    else:

                                        Q1 = (
                                                Q0
                                                * abs(S1) ** (-self._nn)
                                                / (
                                                        abs(S1) ** (-self._nn)
                                                        + abs(S2) ** (-self._nn)
                                                        + abs(S3) ** (-self._nn)
                                                )
                                        )
                                        Q2 = (
                                                Q0
                                                * abs(S2) ** (-self._nn)
                                                / (
                                                        abs(S1) ** (-self._nn)
                                                        + abs(S2) ** (-self._nn)
                                                        + abs(S3) ** (-self._nn)
                                                )
                                        )
                                        Q3 = (
                                                Q0
                                                * abs(S3) ** (-self._nn)
                                                / (
                                                        abs(S1) ** (-self._nn)
                                                        + abs(S2) ** (-self._nn)
                                                        + abs(S3) ** (-self._nn)
                                                )
                                        )

                                        Q1 = np.nan_to_num(Q1)
                                        Q2 = np.nan_to_num(Q2)
                                        Q3 = np.nan_to_num(Q3)

                                        # MaxUpSlope = 0.25  # dam

                                        if S1 > self._MaxUpSlope:
                                            Q1 = 0
                                        else:
                                            Q1 = Q1 * (1 - (abs(S1) / self._MaxUpSlope))

                                        if S2 > self._MaxUpSlope:
                                            Q2 = 0
                                        else:
                                            Q2 = Q2 * (1 - (abs(S2) / self._MaxUpSlope))

                                        if S3 > self._MaxUpSlope:
                                            Q3 = 0
                                        else:
                                            Q3 = Q3 * (1 - (abs(S3) / self._MaxUpSlope))

                                    # ### Reduce Overwash Through Shrub Cells and Save Discharge
                                    # if self._Shrub_ON == 1:
                                    if self._Shrub_ON:
                                        # Cell 1
                                        if i > 0:
                                            if d < ShrubDomainWidth and self._ShrubPercentCover[d, i - 1] > 0:
                                                Q1 = Q1 * ((1 - self._Qshrub_max) * self._ShrubPercentCover[d, i - 1])
                                            elif d < DeadDomainWidth and self._DeadPercentCover[d, i - 1] > 0:
                                                Q1 = Q1 * ((1 - self._Qshrub_max * 0.66) * self._DeadPercentCover[
                                                    d, i - 1])  # Dead shrubs block 2/3 of what living shrubs block
                                            Discharge[TS, d + 1, i - 1] = (
                                                    Discharge[TS, d + 1, i - 1] + Q1
                                            )

                                        # Cell 2
                                        # (ShrubDomainWidth - 1)?
                                        if d < ShrubDomainWidth and self._ShrubPercentCover[d, i] > 0:
                                            Q2 = Q2 * ((1 - self._Qshrub_max) * self._ShrubPercentCover[d, i])
                                        elif d < DeadDomainWidth and self._DeadPercentCover[d, i] > 0:
                                            Q2 = Q2 * ((1 - self._Qshrub_max * 0.66) * self._DeadPercentCover[d, i])
                                        Discharge[TS, d + 1, i] = (
                                                Discharge[TS, d + 1, i] + Q2
                                        )

                                        # Cell 3
                                        if i < (self._BarrierLength - 1):
                                            if d < ShrubDomainWidth and self._ShrubPercentCover[d, i + 1] > 0:
                                                Q3 = Q3 * ((1 - self._Qshrub_max) * self._ShrubPercentCover[d, i + 1])
                                            elif d < DeadDomainWidth and self._DeadPercentCover[d, i + 1] > 0:
                                                Q3 = Q3 * ((1 - self._Qshrub_max * 0.66) * self._DeadPercentCover[
                                                    d, i + 1])
                                            Discharge[TS, d + 1, i + 1] = (
                                                    Discharge[TS, d + 1, i + 1] + Q3
                                            )
                                    else:
                                        # Cell 1
                                        if i > 0:
                                            Discharge[TS, d + 1, i - 1] = (
                                                    Discharge[TS, d + 1, i - 1] + Q1
                                            )

                                        # Cell 2
                                        Discharge[TS, d + 1, i] = (
                                                Discharge[TS, d + 1, i] + Q2
                                        )

                                        # Cell 3
                                        if i < (self._BarrierLength - 1):
                                            Discharge[TS, d + 1, i + 1] = (
                                                    Discharge[TS, d + 1, i + 1] + Q3
                                            )

                                    # ### Calculate Sed Movement
                                    fluxLimit = self._Dmax

                                    # Run-up Regime
                                    if inundation == 0:
                                        if Q1 > self._Qs_min and S1 >= 0:
                                            Qs1 = self._Kr * Q1
                                            if Qs1 > fluxLimit:
                                                Qs1 = fluxLimit
                                        else:
                                            Qs1 = 0

                                        if Q2 > self._Qs_min and S2 >= 0:
                                            Qs2 = self._Kr * Q2
                                            if Qs2 > fluxLimit:
                                                Qs2 = fluxLimit
                                        else:
                                            Qs2 = 0

                                        if Q3 > self._Qs_min and S3 >= 0:
                                            Qs3 = self._Kr * Q3
                                            if Qs3 > fluxLimit:
                                                Qs3 = fluxLimit
                                        else:
                                            Qs3 = 0

                                    # Inundation Regime - Murray and Paola (1994, 1997) Rule 3 with flux limiter
                                    else:
                                        if Q1 > self._Qs_min:
                                            Qs1 = self._Ki * (Q1 * (S1 + C)) ** self._mm
                                            if Qs1 < 0:
                                                Qs1 = 0
                                            elif Qs1 > fluxLimit:
                                                Qs1 = fluxLimit
                                        else:
                                            Qs1 = 0

                                        if Q2 > self._Qs_min:
                                            Qs2 = self._Ki * (Q2 * (S2 + C)) ** self._mm
                                            if Qs2 < 0:
                                                Qs2 = 0
                                            elif Qs2 > fluxLimit:
                                                Qs2 = fluxLimit
                                        else:
                                            Qs2 = 0

                                        if Q3 > self._Qs_min:
                                            Qs3 = self._Ki * (Q3 * (S3 + C)) ** self._mm
                                            if Qs3 < 0:
                                                Qs3 = 0
                                            elif Qs3 > fluxLimit:
                                                Qs3 = fluxLimit
                                        else:
                                            Qs3 = 0

                                    Qs1 = np.nan_to_num(Qs1)
                                    Qs2 = np.nan_to_num(Qs2)
                                    Qs3 = np.nan_to_num(Qs3)

                                    # ### Calculate Net Erosion/Accretion
                                    # KA: note, Ian added bandaid by "or" statement here on 9/28/20
                                    if (
                                            Elevation[TS, d, i] > self._SL or any(z > self._SL for z in Elevation[TS, d + 1:d + 10, i])
                                    ):  # If cell is subaerial, elevation change is determined by difference between
                                        # flux in vs. flux out
                                        if i > 0:
                                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                                        SedFluxIn[TS, d + 1, i] += Qs2

                                        if i < (self._BarrierLength - 1):
                                            SedFluxIn[TS, d + 1, i + 1] += Qs3

                                        Qs_out = Qs1 + Qs2 + Qs3
                                        SedFluxOut[TS, d, i] = Qs_out

                                    else:  # If cell is subaqeous, exponentially decay dep. of remaining sed across bay

                                        if inundation == 0:
                                            Cbb = self._Cbb_r
                                        else:
                                            Cbb = self._Cbb_i

                                        Qs0 = SedFluxIn[TS, d, i] * Cbb

                                        Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                        Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                        Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                        Qs1 = np.nan_to_num(Qs1)
                                        Qs2 = np.nan_to_num(Qs2)
                                        Qs3 = np.nan_to_num(Qs3)

                                        if Qs1 < self._Qs_bb_min:
                                            Qs1 = 0
                                        elif Qs1 > fluxLimit:
                                            Qs1 = fluxLimit
                                        if Qs2 < self._Qs_bb_min:
                                            Qs2 = 0
                                        elif Qs2 > fluxLimit:
                                            Qs2 = fluxLimit
                                        if Qs3 < self._Qs_bb_min:
                                            Qs3 = 0
                                        elif Qs3 > fluxLimit:
                                            Qs3 = fluxLimit

                                        if i > 0:
                                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                                        SedFluxIn[TS, d + 1, i] += Qs2

                                        if i < (self._BarrierLength - 1):
                                            SedFluxIn[TS, d + 1, i + 1] += Qs3

                                        Qs_out = Qs1 + Qs2 + Qs3
                                        SedFluxOut[TS, d, i] = Qs_out

                                    # ### Saline Flooding
                                    (
                                        self._ShrubDomainFemale,
                                        self._ShrubDomainMale,
                                        self._ShrubDomainDead,
                                    ) = self.SalineFlooding(
                                        ShrubDomainWidth,
                                        self._ShrubDomainAll,
                                        self._ShrubDomainFemale,
                                        self._ShrubDomainMale,
                                        self._ShrubDomainDead,
                                        d,
                                        i,
                                        Q0,
                                    )

                        # ### Update Elevation After Every Storm Hour
                        if inundation == 1:
                            ElevationChange = (
                                                      SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                                              ) / substep
                        else:
                            ElevationChange = (
                                                      SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                                              ) / substep
                        Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange

                        # Calculate and save volume of sediment deposited on/behind the island for every hour
                        # (all four methods below should equal the same!)
                        # OWloss = OWloss + np.sum(ElevationChange[1:,:])
                        # OWloss = OWloss + (np.sum(SedFluxIn[TS,1:,:]) - np.sum(SedFluxOut[TS,1:,:])) / substep
                        # OWloss = OWloss + np.sum(SedFluxIn[TS,1,:]) / substep
                        OWloss = OWloss + np.sum(SedFluxOut[TS, 0, :]) / substep

                        # Update amount of burial/erosion for each shrub
                        self._BurialDomain = self.UpdateBurial(
                            self._BurialDomain,
                            ElevationChange,
                            ShrubDomainWidth,
                            self._ShrubDomainAll,
                        )

                    # ### Update Interior Domain After Every Storm
                    InteriorUpdate = Elevation[-1, 1:, :]

                    # Remove all rows of bay without any deposition from the domain
                    check = 1
                    while check == 1:
                        if all(x <= -self._BayDepth for x in InteriorUpdate[-1, :]):
                            InteriorUpdate = np.delete(InteriorUpdate, (-1), axis=0)
                        else:
                            check = 0

                    # Update interior domain
                    self._InteriorDomain = InteriorUpdate

                    # Update Domain widths
                    DomainWidth = np.shape(self._InteriorDomain)[0]
                    (
                        self._ShrubDomainFemale,
                        self._ShrubDomainMale,
                        self._ShrubDomainAll,
                        self._ShrubPercentCover,
                        self._BurialDomain,
                        self._ShrubDomainDead,
                        self._DeadPercentCover,
                    ) = self.UpdateShrubDomains(
                        DomainWidth,
                        ShrubDomainWidth,
                        self._ShrubDomainFemale,
                        self._ShrubDomainMale,
                        self._ShrubDomainAll,
                        self._ShrubPercentCover,
                        self._BurialDomain,
                        self._ShrubDomainDead,
                        self._DeadPercentCover
                    )

        # Record storm data
        self._StormCount.append(numstorm)

        # ###########################################
        # ### Ocean Shoreline Change

        # ### Calculate shoreline/shoreface position change following Lorenzo-Trueba and Ashton (2014)
        (
            self._x_s,
            self._x_t,
            self._x_s_TS,
            self._x_t_TS,
            self._x_b_TS,
            self._s_sf_TS,
            self._QowTS,
            self._QsfTS,
            self._h_b_TS,
        ) = self.LTA_SC(
            self._InteriorDomain,
            OWloss,
            Qdg,
            DuneLoss,
            self._x_s,
            self._x_s_TS,
            self._x_t,
            self._x_t_TS,
            self._x_b_TS,
            self._s_sf_TS,
            InteriorWidth_Avg,
            self._SL,
            self._QowTS,
            self._QsfTS,
            self._time_index,
            self._h_b_TS,
        )

        # NOTE FROM KA: this may break now that I've moved the domain modifications to the beginning

        self._time_index += 1

    @property
    def time_index(self):
        return self._time_index

    @property
    def QowTS(self):
        return self._QowTS

    @property
    def QsfTS(self):
        return self._QsfTS

    @property
    def x_t_TS(self):
        return self._x_t_TS

    @property
    def x_s_TS(self):
        return self._x_s_TS

    @property
    def x_s(self):
        return self._x_s

    @property
    def x_b_TS(self):
        return self._x_b_TS

    @property
    def h_b_TS(self):
        return self._h_b_TS
