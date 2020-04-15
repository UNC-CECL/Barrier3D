import math

import numpy as np

from .parameters import load_params


class Barrier3D:
    def SeaLevel(self, InteriorDomain, DuneDomain, t):
        """Accounts for relative sea level rise, returns updated elevation domains."""

        # from Barrier3D_Parameters import (RSLR, BayDepth)

        # Decrease all elevation this year by RSLR increment
        InteriorDomain = InteriorDomain - self._RSLR
        DuneDomain[t - 1] = DuneDomain[t - 1] - self._RSLR
        InteriorDomain[
            InteriorDomain < -self._BayDepth
        ] = (
            -self._BayDepth
        )  # Bay can't be deeper than BayDepth (assumes equilibrium depth maintained)

        return InteriorDomain, DuneDomain

    def FindWidths(self, InteriorDomain, SL):
        """Finds DomainWidth and InteriorWidth

        DomainWidth is wide enough to fit widest part of island; InteriorWidth
        stores the width of the island at each row of the domain where elevation
        is greater than 0      
        """
        # from Barrier3D_Parameters import (BarrierLength)

        DomainWidth = np.shape(InteriorDomain)[0]  # (dam) analysis:ignore
        InteriorWidth = [0] * self._BarrierLength
        for l in range(self._BarrierLength):
            width = next(
                (
                    index
                    for index, value in enumerate(InteriorDomain[:, l])
                    if value <= SL
                ),
                DomainWidth,
            )
            if width < 0:
                width = 0
            InteriorWidth[l] = width - 1
        # Average width of island
        InteriorWidth_Avg = np.nanmean(InteriorWidth)

        return DomainWidth, InteriorWidth, InteriorWidth_Avg

    def DuneGrowth(self, DuneDomain, t):
        """Grows dune cells"""

        # from Barrier3D_Parameters import (Dmaxel, BermEl, growthparam, DuneWidth, BarrierLength)

        # Set max dune height - depends on beach width according to relationship gathered from VCR data
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
        t,
        ShrubDomainFemale,
        ShrubDomainMale,
        BurialDomain,
        InteriorWidth_Avg,
        DomainWidth,
        ShrubDomainDead,
    ):
        """Main shrub expansion and mortality function"""

        # from Barrier3D_Parameters import (BarrierLength, Dshrub, BermEl, Female, ShrubEl_min, ShrubEl_max, BurialLimit, UprootLimit, TimeFruit, Seedmin, Seedmax, GermRate, disp_mu, disp_sigma)

        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale

        ### Burial
        # Kill all shrubs buried by over depth limit or eroded
        ShrubDomainDead[BurialDomain > self._BurialLimit] = ShrubDomainAll[
            BurialDomain > self._BurialLimit
        ]  # Transfer to dead domain
        ShrubDomainDead[
            BurialDomain >= 0.2
        ] = 0  # If buried by more than 2 m, considered shrub to be completely gone
        ShrubDomainDead[BurialDomain < self._UprootLimit] = ShrubDomainAll[
            BurialDomain < self._UprootLimit
        ]  # Transfer to dead domain

        ShrubDomainFemale[BurialDomain > self._BurialLimit] = 0
        ShrubDomainFemale[BurialDomain < self._UprootLimit] = 0
        ShrubDomainMale[BurialDomain > self._BurialLimit] = 0
        ShrubDomainMale[BurialDomain < self._UprootLimit] = 0

        BurialDomain[BurialDomain > self._BurialLimit] = 0  # Reset burial domain
        BurialDomain[BurialDomain < self._UprootLimit] = 0  # Reset burial domain
        BurialDomain[BurialDomain >= 0.2] = 0  # Reset burial domain
        ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale  # Recalculate

        ### Inundation
        # Remove shrubs that have fallen below Mean High Water (i.e. elevation of 0) (passively inundated by rising back-barrier water elevations)

        TideRange = 0.12  # dam, half range
        for w in range(DomainWidth):
            for l in range(self._BarrierLength):
                if InteriorDomain[w, l] < 0 and ShrubDomainAll[w, l] > 0:
                    if InteriorDomain[w, l] > -TideRange:
                        ShrubDomainDead[w, l] = ShrubDomainAll[
                            w, l
                        ]  # Shrub remains as dead if within tidal range (i.e. MHW - tidal range)
                    ShrubDomainFemale[w, l] = 0
                    ShrubDomainMale[w, l] = 0

        ### Age the shrubs in years
        ShrubDomainFemale[ShrubDomainFemale > 0] += 1
        ShrubDomainMale[ShrubDomainMale > 0] += 1

        ### Randomly drop a seed onto the island each time step
        randX = np.random.randint(0, self._BarrierLength)
        randY = np.random.randint(0, InteriorWidth_Avg)
        if (
            ShrubDomainFemale[randY, randX] == 0
            and ShrubDomainMale[randY, randX] == 0
            and DuneDomainCrest[randX] + self._BermEl >= self._Dshrub
            and InteriorDomain[randY, randX] >= self._ShrubEl_min
            and InteriorDomain[randY, randX] <= self._ShrubEl_max
        ):
            if random.random() > self_Female:
                ShrubDomainFemale[randY, randX] = 1
            else:
                ShrubDomainMale[randY, randX] = 1

        ### Disperse seeds
        for k in range(
            self._BarrierLength
        ):  # Loop through each row of island width (i.e. from ocean to mainland side of island)
            if 0 in ShrubDomainAll:  # <--------------- Working??

                # For all cells with a shrub
                FemaleShrubs = ShrubDomainFemale[:, k]
                I = [
                    index
                    for index, value in enumerate(FemaleShrubs)
                    if value >= self._TimeFruit
                ]
                numShrubCells = len(I)
                # Determine how many seeds in each cell
                Seedspercell = np.random.randint(
                    self._Seedmin, high=self._Seedmax, size=numShrubCells
                )
                # For each shrub cell, determine survival rate for the cell in this year
                SeedSurvYear = self._GermRate * np.random.rand(numShrubCells)

                for i in range(numShrubCells):
                    # For each shrub cell producing seeds, generate a random # of random numbers each representing a single seed
                    randseeds = np.random.rand(Seedspercell[i])
                    # Find how many seeds produced in each cell actually survive
                    Survivors = len(randseeds[randseeds < SeedSurvYear[i]])

                    # Determine distance, rounding to nearest integer
                    if Survivors > 0:
                        DispDist = np.round(
                            np.random.lognormal(
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

                                # Randomly selct one of those points as the target point - this means equal probability of dispersal in every direction (valid assumption for avian dispersal)
                                if len(col) == 0:
                                    targetY = originX
                                    targetX = originY
                                else:
                                    xx = random.randint(0, (len(col)) - 1)
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
                                    targetY >= 0
                                    and targetY < DomainWidth
                                    and targetX >= 0
                                    and targetX < self._BarrierLength
                                    and ShrubDomainFemale[targetY, targetX] == 0
                                    and ShrubDomainMale[targetY, targetX] == 0
                                    and DuneDomainCrest[targetX] + self._BermEl
                                    >= self._Dshrub
                                    and InteriorDomain[targetY, targetX]
                                    >= self._ShrubEl_min
                                    and InteriorDomain[targetY, targetX]
                                    <= self._ShrubEl_max
                                    and ShrubDomainDead[targetY, targetX] < 1
                                ):
                                    # Decide if the tree wll be a female or male
                                    if random.random() > self._Female:
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

    def WaterLevels(self, numstorm, MHW):
        """Returns statistics for a given storm"""

        # from Barrier3D_Parameters import (surge_tide_m, surge_tide_sd, height_mu, height_sigma, period_m, period_sd, beta, duration_mu, duration_sigma)

        ### Generate storm water levels
        Rhigh = []
        Rlow = []
        period = []
        duration = []

        for n in range(numstorm):
            ### Draw statistics for this storm from normal distributions

            surge_tide = np.random.normal(self._surge_tide_m, self._surge_tide_sd)
            if surge_tide < 0:
                surge_tide = 0.01

            height = np.random.lognormal(self._height_mu, self._height_sigma)
            if height < 0:
                height = 0.01

            T = np.random.normal(self._period_m, self._period_sd)
            if T < 0:
                T = 0.01
            L0 = (9.8 * T ** 2) / (2 * math.pi)  # Wavelength
            if L0 < 0:
                L0 = 0.01

            dur = round(np.random.lognormal(self._duration_mu, self._duration_sigma))
            if dur < 8:
                dur = 8

            ### Calculate swash runup

            # Setup
            setup = 0.35 * self._beta * math.sqrt(height * L0)

            # Incident band swash
            Sin = 0.75 * self._beta * math.sqrt(height * L0)

            # Infragravity band swash
            Sig = 0.06 * math.sqrt(height * L0)

            # Swash
            swash = math.sqrt((Sin ** 2) + (Sig ** 2))

            # Rhigh (m NAVD88)
            R2 = 1.1 * (setup + (swash / 2))
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

    def DuneGaps(self, DuneDomain, Dow, bermel, Rhigh):
        """Returns tuple of [gap start index, stop index, avg Rexcess of each gap, alpha: ratio of TWL to dune height]"""
        gaps = []
        start = 0
        stop = 0
        i = start
        while i < (len(Dow) - 1):
            adjacent = Dow[i + 1] - Dow[i]
            if adjacent == 1:
                i = i + 1
            else:
                stop = i

                x = DuneDomain[Dow[start] : (Dow[stop] + 1)]
                Hmean = sum(x) / float(len(x))
                Rexcess = Rhigh - (Hmean + bermel)
                alpha = Rhigh / (Hmean + bermel)
                gaps.append([Dow[start], Dow[stop], Rexcess, alpha])

                start = stop + 1
                i = start
        if i > 0:
            stop = i - 1

            x = DuneDomain[Dow[start] : (Dow[stop] + 1)]
            if len(x) > 0:
                Hmean = sum(x) / float(len(x))
                Rexcess = Rhigh - (Hmean + bermel)
                alpha = Rhigh / (Hmean + bermel)
                gaps.append([Dow[start], Dow[stop], Rexcess, alpha])
        return gaps

    def DiffuseDunes(self, DuneDomain, t):
        """Dune height diffusion in alongshore direction: smoothes dune height by redistributing sand from high dune to neighboring low dune(s)"""

        # from Barrier3D_Parameters import (HdDiffu, BarrierLength, DuneWidth)

        for w in range(self._DuneWidth):
            # Alongshore direction
            for d in range(2, self._BarrierLength):  # Loop L to R
                Ldiff = (
                    DuneDomain[t, d, w] - DuneDomain[t, d - 1, w]
                )  # Calculate height difference
                # Subtract sand from taller dune cell and add to adjacent shorter one (assumes sand within dunes is conserved)
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
        """Kill all immature (< 1 yr-old) shrubs that have been flooded beyond a threshold discharge (Tolliver et al., 1997)"""

        # from Barrier3D_Parameters import (BarrierLength, SalineLimit)

        if d < (ShrubDomainWidth - 1) and i < (self._BarrierLength - 1):
            if Q0 >= self._SalineLimit and ShrubDomainAll[d, i] == 1:
                ShrubDomainDead[d, i] = (
                    ShrubDomainMale[d, i] + ShrubDomainFemale[d, i]
                )  # Transfer to dead shrub domain
                ShrubDomainFemale[d, i] = 0
                ShrubDomainMale[d, i] = 0

        return ShrubDomainFemale, ShrubDomainMale, ShrubDomainDead

    def UpdateBurial(
        self, BurialDomain, ElevationChange, ShrubDomainWidth, ShrubDomainAll
    ):
        """Updates amount of burial/erosion for each shrub"""

        BurialDomain = BurialDomain + ElevationChange[1 : ShrubDomainWidth + 1, :]
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
    ):
        """Updates size of shrub domains"""

        # from Barrier3D_Parameters import (BarrierLength)

        if DomainWidth > ShrubDomainWidth:
            AddRows = np.zeros([DomainWidth - ShrubDomainWidth, self._BarrierLength])
            ShrubDomainFemale = np.vstack([ShrubDomainFemale, AddRows])
            ShrubDomainMale = np.vstack([ShrubDomainMale, AddRows])
            ShrubDomainDead = np.vstack([ShrubDomainDead, AddRows])
            ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
            ShrubPercentCover = np.vstack([ShrubPercentCover, AddRows])
            BurialDomain = np.vstack([BurialDomain, AddRows])
        elif DomainWidth < ShrubDomainWidth:
            RemoveRows = ShrubDomainWidth - DomainWidth
            ShrubDomainFemale = ShrubDomainFemale[0:-RemoveRows, :]
            ShrubDomainMale = ShrubDomainMale[0:-RemoveRows, :]
            ShrubDomainDead = ShrubDomainDead[0:-RemoveRows, :]
            ShrubDomainAll = ShrubDomainFemale + ShrubDomainMale
            ShrubPercentCover = ShrubPercentCover[0:-RemoveRows, :]
            BurialDomain = BurialDomain[0:-RemoveRows, :]

        return (
            ShrubDomainFemale,
            ShrubDomainMale,
            ShrubDomainAll,
            ShrubPercentCover,
            BurialDomain,
            ShrubDomainDead,
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
    ):
        """Finds shoreline change for modeled year following Lorenzo-Trueba and Ashton (2014)"""

        # from Barrier3D_Parameters import (BarrierLength, DShoreface, k_sf, s_sf_eq, RSLR, Qat)

        # Find volume of shoreface/beach/dune sand deposited in island interior and back-barrier
        Qow = (OWloss) / (
            self._BarrierLength
        )  # (dam^3/dam) Volume of sediment lost from shoreface/beach by overwash
        QowTS.append(Qow * 100)  # Save in m^3/m
        if DuneLoss < Qow:
            Qow = Qow - DuneLoss  # Account for dune contribution to overwash volume
        else:
            Qow = 0  # Excess DuneLoss assumed lost offshore

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
            2 * self._RSLR / s_sf
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

        SCR = (
            x_s_dt * -1
        )  # (dam) Positive increase in x-location means shoreline erosion

        return SCR, x_s, x_t, x_s_TS, x_t_TS, x_b_TS, s_sf_TS, QowTS, QsfTS

    def CalcPC(self, ShrubDomainAll, PercentCoverTS, ShrubDomainDead, ShrubArea, t):
        """Calculates percent cover of shrub domain and shrub coverage area (dam^2)"""

        # from Barrier3D_Parameters import PC

        Allshrub_t = ShrubDomainAll.astype("int64")
        Deadshrub_t = ShrubDomainDead.astype("int64")
        ShrubPercentCover = self._PC.take(Allshrub_t)
        DeadPercentCover = self._PC.take(Deadshrub_t)

        #    ShrubPercentCover[ShrubPercentCover == 0] = DeadPercentCover[ShrubPercentCover == 0]
        PercentCoverTS[t] = ShrubPercentCover
        ShrubArea.append(np.count_nonzero(ShrubDomainAll))

        return ShrubPercentCover, PercentCoverTS, ShrubArea

    def __init__(self, **kwds):
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
        self._duration_mu = kwds.pop("duration_mu")
        self._duration_sigma = kwds.pop("duration_sigma")
        self._Female = kwds.pop("Female")
        self._GermRate = kwds.pop("GermRate")
        self._growthparam = kwds.pop("growthparam")
        self._HdDiffu = kwds.pop("HdDiffu")
        self._height_mu = kwds.pop("height_mu")
        self._height_sigma = kwds.pop("height_sigma")
        self._InteriorDomain = kwds.pop("InteriorDomain")
        self._k_sf = kwds.pop("k_sf")
        self._Ki = kwds.pop("Ki")
        self._Kr = kwds.pop("Kr")
        LShoreface = kwds.pop("LShoreface")
        self._mean_storm = kwds.pop("mean_storm")
        self._mm = kwds.pop("mm")
        self._MHW = kwds.pop("MHW")
        self._nn = kwds.pop("nn")
        self._numstorm = kwds.pop("numstorm")
        self.PC = kwds.pop("PC")
        self._period_m = kwds.pop("period_m")
        self._period_sd = kwds.pop("period_sd")
        self._Qat = kwds.pop("Qat")
        self._Qs_min = kwds.pop("Qs_min")
        self._Qshrub_max = kwds.pop("Qshrub_max")
        self._Rin_i = kwds.pop("Rin_i")
        self._Rin_r = kwds.pop("Rin_r")
        self._RSLR = kwds.pop("RSLR")
        self._s_sf_eq = kwds.pop("s_sf_eq")
        self._SalineLimit = kwds.pop("SalineLimit")
        self._SD_storm = kwds.pop("SD_storm")
        self._Seedmax = kwds.pop("Seedmax")
        self._Seedmin = kwds.pop("Seedmin")
        self._Shrub_ON = kwds.pop("Shrub_ON")
        self._ShrubEl_max = kwds.pop("ShrubEL_max")
        self._ShrubEl_min = kwds.pop("ShrubEL_min")
        self._StormStart = kwds.pop("StormStart")
        self._surge_tide_m = kwds.pop("surge_tide_m")
        self._surge_tide_sd = kwds.pop("surge_tide_sd")
        self._threshold_in = kwds.pop("threshold_in")
        self._TimeFruit = kwds.pop("TimeFruit")
        self._TMAX = kwds.pop("TMAX")
        self._UprootLimit = kwds.pop("UprootLimit")

        ### Initialize Shrubs
        # if Shrub_ON ==1: print('Shrubs ON')

        self._PercentCoverTS = [None] * self._TMAX
        self._PercentCoverTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubFemaleTS = [None] * self._TMAX
        self._ShrubFemaleTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubMaleTS = [None] * self._TMAX
        self._ShrubMaleTS[0] = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubDomainFemale = self._ShrubFemaleTS[0]
        self._ShrubDomainMale = self._ShrubMaleTS[0]
        self._ShrubPercentCover = self._PercentCoverTS[0]
        self._BurialDomain = np.zeros([self._DomainWidth, self._BarrierLength])
        self._ShrubArea = [0]

        ### Initialize variables
        self._DomainTS = [
            None
        ] * self._TMAX  # Stores the elevation domain for each timestep
        self._DomainTS[0] = self._InteriorDomain
        # NOTE: SCR may not need to be initialized here
        SCR = 0  # (dam/yr) Change rate of ocean-fronting shoreline: (+) = progradation, (-) = erosion
        self._SCRagg = 0  # Counter variable for shoreline change that is less than 1 cell size per year
        self._ShorelineChangeTS = [0]
        self._ShorelineChange = 0  # (dam) Total amount of shoreline change
        # NOTE: BBShorelineChangeTS and BBShorelineChange are unused
        BBShorelineChangeTS = [0]
        BBShorelineChange = 0  # (dam) Total amount of back-barrier shoreline extension
        self._StormCount = []
        self._InundationCount = 0
        self._RunUpCount = 0
        self._InteriorWidth_AvgTS = [self._DomainWidth]
        self._QowTS = []  # (m^3/m)
        self._x_t = 0  # (dam) Start location of shoreface toe
        self._x_s = self._x_t + LShoreface  # (dam) Start location of shoreline
        self._x_t_TS = [self._x_t]  # (dam) Shoreface toe locations for each time step
        self._x_s_TS = [self._x_s]  # (dam) Shoreline locations for each time step
        self._x_b_TS = [
            (self._x_s + self._DomainWidth)
        ]  # (dam) Bay shoreline locations for each time step
        # NOTE: h_b_TS is unused
        h_b_TS = [self._BermEl]  # (dam) Berm Elevation over time
        self._s_sf_TS = [(self._DShoreface / LShoreface)]
        self._Hd_AverageTS = [Dstart]
        self._QsfTS = [0]  # (m^3/m)

        self._time_index = 1

    @classmethod
    def from_xls(cls, path_to_xls):
        return cls(**load_params(path_to_xls))

    def update(self):
        ### RSLR
        self._InteriorDomain, self._DuneDomain = self.SeaLevel(
            self._InteriorDomain, self._DuneDomain, self._time_index
        )

        ### Find DomainWidth and InteriorWidth
        self._DomainWidth, InteriorWidth, InteriorWidth_Avg = self.FindWidths(
            self._InteriorDomain, self._SL
        )

        ### Check for drowning
        # Barrier drowns if interior width thins to 10 m or fewer or all interior cells are below SL
        if InteriorWidth_Avg < 1:  # if max(InteriorWidth) <= 1:
            print(
                "Barrier has WIDTH DROWNED at t = " + str(self._time_index) + " years!"
            )
            self._TMAX = self._time_index - 1
            return
        if all(j <= self._SL for j in self._InteriorDomain[0, :]):
            print(
                "Barrier has HEIGHT DROWNED at t = " + str(self._time_index) + " years!"
            )
            self._TMAX = self._time_index - 1
            return

        ### Grow Dunes
        self._DuneDomain, Dmax, Qdg = self.DuneGrowth(
            self._DuneDomain, self._time_index
        )

        DuneDomainCrest = self._DuneDomain[self._time_index, :, :].max(
            axis=1
        )  # Maximum height of each row in DuneDomain
        self._Hd_AverageTS.append(
            np.mean(DuneDomainCrest)
        )  # Store average pre-storms dune-heigh for time step

        ###########################################
        ### Shrubs

        if self._Shrub_ON == 1:

            (
                ShrubDomainAll,
                self._ShrubDomainFemale,
                self._ShrubDomainMale,
                self._BurialDomain,
            ) = self.Shrubs(
                self._InteriorDomain,
                DuneDomainCrest,
                self._time_index,
                self._ShrubDomainFemale,
                self._ShrubDomainMale,
                self._BurialDomain,
                InteriorWidth_Avg,
                self._DomainWidth,
            )

        ###########################################
        ### Storms

        OWloss = 0

        if self._time_index >= self._StormStart:
            # Select number of storms for this time step from normal distribution
            self._numstorm = int(self._numstorm)
            self._numstorm = round(
                np.random.normal(self._mean_storm, self._SD_storm)
            )  # analysis:ignore # Comment out this line if using pre-specified static number of storms per year

            if self._numstorm > 0:
                # Draw statistics for each storm
                Rhigh, Rlow, period, dur = self.WaterLevels(self._numstorm, self._MHW)

                ### Individual Storm Impacts
                for n in range(self._numstorm):  # Loop through each individual storm

                    ###########################################
                    ### Dune Erosion

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
                            )  # Amount of dune crest elevation change normalized by pre-storm dune elevation (i.e. a percent change), from Goldstein and Moore (2016)

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

                            # If dune is lowered to essentially zero, allow for chance of regrowth by raising dune height to 5 cm
                            if (
                                self._DuneDomain[self._time_index, Dow[d], w]
                                < self._DuneRestart
                            ):
                                if self._DuneRestart < Dmax:
                                    self._DuneDomain[
                                        self._time_index, Dow[d], w
                                    ] = self._DuneRestart
                                else:
                                    self._DuneDomain[
                                        self._time_index, Dow[d], w
                                    ] = Dmax  # Restart height can't be greater than Dmax

                    # Dune Height Diffusion
                    self._DuneDomain = self.DiffuseDunes(
                        self._DuneDomain, self._time_index
                    )

                    Hd_TSloss = (
                        DuneChange.max(axis=1) / dur[n]
                    )  # Average height of dune loss for each substep during storm

                    ###########################################
                    ### Overwash

                    Iow = 0  # Count of dune gaps in inundation regime
                    Dunes_prestorm = DuneDomainCrest
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2]  # (m)
                        gapwidth = stop - start + 1
                        meandune = (
                            sum(Dunes_prestorm[start : stop + 1]) / gapwidth
                        ) + self._BermEl  # Average elevation of dune gap

                        # Determine number of gaps in inundation regime
                        if Rlow[n] > meandune:
                            Iow += 1

                    # Dertermine Sediment And Water Routing Rules
                    if (
                        len(gaps) > 0 and Iow / len(gaps) >= self._threshold_in
                    ):  # If greater than 25% of dune gaps are in unundation regime, use inundation regime routing
                        inundation = 1
                        #                        #### Calc substep
                        #                        Smax = 0
                        #                        Smin = 0
                        #                        for l in range(self._BarrierLength):
                        #                            for w in range(InteriorWidth[l]):
                        #                                slope = self._InteriorDomain[w,l] - self._InteriorDomain[w+1,l]
                        #                                if slope > Smax:
                        #                                    Smax = slope
                        #                        calc maxQ
                        #                        accrete = f(maxQ, smax)
                        #                        ceil(accrete / (Smax *2)) = substep
                        #                        ####
                        substep = 15  # pre 14Jan20: 5, ~20 needed?
                        self._InundationCount += 1
                    else:
                        inundation = 0
                        surges = (1 / period[n]) * 60 * 60  # Number of surges per hour
                        substep = 1
                        self._RunUpCount += 1

                    # Set Domain
                    add = 25
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
                    inundation = 0

                    # Set Water and Sediment at Dune Crest
                    for q in range(len(gaps)):
                        start = gaps[q][0]
                        stop = gaps[q][1]
                        Rexcess = gaps[q][2]  # (m)
                        gapwidth = stop - start + 1
                        meandune = (
                            sum(Dunes_prestorm[start : stop + 1]) / gapwidth
                        ) + self._BermEl  # Average elevation of dune gap

                        # Calculate discharge through each dune cell
                        Vdune = math.sqrt(2 * 9.8 * (Rexcess * 10)) / 10  # (dam/surge)
                        Qdune = Vdune * Rexcess * 3600  # (dam^3/hr)

                        # Set discharge at dune gap
                        Discharge[:, 0, start:stop] = Qdune

                        if inundation == 1:  # Inundation regime
                            Rin = self._Rin_i
                            slopes = []
                            # Find average slope of interior
                            for u in range(1, width - 1):
                                for v in range(self._BarrierLength):
                                    if (
                                        Elevation[0, u, v] > self._SL
                                        and Elevation[0, u + 1, v] > self._SL
                                    ):
                                        SS = Elevation[0, u, v] - Elevation[0, u + 1, v]
                                        slopes.append(SS)
                            AvgSlope = sum(slopes) / len(slopes)
                            C = 2 * AvgSlope
                            # Set sediment at dune gap
                            SedFluxIn[:, 0, start:stop] = (
                                self._Ki * (Qdune * (AvgSlope + C)) ** self._mm
                            )
                        else:  # Run-up regime
                            Rin = self._Rin_r
                            # Set sediment at dune gap
                            SedFluxIn[:, 0, start:stop] = self._Kr * Qdune

                    ### Run Flow Routing Algorithm
                    for TS in range(duration):

                        ShrubDomainWidth = np.shape(self._ShrubDomainFemale)[0]

                        if TS > 0:
                            Elevation[TS, 1:, :] = Elevation[
                                TS - 1, 1:, :
                            ]  # Begin timestep with elevation from end of last
                            Elevation[TS, 0, :] = (
                                Dunes - Hd_TSloss * TS
                            )  # Reduce dune in height linearly over course of storm

                        for d in range(width - 1):
                            # Reduce discharge across row via infiltration
                            if d > 0:
                                Discharge[TS, d, :][Discharge[TS, d, :] > 0] -= Rin
                            Discharge[TS, d, :][Discharge[TS, d, :] < 0] = 0

                            for i in range(self._BarrierLength):
                                if Discharge[TS, d, i] > 0:

                                    Q0 = Discharge[TS, d, i]

                                    ### Calculate Slopes
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

                                    ### Calculate Discharge To Downflow Neighbors
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

                                    ### Reduce Overwash Through Shrub Cells and Save Discharge
                                    if self._Shrub_ON == 1:
                                        # Cell 1
                                        if i > 0:
                                            if d < (ShrubDomainWidth):
                                                if (
                                                    self._ShrubPercentCover[d, i - 1]
                                                    > 0
                                                ):
                                                    Q1 = Q1 * (
                                                        self._Qshrub_max
                                                        * self._ShrubPercentCover[
                                                            d, i - 1
                                                        ]
                                                    )
                                            Discharge[TS, d + 1, i - 1] = (
                                                Discharge[TS, d + 1, i - 1] + Q1
                                            )

                                        # Cell 2
                                        if d < (ShrubDomainWidth - 1):
                                            if self._ShrubPercentCover[d, i] > 0:
                                                Q2 = Q2 * (
                                                    self._Qshrub_max
                                                    * self._ShrubPercentCover[d, i]
                                                )
                                        Discharge[TS, d + 1, i] = (
                                            Discharge[TS, d + 1, i] + Q2
                                        )

                                        # Cell 3
                                        if i < (self._BarrierLength - 1):
                                            if d < (ShrubDomainWidth):
                                                if (
                                                    self._ShrubPercentCover[d, i + 1]
                                                    > 0
                                                ):
                                                    Q3 = Q3 * (
                                                        self._Qshrub_max
                                                        * self._ShrubPercentCover[
                                                            d, i + 1
                                                        ]
                                                    )
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

                                    ### Calculate Sed Movement
                                    # Run-up regime
                                    if inundation == 0:
                                        if Q1 > self._Qs_min and S1 >= 0:
                                            Qs1 = self._Kr * Q1
                                        else:
                                            Qs1 = 0

                                        if Q2 > self._Qs_min and S2 >= 0:
                                            Qs2 = self._Kr * Q2
                                        else:
                                            Qs2 = 0

                                        if Q3 > self._Qs_min and S3 >= 0:
                                            Qs3 = self._Kr * Q3
                                        else:
                                            Qs3 = 0

                                    # Inundation Regime (Murray and Paola (1994, 1997) Rule 3)
                                    else:
                                        if Q1 > self._Qs_min:
                                            Qs1 = self._Ki * (Q1 * (S1 + C)) ** self._mm
                                        else:
                                            Qs1 = 0

                                        if Q2 > self._Qs_min:
                                            Qs2 = self._Ki * (Q2 * (S2 + C)) ** self._mm
                                        else:
                                            Qs2 = 0

                                        if Q3 > self._Qs_min:
                                            Qs3 = self._Ki * (Q3 * (S3 + C)) ** self._mm
                                        else:
                                            Qs3 = 0

                                    #                                    # Inundation Regime (Murray and Paola (1994, 1997) Rule 5)
                                    #                                    else:
                                    #                                        Qt = np.mean(Discharge[0,0,:]) # Typical discharge
                                    #                                        St = AvgSlope # Typical slope
                                    #                                        ThD = 2 # Threshold denomenator, threshold therefore around half of the typical stream power
                                    #                                        Th = Qt * St / ThD
                                    #                                        if Q1 > self._Qs_min:
                                    #                                            Qs1 = self._Ki * (Q1 * (S1 + C) - Th)**self._mm
                                    #                                        else:
                                    #                                            Qs1 = 0
                                    #
                                    #                                        if Q2 > self._Qs_min:
                                    #                                            Qs2 = self._Ki * (Q2 * (S2 + C) - Th)**self._mm
                                    #                                        else:
                                    #                                            Qs2 = 0
                                    #
                                    #                                        if Q3 > self._Qs_min:
                                    #                                            Qs3 = self._Ki * (Q3 * (S3 + C) - Th)**self._mm
                                    #                                        else:
                                    #                                            Qs3 = 0

                                    ### Calculate Net Erosion/Accretion
                                    if (
                                        Elevation[TS, d, i] > self._SL
                                    ):  # If cell is subaerial, elevation change is determined by difference between flux in vs. flux out
                                        if i > 0:
                                            SedFluxIn[TS, d + 1, i - 1] += Qs1

                                        SedFluxIn[TS, d + 1, i] += Qs2

                                        if i < (self._BarrierLength - 1):
                                            SedFluxIn[TS, d + 1, i + 1] += Qs3

                                        Qs_out = Qs1 + Qs2 + Qs3
                                        SedFluxOut[TS, d, i] = Qs_out

                                    else:  # If cell is subaqeous, exponentially decay deposition of remaining sediment across bay
                                        Qs_bb_min = 0.0001  # Assumes sediment flux less than this threshold can be ignored
                                        if inundation == 1:  # Inundation regime

                                            SedFluxOut[TS, d, i] = (
                                                SedFluxIn[TS, d, i] * self._Cbb_i
                                            )
                                            Qs0 = SedFluxOut[TS, d, i]

                                            Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                            Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                            Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                            if i > 0 and Qs1 > Qs_bb_min:
                                                SedFluxIn[TS, d + 1, i - 1] += Qs1

                                            if Qs1 > Qs_bb_min:
                                                SedFluxIn[TS, d + 1, i] += Qs2

                                            if (
                                                i < (self._BarrierLength - 1)
                                                and Qs1 > Qs_bb_min
                                            ):
                                                SedFluxIn[TS, d + 1, i + 1] += Qs3

                                        else:  # Run-up regime

                                            SedFluxOut[TS, d, i] = (
                                                SedFluxIn[TS, d, i] * self._Cbb_r
                                            )
                                            Qs0 = SedFluxOut[TS, d, i]

                                            Qs1 = Qs0 * Q1 / (Q1 + Q2 + Q3)
                                            Qs2 = Qs0 * Q2 / (Q1 + Q2 + Q3)
                                            Qs3 = Qs0 * Q3 / (Q1 + Q2 + Q3)

                                            if i > 0 and Qs1 > Qs_bb_min:
                                                SedFluxIn[TS, d + 1, i - 1] += Qs1

                                            if Qs1 > Qs_bb_min:
                                                SedFluxIn[TS, d + 1, i] += Qs2

                                            if (
                                                i < (self._BarrierLength - 1)
                                                and Qs1 > Qs_bb_min
                                            ):
                                                SedFluxIn[TS, d + 1, i + 1] += Qs3

                                    ### Saline Flooding
                                    (
                                        self._ShrubDomainFemale,
                                        self._ShrubDomainMale,
                                    ) = self.SalineFlooding(
                                        ShrubDomainWidth,
                                        ShrubDomainAll,
                                        self._ShrubDomainFemale,
                                        self._ShrubDomainMale,
                                        d,
                                        i,
                                        Q0,
                                    )

                        ### Update Elevation After Every Storm Hour
                        if inundation == 1:
                            ElevationChange = (
                                SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                            ) / substep
                        else:
                            ElevationChange = (
                                SedFluxIn[TS, :, :] - SedFluxOut[TS, :, :]
                            ) / substep  # / surges #?
                        Elevation[TS, :, :] = Elevation[TS, :, :] + ElevationChange

                        # Calculate and save volume of sedimen deposited on/behind the island for every hour
                        OWloss = OWloss + np.sum(ElevationChange[1:, :])

                        # Update amount of burial/erosion for each shrub
                        self._BurialDomain = self.UpdateBurial(
                            self._BurialDomain,
                            ElevationChange,
                            ShrubDomainWidth,
                            ShrubDomainAll,
                        )

                    ### Update Interior Domain After Every Storm
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
                    self._DomainWidth = np.shape(self._InteriorDomain)[0]
                    (
                        self._ShrubDomainFemale,
                        self._ShrubDomainMale,
                        ShrubDomainAll,
                        self._ShrubPercentCover,
                        self._BurialDomain,
                    ) = self.UpdateShrubDomains(
                        self._DomainWidth,
                        ShrubDomainWidth,
                        self._ShrubDomainFemale,
                        self._ShrubDomainMale,
                        ShrubDomainAll,
                        self._ShrubPercentCover,
                        self._BurialDomain,
                    )

        # Record storm data
        self._StormCount.append(self._numstorm)

        ###########################################
        ### Ocean Shoreline Change

        ### Calculate shoreline/shoreface position change following Lorenzo-Trueba and Ashton (2014)
        (
            SCR,
            self._x_s,
            self._x_t,
            self._x_s_TS,
            self._x_t_TS,
            self._x_b_TS,
            self._s_sf_TS,
            self._QowTS,
            self._QsfTS,
        ) = self.LTA_SC(
            self._InteriorDomain,
            OWloss,
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
        )

        ### Erode/Prograde Ocean Shoreline
        drown_break = 0

        self._SCRagg = (
            self._SCRagg + SCR
        )  # Account for any residual shoreline change (less than cell size) from previous time step

        if abs(self._SCRagg) >= 1:
            sc = math.floor(abs(self._SCRagg))
            if (
                self._SCRagg > 0
            ):  # Positive = prograde, add row(s) to front of interior domain
                for d in range(sc):
                    self._InteriorDomain = np.vstack(
                        [
                            self._DuneDomain[self._time_index, :, -1] + self._BermEl,
                            self._InteriorDomain,
                        ]
                    )  # New interior row added with elevation of previous dune field (i.e. previous dune row now part of interior)

                    newDuneHeight = np.ones([self._BarrierLength]) * (
                        0.01
                        + (
                            -0.005
                            + (0.005 - (-0.005)) * np.random.rand(self._BarrierLength)
                        )
                    )
                    self._DuneDomain[self._time_index, :, :] = np.roll(
                        self._DuneDomain[self._time_index, :, :], 1, axis=1
                    )
                    self._DuneDomain[self._time_index, :, 0] = newDuneHeight

                    # Update shrub domains too
                    self._ShrubDomainFemale = np.concatenate(
                        (np.zeros([1, self._BarrierLength]), self._ShrubDomainFemale)
                    )
                    self._ShrubDomainMale = np.concatenate(
                        (np.zeros([1, self._BarrierLength]), self._ShrubDomainMale)
                    )
                    self._ShrubPercentCover = np.concatenate(
                        (np.zeros([1, self._BarrierLength]), self._ShrubPercentCover)
                    )
                    self._BurialDomain = np.concatenate(
                        (np.zeros([1, self._BarrierLength]), self._BurialDomain)
                    )
                self._DomainWidth = (
                    self._DomainWidth + sc
                )  # Update width of interior domain
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
                    self._InteriorDomain = np.delete(self._InteriorDomain, (0), axis=0)
                    if np.shape(self._InteriorDomain)[0] <= 0:
                        drown_break = 1
                        return
                    self._DuneDomain[self._time_index, :, :] = np.roll(
                        self._DuneDomain[self._time_index, :, :], -1, axis=1
                    )
                    self._DuneDomain[self._time_index, :, -1] = newDuneHeight

                    # Update shrub domains too
                    self._ShrubDomainFemale = np.delete(
                        self._ShrubDomainFemale, (0), axis=0
                    )
                    self._ShrubDomainMale = np.delete(
                        self._ShrubDomainMale, (0), axis=0
                    )
                    self._ShrubPercentCover = np.delete(
                        self._ShrubPercentCover, (0), axis=0
                    )
                    self._BurialDomain = np.delete(self._BurialDomain, (0), axis=0)
                self._DomainWidth = (
                    self._DomainWidth - sc
                )  # Update width of interior domain
                self._ShorelineChange = self._ShorelineChange - sc
                self._ShorelineChangeTS.append(-sc)
                self._SCRagg = self._SCRagg + sc  # Reset, leaving residual
        else:
            self._ShorelineChangeTS.append(0)

        ### Check for drowning
        if drown_break == 1:
            print(
                "Barrier has WIDTH DROWNED at t = " + str(self._time_index) + " years"
            )
            self._TMAX = self._time_index - 1
            return
        elif all(j <= self._SL for j in self._InteriorDomain[0, :]):
            print(
                "Barrier has HEIGHT DROWNED at t = " + str(self._time_index) + " years"
            )
            self._TMAX = self._time_index - 1
            return

        ### Recalculate and save DomainWidth and InteriorWidth
        self._DomainWidth, InteriorWidth, InteriorWidth_Avg = self.FindWidths(
            self._InteriorDomain, self._SL
        )
        self._InteriorWidth_AvgTS.append(InteriorWidth_Avg)

        ###########################################
        ### Save domains of this timestep
        self._DomainTS[self._time_index] = self._InteriorDomain

        self._ShrubFemaleTS[self._time_index] = self._ShrubDomainFemale
        self._ShrubMaleTS[self._time_index] = self._ShrubDomainMale
        ShrubDomainAll = self._ShrubDomainMale + self._ShrubDomainFemale

        # Calculate Percent Cover and Area
        self._ShrubPercentCover, self._PercentCoverTS, self._ShrubArea = self.CalcPC(
            ShrubDomainAll, self._PercentCoverTS, self._ShrubArea, self._time_index
        )
