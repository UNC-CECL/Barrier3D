# Plotting functions for

# ~ Barrier3D ~
# A spatially explicit exploratory model of barrier island evolution in three dimensions


"""----------------------------------------------------
Copyright (C) 2020 Ian R.B. Reeves
Full copyright notice located in main Barrier3D.py file
----------------------------------------------------"""

# Version Number: 5
# Updated: 12 January 2021 by K. Anarde

import math
import os

import imageio  # analysis:ignore
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # analysis:ignore
from scipy import signal

# ==============================================================================================================================================
# PLOTTING & CALCULATION FUNCTIONS
# ==============================================================================================================================================


# ===================================================
# 2: Elevation Domain For Last Time Step


def plot_ElevTMAX(
    TMAX,
    DuneDomain,
    DomainTS,
    BermEl,
    Shrub_ON,
    PercentCoverTS,
    DeadPercentCoverTS,
    DuneWidth,
):

    TMAX = TMAX - 1
    Dunes = (DuneDomain[TMAX, :, :] + BermEl) * 10
    Dunes = np.rot90(Dunes)
    Dunes = np.flipud(Dunes)
    Domain = DomainTS[TMAX] * 10
    Domain = np.vstack([Dunes, Domain])
    if Shrub_ON == 1:
        Shrubs = PercentCoverTS[TMAX]
        Dead = DeadPercentCoverTS[TMAX]
        Sy, Sx = np.argwhere(Shrubs > 0).T
        Sz = Shrubs[Sy, Sx] * 80
        Dy, Dx = np.argwhere(Dead > 0).T
        Dz = Dead[Dy, Dx] * 80
    elevFig1 = plt.figure(figsize=(14, 5))
    ax = elevFig1.add_subplot(111)
    cax = ax.matshow(
        Domain, origin="lower", cmap="terrain", vmin=-1.1, vmax=4.0
    )  # , interpolation='gaussian') # analysis:ignore
    if Shrub_ON == 1:
        ax.scatter(
            Sx,
            Sy + DuneWidth,
            marker="$*$",
            s=Sz,
            c="black",
            alpha=0.7,
            edgecolors="none",
        )
        ax.scatter(
            Dx,
            Dy + DuneWidth,
            marker="$*$",
            s=Dz,
            c="red",
            alpha=0.7,
            edgecolors="none",
        )
    ax.xaxis.set_ticks_position("bottom")
    # cbar = elevFig1.colorbar(cax)
    # cbar.set_label('Elevation (m)', rotation=270)
    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Cross-Shore Distance (dam)")
    plt.title("Interior Elevation (m)")
    timestr = "Time = " + str(TMAX) + " yrs"
    plt.text(1, 1, timestr)
    plt.tight_layout()
    plt.show()
    name = "Output/FinalElevation"
    # elevFig1.savefig(name)
    plt.show()


# ===================================================
# 3: Elevation Domain Frames


def plot_ElevFrames(TMAX, DomainTS):

    for t in range(TMAX):
        elevFig1 = plt.figure(figsize=(14, 5))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            DomainTS[t], origin="lower", cmap="terrain", vmin=-1.1, vmax=4.0
        )
        ax.xaxis.set_ticks_position("bottom")
        cbar = elevFig1.colorbar(cax)  # analysis:ignore
        # cbar.set_label('Elevation (m)', rotation=270)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Diatance (dam)")
        plt.title("Interior Elevation")
        timestr = "Time = " + str(t) + " yrs"
        plt.text(1, 1, timestr)
        # plt.show()
        name = "Output/SimFrames/elev_" + str(t)
        # elevFig1.savefig(name)


# ===================================================
# 4: Animation Frames of Barrier and Dune Elevation


def plot_ElevAnimation(
    InteriorWidth_AvgTS,
    ShorelineChange,
    DomainTS,
    DuneDomain,
    SL,
    x_s_TS,
    Shrub_ON,
    PercentCoverTS,
    TMAX,
    DeadPercentCoverTS,
):

    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(
        max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange) + OriginY + 35
    )  # was +15

    for t in range(TMAX):
        # Build beach elevation domain
        BeachDomain = np.zeros([BeachWidth, BarrierLength])
        berm = math.ceil(BeachWidth * 0.65)
        BeachDomain[berm : BeachWidth + 1, :] = BermEl
        add = (BermEl - SL) / berm
        for i in range(berm):
            BeachDomain[i, :] = SL + add * i

        # Make animation frame domain
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t, :, :] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength]) * -1
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
            Sz = Shrubs[Sy, Sx] * 80  # 30
            Dy, Dx = np.argwhere(Dead > 0).T
            Dz = Dead[Dy, Dx] * 80  # 22

        # Plot and save
        elevFig1 = plt.figure(figsize=(7, 12))
        ax = elevFig1.add_subplot(111)
        cax = ax.matshow(
            AnimateDomain, origin="lower", cmap="terrain", vmin=-1.1, vmax=4.0
        )  # , interpolation='gaussian') # analysis:ignore
        if Shrub_ON == 1:
            # ax.scatter(Sx, Sy, marker='$*$', s=Sz, c='black', alpha=0.35, edgecolors='none')
            # ax.scatter(Dx, Dy, marker='$*$', s=Dz, facecolors='none', edgecolors='maroon', alpha=0.35)
            ax.scatter(
                Sx, Sy, marker="$*$", s=Sz, c="black", alpha=0.6, edgecolors="none"
            )
            ax.scatter(
                Dx, Dy, marker="$*$", s=Dz, c="maroon", alpha=0.5, edgecolors="none"
            )
        ax.xaxis.set_ticks_position("bottom")
        # cbar = elevFig1.colorbar(cax)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Diatance (dam)")
        plt.title("Interior Elevation")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        newpath = "Output/SimFrames/"
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        plt.text(1, 1, timestr)
        name = "Output/SimFrames/elev_" + str(t)
        elevFig1.savefig(name)  # dpi=200
        plt.close(elevFig1)

    frames = []
    for filenum in range(TMAX):
        filename = "Output/SimFrames/elev_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("Output/SimFrames/elev.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# ===================================================
# 5: Cross-Shore Transect Every 100 m Alongshore For Last Time Step


def plot_XShoreTransects(InteriorDomain, DuneDomain, SL, TMAX):

    # Build beach elevation
    BW = 6  # beach width (dam) for illustration purposes
    BeachX = np.zeros(BW)
    berm = math.ceil(BW * 0.5)
    BeachX[berm : BW + 1] = BermEl
    add = (BermEl - SL) / berm
    for i in range(berm):
        BeachX[i] = SL + add * i
    # Plot full tranects
    plt.figure()
    for v in range(0, BarrierLength, 20):
        CrossElev = InteriorDomain[:, v]
        Dunes = DuneDomain[TMAX - 1, v, :] + BermEl
        CrossElev1 = np.insert(CrossElev, 0, Dunes)
        CrossElev2 = np.insert(CrossElev1, 0, BeachX)
        CrossElev = CrossElev2 * 10  # Convert to meters
        plt.plot(CrossElev)
    fig = plt.gcf()
    fig.set_size_inches(14, 6)
    plt.hlines(SL, -1, len(CrossElev + 1), colors="dodgerblue")
    plt.xlabel("Cross-Shore Distance (dam)")
    plt.ylabel("Elevation (m)")
    plt.title("Cross-shore Topo Transects")
    plt.show()
    name = "Output/Profiles"
    fig.savefig(name)


# ===================================================
# 6: Shoreline Positions Over Time


def plot_ShorelinePositions(x_s_TS, x_b_TS):

    scts = [(x - x_s_TS[0]) * -10 for x in x_s_TS]
    bscts = [(x - x_s_TS[0]) * -10 for x in x_b_TS]
    shorelinefig = plt.figure()
    plt.plot(scts, "b")
    plt.plot(bscts, "g")
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    plt.ylabel("Shoreline Position (m)")
    plt.xlabel("Year")
    plt.show()
    name = "Output/Shorelines"
    shorelinefig.savefig(name)


# ===================================================
# 7: Shoreline Change Rate Over Time


def plot_ShorelineChangeRate(x_s_TS):

    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    ratefig = plt.figure()
    plt.plot(rate)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Year")
    plt.ylabel("Shoreline Erosion Rate(m/yr)")
    plt.show()
    name = "Output/ShorelineRate"
    ratefig.savefig(name)


# ===================================================
# 8: Run-up vs Inundation count


def plot_RuInCount(RunUpCount, InundationCount):

    objects = ("Run Up", "Inundation")
    y_pos = np.arange(len(objects))
    plt.figure()
    plt.bar(y_pos, [RunUpCount, InundationCount])
    plt.xticks(y_pos, objects)
    plt.ylabel("Count")


# ===================================================
# 9: Shoreface LTA14 transects over time


def plot_LTATransects(
    SL, TMAX, x_b_TS, x_t_TS, x_s_TS, RSLR, DShoreface, BayDepth, BermEl
):

    xmax = x_b_TS[TMAX - 1] + 20

    SFfig = plt.figure(figsize=(20, 5))
    colors = plt.cm.jet(np.linspace(0, 1, TMAX))

    for t in range(0, TMAX, 5):  # Plots one transect every 25 years
        # Create data points
        Tx = x_t_TS[t]
        Ty = ((SL + (t * RSLR[t])) - DShoreface) * 10
        Sx = x_s_TS[t]
        Sy = (SL + (t * RSLR[t])) * 10
        Bx = x_b_TS[t]
        By = ((SL + (t * RSLR[t])) - BayDepth) * 10
        Hx1 = Sx
        Hy1 = ((t * RSLR[t]) + BermEl) * 10
        Hx2 = Bx
        Hy2 = Hy1
        Mx = xmax
        My = By

        x = [Tx, Sx, Hx1, Hx2, Bx, Mx]
        y = [Ty, Sy, Hy1, Hy2, By, My]

        # Plot
        plt.plot(x, y, color=colors[t])

    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Elevation (m)")
    plt.title("Shoreface Evolution")
    plt.show()

    # Save
    name = "Output/Shoreface"
    SFfig.savefig(name)


# ===================================================
# 10: Average Island Elevation Over Time


def plot_AvgIslandElev(h_b_TS):

    beTS = h_b_TS
    plt.figure()
    plt.plot(beTS)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Time (yrs)")
    plt.ylabel("Average Island Elevation (m)")
    plt.show()


# ===================================================
# 11: Shoreface Slope Over Time


def plot_ShorefaceSlope(s_sf_TS):

    ssfTS = s_sf_TS
    plt.figure()
    plt.plot(ssfTS)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Time (yrs)")
    plt.ylabel("Shoreface Slope")
    plt.title("Shoreface Slope Over Time")
    plt.show()


# ===================================================
# 12: Average Interior Width Over Time


def plot_AvgInteriorWidth(InteriorWidth_AvgTS):

    aiw = InteriorWidth_AvgTS
    plt.figure()
    plt.plot(aiw)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Time (yrs)")
    plt.ylabel("Average Interior Width (dam)")
    plt.title("Average Interior Width Over Time")
    plt.show()


# ===================================================
# 13: Shoreface Overwash Flux Over Time


def plot_OverwashFlux(QowTS):

    plt.figure()
    plt.plot(QowTS)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Time (yrs)")
    plt.ylabel("Qow (m^3)")
    plt.title("Overwash Flux")
    # plt.ylim(0,210)
    plt.show()


# ===================================================
# 14: Width, Berm Elevation, SF Slope, Shoreline Change, and Overwash Flux Over Time (all in one)


def plot_StatsSummary(
    s_sf_TS, x_s_TS, TMAX, InteriorWidth_AvgTS, QowTS, QsfTS, Hd_AverageTS
):

    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(14, 18)
    plt.rcParams.update({"font.size": 17})

    # Shoreface Slope
    plt.subplot(6, 1, 1)
    ssfTS = s_sf_TS
    plt.plot(ssfTS)
    plt.hlines(s_sf_eq, 0, TMAX - 1, colors="black", linestyles="dashed")

    plt.ylabel("Shoreface Slope")

    # Interior Width
    plt.subplot(6, 1, 2)
    aiw = [a * 10 for a in InteriorWidth_AvgTS]
    plt.plot(aiw)
    plt.ylabel("Avg. Width (m)")  # Avergae Interior Width

    # Shoreline Change
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]
    plt.subplot(6, 1, 3)
    plt.plot(scts)
    plt.ylabel("Shoreline Position (m)")

    # Overwash Flux
    plt.subplot(6, 1, 4)
    plt.plot(QowTS)
    #    movingavg = np.convolve(QowTS, np.ones((50,))/50, mode='valid')
    #    movingavg = [i * 10 for i in movingavg]
    #    plt.plot(movingavg, 'r--')
    plt.ylabel("Qow (m^3/m)")

    # Shoreface Flux
    plt.subplot(6, 1, 5)
    plt.plot(QsfTS)
    plt.ylabel("Qsf (m^3/m)")

    # Dune Height
    aHd = [a * 10 for a in Hd_AverageTS]
    plt.subplot(6, 1, 6)
    plt.plot(aHd)
    plt.xlabel("Year")
    plt.ylabel("Avg. Dune Height (m)")  # Average Dune Height

    plt.show()
    name = "Output/Stats"
    fig.savefig(name)


# ===================================================
# 15: 3D Plot of Island Domain For Last Time Step


def plot_3DElevTMAX(TMAX, t, SL, DuneDomain, DomainTS):

    if TMAX > t:
        TMAX = t
    # Build beach elevation domain
    BW = 6
    Beach = np.zeros([BW, BarrierLength])
    berm = math.ceil(BW * 0.65)
    Beach[berm : BW + 1, :] = BermEl
    add = (BermEl - SL) / berm
    for i in range(berm):
        Beach[i, :] = SL + add * i

    # Construct frame
    Dunes = [(DuneDomain[TMAX] + BermEl) * 10] * DuneWidth
    Water = np.zeros([3, BarrierLength])
    Domain = DomainTS[TMAX] * 10
    Domain[Domain < 0] = 0
    Domain = np.vstack([Water, Beach, Dunes, Domain, Water])
    Dlen = np.shape(Domain)[1]
    Dwid = np.shape(Domain)[0]
    fig = plt.figure(figsize=(12, 9))
    # fig.set_size_inches(12,7)
    ax = fig.add_subplot(111, projection="3d")
    scale_x = 1
    scale_y = Dwid / Dlen
    scale_z = 4 / Dlen * 3
    ax.get_proj = lambda: np.dot(
        Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1])
    )
    X = np.arange(Dlen)
    Y = np.arange(Dwid)
    X, Y = np.meshgrid(X, Y)
    Z = Domain
    ax.plot_surface(
        X, Y, Z, cmap="terrain", alpha=1, vmin=-1.1, vmax=4.0, linewidth=0, shade=True
    )
    ax.set_zlim(0, 4)

    # Plot shrubs
    # Shrubs = PercentCoverTS[TMAX]
    # Shrubs[Shrubs>0] = 1
    # Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
    # Shrubs = Shrubs * Domain
    # Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1
    # Shrubs[Shrubs<1] = None
    # ax.scatter(X, Y+1, Shrubs, s=30, c='black')

    ax.view_init(10, 155)
    plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3)  # mostly centered
    plt.show()
    name = "Output/Domain3D"
    fig.savefig(name, dpi=200)


# ===================================================
# 16: 3D Animation Frames of Island Elevation and Shrubs (no translation)


def plot_3DElevFrames(DomainTS, SL, TMAX, DuneDomain):

    for t in range(0, len(DomainTS)):
        # Build beach elevation domain
        BW = 6
        Beach = np.zeros([BW, BarrierLength])
        berm = math.ceil(BW * 0.65)
        Beach[berm : BW + 1, :] = BermEl
        add = (BermEl - SL) / berm
        for i in range(berm):
            Beach[i, :] = SL + add * i

        # Construct frame
        Dunes = [(DuneDomain[TMAX] + BermEl) * 10] * DuneWidth
        Water = np.zeros([3, BarrierLength])
        Domain = DomainTS[TMAX] * 10
        Domain = np.vstack([Water, Beach, Dunes, Domain, Water])
        Dlen = np.shape(Domain)[1]
        Dwid = np.shape(Domain)[0]
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection="3d")
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 4
        ax.get_proj = lambda: np.dot(
            Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1])
        )
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = Domain
        ax.plot_surface(
            X,
            Y,
            Z,
            cmap="terrain",
            alpha=1,
            vmin=-1.1,
            vmax=4.0,
            linewidth=0,
            shade=True,
        )
        ax.set_zlim(0, 4)

        # Plot shrubs
        #    Shrubs = PercentCoverTS[t]
        #    Shrubs[Shrubs>0] = 1
        #    Shrubs = np.vstack([np.zeros([DuneWidth,BarrierLength]), Shrubs])
        #    Shrubs = Shrubs * Domain
        #    Shrubs[Shrubs>0] = Shrubs[Shrubs>0] + 0.1
        #    Shrubs[Shrubs<1] = None
        #    ax.scatter(X, Y+1, Shrubs, s=30, c='black')

        timestr = "Time = " + str(t) + " yrs"
        ax.set_ylabel(timestr)
        ax.view_init(20, 155)
        plt.subplots_adjust(
            left=-1.2, right=1.3, top=2.2, bottom=-0.3
        )  # mostly centered
        plt.show()
        name = "Output/SimFrames/3D_" + str(t)
        fig.savefig(name, dpi=150)


# ===================================================
# 17: 3D Animation Frames of Island Evolution (with barrier translation)


def plot_3DElevAnimation(
    DomainTS, SL, TMAX, DuneDomain, DomainWidth, x_s_TS, ShorelineChange
):

    BW = 6
    AniDomainWidth = DomainWidth + round(BW) + 12 + abs(ShorelineChange)
    OriginY = 5
    for t in range(0, TMAX):

        # Build beach elevation domain
        BeachDomain = np.zeros([BW, BarrierLength])
        berm = math.ceil(BW * 0.65)
        BeachDomain[berm : BW + 1, :] = BermEl
        add = (BermEl - SL) / berm
        for i in range(berm):
            BeachDomain[i, :] = SL + add * i

        # Make animation frame
        Domain = DomainTS[t] * 10
        Domain[Domain < 0] = 0
        Dunes = (DuneDomain[t, :, :] + BermEl) * 10
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
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection="3d")
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 1  # 4 / Dlen * 4
        ax.get_proj = lambda: np.dot(
            Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1])
        )
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = AnimateDomain
        ax.plot_surface(
            X,
            Y,
            Z,
            cmap="terrain",
            alpha=1,
            vmin=-1.1,
            vmax=4.0,
            linewidth=0,
            shade=True,
        )
        ax.set_zlim(0, 4)

        timestr = "Time = " + str(t) + " yrs"
        ax.set_ylabel(timestr, labelpad=50)
        ax.view_init(20, 150)  # ax.view_init(20,155)
        #    plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        # plt.subplots_adjust(left=-0.7, right=1.3, top=2.2, bottom=-0.3) # mostly centered previous
        plt.subplots_adjust(left=-0.25, right=1.05, top=2.3, bottom=-0.2)
        # plt.show()
        name = "Output/SimFrames/3DAni_" + str(t)
        fig.savefig(name, dpi=150)
        plt.close(fig)

    frames = []
    for filenum in range(TMAX):
        filename = "Output/SimFrames/3DAni_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("Output/SimFrames/3DAni.gif", frames, "GIF-FI")


# ===================================================
# 18: Shrub Age Domain at Simulation End


def plot_ShrubAgeTMAX(ShrubDomainAll, ShrubDeadTS):

    if Shrub_ON == 1:
        ShrubAll = ShrubDomainAll
        shrubFig1 = plt.figure(figsize=(14, 5))
        ax = shrubFig1.add_subplot(111)
        cax = ax.matshow(
            ShrubAll, origin="lower", cmap="afmhot_r", vmin=0, vmax=10
        )  # analysis:ignore
        ax.xaxis.set_ticks_position("bottom")
        cbar = shrubFig1.colorbar(cax)
        cbar.set_label("Shrub Age", rotation=270, labelpad=20)

        Dead = ShrubDeadTS[-1]
        Dy, Dx = np.where(Dead > 0)
        ax.scatter(Dx, Dy, marker="x", s=20, color="black", alpha=0.25)

        plt.xlabel("Alongshore Distance (dm)")
        plt.ylabel("Cross-Shore Diatance (dm)")
        plt.title("Final Shrub Age")
        plt.show()


# ===================================================
# 19: Percent Cover Domain at Simulation End


def plot_ShrubPercentCoverTMAX(PercentCoverTS, TMAX, DeadPercentCoverTS):

    ShrubPC = PercentCoverTS[TMAX - 1]
    shrubFig2 = plt.figure(figsize=(10, 5))
    ax = shrubFig2.add_subplot(111)
    cax = ax.matshow(
        ShrubPC, origin="lower", cmap="YlGn", vmin=0, vmax=1
    )  # analysis:ignore
    ax.xaxis.set_ticks_position("bottom")
    cbar = shrubFig2.colorbar(cax)
    cbar.set_label("Shrub Percent Cover", rotation=270, labelpad=20)

    Dead = DeadPercentCoverTS[TMAX - 1]
    Dy, Dx = np.argwhere(Dead > 0).T
    Dz = Dead[Dy, Dx] * 20
    ax.scatter(Dx, Dy, marker="x", s=Dz, color="k", alpha=0.25)

    plt.xlabel("Alongshore Distance (dm)")
    plt.ylabel("Cross-Shore Diatance (dm)")
    plt.title("Final Shrub Percent Cover")
    plt.show()
    name = "Output/PercentCover"
    # shrubFig2.savefig(name)
    shrubFig2.show()


# ===================================================
# 20: Shrub Area Over Time


def plot_ShrubArea(ShrubArea):

    if Shrub_ON == 1:
        area = ShrubArea
        plt.figure()
        plt.plot(area)
        fig = plt.gcf()
        fig.set_size_inches(14, 3)
        plt.xlabel("Year")
        plt.ylabel("Shrub Area (dam^2)")
        plt.show()
        name = "Output/ShrubArea"
        fig.savefig(name)


# ===================================================
# 21: Storm count over time


def plot_StormCount(StormCount):

    plt.figure()
    plt.plot(StormCount)
    fig = plt.gcf()
    fig.set_size_inches(14, 5)
    plt.xlabel("Year")
    plt.ylabel("Number of Storms")
    plt.title("Storm Count")
    plt.show()


# ===================================================
# 22: Alongshore Dune Height Over Time


def plot_AlongshoreDuneHeight(DuneDomain):

    Dunes = DuneDomain.max(axis=2)

    # Plot full tranects
    plt.figure()
    for x in range(0, BarrierLength, 10):
        Hd_TS = Dunes[:, x]
        Hd_TS = Hd_TS * 10  # Convert to meters
        plt.plot(Hd_TS)
    fig = plt.gcf()
    fig.set_size_inches(14, 6)
    plt.xlabel("Year")
    plt.ylabel("Dune Height (m)")
    plt.title("Dune Height Alongshore")
    plt.show()
    name = "Output/Dunes_Alongshore"
    fig.savefig(name)


# ===================================================
# 23: Calculate shoreline change periodicity


def calc_ShorelinePeriodicity(TMAX, x_s_TS):

    # Shoreline Change & Change Rate Over Time
    scts = [(x - x_s_TS[0]) * 10 for x in x_s_TS]

    # Filter
    win = 15
    poly = 3
    der1 = signal.savgol_filter(scts, win, poly, deriv=1)

    HitDown = []  # Slow-downs
    HitUp = []  # Speed-ups

    window1 = 3  # Maximum allowed length for gaps in slow periods
    window2 = 15  # Minimum length required for slow periods, including gaps
    buffer = 3
    thresh1 = 0.5  # Max slope for slow periods
    thresh2 = 1.5

    # Find slow periods
    der_under = np.where(der1 < thresh1)[0]

    if len(der_under) > 0:

        gaps = np.diff(der_under) > window1
        peak_start = np.insert(der_under[1:][gaps], 0, der_under[0])
        peak_stop = np.append(der_under[:-1][gaps], der_under[-1])

        for n in range(len(peak_stop)):
            if peak_stop[n] - peak_start[n] > window2:
                if len(HitDown) == 0:
                    if peak_start[n] > buffer:
                        HitDown.append(peak_start[n])
                    if peak_stop[n] < len(scts) - buffer:
                        HitUp.append(peak_stop[n])
                else:
                    gap_length = peak_start[n] - HitUp[-1]
                    gap_slope = (scts[peak_start[n]] - scts[HitUp[-1]]) / gap_length
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
                FastDur.append(HitDown[n + 1] - HitUp[n])
        else:
            for n in range(Slows - 1):
                SlowDur.append(HitUp[n + 1] - HitDown[n])
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


# ===================================================
# 24: Average Dune Height and Storm TWL


def plot_DuneStorm(Hd_AverageTS, StormSeries, TMAX):

    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(14, 6)

    # Dune Height
    aHd = [a * 10 for a in Hd_AverageTS]  # Use this for dune height
    # aHd = [(a + BermEl) * 10 for a in Hd_AverageTS] # Use this for dune elevation
    plt.plot(aHd, color="teal")
    plt.xlabel("Year")
    plt.ylabel("Avg. Dune Height (m)")  # Use this for dune height
    # plt.ylabel('Avg. Dune Elevation (m)') # Use this for dune elevation

    # Storms
    if TMAX >= 1000:
        stop = len(StormSeries)
    else:
        stop = np.where(StormSeries[:, 0] >= TMAX)[0][0]

    stormX = StormSeries[0:stop, 0]  # Yr
    stormY = StormSeries[0:stop, 1]  # Rhigh

    stormX = stormX[stormY > BermEl]  # Remove storms with TWL lower than berm
    stormY = stormY[stormY > BermEl]
    stormY = [
        (a - BermEl) * 10 for a in stormY
    ]  # Height relative to berm # Use this for dune height
    # stormY = [(a) * 10 for a in stormY] # Use this for dune elevation

    plt.scatter(stormX, stormY, c="r", marker="*")

    plt.show()
    name = "Output/DuneStorm"
    fig.savefig(name)


# ===================================================
# 25: Seabed Profile


def plot_SeabedProfile(SL, TMAX, x_t_TS):

    Tx = []
    Ty = []
    for t in range(TMAX):
        Tx.append(x_t_TS[t] * 10)
        Ty.append(((SL + np.sum(RSLR[0:t])) - DShoreface) * 10)

    Sby = np.linspace(Ty[0], Ty[-1], num=TMAX)
    Sbx = np.linspace(0, Tx[-1], num=TMAX)

    # Plot
    plt.figure(figsize=(14, 5))
    plt.plot(Sbx, Sby, "slategray", linestyle="--")
    plt.plot(Tx, Ty, "chocolate")
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    plt.title("Seabed Profile")
    plt.show()


# ===================================================
# 26: Animation Frames of Shrub Percent Cover and Barrier Migration


def plot_ShrubAnimation(
    InteriorWidth_AvgTS,
    ShorelineChange,
    DomainTS,
    DuneDomain,
    SL,
    x_s_TS,
    Shrub_ON,
    PercentCoverTS,
    TMAX,
    DeadPercentCoverTS,
):

    BeachWidth = 6
    OriginY = 10
    AniDomainWidth = int(
        max(InteriorWidth_AvgTS) + BeachWidth + abs(ShorelineChange) + OriginY + 35
    )  # was +15

    for t in range(TMAX):
        # Build beach elevation domain
        BeachDomain = np.zeros([BeachWidth, BarrierLength])
        berm = math.ceil(BeachWidth * 0.65)
        BeachDomain[berm : BeachWidth + 1, :] = BermEl
        add = (BermEl - SL) / berm
        for i in range(berm):
            BeachDomain[i, :] = SL + add * i

        # Make animation frame domain
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t, :, :] + BermEl) * 10
        Dunes = np.rot90(Dunes)
        Dunes = np.flipud(Dunes)
        Beach = BeachDomain * 10
        Domain = np.vstack([Beach, Dunes, Domain])
        Domain[Domain < 0] = -1
        Domain[Domain > 0] = 1
        AnimateDomain = np.ones([AniDomainWidth + 1, BarrierLength]) * -1
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
            Sz = Shrubs[Sy, Sx]
            Dy, Dx = np.argwhere(Dead > 0).T

        # Plot and save
        shrubfig = plt.figure(figsize=(10, 12))
        ax = shrubfig.add_subplot(111)
        cax = ax.matshow(AnimateDomain, origin="lower", cmap="Blues_r")
        im = ax.scatter(
            Sx,
            Sy,
            marker="o",
            s=32,
            c=Sz,
            cmap="YlGn",
            vmin=0,
            vmax=1,
            alpha=1,
            edgecolors="none",
        )
        ax.scatter(
            Dx, Dy, marker="o", s=22, facecolors="none", edgecolors="maroon", alpha=0.4
        )
        ax.xaxis.set_ticks_position("bottom")
        cbar = shrubfig.colorbar(im, ax=ax)
        cbar.set_label("Shrub Percent Cover", rotation=270, labelpad=20)
        plt.xlabel("Alongshore Distance (dam)")
        plt.ylabel("Cross-Shore Diatance (dam)")
        plt.title("Interior Elevation")
        plt.tight_layout()
        timestr = "Time = " + str(t) + " yrs"
        newpath = "Output/SimFrames/"
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        plt.text(1, 1, timestr)
        name = "Output/SimFrames/shrubani_" + str(t)
        shrubfig.savefig(name)  # dpi=200
        plt.close(shrubfig)

    frames = []
    for filenum in range(TMAX):
        filename = "Output/SimFrames/shrubani_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("Output/SimFrames/shrubani.gif", frames, "GIF-FI")
    print()
    print("[ * GIF successfully generated * ]")


# ===================================================
# 27: 3D Animation of Island Evolution With Shoreface and Bay


def plot_3DElevAnimation_Super(
    DomainTS, SL, TMAX, DuneDomain, DomainWidth, x_s_TS, x_t_TS, ShorelineChange
):  # <---------------------Does not work yet!

    BW = 6
    AniDomainWidth = int(
        DomainWidth + round(BW) + 12 + abs(ShorelineChange) + LShoreface
    )
    OriginY = 5
    for t in range(TMAX):

        # Build beach elevation domain
        BeachDomain = np.zeros([BW, BarrierLength])
        berm = math.ceil(BW * 0.65)
        BeachDomain[berm : BW + 1, :] = BermEl
        add = (BermEl - SL) / berm
        for i in range(berm):
            BeachDomain[i, :] = SL + add * i

        # Build shoreface elevation domain
        l_sf = int(x_s_TS[t] - x_t_TS[t])
        SFDomain = np.zeros([l_sf, BarrierLength])
        add = DShoreface / l_sf
        for i in range(l_sf):
            SFDomain[i, :] = -DShoreface + add * i

        # Make animation frame
        Domain = DomainTS[t] * 10
        Dunes = (DuneDomain[t, :, :] + BermEl) * 10
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
        AnimateDomain[0, :] = DShoreface * 10
        AnimateDomain[0:OriginTstart, :] = DShoreface * 10
        AnimateDomain[OriginTstop:-1, :] = BayDepth * 10
        AnimateDomain[OriginTstart:OriginTstop, :] = Domain

        Dlen = np.shape(AnimateDomain)[1]
        Dwid = np.shape(AnimateDomain)[0]
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection="3d")
        scale_x = 1
        scale_y = Dwid / Dlen
        scale_z = 4 / Dlen * 1  # 4 / Dlen * 4
        ax.get_proj = lambda: np.dot(
            Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1])
        )
        X = np.arange(Dlen)
        Y = np.arange(Dwid)
        X, Y = np.meshgrid(X, Y)
        Z = AnimateDomain
        ax.plot_surface(
            X,
            Y,
            Z,
            cmap="terrain",
            alpha=1,
            vmin=-1.1,
            vmax=4.0,
            linewidth=0,
            shade=True,
        )
        ax.set_zlim(0, 4)

        timestr = "Time = " + str(t) + " yrs"
        ax.set_ylabel(timestr, labelpad=50)
        ax.view_init(20, 150)  # ax.view_init(20,155)
        # plt.subplots_adjust(left=-1.2, right=1.3, top=2.2, bottom=-0.3) # mostly centered
        # plt.subplots_adjust(left=-0.7, right=1.3, top=2.2, bottom=-0.3) # mostly centered previous
        # plt.subplots_adjust(left=-0.25, right=1.05, top=2.3, bottom=-0.2)
        # plt.show()
        name = "Output/SimFrames/3DAni_" + str(t)
        fig.savefig(name, dpi=150)
        plt.close(fig)

    frames = []
    for filenum in range(TMAX):
        filename = "Output/SimFrames/3DAni_" + str(filenum) + ".png"
        frames.append(imageio.imread(filename))
    imageio.mimsave("Output/SimFrames/3DAni.gif", frames, "GIF-FI")
