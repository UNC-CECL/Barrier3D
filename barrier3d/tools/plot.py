import matplotlib.pyplot as plt
import numpy as np


def plot_dune_height(dune_height, max_dune_height):
    dune_crest = dune_height.max(axis=2)

    fig = plt.figure(figsize=(14, 8))
    plt.rcParams.update({"font.size": 13})
    ax = fig.add_subplot(111)
    ax.matshow(
        (dune_crest) * 10,
        origin="lower",
        cmap="bwr",
        aspect="auto",
        vmin=0,
        vmax=max_dune_height * 10,
    )
    ax.xaxis.set_ticks_position("bottom")  # analysis:ignore

    plt.xlabel("Alongshore Distance (dam)")
    plt.ylabel("Year")
    plt.title("Dune Height (m)")

    fig.show()


def plot_shrub_percent_cover_tmax(PercentCoverTS, TMAX, DeadPercentCoverTS):
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

    shrubFig2.show()


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
    ax.matshow(
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
    # name = "Output/FinalElevation"
    # elevFig1.savefig(name)
    plt.show()
