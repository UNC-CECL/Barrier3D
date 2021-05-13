import matplotlib.pyplot as plt


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
