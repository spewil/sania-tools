import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from utils import plot
import numpy as np
matplotlib.use("MacOSX")

from sequence import sc
from pca import bin_projections, bin_components, num_projections, num_components, mut_projections
from frequency import unique_mutations, counts

base_pairs = ["-", "N", "W", "A", "T", "C", "G"]
colors = [
    "tab:green", "tab:pink", "tab:red", "tab:blue", "tab:purple", "tab:orange",
    "tab:cyan", "tab:gray", "tab:gray", "white"
]
listed_cmap = ListedColormap(colors)


def make_num_colorbar():
    pass


def plot_examples(colormaps):
    """
    Helper function to plot data with associated colormap.
    """
    np.random.seed(19680801)
    data = np.random.uniform(low=0, high=9, size=(30, 30))
    n = len(colormaps)
    fig, axs = plt.subplots(1,
                            n,
                            figsize=(n * 2 + 2, 3),
                            constrained_layout=True,
                            squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=0, vmax=10)
        cbar = fig.colorbar(psm, ax=ax, ticks=range(10))
        cbar.ax.set_yticklabels(base_pairs + [
            " ",
            " ",
            "Match",
        ])
    plt.show()


def plot_sequence_stack(stack):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    img = ax.imshow(stack,
                    aspect='auto',
                    cmap=listed_cmap,
                    vmin=0,
                    vmax=len(colors))
    cbar = fig.colorbar(img, ax=ax)
    labels = base_pairs + [" ", " ", "Match"]
    cbar.set_ticks(np.linspace(0.5, len(colors) - 0.5, len(colors)))
    cbar.set_ticklabels(labels)
    ax.set_yticks([i - 0.5 for i in list(range(stack.shape[0]))])
    ax.set_yticklabels(list(range(1, stack.shape[0] + 1)))
    ax.yaxis.grid(True, which='major', color="r")
    ax.set_ylabel("sequence")
    ax.set_xlabel("nucleotide")
    return fig, ax


def plot_PC_projections(projections, title=None):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.plot(projections[0], projections[1], 'ok', alpha=0.2)
    ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")
    if title is not None:
        ax.set_title(title)
    return fig, ax


fig, ax = plot_sequence_stack(sc.mutation_stack)
fig, ax = plot_sequence_stack(unique_mutations)

plot_PC_projections(bin_projections, title="binary")
plot_PC_projections(num_projections, title="numerical")
plot_PC_projections(mut_projections, title="mutations")

# fig, axes = plt.subplots(2, 1, figsize=(8, 8))
# axes[0].imshow(unique_sequences, aspect="auto")
# axes[1].bar(range(len(counts)), counts)
plt.show()
