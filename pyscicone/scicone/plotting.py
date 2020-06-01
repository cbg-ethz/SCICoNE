from scicone.constants import *
from scicone.utils import cluster_clones
import numpy as np
import string
from copy import deepcopy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
import seaborn as sns

sns.set_style("ticks", {"axes.grid": True})

datacmap = matplotlib.colors.LinearSegmentedColormap.from_list("cmap", BLUE_WHITE_RED)
cnvcmap = matplotlib.colors.LinearSegmentedColormap.from_list("cnvcmap", BLUE_WHITE_RED, 5)

def plot_matrix(data, cbar_title="", mode='data', chr_stops_dict=None,
                labels=None, textfontsize=24, tickfontsize=22, bps=None,
                figsize=(24,8), dpi=100, vmax=None, output_path=None):
    if mode == 'data':
        cmap = datacmap
    elif mode == 'cnv':
        cmap = cnvcmap
    else:
        raise AttributeError('mode argument must be one of \'data\' or \'cnv\'')

    data_ = np.array(data, copy=True)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    if labels is not None:
        labels = np.array(labels).ravel()
        labels_ = np.array(labels, copy=True)

        if mode == 'cnv':
            data_, labels_ = cluster_clones(data_, labels_, within_clone=False)
        else:
            data_, labels_ = cluster_clones(data_, labels_, within_clone=True)

        ticks = dict()
        unique_labels = np.unique(labels_)
        n_unique_labels = len(unique_labels)
        for label in unique_labels: # sorted
            # Get last pos
            t = np.where(labels_ == label)[0]
            if len(t) > 1:
                t = t[-1]
            ticks[label] = t
        gs = GridSpec(1, 2, wspace=0.05, width_ratios=[1, 40])
        ax = fig.add_subplot(gs[0])
        bounds = [0] + list(ticks.values())
        subnorm = matplotlib.colors.BoundaryNorm(bounds, n_unique_labels)
        clonecmap = matplotlib.colors.ListedColormap(list(LABEL_COLORS_DICT.values())[:n_unique_labels])
        cb = matplotlib.colorbar.ColorbarBase(
            ax,
            cmap=clonecmap,
            norm=subnorm,
            boundaries=bounds,
            spacing="proportional",
            orientation="vertical",
        )
        cb.outline.set_visible(False)
        cb.ax.set_ylabel("Clones", fontsize=textfontsize)
        ax.yaxis.set_label_position("left")
        for j, lab in enumerate(ticks.keys()):
            cb.ax.text(
                0.5,
                ((bounds[j + 1] + bounds[j]) / 2) / bounds[-1],
                lab,
                ha="center",
                va="center",
                rotation=90,
                color="w",
                fontsize=tickfontsize,
            )
        cb.set_ticks([])

        ax = fig.add_subplot(gs[1])
    else:
        ax = plt.gca()
    im = plt.pcolormesh(data_, cmap=cmap, rasterized=True)
    ax = plt.gca()
    plt.ylabel('Cells', fontsize=textfontsize)
    plt.xlabel('Bins', fontsize=textfontsize)
    if bps is not None:
        ax.vlines(bps, *ax.get_ylim(), colors="k", linestyles="dashed", linewidth=2.)
    if chr_stops_dict is not None:
        chr_stops_chrs = list(chr_stops_dict.keys())
        chr_stops_bins = list(chr_stops_dict.values())
        chr_means = []
        chr_means.append(chr_stops_bins[0] / 2)
        for c in range(1, len(chr_stops_bins)):
            aux = (chr_stops_bins[c] + chr_stops_bins[c - 1]) / 2
            chr_means.append(aux)
        chrs_ = deepcopy(chr_stops_chrs)
        chrs_[12] = f'{chr_stops_chrs[12]} '
        chrs_[12]
        chrs_[20] = ""
        chrs_[21] = f'  {chr_stops_chrs[21]}'
        plt.xticks(chr_means, np.array(chrs_), rotation=0, fontsize=tickfontsize)
        ax.vlines(list(chr_stops_dict.values())[:-1], *ax.get_ylim(), ls='--', linewidth=2.5)
        plt.xlabel('Chromosomes', fontsize=textfontsize)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(tickfontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(tickfontsize)
    if labels is not None:
        plt.yticks([])

    axins = inset_axes(
            ax,
            width="2%",  # width = 5% of parent_bbox width
            height="85%",
            loc="lower left",
            bbox_to_anchor=(1.045, 0.0, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
    cb = plt.colorbar(im, cax=axins)
    if vmax is not None:
        im.set_clim(vmin=0, vmax=vmax)
    cb.ax.tick_params(labelsize=tickfontsize)
    cb.outline.set_visible(False)
    cb.ax.set_title(cbar_title, y=1.05, fontsize=textfontsize)

    if mode == 'cnv':
        if vmax is None or vmax == 4 :
            im.set_clim(vmin=0, vmax=4)
            cb.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
            cb.set_ticklabels(["0", "1", "2", "3", "4+"])

    if output_path is not None:
        print("Creating {}...".format(output_path))
        plt.savefig(output_path, bbox_inches="tight")
        plt.close()
        print("Done.")
    else:
        plt.show()
