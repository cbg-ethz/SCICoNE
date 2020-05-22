import numpy as np
from copy import deepcopy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

sns.set_style("ticks", {"axes.grid": True})

colors = ["#2040C8", "white", "#EE241D"]
datacmap = matplotlib.colors.LinearSegmentedColormap.from_list("cmap", colors)
cnvcmap = matplotlib.colors.LinearSegmentedColormap.from_list("cnvcmap", colors, 5)

def plot_matrix(data, cbar_title="", mode='data', chr_stops_dict=None,
                textfontsize=24, tickfontsize=22, bps=None,
                figsize=(24,8), dpi=100, vmax=None, output_path=None):
    if mode == 'data':
        cmap = datacmap
    elif mode == 'cnv':
        cmap = cnvcmap
    else:
        raise AttributeError('mode argument must be one of \'data\' or \'cnv\'')

    fig = plt.figure(figsize=figsize, dpi=dpi)
    im = plt.pcolormesh(data, cmap=cmap, rasterized=True)
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
