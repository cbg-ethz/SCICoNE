import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle, Polygon
from tqdm import tqdm
import os
import re
import subprocess

sns.set_style("ticks", {"axes.grid": True})
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

def plot_bins(
    bins,
    chr_stops_dict=None,
    cbar_title="",
    ylabel="Cells",
    annotations=None,
    bps=None,
    cluster=False,
    figsize=(14, 4),
    vlines=True,
    clone_cmap=None,
    vmin=0,
    vmax=2,
    title="",
    final_normal_cell=None,
    mid_normal_cell=None,
    mid_malignant_cell=None,
    norm=None,
    output_path=None,
    dpi=300,
):
    if cluster:
        if annotations is not None:
            # Within each clone, sort the cells's normalised regions via hierarchical clustering
            annotations = np.array(annotations)
            for c_id in np.unique(annotations):
                Z = ward(pdist(bins[np.where(annotations == c_id)[0]]))
                hclust_index = leaves_list(Z)
                bins[np.where(annotations == c_id)[0]] = bins[
                    np.where(annotations == c_id)[0]
                ][hclust_index]
        else:
            Z = ward(pdist(bins))
            hclust_index = leaves_list(Z)
            bins = bins[hclust_index]

    if annotations is not None:
        ticks = dict()
        annotations = np.array(annotations).ravel()
        for c_id in np.unique(annotations):
            # Get last pos
            t = np.where(annotations == c_id)[0]
            if len(t) > 1:
                t = t[-1]
            ticks[c_id] = t
    if clone_cmap is None:
        clone_cmap = matplotlib.cm.get_cmap("Dark2")
        subset_colors = clone_cmap(np.arange(0, len(np.unique(annotations)), 1))
        clone_cmap = matplotlib.colors.ListedColormap(subset_colors)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    if annotations is not None:
        gs = GridSpec(1, 2, wspace=0.05, width_ratios=[1, 40])
        ax1 = fig.add_subplot(gs[0])
        bounds = [0] + list(ticks.values())
        subnorm = matplotlib.colors.BoundaryNorm(bounds, clone_cmap.N)
        cb = matplotlib.colorbar.ColorbarBase(
            ax1,
            cmap=clone_cmap,
            norm=subnorm,
            boundaries=bounds,
            spacing="proportional",
            orientation="vertical",
        )
        cb.outline.set_visible(False)
        cb.ax.set_ylabel("Clones")
        ax1.yaxis.set_label_position("left")
        for j, lab in enumerate(ticks.keys()):
            cb.ax.text(
                0.5,
                ((bounds[j + 1] + bounds[j]) / 2) / bounds[-1],
                int(lab),
                ha="center",
                va="center",
                rotation=90,
                color="w",
            )
        cb.set_ticks([])

        ax2 = fig.add_subplot(gs[1])
    else:
        ax2 = plt.gca()

    if len(np.unique(bins)) < 6:
        vmax = 4
        cmap = matplotlib.colors.ListedColormap(sns.diverging_palette(220, 10, n=5))
    else:
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
    im = plt.pcolormesh(bins, cmap=cmap, clim=(vmin, vmax), norm=norm)
    plt.ylabel(ylabel)
    plt.yticks([])
    plt.title(title)

    if bps is not None:
        ax2.vlines(bps, *ax2.get_ylim(), colors="b", linestyles="dashed")

    plt.xlabel("Bins")
    if chr_stops_dict is not None:
        plt.xlabel("Chromosomes")
        chr_stops_chrs = list(chr_stops_dict.keys())
        chr_stops_bins = list(chr_stops_dict.values())

        chr_means = []
        chr_means.append(chr_stops_bins[0] / 2)
        for c in range(1, len(chr_stops_bins)):
            aux = (chr_stops_bins[c] + chr_stops_bins[c - 1]) / 2
            chr_means.append(aux)

        if vlines:
            ax2.vlines(chr_stops_bins[:-1], *ax2.get_ylim(), ls="--")
            plt.xticks(chr_means, chr_stops_chrs, rotation=0)
        else:
            xticks = [0] + chr_stops_bins
            xticks_labels = chr_stops_chrs
            ax2.xaxis.set_major_locator(
                ticker.FixedLocator(xticks)
            )  # locate major ticks
            ax2.xaxis.set_minor_locator(
                ticker.FixedLocator(chr_means)
            )  # locate minor ticks

            ax2.xaxis.set_major_formatter(
                ticker.NullFormatter()
            )  # hide major tick labels
            ax2.xaxis.set_minor_formatter(
                ticker.FixedFormatter(xticks_labels)
            )  # show minor labels

            for tick in ax2.xaxis.get_minor_ticks():
                tick.tick1line.set_markersize(0)
                tick.tick2line.set_markersize(0)
                tick.label1.set_horizontalalignment("center")

    if final_normal_cell is not None:
        ax2.hlines(final_normal_cell, *ax2.get_xlim(), ls="-")
    if mid_normal_cell is not None and mid_malignant_cell is not None:
        plt.yticks(
            [mid_malignant_cell, mid_normal_cell],
            ["tumor", "immune"],
            va="center",
            rotation=90,
        )

    axins = inset_axes(
        ax2,
        width="2%",  # width = 5% of parent_bbox width
        height="85%",
        loc="lower left",
        bbox_to_anchor=(1.05, 0.0, 1, 1),
        bbox_transform=ax2.transAxes,
        borderpad=0,
    )
    cb = plt.colorbar(im, cax=axins)
    if len(np.unique(bins)) < 6:
        im.set_clim(vmin=0, vmax=4)
        cb.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
        cb.set_ticklabels(["0", "1", "2", "3", "4+"])
    else:
        if vmax is not None and vmin is not None:
            im.set_clim(vmin=vmin, vmax=vmax)
        else:
            cb.set_ticks([0, 1, 2])
            cb.set_ticklabels(["0", "1", "2+"])
    cb.outline.set_visible(False)
    cb.ax.set_title(cbar_title, y=1.01)

    sns.despine(left=True, right=True, top=True, bottom=True)

    if output_path is not None:
        print("Creating {}...".format(output_path))
        plt.savefig(output_path, bbox_inches="tight")
        plt.close()
        print("Done.")
    else:
        plt.show()


def _plot_bins(
    bins,
    chr_stops_dict,
    bps=None,
    cluster=False,
    figsize=(14, 4),
    vlines=True,
    cmap=None,
    vmax=2,
    title="Normalised counts per filtered bin",
    output_path=None,
):
    chr_stops_chrs = list(chr_stops_dict.keys())
    chr_stops_bins = list(chr_stops_dict.values())

    if cluster:
        Z = ward(pdist(bins))
        hclust_index = leaves_list(Z)
        bins = bins[hclust_index]

    if cmap is None:
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        if len(np.unique(bins)) < 6:
            cmap = matplotlib.colors.ListedColormap(sns.diverging_palette(220, 10, n=5))
            vmax = 4

    chr_means = []
    chr_means.append(chr_stops_bins[0] / 2)
    for c in range(1, len(chr_stops_bins)):
        aux = (chr_stops_bins[c] + chr_stops_bins[c - 1]) / 2
        chr_means.append(aux)

    fig = plt.figure(figsize=figsize, dpi=300)
    im = plt.pcolormesh(bins, cmap=cmap, vmax=vmax, vmin=0, clip_on=False)
    plt.xlabel("chromosomes")
    plt.ylabel("cells")
    cb = fig.colorbar(im)
    if len(np.unique(bins)) < 6:
        if vmax is not None:
            im.set_clim(vmin=0, vmax=vmax)
        cb.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
        cb.set_ticklabels(["0", "1", "2", "3", "4+"])
    cb.outline.set_visible(False)
    sns.despine(left=True, bottom=True, right=True, top=True)
    plt.title(title)

    ax = plt.gca()
    if bps is not None:
        ax.vlines(bps, *ax.get_ylim(), colors="b", linestyles="dashed")

    if vlines:
        ax.vlines(chr_stops_bins[:-1], *ax.get_ylim(), ls="--")
        plt.xticks(chr_means, chr_stops_chrs, rotation=0)
    else:
        xticks = [0] + chr_stops_bins
        xticks_labels = chr_stops_chrs
        ax.xaxis.set_major_locator(ticker.FixedLocator(xticks))  # locate major ticks
        ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_means))  # locate minor ticks

        ax.xaxis.set_major_formatter(ticker.NullFormatter())  # hide major tick labels
        ax.xaxis.set_minor_formatter(
            ticker.FixedFormatter(xticks_labels)
        )  # show minor labels

        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)
            tick.label1.set_horizontalalignment("center")

    if output_path is not None:
        print("Creating {}...".format(output_path))
        plt.savefig(output_path, bbox_inches="tight")
        plt.close()
        print("Done.")
    else:
        plt.show()


def convert_event_region_to_gene(
    region_event_str,
    bin_gene_region_df,
    priority_only=False,
    genes_to_highlight=None,
    highlight_color="red",
    genes_only=False,
):
    """
        Returns a string indicating gene-wise events in affected region
        Examples:
                "+2R174"     -> ["+2BRAF", "+2MALAT1"]
                "-1R656:658" -> ["-1JAK2", "-1MLNA", "-1CDK4"]
        :param region_event_str: str
        :param bin_gene_region_df: DataFrame
            with (gene_name, region, original_bin) fields
        :param priority_only: boolean
            indicating if only priority genes should be returned
        :param genes_to_highlight: list
            genes that should be displayed in a different color
        :param highlight_color: str
            color to use in genes to highlight
        :return: list of str
    """
    # Get event (-2, -1, +1, +2, etc)
    event_str = region_event_str[:2]
    region_str = region_event_str[3:]
    if ":" in region_str:  # multiple regions: "-1R656:658"
        aux = [int(region) for region in region_str.split(":")]
        region_list = np.arange(aux[0], aux[1] + 1)
    else:
        region_list = [int(region_str)]

    gene_list = []
    for region in region_list:
        genes_in_region = get_genes_in_region(
            region, bin_gene_region_df, priority_only=priority_only
        )
        gene_list.append(genes_in_region)

    gene_list = [item for sublist in gene_list for item in sublist]

    # Highlight some genes
    if genes_to_highlight is not None:
        for index, gene in enumerate(gene_list):
            if gene in genes_to_highlight:
                gene_list[index] = (
                    "<font color="
                    + "'"
                    + highlight_color
                    + "'"
                    + ">"
                    + gene
                    + "</font>"
                )

    gene_string = "[" + ",".join(gene_list) + "]"
    if len(gene_list) == 0:
        gene_event_str = ""
    else:
        if not genes_only:
            gene_event_str = event_str + gene_string
        else:
            gene_event_str = gene_string

    return gene_event_str


def convert_node_regions_to_genes(
    node_str,
    bin_gene_region_df,
    priority_only=False,
    genes_to_highlight=None,
    highlight_color="red",
    max_genes_per_line=10,
):
    """
            Returns a string indicating gene events and total number of
            amplifications and deletions in node
            Examples:
                    "+2R174 +2R291:293 -1R656:658" -> "+2[BRAF,MALAT1] -1[JAK2,MLNA,CDK4]\n(3+, 1-)"
            :param node_str: str
            :param bin_gene_region_df: DataFrame
                with (gene_name, region, original_bin) fields
            :param priority_only: boolean
                indicating if only priority genes should be returned
            :param genes_to_highlight: list
                genes that should be displayed in a different color
            :return: str
        """
    region_event_strs = node_str.split(" ")
    num_events = len(region_event_strs)

    num_amplifications = 0

    str_dict = dict()
    possible_events = ["+4", "+3", "+2", "+1", "-1", "-2", "-3", "+4"]
    for event in possible_events:
        str_dict[event] = []

    for region_event_str in region_event_strs:
        # Get event (-2, -1, +1, +2, etc)
        event_str = region_event_str[:2]
        region_str = region_event_str[3:]
        if ":" in region_str:  # multiple regions: "-1R656:658"
            aux = [int(region) for region in region_str.split(":")]
            region_list = np.arange(aux[0], aux[1] + 1)
        else:
            region_list = [int(region_str)]

        gene_list = []
        for region in region_list:
            genes_in_region = get_genes_in_region(
                region, bin_gene_region_df, priority_only=priority_only
            )
            gene_list.append(genes_in_region)

        gene_list = [item for sublist in gene_list for item in sublist]

        str_dict[event_str].append(gene_list)

        if region_event_str[0] == "+":
            num_amplifications += 1

    for key in str_dict:
        str_dict[key] = [item for sublist in str_dict[key] for item in sublist]

    node_str = ""
    newline = "<br/>"

    num_deletions = num_events - num_amplifications
    num_events, num_amplifications, num_deletions
    num_events_str = "{}+, {}-".format(num_amplifications, num_deletions)

    node_str = (
        "<i> </i><font point-size='30'>"
        + "("
        + num_events_str
        + ")"
        + "</font><i> </i>"
        + newline
        + " "
        + newline
    )

    i = 0
    for key in str_dict:
        i += 1
        if len(str_dict[key]) != 0:
            # Show only one occurence of each gene
            str_dict[key] = np.unique(np.array(str_dict[key])).tolist()
            color = "blue"
            if key[0] == "+":
                color = "red"
            node_str = (
                node_str
                + "<font color='{}' point-size='26'><b>".format(color)
                + key
                + "</b></font>"
            )
            node_str = (
                node_str
                + "["
                + ",".join(
                    f"{x}" + newline if (i + 1) % max_genes_per_line == 0 else str(x)
                    for i, x in enumerate(str_dict[key])
                )
                + "]"
            )
            node_str = "<i> </i>" + node_str + "<i> </i>"
            if i < len(str_dict.keys()):
                node_str = node_str + " " + newline + " " + newline

    # If newline followed by ']', replace with newline after ']'
    node_str = node_str.replace(newline + "]", "]" + newline)

    # If newline followed by ',', replace with newline after ','
    node_str = node_str.replace(newline + ",", "," + newline)

    # If newline followed by newline, remove one
    for m in re.finditer(newline, node_str):
        index = m.start()
        if node_str[index + len(newline) : index + 2 * len(newline)] == newline:
            node_str = "".join(
                (
                    node_str[: index + len(newline)],
                    "",
                    node_str[index + 2 * len(newline) :],
                )
            )

    # highlight genes
    if genes_to_highlight is not None:
        for gene in genes_to_highlight:
            node_str = node_str.replace(
                gene,
                "<b><font color="
                + "'"
                + highlight_color
                + "'"
                + ">"
                + gene
                + "</font></b>",
            )

    if node_str == "":
        node_str = num_events_str

    return node_str


def tree_to_graphviz(
    tree_path,
    node_sizes=None,
    gene_labels=False,
    bin_gene_region_df=None,
    genes_to_highlight=None,
    highlight_color="yellow",
    max_genes_per_line=6,
    output_path=None,
):
    """
        reads the file containing trees converts it to graphviz format
        :param tree_path: path to the tree file.
        :param node_sizes: dictionary containing the size of each node
        :param gene_labels: whether to label nodes with genes
        :param bin_gene_region_df: Pandas DataFrame with gene-region correspondence
        :param genes_to_highlight: List containing genes to highlight
        :return: string object containing the graphviz formatted tree
    """
    with open(tree_path) as f:
        list_tree_file = list(f)

    graphviz_header = [
        "digraph { \n",
        'node [style=filled,color="#D4C0D6",fontsize=20,margin=0,shape=oval]'
        'edge [arrowhead=none, color="#602A86"]',
    ]

    graphviz_labels = []
    graphviz_links = []

    graphviz_labels.append("0[label=<<font point-size='30'> Neutral </font>>]")  # root

    for line in list_tree_file:
        if line.startswith("node 0:"):
            continue
        elif line.startswith("node"):
            comma_splits = line.split(",")

            comma_first = re.split(" |:", comma_splits[0])
            node_id = comma_first[1]
            p_id = comma_first[4]
            comma_rest = comma_splits[1:]
            comma_rest[0] = comma_rest[0].lstrip("[")
            comma_rest[-1] = comma_rest[-1].rstrip("]\n")
            merged_labels = []
            [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
            k_end = k_begin
            for term in comma_rest[1:]:  # events vector
                [k, v] = (int(x) for x in term.split(":"))
                if k == k_end + 1 and v == previous_v:
                    k_end = k  # update the end
                else:
                    if k_begin == k_end:
                        merged_labels.append(f"{previous_v:+}R{k_begin}")
                    else:
                        merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
                    k_begin = k_end = k
                previous_v = v
            # print the last one
            if k_begin == k_end:
                merged_labels.append(f"{previous_v:+}R{k_begin}")
            else:
                merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")

            str_merged_labels = (
                "<i> </i>"
                + " ".join(
                    f"{x}<i> </i><br/>" if i % 10 == 0 and i > 0 else str(x)
                    for i, x in enumerate(merged_labels)
                )
                + "<i> </i>"
            )
            if gene_labels and bin_gene_region_df is not None:
                node_str = " ".join(merged_labels)  # "+1R75 +1R218:219 +1R221:223"
                str_merged_labels = convert_node_regions_to_genes(
                    node_str,
                    bin_gene_region_df,
                    priority_only=True,
                    genes_to_highlight=genes_to_highlight,
                    highlight_color=highlight_color,
                    max_genes_per_line=max_genes_per_line,
                )

            newline = "<br/>"
            endline = "<i> </i>"

            # If there are newlines followed by endlines, switch
            str_merged_labels = str_merged_labels.replace(
                newline + endline, endline + newline
            )

            # Remove whatever comes after the last endline position
            new_end_pos = (
                [m.start() for m in re.finditer(endline, str_merged_labels)][-1]
                + len(endline)
                - 1
            )
            if len(str_merged_labels) > new_end_pos + 1:
                str_merged_labels = str_merged_labels[: new_end_pos + 1]

            # Replace multiple endlines with one
            while endline * 2 in str_merged_labels:
                str_merged_labels = str_merged_labels.replace(endline * 2, endline)

            # # If node string ends with newlines, remove them
            # par = " " + newline + " " + newline
            # l = list(str_merged_labels)
            # if l[-len(par):] == list(par):
            #     str_merged_labels = ''.join(l[:-len(par)])

            # Add node size
            if node_sizes is not None:
                try:
                    node_size = node_sizes[node_id]
                except:
                    node_size = 0
                str_merged_labels = str_merged_labels + " " + newline + " " + newline
                str_merged_labels = (
                    str_merged_labels
                    + "<font point-size='26'>"
                    + str(int(node_size))
                    + " cell"
                )
                if int(node_size) > 1 or int(node_size) == 0:
                    str_merged_labels = str_merged_labels + "s"
                str_merged_labels = str_merged_labels + "</font>"

            graphviz_labels.append(
                f"{node_id}[label=<{str_merged_labels}>]"
            )  # use < > to allow HTML
            graphviz_links.append(f"{p_id} -> {node_id}")

    txt = graphviz_header + graphviz_labels + graphviz_links + ["}"]

    if output_path is not None:
        with open(output_path, "w") as file:
            for line in txt:
                file.write(f"{line}\n")

    return txt


def plot_tree_graphviz(tree_graphviz_path, output_path):
    try:
        format = output_path.split(".")[-1]
    except:
        format = "png"

    try:
        cmd_output = subprocess.run(
            [
                "dot",
                f"-T{format}:cairo",
                f"{tree_graphviz_path}",
                "-o",
                f"{output_path}",
            ]
        )
    except subprocess.SubprocessError as e:
        print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
    else:
        print(f"subprocess out: {cmd_output}")
        print(f"stdout: {cmd_output.stdout}\n stderr: {cmd_output.stderr}")


def plot_heatmap(gene_cn_df, output_path=None):
    if "is_imputed" in gene_cn_df.columns:
        is_imputed = gene_cn_df["is_imputed"]
        gene_cn_df = gene_cn_df.drop(columns=["is_imputed"])

    if np.all(~np.isnan(gene_cn_df)):
        annot = np.array(gene_cn_df.astype(int).astype(str))
    else:
        annot = np.array(gene_cn_df.astype(str))
    annot[np.where(gene_cn_df == 4)] = ["4+"] * len(np.where(gene_cn_df == 4)[0])
    annot = annot.astype(str)

    figure_width = gene_cn_df.shape[0] / 2 + 1.5
    plt.figure(figsize=(8, figure_width), dpi=300)
    cmap = matplotlib.colors.ListedColormap(sns.diverging_palette(220, 10, n=5))
    heatmap = sns.heatmap(
        gene_cn_df,
        annot=annot,
        cmap=cmap,
        vmin=0,
        vmax=4,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"ticks": [0, 1, 2, 3, 4]},
        fmt="",
    )
    colorbar = heatmap.collections[0].colorbar
    colorbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
    colorbar.set_ticklabels(["0", "1", "2", "3", "4+"])
    heatmap.set_title("Copy number values of genes per subclone")
    heatmap.set_facecolor("#656565")
    plt.xlabel("Subclone")
    # b, t = plt.ylim()  # discover the values for bottom and top
    # b += 0.5  # Add 0.5 to the bottom
    # t -= 0.5  # Subtract 0.5 from the top
    # plt.ylim(b, t)  # update the ylim(bottom, top) values

    # ax = heatmap.ax_heatmap
    if is_imputed is not None:
        for i in range(is_imputed.shape[0]):
            if is_imputed.iloc[i]:
                # heatmap.add_patch(Rectangle((j, is_imputed.shape[0]-1-i), closed=True, fill=False, edgecolor='gray', lw=3))
                box = np.array(
                    [
                        [0, i],
                        [gene_cn_df.shape[1], i],
                        [gene_cn_df.shape[1], i + 1],
                        [0, i + 1],
                    ]
                )
                heatmap.add_patch(
                    Polygon(
                        box,
                        closed=True,
                        fill=False,
                        edgecolor="gray",
                        lw=1.5,
                        ls="--",
                        clip_on=False,
                    )
                )

    if output_path is not None:
        plt.savefig(output_path, bbox_inches="tight")
    else:
        plt.show()


def plot_profile(
    cnv_arr,
    chr_stops,
    subclone_ids=None,
    offset_sizes=0.1,
    ymax=5,
    s=1,
    cmap=None,
    colors_idx=None,
    chr_mean_pos=None,
    figsize=(14, 4),
    yticksize=11,
    ncol=None,
    output_path=None,
):
    """
    Plots CNV profiles
    :param cnv_arr: (n_subclones, n_bins) array of CNVs
    :param chr_stops: (chr, pos) DataFrame with final bin of each chromosome
    (optional) :param subclone_ids: list of string identifiers for each clone
    :return: axis
    """
    if len(cnv_arr.shape) == 1:
        cnv_arr = cnv_arr.reshape(1, -1)

    if chr_mean_pos is None:
        # Customize minor tick labels
        chr_mean_pos = []
        chr_mean_pos.append(chr_stops["Unnamed: 0"][0] / 2)
        for c in range(1, 24):
            aux = (chr_stops["Unnamed: 0"][c] + chr_stops["Unnamed: 0"][c - 1]) / 2
            chr_mean_pos.append(aux)

    n_subclones = cnv_arr.shape[0]
    n_bins = cnv_arr.shape[1]

    if subclone_ids is not None:
        if not (isinstance(subclone_ids, list)):
            subclone_ids = [subclone_ids]
        assert len(subclone_ids) == n_subclones
    else:
        subclone_ids = np.arange(n_subclones).astype(str).tolist()

    if cmap is None:
        cmap = matplotlib.cm.get_cmap("Dark2")
    if colors_idx is None:
        colors_idx = np.arange(n_subclones)

    offsets = [0.0] * n_subclones
    if n_subclones > 1:
        offset = offset_sizes
        offsets = [offset * c for c in range(n_subclones)]
        offsets = offsets - np.mean(offsets)

    fig = plt.figure(figsize=figsize, dpi=300)

    for c in range(n_subclones):
        plt.scatter(
            range(n_bins),
            cnv_arr[c].astype(int) + offsets[c],
            s=s,
            label="subclone {}".format(subclone_ids[c]),
            c=np.array(cmap(colors_idx[c])).reshape(1, -1),
            alpha=1,
        )
        plt.ylim([0 - 0.2, ymax])
        plt.xlim([0 - 0.01 * n_bins, n_bins + 0.01 * n_bins])
        plt.xlabel("chromosome")
        plt.ylabel("copy number")

    plt.tick_params(axis="y", which="major", labelsize=yticksize)

    ax = plt.gca()
    xticks = [0] + list(chr_stops["Unnamed: 0"])
    xticks_labels = list(chr_stops["chr"])
    ax.xaxis.set_major_locator(ticker.FixedLocator(xticks))  # locate major ticks
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_mean_pos))  # locate minor ticks

    ax.xaxis.set_major_formatter(ticker.NullFormatter())  # hide major tick labels
    ax.xaxis.set_minor_formatter(
        ticker.FixedFormatter(xticks_labels)
    )  # show minor labels

    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment("center")

    if n_subclones > 1:
        if ncol is None:
            ncol = n_subclones
        lgnd = plt.legend(
            bbox_to_anchor=(0, 1.02, 1, 0.2),
            loc="lower center",
            borderaxespad=0,
            ncol=ncol,
        )
        for c in range(n_subclones):
            lgnd.legendHandles[c]._sizes = [40]
            lgnd.legendHandles[c].set_alpha(1)
    else:
        plt.title("Subclone {}".format(subclone_ids[0]))

    if output_path is not None:
        print("Creating {}...".format(output_path))
        plt.savefig(output_path, bbox_inches="tight")
        print("Done.")
    else:
        plt.show()


def plot_cluster_cnvs(
    cnvs_arr,
    chr_stops,
    subclone_ids=None,
    offset_sizes=0.0,
    s=1,
    ymax=None,
    ncol=None,
    output_path=None,
):
    if len(cnvs_arr.shape) == 1:
        cnvs_arr = cnvs_arr.reshape(1, -1)

    n_subclones = cnvs_arr.shape[0]

    if subclone_ids is not None:
        if not (isinstance(subclone_ids, list)):
            subclone_ids = [subclone_ids]
        assert len(subclone_ids) == n_subclones
    else:
        subclone_ids = np.arange(n_subclones).astype(str).tolist()

    if ymax is None:
        ymax = np.nanmax(cnvs_arr) + 1

    # Customize minor tick labels
    chr_mean_pos = []
    chr_mean_pos.append(chr_stops["Unnamed: 0"][0] / 2)
    for c in range(1, 24):
        aux = (chr_stops["Unnamed: 0"][c] + chr_stops["Unnamed: 0"][c - 1]) / 2
        chr_mean_pos.append(aux)

    cnvs_arr[np.where(cnvs_arr != np.nan)] = cnvs_arr[
        np.where(cnvs_arr != np.nan)
    ].astype(int)

    # Plot each profile separately
    for c in range(n_subclones):
        cluster_path = (
            output_path + "__cluster_profile_" + str(c) + ".png"
            if output_path is not None
            else None
        )
        plot_profile(
            cnvs_arr[c],
            chr_stops,
            subclone_ids=subclone_ids[c],
            s=s,
            ymax=ymax,
            colors_idx=[c],
            chr_mean_pos=chr_mean_pos,
            figsize=(14, 4),
            output_path=cluster_path,
        )

    # Plot all profiles overlapping
    overlapping_path = (
        output_path + "__cluster_profile_overlapping.png"
        if output_path is not None
        else None
    )
    plot_profile(
        cnvs_arr,
        chr_stops,
        subclone_ids=subclone_ids,
        s=s,
        ymax=ymax,
        offset_sizes=offset_sizes,
        colors_idx=np.arange(n_subclones),
        chr_mean_pos=chr_mean_pos,
        figsize=(14, 8),
        yticksize=13,
        ncol=ncol,
        output_path=overlapping_path,
    )
