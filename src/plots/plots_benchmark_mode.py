import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.switch_backend('agg')
sns.set(style="white")
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from statannot import add_stat_annotation
from src.plots.plots_utils import *


def plot_density_by_class(data: pd.DataFrame,
                          thresholds: list,
                          fname: str,
                          min_predicted: int = 20):
    """
    Plot class distribution for a given tool taking
    into account the reference threshold
    :param pd.DataFrame data: Scores and class
        for a single tool
    :param list thresholds: Reference thresholds
    :param str fname: Output file name
    :param int min_predicted: Minimum number of variants
    predicted to draw plots. Default: `20`
    """

    fraction_nan = round(np.sum(np.sum(data.iloc[:, 0].isnull())) / data.shape[0] * 100, 2)
    df_ = data[~data.iloc[:, 0].isnull()]

    if df_.shape[0] > min_predicted:
        pathogenic = df_.loc[df_['label']]
        benign = df_.loc[~df_['label']]

        sns.kdeplot(benign.iloc[:, 0],
                    shade=True,
                    linewidth=3,
                    label="Benign")

        sns.kdeplot(pathogenic.iloc[:, 0],
                    shade=True,
                    linewidth=3,
                    label="Pathogenic")

        plt.axvline([x[2] for x in thresholds if x[0] == list(df_)[0]], color='r', linestyle="--")

        thresh_line = Line2D([], [], color='r', linestyle='--',
                             markersize=10,
                             markeredgewidth=1.5,
                             label='Reference threshold')
        plt.legend(handles=[thresh_line])

        plt.xlabel(list(df_)[0] + " ({}% missing data)".format(fraction_nan))
        ends = (max(df_.iloc[:, 0]) - min(df_.iloc[:, 0])) * 0.05
        plt.xlim(min(df_.iloc[:, 0]) - ends, max(df_.iloc[:, 0]) + ends)
        plt.savefig(fname + ".pdf")
        plt.close()


def plot_allele_frequency(df: pd.DataFrame,
                          fname: str,
                          gnomad_col: str = "gnomAD_genomes"):
    """
    Plots allele frequencies for each class

    :param pd.DataFrame df: Input df
    :param str fname: Output basemame
    :param str gnomad_col: Column name that accounts
        for allele frequencies. Default: `gnomAD_genomes`.
        If column does not exist, analysis will be skipped.
    """
    if gnomad_col not in df.columns:
        return

    df[gnomad_col] = pd.to_numeric(df[gnomad_col])
    df['grouper'] = df['outcome'].astype(str) + '\nN = ' + df['count_class'].astype(str)
    order = sorted(list(df['grouper'].unique()))
    ax = sns.boxplot(data=df, x="grouper", order=order, y=gnomad_col)
    try:
        add_stat_annotation(ax, data=df, x="grouper", y=gnomad_col,
                            order=order,
                            box_pairs=[tuple(order)],
                            test='Mann-Whitney',
                            text_format='star',
                            loc='inside',
                            verbose=0,
                            pvalue_format_string='{:.4f}')
        plt.xlabel("")
        plt.ylabel("Allele frequency")
        plt.tight_layout()
        out = fname + '.pdf'
        plt.savefig(out)
        plt.close()
    except ValueError:
        plt.close()
        pass


def plot_unscored(data: pd.DataFrame, fname: str):
    """
    Plot fraction of unscored variants by each
    tool in `data`.

    :param pd.DataFrame data: Df with stats
        for each tool
    :param str fname: Output file
    """

    ax = sns.barplot(x="fraction_nan",
                     y="tool",
                     color="skyblue",
                     data=data)

    ax.set(xlabel='Fraction of unscored variants', ylabel='')
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_metrics(data: pd.DataFrame, fname: str, metric: str):
    """
    Plot tools ranked by a given metric

    :param pd.DataFrame data: Df with stats
        for each tool
    :param str fname: Output name
    :param str metric: Metric used
        to rank the tools
    """
    my_range_coverage = list(range(1, len(data.index) + 1))
    my_range_specificity = np.arange(1 - 0.2, len(data.index)).tolist()
    my_range_sensitivity = np.arange(1 + 0.2, len(data.index) + 0.5).tolist()

    data = data.sort_values(metric)

    fig, ax = plt.subplots(figsize=(10, 8)) if data.shape[0] > 25 else plt.subplots(figsize=(8, 6))

    _target_col = "specificity" if "accuracy" in metric else "precision"
    plt.scatter(data[_target_col], my_range_specificity,
                color='skyblue',
                alpha=1,
                marker='s',
                edgecolors='black',
                linewidths=0.5,
                label=_target_col.capitalize())

    plt.scatter(data['sensitivity'], my_range_sensitivity,
                color='yellow',
                alpha=0.75,
                marker='v',
                edgecolors='black',
                linewidths=0.5,
                label='Sensitivity')

    plt.scatter(data['coverage'],
                my_range_coverage,
                color='brown',
                alpha=0.75,
                marker='o',
                edgecolors='black',
                linewidths=0.5,
                label='Fraction_predictable')

    rectangule_attr = [_target_col, 'sensitivity', 'coverage']
    i = 0
    for idx, rows in data.iterrows():
        val_range = list(rows[rectangule_attr])
        width = [min(val_range), max(val_range)]
        height = 1
        ax.add_patch(
            plt.Rectangle(xy=(width[0] - 0.01, my_range_coverage[i] - 0.4),
                          width=width[1] - width[0] + 0.02,
                          height=height,
                          linewidth=0.75,
                          color='gray',
                          fill=False))
        i += 1

    ax.grid(axis='x', linestyle='dashed')
    plt.legend(loc="upper right",
               bbox_to_anchor=(1.4, 1),
               borderaxespad=0,
               prop=dict(size=8))
    plt.subplots_adjust(left=0.35)

    plt.yticks(my_range_coverage, data['tool'] + " (" + data[metric].astype(str) + ")")

    if all(col in data.columns for col in ['total', 'total_p', 'total_n']):
        plt.title("#variants: {} ({} pos, {} neg)".format(data["total"].iloc[0],
                                                          data["total_p"].iloc[0],
                                                          data["total_n"].iloc[0]))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools_total_correct(data: pd.DataFrame, fname: str):
    """
    Plots the total number of correct
    classification and ranks the display
    by the top tool

    :param pd.DataFrame data: Df with stats
        for each tool
    :param str fname: Output file
    """
    df = data.iloc[::-1]
    set_style()
    ind = np.arange(df.shape[0])
    df.sort_values(by=['correct'], inplace=True)
    w = 0.8

    fig, ax = plt.subplots()
    plt.xlim(df['total'][0])
    ax.barh(ind, df['correct'],
            align='center',
            color='green',
            zorder=10,
            alpha=0.8,
            height=w,
            edgecolor="black",
            linewidth=0.75)

    ax.set(yticks=ind, yticklabels=df['tool'])
    ax.grid(True)
    ax.margins(0.00)
    ax.invert_xaxis()
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xlabel('# Number of correctly predicted variants')
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools_barplot_only_correct(data: pd.DataFrame, fname: str):
    """
    Plot tools barplot for correct
    prediction alone separated by
    class labels where tools are
    sorted by the correct number of
    pathogenic variants identified

    :param pd.DataFrame data: Df with stats
        for each tool
    :param str fname: Output name
    """
    df = data.iloc[::-1]
    set_style()

    w = 0.5
    fig, axes = plt.subplots(ncols=2, sharey='all')
    df.sort_values('tp', inplace=True, ascending=True)
    ind = np.arange(df.shape[0])

    axes[0].barh(ind, df['tp'],
                 align='center',
                 color='green',
                 zorder=10,
                 alpha=0.8,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    axes[1].barh(ind, df['tn'],
                 align='center',
                 color='green',
                 zorder=10,
                 alpha=0.8,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    try:
        axes[0].set_xlim(0, df["tp"][0] + df["fn"][0])
        axes[1].set_xlim(0, df["tn"][0] + df["fp"][0])
    except KeyError:
        pass

    axes[0].set(yticks=ind, yticklabels=df['tool'])
    axes[0].set_xlabel('# Pathogenic Variants')
    axes[1].set_xlabel('# Benign Variants')

    for ax in axes.flat:
        ax.grid(False)
        ax.margins(0.00)
        ax.tick_params(axis='both', which='major', labelsize=8)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools_barplot(data: pd.DataFrame, fname: str, metric: str):
    """
    Tools performance barplot with two
    axis (pathogenic and benign).

    :param pd.DataFrame data: Input df
    :param str fname: Output basename
    :param str metric: Metric to evaluate
    """

    # revert order
    df = data.iloc[::-1]

    set_style()
    w = 0.8

    fig, axes = plt.subplots(ncols=2,
                             sharey='all',
                             figsize=(15, 15))
    ind = np.arange(df.shape[0])

    # Pathogenic variants
    axes[0].barh(ind, df['tp'],
                 align='center',
                 color='darkblue',
                 alpha=0.7,
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    axes[0].barh(ind, df['fn'],
                 left=df['tp'],
                 align='center',
                 alpha=0.8,
                 color='darkred',
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    axes[0].barh(ind, df['mp'],
                 left=df['scored_p'],
                 align='center',
                 color='lightgrey',
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 alpha=0.9,
                 linewidth=0.75)

    # Benign variamts
    axes[1].barh(ind, df['tn'],
                 align='center',
                 color='darkblue',
                 alpha=0.7,
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    axes[1].barh(ind, df['fp'],
                 left=df['tn'],
                 align='center',
                 alpha=0.8,
                 color='darkred',
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 linewidth=0.75)

    axes[1].barh(ind, df['mn'],
                 left=df['scored_n'],
                 align='center',
                 color='lightgrey',
                 zorder=10,
                 height=w,
                 edgecolor='black',
                 alpha=0.9,
                 linewidth=0.75)

    axes[0].invert_xaxis()
    axes[0].set(yticks=ind,
                yticklabels=["{} ({:.2f})".format(t, p) for
                             t, p in df[['tool', metric]].values.tolist()])

    axes[0].set_xlabel('# Pathogenic Variants')
    axes[1].set_xlabel('# Benign Variants')
    for ax in axes.flat:
        ax.grid(False)
        ax.margins(0.00)
        ax.tick_params(axis='both', which='major', labelsize=8)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)

    colors = ["darkblue", "darkred", "lightgrey"]
    patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
    axes[1].legend(patches, ["Correct", "Incorrect", "Unpredictable"],
                   bbox_to_anchor=(1.7, 1),
                   borderaxespad=0,
                   loc='upper right',
                   prop=dict(size=8))

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


###############################
## Heatmap related functions ##
###############################
GT_LABEL = '(*)True labels'


def prepare_dataset_for_heatmap(df):
    """
    Convert boolean dataframe into numbers
    It discards info about each variant

    :param pd.DataFrame df: Boolean dataframe to process
    :return pd.DataFrame: Dataframe ready for plotting
    """
    df['blank'] = pd.Series(np.nan, index=np.arange(df.shape[0]))
    df = df[['label'] + [col for col in df.columns if '_prediction' in col]].copy()
    df.columns = [GT_LABEL] + [col.replace("_prediction", "") for col in df.columns if '_prediction' in col]
    f = lambda x: pd.isna(x) and np.nan or (1 - int(x))
    for col in df.columns:
        df[col] = df[col].apply(f)
        if df[col].isnull().all():
            del df[col]
    df = df.fillna(-1)
    df = df.reset_index(drop=True)
    df = df.sort_values(by=GT_LABEL)
    return df


def plot_heatmap(df, location, output_dir,
                 cluster_rows: bool = False,
                 skip_preparation: bool = False,
                 prefix: str = None):
    """
    Plot a binary heatmap of tools performance

    :param pd.DataFrame df: Boolean dataframe with predictions after
        running `apply_thresholds` method
    :param str location: Location of the filtered dataframe to write
        the output file
    :param str output_dir: Output directory
    :param bool cluster_rows: Cluster by rows. Default: `False`
    :param bool skip_preparation: Whether input dataset is ready. Default: `False`
    :param str prefix: Extra string to add as prefix
    """

    _name = "heatmap_" if prefix is None else "{}_heatpmap_".format(prefix)
    outfile = os.path.join(output_dir, "heatmap_" + location + '.pdf')
    if len(df.columns) < 2:
        return

    if skip_preparation is False:
        df = prepare_dataset_for_heatmap(df)

    # sns.heatmap(df, cmap=['ivory', 'lightsteelblue', 'darksalmon'],
    #             linecolor='black',
    #             linewidths=0.0,
    #             cbar=False)

    cm = sns.clustermap(df, row_cluster=cluster_rows, cmap=['ivory', 'lightsteelblue', 'darksalmon'],
                        linecolor='black',
                        linewidths=0.0,
                        cbar=False,
                        annot=False,
                        xticklabels=1,
                        yticklabels=False)

    # Aesthetics of the clustergrid
    cm.ax_row_dendrogram.set_visible(False)
    cm.ax_col_dendrogram.set_visible(True)
    cm.cax.set_visible(False)
    ax = cm.ax_heatmap
    ax.set_xticklabels(ax.get_xticklabels(), rotation=80)

    # Add legend
    legend_patch = [mpatches.Patch(facecolor=c, edgecolor='black', linewidth=0.5) for c in
                    ['darksalmon', 'lightsteelblue', 'ivory']]

    legend = ax.legend(legend_patch,
                       ['Pathogenic', 'Benign', 'Unpredictable'],
                       loc='upper right',
                       bbox_to_anchor=(1.25, 1),
                       frameon=True,
                       handlelength=1,
                       handleheight=1)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.0)

    # Add vertical lines
    ax.vlines(range(0, df.shape[1]), ymin=0, ymax=df.shape[0], colors='black', linewidths=0.5)
    ax.axhline(y=0, color='k', linewidth=5)
    ax.axhline(y=df.shape[0], color='k', linewidth=5)
    ax.axvline(x=0, color='k', linewidth=5)
    ax.axvline(x=df.shape[1], color='k', linewidth=5)

    plt.savefig(outfile, bbox_inches='tight', pad_inches=0)
    plt.close()
