import matplotlib.pyplot as plt
import os
import seaborn as sns
import matplotlib.patches as mpatches
import pandas as pd
from typing import List


def plot_area(df: pd.DataFrame, outdir: str):

    """
    Area plot that gives an overview of the
    tools predictions for a given set of
    input variants
    :param pd.DataFrame df: Input df with
    ratios of tools that predict a given
    outcome.
    :param str outdir: Output directory to
        save the fig
    """

    fig, ax = plt.subplots(figsize=(12, 12))
    map_colors = {"unpredictable": "grey",
                  "is_benign": "darkblue",
                  "is_pathogenic": "darkred"}

    ax = df.plot.area(color=[map_colors.get(x) for x in df.columns], alpha=0.65)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_ylabel('Fraction of tools')
    ax.set_xlabel('Variants')
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'areaplot.pdf'))
    plt.close()


def plot_heatmap(df: pd.DataFrame, outdir: str, display_annot: bool = False):
    """

    Draw heatmaps that color the ratio of tools
    that predict pathogenicity or unpredictability.
    Each row is a variant, the more tools predict
    a given outcome, the redder the color

    :param pd.DataFrame df: Input df with
    ratios of tools that predict a given
    outcome.
    :param str outdir: Output directory to
        save the fig
    :param bool display_annot: Display row
        annotations in the heatmap (variant ID).
        Default: `False`. Useful when few
        variants are being plotted.
    :return:
    """
    fig, ax = plt.subplots(1, 2, figsize=(5,7))

    sns.heatmap(df[["is_pathogenic"]], ax=ax[0], vmax=1, vmin=0, cmap="OrRd", yticklabels=display_annot)
    sns.heatmap(df[["unpredictable"]], ax=ax[1], cbar_kws={'label': '%'}, vmax=100, vmin=0,  cmap="OrRd",
                yticklabels= False)
    ax[0].set_ylabel('')
    ax[1].set_ylabel('')
    plt.tight_layout()

    if display_annot:
        plt.savefig(os.path.join(outdir, 'heatmap_with_annot.pdf'))
    else:
        plt.savefig(os.path.join(outdir, 'heatmap_all.pdf'))
    plt.close()


def plot_heatmap_toptools(df: pd.DataFrame, filters, outdir):
    """
    Draw heatmap for the top tools

   :param pd.DataFrame df: Df with predictions
        for each variant
    :param List filters: Location filters to employ
        so that variants at each location are processed
        independently
    :param str outdir: Output directory
    """

    for filter_name, _func in filters:
        df_f = _func(df).drop(["location"], axis=1).copy()

        if df_f.shape[0] > 0:
            plt.clf()
            f, ax = plt.subplots(figsize=(0.5 * len(df_f.columns), 6))
            colors = ["slateblue", "w", "indianred"]
            sns.heatmap(df_f, cmap=colors, ax=ax, cbar=False)
            if df_f.shape[0] > 40:
                plt.yticks([])
                ax.set_ylabel("{} variants".format(df_f.shape[0]))
            else:
                ax.set_ylabel('')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
            legend_ax = f.add_axes([.7, .5, 1, .1])
            legend_ax.axis('off')

            patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
            legend_ax.legend(patches, ["Benign", "Unpredictable", "Pathogenic"], handlelength=0.8, loc='lower left')
            plt.savefig(os.path.join(outdir, "top_tools_{}.pdf".format(filter_name)),
                        bbox_inches='tight', pad_inches=0)
            plt.close()


def plot_tool_score_distribution(df: pd.DataFrame,
                                 tool: str,
                                 thresholds: List,
                                 outdir: str):
    """
    Plot score distribution for a given tool
    :param pd.DataFrame df:
    :param str tool:
    :param List thresholds: List of tools with the
        reference thresholds
    :param outdir:
    :return:
    """

    df[tool].dropna(inplace=True)
    sns.distplot(df[tool], hist=True, kde=False,
                 bins=30, color='slateblue',
                 hist_kws={'edgecolor': 'black'})

    plt.axvline([x[2] for x in thresholds if x[0] == tool], color='r', linestyle="--")
    plt.savefig(os.path.join(outdir, "score_dist_{}.pdf".format(tool)))
    plt.close()
