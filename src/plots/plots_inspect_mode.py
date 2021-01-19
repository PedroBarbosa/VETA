import os
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D


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
                  np.nan: "grey",
                  "is_benign": "darkblue",
                  "is_pathogenic": "darkred"}

    ax = df.plot.area(color=[map_colors.get(x) for x in df.columns], alpha=0.65)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_ylabel('Fraction of tools')
    ax.set_xlabel('Variants (N={})'.format(df.shape[0]))
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'overall_predicted_ratios.pdf'))
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
    fig, ax = plt.subplots(1, 2, figsize=(5, 7))
    if display_annot:
        if df.shape[0] > 0:
            df = df.sort_values(['is_pathogenic']).head(50)
            outfile = os.path.join(outdir, 'predictability_trade_off_top_candidates.pdf')

        else:
            plt.close()
            return

    else:
        outfile = os.path.join(outdir, 'predictability_trade_off.pdf')

    sns.heatmap(df[["is_pathogenic"]], ax=ax[0], cbar_kws={'label': 'Fraction of tools'},
                vmax=1,
                vmin=0,
                cmap="OrRd",
                yticklabels=display_annot)

    sns.heatmap(df[["unpredictable"]], ax=ax[1], cbar_kws={'label': 'Percentage (%)'},
                vmax=100,
                vmin=0,
                cmap="OrRd",
                yticklabels=False)

    ax[0].set_ylabel('All variants ({})'.format(df.shape[0]))
    ax[1].set_ylabel('')
    plt.tight_layout()

    plt.savefig(outfile)
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
    from src.plots.plots_benchmark_mode import plot_heatmap

    if 'label' in df.columns:
        df.rename(columns={"label": "Ground Truth (*)"}, inplace=True)

    for filter_name, _func in filters:

        df_f = _func(df)._get_numeric_data().copy()

        if df_f.shape[0] < 5:
            continue

        plot_heatmap(df_f, filter_name, outdir,
                     cluster_rows=True,
                     skip_preparation=True,
                     prefix="top_tools")


def plot_accuracy(stats_df: pd.DataFrame,
                  metric: str,
                  location: str,
                  out_dir: str):
    """
    Plots weighted_accuracy bars

    :param pd.DataFrame stats_df: Stats dataframe
    :param str metric: Metric to rank
    :param str location: Variants location
    :param str out_dir: Output directory
    """

    fig, ax = plt.subplots()
    stats_df = stats_df.sort_values([metric])
    plt.barh(range(stats_df.shape[0]), stats_df[metric],
             color='darkgrey',
             edgecolor='k',
             linewidth=1)
    plt.xlabel(metric)

    plt.xlim(left=0)
    plt.yticks(range(stats_df.shape[0]), stats_df['tool'])
    plt.ylim(-1, stats_df.shape[0])
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, "performance_{}.pdf".format(location)))
    plt.close()


def plot_tool_score_distribution(_df: pd.DataFrame,
                                 tool: str,
                                 thresholds: List,
                                 outdir: str):
    """
    Plot score distribution for a given tool
    :param pd.DataFrame _df: Input df
    :param str tool: Specific tool to analyse
    :param List thresholds: List of tools with the
        reference thresholds
    :param outdir: Output directory
    """
    plt.subplots(figsize=(7, 5))
    nas = sum(pd.isnull(_df[tool]))

    df = _df.copy()
    df[tool].dropna(inplace=True)
    p = sns.histplot(df[tool], bins=30, color='slateblue', edgecolor='black')
    p.tick_params(labelsize=13)
    plt.xlabel('{} ({} % missing data)'.format(tool, round(nas / _df.shape[0] * 100, 2)), fontsize=13)
    p.yaxis.get_label().set_fontsize(13)
    plt.axvline([x[2] for x in thresholds if x[0] == tool],
                color='r',
                linestyle="--")

    legend_element = [Line2D([0], [0], color='r', lw=4, label='Reference threshold')]
    plt.legend(handles=legend_element, loc='best')
    plt.savefig(os.path.join(outdir, "score_dist_{}.pdf".format(tool)))
    plt.close()
