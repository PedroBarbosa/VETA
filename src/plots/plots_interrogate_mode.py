import os
from typing import List
import logging
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

    map_colors = {"Unpredictable": "grey",
                  np.nan: "grey",
                  "Is benign": "darkblue",
                  "Is pathogenic": "darkred"}
    
    ax = df.sort_values(['Is pathogenic', 'Is benign'], ascending=False).plot.area(color=[map_colors.get(x)
                      for x in df.columns], alpha=0.65)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_ylabel('Fraction of tools')
    ax.set_xlabel('Variants (N={})'.format(df.shape[0]))
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'overall_predicted_ratios.pdf'))
    plt.close()


def plot_heatmap(df: pd.DataFrame,
                 outdir: str,
                 display_annot: bool = False,
                 benign_too: pd.DataFrame = None):

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
    :param pd.DataFrame benign_too: Display 
    also variants for which most tools predict 
    to be benign. Default: None, display
    unpredictable variants
    :return:
    """

    plot_benign = False
    fig, ax = plt.subplots(1, 2, figsize=(5, 7))

    if display_annot:

        if isinstance(benign_too, pd.DataFrame):
            if benign_too.shape[0] > 0 and df.shape[0] > 0:
                all_path = df[df['Is pathogenic'] == 1].copy()
                if all_path.shape[0] > 50:
                    df = all_path
                else:
                    df = df.sort_values(['Is pathogenic']).head(50)
                benign_too = benign_too.sort_values(['Is benign']).head(50)
                outfile = os.path.join(
                    outdir, 'predictability_trade_off_top_patho_and_benign.pdf')
                plot_benign = True
            else:
                plt.close()
                return

        elif df.shape[0] > 0:
            df = df.sort_values(['Is pathogenic']).head(50)
            outfile = os.path.join(
                outdir, 'predictability_trade_off_top_candidates.pdf')
        else:
            plt.close()
            return
    else:
        outfile = os.path.join(outdir, 'predictability_trade_off.pdf')

    if plot_benign:

        sns.heatmap(df[["Is pathogenic"]][::-1], ax=ax[0],
                    cbar=False,
                    vmax=1,
                    vmin=0,
                    cmap="OrRd",
                    linewidths=.5,
                    linecolor='k',
                    yticklabels=display_annot)

        sns.heatmap(benign_too[["Is benign"]][::-1], ax=ax[1],
                    cbar_kws={'label': 'Fraction of tools'},
                    vmax=1,
                    vmin=0,
                    cmap="OrRd",
                    linewidths=.5,
                    linecolor='k',
                    annot_kws={'size': 1},
                    yticklabels=display_annot)

        ax[0].set_yticklabels(ax[0].get_ymajorticklabels(), fontsize=10)
        ax[1].set_yticklabels(ax[1].get_ymajorticklabels(), fontsize=10)

        ax[0].set_ylabel('')
        ax[1].set_ylabel('')

    else:
        _linewidth = 0.001 if df.shape[0] > 2000 else .5
        sns.heatmap(df[["Is pathogenic"]][::-1], ax=ax[0], cbar_kws={'label': 'Fraction of tools'},
                    vmax=1,
                    vmin=0,
                    cmap="OrRd",
                    linewidths=_linewidth,
                    linecolor='k',
                    yticklabels=display_annot)

        sns.heatmap(df[["Unpredictable"]], ax=ax[1], cbar_kws={'label': 'Percentage (%)'},
                    vmax=100,
                    vmin=0,
                    cmap="OrRd",
                    linewidths=_linewidth,
                    linecolor='k',
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
    from plots.plots_benchmark_mode import plot_heatmap

    if 'label' in df.columns:
        df.rename(columns={"label": "Ground Truth (*)"}, inplace=True)
    
    os.makedirs(os.path.join(outdir, 'out_heatmaps'), exist_ok=True)
    outdir = os.path.join(outdir, 'out_heatmaps')

    for filter_name in filters:
        
        df_f = df[df.location == filter_name].copy()
        df_f = df_f.drop(columns=['variant_id', 'HGVSc', 'HGVSg', 'location', "SYMBOL"])
        if df_f.shape[0] < 5:
            continue

        fname = os.path.join(
            outdir, "top_tools_heatmap_{}.pdf".format(filter_name))

        try:
            plot_heatmap(df_f, fname,
            cluster_rows=True, 
            skip_preparation=True)
        except RecursionError:
            logging.info("WARN: Error when generating heatmap for '{}' variants".format(filter_name))
            pass


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
    stats_df = stats_df.sort_values([metric])
    n_tools = len(stats_df)
    if n_tools < 5:
        figsize = (3, 3)
    elif 5 <= n_tools <= 30:
        figsize = (5, 5)
    else:
        figsize = (7, 7)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.barh(range(stats_df.shape[0]), stats_df[metric],
            color='darkgrey',
            edgecolor='k',
            linewidth=1)

    plt.axvline(0.9, color='b', linestyle='--')
    xlabel = (metric[:1].upper() + metric[1:]).replace("_", " ")
    _title = (location[:1].upper() + location[1:]).replace("_", " ")
    plt.xlabel(xlabel)
    plt.title("{} variants (N={})".format(
        _title, str(stats_df.iloc[0, :].total)))

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
    os.makedirs(os.path.join(outdir, 'out_score_distribution'), exist_ok=True)
    outdir = os.path.join(outdir, 'out_score_distribution')
    plt.subplots(figsize=(7, 5))
    nas = sum(pd.isnull(_df[tool]))

    df = _df.copy()
    df[tool].dropna(inplace=True)
    p = sns.histplot(df[tool], bins=30, color='slateblue', edgecolor='black')
    p.tick_params(labelsize=13)
    plt.xlabel('{} ({} % missing data)'.format(
        tool, round(nas / _df.shape[0] * 100, 3)), fontsize=13)

    p.yaxis.get_label().set_fontsize(13)
    plt.axvline([x[2] for x in thresholds if x[0] == tool],
                color='r',
                linestyle="--")

    legend_element = [Line2D([0], [0], color='r', lw=4,
                             label='Reference threshold')]

    plt.legend(handles=legend_element, loc='best')
    plt.savefig(os.path.join(outdir, "score_dist_{}.pdf".format(tool)))
    plt.close()
