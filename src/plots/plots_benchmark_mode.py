import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import logging
plt.switch_backend('agg')
sns.set(style="white")
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from matplotlib.lines import Line2D
import matplotlib.colors as c
import matplotlib.patches as mpatches
from statannotations.Annotator import Annotator
from plots.plots_utils import *


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
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    
    fraction_nan = round(np.sum(np.sum(data.iloc[:, 0].isnull())) / data.shape[0] * 100, 3)
    df_ = data[~data.iloc[:, 0].isnull()]

    if df_.shape[0] > min_predicted:

        sns.kdeplot(x=df_.columns[0],
                    data=df_,
                    palette={True: 'firebrick', False: 'royalblue'},
                    shade=True,
                    linewidth=1.5,
                    hue='label',
                    legend=False,
                    warn_singular=False)

        plt.axvline(x=[x[2] for x in thresholds if x[0] == list(df_)[0]][0], color='k', linestyle="--")
        legend_d = {'Pathogenic': 'firebrick',
                    'Benign': 'royalblue',
                    'Reference threshold': 'black'}

        handles = []
        for _l, _c in legend_d.items():
            style = '--' if _l not in ['Benign', 'Pathogenic'] else '-'

            handles.append(Line2D([], [], color=_c, linestyle=style,
                                  markersize=10,
                                  markeredgewidth=1.5,
                                  label=_l))
        plt.legend(handles=handles)

        plt.xlabel(list(df_)[0] + " ({}% missing data)".format(fraction_nan))
        ends = (max(df_.iloc[:, 0]) - min(df_.iloc[:, 0])) * 0.05
        plt.xlim(min(df_.iloc[:, 0]) - ends, max(df_.iloc[:, 0]) + ends)
        plt.savefig(fname + ".pdf")
        plt.close()


def plot_allele_frequency(df: pd.DataFrame,
                          fname: str,
                          af_col: str = "gnomADg_AF"):
    """
    Plots allele frequencies for each class

    :param pd.DataFrame df: Input df
    :param str fname: Output basemame
    :param str af_col: Column name that accounts
        for allele frequencies. Default: `gnomAD_genomes`.
        If column does not exist, analysis will be skipped.
    """

    if af_col not in df.columns:
        return
    else:
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        
    df['grouper'] = df['outcome'].astype(str) + '\nN = ' + df['count_class'].astype(str)
    order = sorted(list(df['grouper'].unique()))

    ax = sns.boxplot(data=df, x="grouper", order=order, y=af_col)
    try:
        if len(order) == 2:
            annotator = Annotator(ax, 
                                pairs=[(order[0], order[1])], 
                                data=df,
                                x='grouper', 
                                y=af_col, 
                                order=order)
            annotator.configure(test='Mann-Whitney',
                                text_format='star',
                                loc='inside',
                                verbose=0,
                                pvalue_format_string='{:.4f}')
            annotator.apply_and_annotate()

        plt.xlabel("")
        plt.ylabel("Allele frequency")
        plt.tight_layout()
        plt.savefig(fname)
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
    _data = data[['tool', 'mp', 'mn', 'fraction_nan', 'total']].copy()

    _data.loc[:, 'Pathogenic'] = round(data.mp / data.total, 3)
    _data.loc[:, 'Benign'] = round(data.mn / data.total, 3)
    _data = _data[['tool', 'Pathogenic', 'Benign', 'fraction_nan']].set_index('tool')

    os.makedirs(os.path.dirname(fname), exist_ok=True)

    plt.figure()
    plt.rcParams.update({'font.size': 10}) 
    ax = _data.sort_values('fraction_nan').drop('fraction_nan', axis=1).plot.barh(stacked=True, 
                                                                                  width=1,
                                                                                  color={'Pathogenic': 'darkred', 'Benign': 'skyblue'},
                                                                                  edgecolor='k',
                                                                                  alpha=0.7,
                                                                                  fontsize=8,
                                                                                  figsize=(4.25, 4.5))
    ax.legend(fontsize=6)
    ax.set(xlabel='Fraction of unscored variants', ylabel='')
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(fname)
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
    if data.empty:
        return 
    os.makedirs(os.path.dirname(fname), exist_ok=True)

    try:
        data = data[data.fraction_nan < 0.95]
    except AttributeError:
        pass

    data = data.sort_values(metric)

    my_range_coverage = list(range(1, len(data.index) + 1))
    my_range_specificity = np.arange(1 - 0.2, len(data.index)).tolist()
    my_range_sensitivity = np.arange(1 + 0.2, len(data.index) + 0.5).tolist()
    
    n_tools = data.shape[0]
    if n_tools <= 10:
        figsize = (6 , 3.4)
        _bbox = 1.6
        _left = 0.25
        _right = 0.75
        m_title_l = 0.03
        
    elif n_tools < 20:
        figsize = (6.3 , 4.4)
        _bbox = 1.4
        _left = 0.2 
        _right = 0.8
        m_title_l = 0.03
    elif 20 <= n_tools <= 30:
        figsize = (7.5 , 6.5)
        _bbox = 1.3
        _left = 0.3
        _right = 0.8
        m_title_l = 0.04
    else:
        figsize = (8.2, 7)
        _bbox = 1.35
        _left = 0.3
        _right = 0.8
        m_title_l = 0.05
    fig, ax = plt.subplots(figsize=figsize)
    
    if n_tools > 10:
        plt.text(m_title_l, 0.9, "Tool({})".format(metric), fontsize=10, transform=plt.gcf().transFigure)

    _target_col = "precision" if metric in ['F1', 'weighted_F1'] else "specificity"

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
    plt.subplots_adjust(left=_left, right=_right)
    plt.legend(loc="upper right",
               bbox_to_anchor=(_bbox, 1),
               borderaxespad=0,
               prop=dict(size=8))

    if n_tools <= 10:
        plt.yticks(my_range_coverage, data['tool'] + " (" + data[metric].astype(str) + ")", fontsize=10)
    else:
        plt.yticks(my_range_coverage, data['tool'] + " (" + data[metric].astype(str) + ")")

    if all(col in data.columns for col in ['total', 'total_p', 'total_n']):
        plt.title("#variants: {} ({} pos, {} neg)".format(data["total"].iloc[0],
                                                          data["total_p"].iloc[0],
                                                          data["total_n"].iloc[0]))

    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
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
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    
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
    plt.savefig(fname)
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
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    df = data.iloc[::-1]
    df = df[df.fraction_nan < 0.95]
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
    plt.savefig(fname)
    plt.close()


def plot_tools_barplot(data: pd.DataFrame, fname: str, metric: str):
    """
    Tools performance barplot with two
    axis (pathogenic and benign).

    :param pd.DataFrame data: Input df
    :param str fname: Output basename
    :param str metric: Metric to evaluate
    """
    if data.empty:
        return
    
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    
    # revert order
    df = data.iloc[::-1]
    df = df[df.fraction_nan < 0.95]
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
    plt.savefig(fname)
    plt.close()
    sns.reset_defaults()


def plot_curve(data: list,
               fname: str,
               class_counts: tuple,
               is_roc: bool = True,
               max_nan_allowed: float = 0.5):
    """
    Plot ROC or pr curves for all tools in data

    :param list data: ROC analysis results for
    all the tools
    :param str fname: Output basename
    :param tuple class_counts: Number of positive and
    negative variants
    :param bool is_roc: Whether analysis refers to
    ROC curve. If `False`, precision-recall curves are
    drawn. Default: `True`
    :param float max_nan_allowed: Maximum
    fraction of missing variants per tool allowed
    for the curve to be drawn. Default: `0.5`
    """
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    sns.set_style("white")

    data = [x for x in data if len(x) > 0]
    
    if data:
        if is_roc:
            colnames = ['tool', 'fraction_nan', 'label', 'thresholds', 'True Positive Rate (TPR)',
                        'False Positive Rate (FPR)', 'roc_auc']
            to_explode = ['thresholds', 'True Positive Rate (TPR)', 'False Positive Rate (FPR)']
        else:
            colnames = ['tool', 'fraction_nan', 'label', 'thresholds', 'Recall', 'Precision', 'ap_score']
            to_explode = ['thresholds', 'Recall', 'Precision']

        df_metrics = pd.DataFrame.from_records(data, columns=colnames)
        df_metrics = df_metrics[df_metrics.thresholds.notna()]
        df_metrics = df_metrics.explode(to_explode).reset_index()

        if is_roc:
            df_metrics['True Positive Rate (TPR)'] = pd.to_numeric(df_metrics['True Positive Rate (TPR)'])
            df_metrics['False Positive Rate (FPR)'] = pd.to_numeric(df_metrics['False Positive Rate (FPR)'])
            df_metrics["tool_with_roc_auc"] = df_metrics["label"] + " auROC=" + \
                                            df_metrics["roc_auc"].round(3).map(str) + ")"
            hue = "tool_with_roc_auc"
            x = "False Positive Rate (FPR)"
            y = "True Positive Rate (TPR)"
            df_metrics = df_metrics.sort_values('roc_auc', ascending=False)
            
            no_min_scored = df_metrics[df_metrics['fraction_nan'] > max_nan_allowed].tool.unique().tolist()
            if no_min_scored:
                logging.info('ROC curves will not be drawn for the following tool(s) (more than {} fraction of missing predictions): {}'.format(max_nan_allowed, ','.join(no_min_scored)))
        else:
            df_metrics['Recall'] = pd.to_numeric(df_metrics['Recall'])
            df_metrics['Precision'] = pd.to_numeric(df_metrics['Precision'])
            df_metrics["tool_with_ap_score"] = df_metrics["label"] + " AP=" + \
                                            df_metrics["ap_score"].round(3).map(str) + ")"
            hue = "tool_with_ap_score"
            x = "Recall"
            y = "Precision"
            df_metrics = df_metrics.sort_values('ap_score', ascending=False)

        df_metrics = df_metrics[df_metrics['fraction_nan'] <= max_nan_allowed]

        # Since S-CAP has several different reference
        # threshold, S-CAP is removed from these analyses
        df_metrics = df_metrics[~df_metrics.tool.str.contains("S-CAP")]

        # If many tools to plot, change color pallette
        if df_metrics.tool.unique().size > 12:
            sns.set(rc={'figure.figsize':(6, 4)})
            sns.set_palette(sns.mpl_palette("magma", df_metrics.tool.unique().size))
            fontsize="small"
        else:
            sns.set(rc={'figure.figsize':(6, 4)})
            sns.set_palette(sns.color_palette("Paired"))
            fontsize="medium"
            
        ax = sns.lineplot(x=x, y=y,
                        data=df_metrics,
                        hue=hue)
        ax.set_aspect(1)

        plt.title("N pos = {}; N neg = {}".format(class_counts[0], class_counts[1]))
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False, fontsize=fontsize)
        plt.ylim(0, 1.05)
        plt.tight_layout()
        plt.savefig(fname, bbox_inches='tight')
        plt.close()
        sns.reset_defaults()
        
        _metric = 'roc_auc' if is_roc else 'ap_score'   
        df_metrics = df_metrics[['tool', _metric]].drop_duplicates()
        return dict(zip(df_metrics.tool, df_metrics[_metric])) 

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
    f = lambda x: pd.isna(x) and np.nan or int(x)
    for col in df.columns:
        df[col] = df[col].apply(f)
        if df[col].isnull().all():
            del df[col]
    df = df.fillna(-1)
    df = df.reset_index(drop=True)
    df = df.sort_values(by=GT_LABEL)
    return df


def plot_heatmap(df,
                 fname: str,
                 cluster_rows: bool = False,
                 skip_preparation: bool = False):
    """
    Plot a binary heatmap of tools performance

    :param pd.DataFrame df: Boolean dataframe with predictions after
        running `apply_thresholds` method
    :param str output_dir: Output directory
    :param bool cluster_rows: Cluster by rows. Default: `False`
    :param bool skip_preparation: Whether input dataset is ready. Default: `False`
    """
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    
    if len(df.columns) < 2:
        return

    if skip_preparation is False:
        df = prepare_dataset_for_heatmap(df)

    colors = {"ivory": -1, "lightsteelblue": 0, "darksalmon": 1}
    l_colors = sorted(colors, key=colors.get)
    cmap = c.ListedColormap(l_colors)

    cm = sns.clustermap(df, row_cluster=cluster_rows, cmap=cmap,
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

    plt.savefig(fname, bbox_inches='tight', pad_inches=0)
    plt.close()
