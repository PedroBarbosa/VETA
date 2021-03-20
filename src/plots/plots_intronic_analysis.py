import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.switch_backend('agg')
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from src.plots.plots_utils import *
from src.predictions.filters import filter_intronic_bins


def plot_general_bin_info(df_i: pd.DataFrame, outdir: str, af_column: str):
    """
    Barplot with variant counts on each
    intronic bin
    
    :param pd.DataFrame df_i: Input df
    :param str outdir: Output directory
    :param str af_column: Allele frequency column
    """

    df = df_i.copy()
    df.groupby(['intron_bin', 'outcome']).size().to_csv(os.path.join(outdir,
                                                                     'intronic_bin_counts.tsv'),
                                                        sep="\t")

    # if there are no two labels
    if df[df['label'].isin([True])].empty or df[df['label'].isin([False])].empty:
        dic = {}

        for _bin in filter_intronic_bins:
            if _bin[0] not in {"all_intronic", "all_except_0-2", "all_except_0-10"}:
                dic[_bin[0]] = np.sum(df.intron_bin == _bin[0])

        plt.bar(dic.keys(), dic.values(), color="silver")

    else:
        ax = sns.countplot(x="intron_bin",
                           hue="outcome",
                           order=[x[0] for x in filter_intronic_bins if
                                  x[0] not in {"all_intronic", "all_except_0-2", "all_except_0-10"}],
                           data=df,
                           linewidth=1,
                           edgecolor='k',
                           palette={"Benign": "skyblue",
                                    "Pathogenic": "chocolate"
                                    }
                           )
        ax.get_legend().set_title('')

    plt.xlabel("Intron base-pair bins")
    plt.ylabel("Variant counts")
    out = os.path.join(outdir, 'intronic_bin_counts.pdf')
    plt.savefig(out)

    df_zoom = df[(~df['intron_bin'].str.match('0-2')) &
                 (~df['intron_bin'].str.match('3-10')) &
                 (df['outcome'] == "Pathogenic")]

    ylim = df_zoom['intron_bin'].value_counts().max() + (df_zoom['intron_bin'].value_counts().max() * 0.05)

    out_zoomed = os.path.join(outdir, 'intronic_bin_counts_zoomed.pdf')
    plt.ylim(0, ylim)
    plt.savefig(out_zoomed)
    plt.close()

    if af_column in df.columns:
        df[af_column] = pd.to_numeric(df[af_column], downcast='float')

        # df[af_column] = df[af_column].astype('float64')
        # df[af_column].apply(lambda x: '%.10f' % float(x))

        if df[df['label'].isin([True])].empty or df[df['label'].isin([False])].empty:

            ax = sns.scatterplot(x="intron_offset",
                                 y=af_column,
                                 data=df)

        else:
            ax = sns.scatterplot(x="intron_offset",
                                 y=af_column,
                                 data=df,
                                 hue="outcome",
                                 s=10,
                                 palette={"Benign": "skyblue",
                                          "Pathogenic": "chocolate"})
            ax.get_legend().set_title('')

        ax.set(xlabel='Variant intronic offset', ylabel="{} allele frequency".format(af_column))
        plt.tight_layout()
        out = os.path.join(outdir, 'intronic_offset_distribution.pdf')
        plt.savefig(out)
        plt.close()


def plot_curve(data: list,
               fname: str,
               class_counts: tuple,
               is_roc: bool = True,
               min_score_fraction: float = 0.5):
    """
    Plot ROC or pr curves for all tools at a given bin

    :param list data: ROC analysis results for
    all the tools at the given intronic bin
    :param str fname: Output basename
    :param tuple class_counts: Number of positive and
    negative variants at the given intronic bin
    :param bool is_roc: Whether analysis refers to
    ROC curve. If `False`, precision-recall curves are
    drawn. Default: `True`
    :param float min_score_fraction: Minimum
    fraction of predictive power of a given
    tool for the curve to be drawn. Default: `0.5`
    """
    if is_roc:
        colnames = ['tool', 'fraction_nan', 'label', 'thresholds', 'True Positive Rate (TPR)',
                    'False Positive Rate (FPR)', 'roc_auc']
        to_explode = ['thresholds', 'True Positive Rate (TPR)', 'False Positive Rate (FPR)']
    else:
        colnames = ['tool', 'fraction_nan', 'label', 'thresholds', 'Recall', 'Precision', 'ap_score']
        to_explode = ['thresholds', 'Recall', 'Precision']

    df_metrics = pd.DataFrame.from_records(data, columns=colnames)
    df_metrics = df_metrics.reset_index().apply(lambda x: x.explode() if x.name in to_explode else x)

    if is_roc:
        df_metrics['True Positive Rate (TPR)'] = pd.to_numeric(df_metrics['True Positive Rate (TPR)'])
        df_metrics['False Positive Rate (FPR)'] = pd.to_numeric(df_metrics['False Positive Rate (FPR)'])
        df_metrics["tool_with_roc_auc"] = df_metrics["label"] + " auROC=" + \
                                          df_metrics["roc_auc"].round(2).map(str) + ")"
        hue = "tool_with_roc_auc"
        x = "False Positive Rate (FPR)"
        y = "True Positive Rate (TPR)"
        df_metrics = df_metrics.sort_values('roc_auc', ascending=False)
    else:
        df_metrics['Recall'] = pd.to_numeric(df_metrics['Recall'])
        df_metrics['Precision'] = pd.to_numeric(df_metrics['Precision'])
        df_metrics["tool_with_ap_score"] = df_metrics["label"] + " AP=" + \
                                           df_metrics["ap_score"].round(2).map(str) + ")"
        hue = "tool_with_ap_score"
        x = "Recall"
        y = "Precision"
        df_metrics = df_metrics.sort_values('ap_score', ascending=False)

    df_metrics = df_metrics[df_metrics['fraction_nan'] <= min_score_fraction]

    # Since S-CAP has several different reference
    # threshold, S-CAP is removed from these analyses
    df_metrics = df_metrics[~df_metrics.tool.str.contains("S-CAP")]

    # If many tools to plot, change color pallette
    if df_metrics.tool.unique().size > 12:
        sns.set_palette(sns.mpl_palette("magma", df_metrics.tool.unique().size))
    else:
        sns.set_palette(sns.color_palette("Paired"))

    ax = sns.lineplot(x=x, y=y,
                      data=df_metrics,
                      hue=hue)
    ax.set_aspect(1.15)
    plt.title("N pos = {}; N neg = {}".format(class_counts[0], class_counts[1]))
    plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
    plt.ylim(0, 1.05)
    plt.tight_layout()
    out = fname + '.pdf'
    plt.savefig(out)
    plt.close()
    sns.reset_defaults()


def plot_metrics_by_bin(df: pd.DataFrame, fname: str):
    """
    :param pd.DataFrame df: Df with a list of
        metrics per each intronic bin
    :param str fname: Output basename
    """

    metrics = {"F1": "F1_Score",
               "weighted_F1": "F1 score (weighted)",
               "fraction_nan": "Fraction unscored"}

    for metric, description in metrics.items():

        # Ploting auROCs and auROCpr requires removal of S-CAP
        # if metric == "auROC":
        #     df = df[~df["tool"].str.contains("S-CAP")]
        df[metric] = pd.to_numeric(df[metric])
        sns.color_palette("dark")
        sns.catplot(x="bin",
                    y=metric,
                    kind='point',
                    order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-2",
                                                                              "all_except_0-10"}],
                    data=df, linestyles="--", linewidth=0.2, scale=0.5, aspect=1.2,
                    legend=False,
                    dodge=True,
                    hue="tool")

        plt.legend(loc="upper right",
                   bbox_to_anchor=(1.2, 1),
                   borderaxespad=0,
                   prop=dict(size=8))
        plt.yticks(np.arange(0, 1.05, 0.1))
        plt.subplots_adjust(left=0.35)

        plt.xlabel("Intron bin (bp)")
        plt.ylabel(description)
        plt.tight_layout()
        out = fname + '_' + metric + '.pdf'
        plt.savefig(out)
        plt.close()
