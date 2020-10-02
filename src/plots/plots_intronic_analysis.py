import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.switch_backend('agg')
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from src.plots.plots_utils import *
from src.predictions.filters import filter_intronic_bins


def plot_general_bin_info(df: pd.DataFrame, outdir: str):
    """
    Barplot with variant counts on each
    intronic bin
    
    :param pd.DataFrame df: Input df
    :param str outdir: Output directory
    """

    # if there are no two labels
    if df[df['label'].isin([True])].empty or df[df['label'].isin([False])].empty:
        dic = {}

        for _bin in filter_intronic_bins:
            if _bin[0] not in {"all_intronic", "all_except_0-10"}:
                dic[_bin[0]] = np.sum(df.intron_bin == _bin[0])

        plt.bar(dic.keys(), dic.values(), color="silver")

    else:
        ax = sns.countplot(x="intron_bin",
                           hue="outcome",
                           order=[x[0] for x in filter_intronic_bins if
                                  x[0] not in {"all_intronic", "all_except_0-10"}],
                           data=df,
                           palette={"Benign": "skyblue",
                                    "Pathogenic": "chocolate"
                                    }
                           )
        ax.get_legend().set_title('')

    plt.xlabel("Intron base-pair bins")
    plt.ylabel("Variant counts")
    out = os.path.join(outdir, 'intronic_bin_counts.pdf')
    plt.savefig(out)
    plt.close()

    df_zoom = df[(~df['intron_bin'].str.match('0-10')) & (df['outcome'] == "Pathogenic")]
    ylim = df_zoom['intron_bin'].value_counts().max() + (df_zoom['intron_bin'].value_counts().max() * 0.05)
    out_zoomed = os.path.join(outdir, 'intronic_bin_counts_zoomed.pdf')
    plt.ylim(0, ylim)
    plt.savefig(out_zoomed)
    plt.close()

    df['gnomAD_genomes'] = pd.to_numeric(df['gnomAD_genomes'], downcast='float')
    # df['gnomAD_genomes'] = df['gnomAD_genomes'].astype('float64')
    # df['gnomAD_genomes'].apply(lambda x: '%.10f' % float(x))

    if df[df['label'].isin([True])].empty or df[df['label'].isin([False])].empty:
        ax = sns.scatterplot(x="intron_offset",
                             y="gnomAD_genomes",
                             data=df)

    else:
        ax = sns.scatterplot(x="intron_offset",
                             y="gnomAD_genomes",
                             data=df,
                             hue="outcome",
                             palette={"Benign": "skyblue",
                                      "Pathogenic": "chocolate"
                                      },
                             s=20
                             )
        ax.get_legend().set_title('')

    ax.set(xlabel='Variant intronic offset', ylabel="gnomAD allele frequency")
    plt.tight_layout()
    out = os.path.join(outdir, 'intronic_offset_distribution.pdf')
    plt.savefig(out)
    plt.close()


def plot_ROCs(df: pd.DataFrame, fname: str,
              n_positive_class: int,
              min_score_fraction: float = 0.3):
    """
    Plot ROC curves for all tools at a given bin
    :param pd.DataFrame df: Variants at a given
        intronic bin
    :param str fname: Output basename
    :param int n_positive_class: Number of
        instances of positive class
    :param float min_score_fraction: Minimum
    fraction of predictive power of a given
    tool for the curve to be drawn. Default: `0.3`
    """

    df_metrics = df[df['fraction_nan'] <= min_score_fraction].copy()
    df_metrics["tool_with_roc_auc"] = df_metrics["tool"] + " auROC=" + df_metrics["ROC-auc"].round(2).map(str) + ")"
    df_metrics["tool_with_pr_auc"] = df_metrics["tool"] + " prAUC=" + df_metrics["PR-auc"].round(2).map(str) + ")"
    df_metrics["tool_with_f1"] = df_metrics["tool"] + " F1=" + df_metrics["F1"].round(2).map(str) + ")"

    # Since S-CAP has several different reference
    # threshold, S-CAP is removed from these analyses
    df_metrics = df_metrics[~df_metrics.tool.str.contains("S-CAP")]

    # If many tools to plot, change color pallette
    if df_metrics.tool.unique().size > 16:
        sns.set_palette(sns.mpl_palette("GnBu_d", df_metrics.tool.unique().size))

    #############
    ### prAUC ###
    #############
    ax = sns.lineplot(x="recall",
                      y="precision",
                      data=df_metrics,
                      hue="tool_with_pr_auc")

    ax.get_legend().set_title('')
    plt.title("N positive = {}".format(n_positive_class))
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc=2,
               borderaxespad=0.)

    # plt.rcParams.update({'font.size': 30})
    # ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.ylim(0, 1.05)
    out = fname + '_pr.pdf'
    plt.savefig(out)
    plt.close()

    # grouped = df_metrics.groupby('tool_with_pr_auc')
    # for name, group in grouped:
    #     tool = name.split("(")[0]
    #     ax = sns.lineplot(x="recall", y="precision", data=group, hue="tool_with_pr_auc")
    #     l=[0, 40, 80, 99]
    #     i=0
    #     for x, y, s in zip(group['recall'], group['precision'], group['threshold']):
    #         if i in l:
    #             plt.text(x - 0.004, y + 0.002, round(s, 2), fontdict={'size': 6})
    #         i+=1
    #     out = fname + "_{}".format(tool) + '_pr.pdf'
    #     plt.savefig(out)
    #     plt.close()

    ###########
    ### AUC ###
    ###########
    ax = sns.lineplot(x="FPR",
                      y="recall",
                      data=df_metrics,
                      hue="tool_with_roc_auc")

    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")
    ax.get_legend().set_title('')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
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

    # sns.set_style("darkgrid")

    metrics = {"F1": "F1_Score",
               "weighted_F1": "F1 score (weighted)", 
               "fraction_nan": "_fraction_unscored.pdf",
               "auROC": "auROC",
               "prROC": "prROC"}

    for metric, description in metrics.items():

        # Ploting auROCs and auROCpr requires removal of S-CAP
        if metric == "auROC":
            df = df[~df["tool"].str.contains("S-CAP")]

        sns.catplot(x="bin",
                    y=metric,
                    kind='point',
                    order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                    data=df, linestyles="--", scale=0.7, aspect=0.9,
                    hue="tool")

        plt.xlabel("Intron bin (bp)")
        plt.ylabel(description)
        plt.tight_layout()
        out = fname + '_' + metric + '.pdf'
        plt.savefig(out)
        plt.close()