import numpy as np
import seaborn as sns
import pandas as pd
from preprocessing.utils import ratio
import os
from collections import defaultdict
from osutils import ensure_folder_exists
from matplotlib import lines
from matplotlib.legend import Legend
from filters import filter_intronic_bins
from statannot import add_stat_annotation
from scipy import stats
from plots.plots_utils import *
import matplotlib.pyplot as plt
plt.switch_backend('agg')
cmap = sns.diverging_palette(220, 10, as_cmap=True)


def plot_general_bin_info(df, bins, fname):

    if df[df['class'].isin([True])].empty or df[df['class'].isin([False])].empty:
        dic = {}
        for bin in bins:
            if bin[0] not in {"all_intronic", "all_except_0-10"}:
                dic[bin[0]] = (df.intron_bin == bin[0]).sum()

        plt.bar(dic.keys(), dic.values(), color="silver")
    
    else:
        ax = sns.countplot(x="intron_bin", hue="outcome",
                           order=[x[0] for x in bins if x[0] not in {"all_intronic", "all_except_0-10"}],
                           data=df,
                           palette={"Benign": "skyblue", "Pathogenic": "chocolate"})
        ax.get_legend().set_title('')

    plt.xlabel("Intron base-pair bins")
    plt.ylabel("Variant counts")
    out = fname + '_bin_counts.pdf'
    plt.savefig(out)
    plt.close()

    df['gnomAD_genomes'] = pd.to_numeric(df['gnomAD_genomes'], downcast='float')
    # df['gnomAD_genomes'] = df['gnomAD_genomes'].astype('float64')
    # df['gnomAD_genomes'].apply(lambda x: '%.10f' % float(x))
    if df[df['class'].isin([True])].empty or df[df['class'].isin([False])].empty:
        ax = sns.scatterplot(x="intron_offset", y="gnomAD_genomes", data=df)

    else:
        ax = sns.scatterplot(x="intron_offset", y="gnomAD_genomes", data=df, hue="outcome",
                         palette={"Benign": "skyblue", "Pathogenic": "chocolate"},
                         s=20)
        ax.get_legend().set_title('')

    ax.set(xlabel='Variant intronic offset', ylabel="gnomAD allele frequency")
    plt.tight_layout()
    out = fname + '_offset_distribution.pdf'
    plt.savefig(out)
    plt.close()


def plot_ROCs(df_metrics, fname, n_positive_class, min_score_fraction=0.3):
    df_metrics = df_metrics[df_metrics['fraction_nan'] <= min_score_fraction]
    df_metrics["tool_with_roc_auc"] = df_metrics["tool"] + " auROC=" + df_metrics["ROC-auc"].round(2).map(str) + ")"
    df_metrics["tool_with_pr_auc"] = df_metrics["tool"] + " prAUC=" + df_metrics["PR-auc"].round(2).map(str) + ")"
    df_metrics["tool_with_f1"] = df_metrics["tool"] + " F1=" + df_metrics["F1"].round(2).map(str) + ")"

    #plt.figure(figsize=(15, 11))
    ax = sns.lineplot(x="recall", y="precision", data=df_metrics, hue="tool_with_pr_auc")
    ax.get_legend().set_title('')
    plt.title("N positive = {}".format(n_positive_class))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.rcParams.update({'font.size': 30})
    #ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.ylim(0, 1.05)
    out = fname + '_pr.pdf'
    plt.savefig(out)
    plt.close()

    #grouped = df_metrics.groupby('tool_with_pr_auc')
    #for name, group in grouped:
    #    tool = name.split("(")[0]
    #    ax = sns.lineplot(x="recall", y="precision", data=group, hue="tool_with_pr_auc")
    #    l=[0, 40, 80, 99]
    #    i=0
    #    for x, y, s in zip(group['recall'], group['precision'], group['threshold']):
    #        if i in l:
    #            plt.text(x - 0.004, y + 0.002, round(s, 2), fontdict={'size': 6})
    #        i+=1
    #    out = fname + "_{}".format(tool) + '_pr.pdf'
    #    plt.savefig(out)
    #    plt.close()


    ax = sns.lineplot(x="FPR", y="recall", data=df_metrics, hue="tool_with_roc_auc")
    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")
    ax.get_legend().set_title('')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    out = fname + '.pdf'
    plt.savefig(out)
    plt.close()


def plot_auROC_by_bin(df, fname):

    #sns.set_style("darkgrid")
    sns.catplot(x="bin", y="auROC", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                data=df, linestyles="--", scale=0.7, aspect=0.9,
                hue="tool")
    plt.xlabel("Intron bin (bp)")
    plt.ylabel("auROC")
    plt.tight_layout()
    out = fname + '_auROC.pdf'
    plt.savefig(out)
    plt.close()

    sns.catplot(x="bin", y="prROC", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                data=df,
                linestyles="--", scale=0.7, aspect=0.9,
                hue="tool")
    plt.xlabel("Intron bin (bp)")
    plt.ylabel("prAUC")
    plt.tight_layout()
    out = fname + '_prROC.pdf'
    plt.savefig(out)
    plt.close()


    sns.catplot(x="bin", y="F1", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                data=df, linestyles="--", scale=0.7, aspect=0.9,
                hue="tool",
                )

    plt.xlabel("Intron bin (bp)")
    plt.ylabel("F1 score")
    plt.tight_layout()
    out = fname + '_F1.pdf'
    plt.savefig(out)
    plt.close()

    sns.catplot(x="bin", y="weighted_F1", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                data=df, linestyles="--", scale=0.7, aspect=0.9,
                hue="tool",
                )

    plt.xlabel("Intron bin (bp)")
    plt.ylabel("F1 score (weighted)")
    plt.tight_layout()
    out = fname + '_weighted_F1.pdf'
    plt.savefig(out)
    plt.close()

    sns.catplot(x="bin", y="fraction_nan", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] not in {"all_intronic", "all_except_0-10"}],
                data=df, linestyles="--", scale=0.7, aspect=0.9, dodge=True,
                hue="tool")

    plt.xlabel("Intron bin (bp)")
    plt.ylabel("Fraction of NaN")
    plt.ylim(-0.01, 1)
    plt.tight_layout()
    out = fname + '_fraction_unscored.pdf'
    plt.savefig(out)
    plt.close()


def plot_allele_frequency(df, fname):
    try:
        df['grouper'] = df['outcome'].astype(str) + '\nN = ' + df['count_class'].astype(str)
        order = sorted(list(df['grouper'].unique()))
        ax = sns.boxplot(data=df, x="grouper", order=order,  y="gnomAD_genomes")
        add_stat_annotation(ax, data=df, x="grouper", y="gnomAD_genomes",
                            order=order,
                            box_pairs=[tuple(order)],
                            test='Mann-Whitney',
                            text_format='star',
                            loc='inside',
                            verbose=0,
                            pvalue_format_string='{:.4f}')
        plt.xlabel("")
        plt.ylabel("gnomAD frequency")
        plt.tight_layout()
        out = fname + '.pdf'
        plt.savefig(out)
        plt.close()

    except ValueError:
        pass


