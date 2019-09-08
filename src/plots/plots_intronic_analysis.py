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

    print(df['class'].value_counts())
    if df[df['class'].isin([True])].empty or df[df['class'].isin([False])].empty:
        dic = {}
        for bin in bins:
            if bin[0] != "all_intronic":
                dic[bin[0]] = (df.intron_bin == bin[0]).sum()

        plt.bar(dic.keys(), dic.values(), color="silver")
    
    else:
        ax = sns.countplot(x="intron_bin", hue="outcome",
                           order=[x[0] for x in bins if x[0] != "all_intronic"],
                           data=df,
                           palette={"Benign": "skyblue", "Pathogenic": "chocolate"})
        ax.get_legend().set_title('')

    plt.xlabel("Intron base-pair bins")
    plt.ylabel("Variant counts")
    out = fname + '_bin_counts.pdf'
    plt.savefig(out)
    plt.close()

    df['gnomAD_genomes'] = df['gnomAD_genomes'].astype('float64')
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


def plot_ROCs(df_metrics, fname, n_positive_class):

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

    sns.set_style("darkgrid")
    sns.catplot(x="bin", y="auROC", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] != "all_intronic"],
                data=df,
                hue="tool")
    plt.xlabel("Intron bin")
    plt.ylabel("auROC")
    plt.tight_layout()
    out = fname + '_auROC.pdf'
    plt.savefig(out)
    plt.close()

    sns.catplot(x="bin", y="prROC", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] != "all_intronic"],
                data=df,
                hue="tool")
    plt.xlabel("Intron bin")
    plt.ylabel("prAUC")
    plt.tight_layout()
    out = fname + '_prROC.pdf'
    plt.savefig(out)
    plt.close()


    sns.catplot(x="bin", y="F1", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] != "all_intronic"],
                data=df,
                hue="tool",
                color="grey",
                )

    plt.xlabel("Intron bin")
    plt.ylabel("F1 score")
    plt.tight_layout()
    out = fname + 'F1.pdf'
    plt.savefig(out)
    plt.close()

    sns.catplot(x="bin", y="fraction_nan", kind='point',
                order=[i[0] for i in filter_intronic_bins if i[0] != "all_intronic"],
                data=df,
                hue="tool")

    plt.xlabel("Intron bin")
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
        print(df["outcome"].value_counts())
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


