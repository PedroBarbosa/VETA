import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from preprocessing.osutils import ensure_folder_exists
import seaborn as sns

plt.switch_backend('agg')
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from plots.plots_utils import *
from predictions.filters import filter_intronic_bins


def plot_general_bin_info(_df: pd.DataFrame, 
                          outdir: str, 
                          variant_class: str,
                          aggregate_classes: bool,
                          af_column: str):
    """
    Barplot with variant counts on each
    intronic bin
    
    :param pd.DataFrame _df: Input df
    :param str outdir: Output directory
    :param str variant_class: Variant class
    :oaram bool aggregate_classes: Whether splice site are analyzed together with splice region variants
    :param str af_column: Allele frequency column
    """
    
    def _generate_plot(df, filt_to_exclude):
        if df[df['label'].isin([True])].empty or df[df['label'].isin([False])].empty:
            dic = {}

            for _bin in filter_intronic_bins:
                if _bin[0] not in filt_to_exclude:
                    dic[_bin[0]] = np.sum(df.intron_bin == _bin[0])

            plt.bar(dic.keys(), dic.values(), color="silver")

        else:
            ax = sns.countplot(x="intron_bin",
                            hue="outcome",
                            order=[x[0] for x in filter_intronic_bins if
                                    x[0] not in filt_to_exclude],
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
        
    outdir = os.path.join(outdir, "variant_counts")
    ensure_folder_exists(outdir)
    
    fname = os.path.join(outdir, 'counts_per_intronic_bin_{}'.format(variant_class))
    df = _df.copy()
    df.groupby(['intron_bin', 'outcome']).size().to_csv(fname + '.tsv',
                                                        sep="\t")
    
   
    filt_to_exclude = ["all_intronic", "all_except_1-2", "all_except_1-10"]
    if aggregate_classes:
        agg_bins = ['1-2', '3-10']
        opposite = ['1-10']
    else:
        agg_bins = ['1-10']
        opposite = ['1-2', '3-10']
        
    filt_to_exclude.extend(agg_bins)
    
    _generate_plot(df, filt_to_exclude)
    plt.savefig(fname + '.pdf')
    plt.close()
    
    df_zoom = df[(~df['intron_bin'].str.match('1-2')) &
                 (~df['intron_bin'].str.match('3-10')) &
                 (~df['intron_bin'].str.match('1-10'))]
   
    filt_to_exclude.extend(opposite)
    _generate_plot(df_zoom, filt_to_exclude)

    ylim = df_zoom['intron_bin'].value_counts().max() + (df_zoom['intron_bin'].value_counts().max() * 0.05)
    out_zoomed = fname + '_zoomed.pdf'
    plt.ylim(0, ylim)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper right")
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
        out = os.path.join(outdir, 'AF_distribution.pdf')
        plt.savefig(out)
        plt.close()


def plot_metrics_by_bin(df: pd.DataFrame, fname: str, aggregate_classes: bool):
    """
    :param pd.DataFrame df: Df with a list of
        metrics per each intronic bin
    :param str fname: Output basename
    :param bool aggregate_classes: Whether VETA run is meant to 
    aggregate classes into higher level concepts.
    """
    
    bins_with_scores = df.bin.unique().tolist()

    bins_to_exclude = ['1-2', '3-10'] if aggregate_classes else ['1-10']
    bins_to_exclude.extend(["all_intronic", "all_except_1-2", "all_except_1-10"])
    variant_class = df.name
    metrics = {"F1": "F1_Score",
               "weighted_F1": "F1 score (weighted)",
               "fraction_nan": "Fraction unscored"}
    
    avg = df.groupby('tool').agg({'F1':'mean',
                             'weighted_F1': 'mean',
                             'fraction_nan': 'mean'})
    
    avg = avg.round(2)
    avg.columns = ["avg_{}".format(x) for x in avg.columns]
    df = pd.merge(df, avg, on="tool", how='left')
   
    n_tools = df.tool.unique().size 
    if n_tools > 10:
        a = sns.diverging_palette(250, 30, l=65, n=n_tools, center="dark")
        fontsize="small"
    else:
        a = sns.color_palette("Paired", n_colors=n_tools)
        fontsize="medium"

    n_tools = df.tool.unique().size 
    # if n_tools > 10:
    #     sns.set_palette(sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True), n_colors=n_tools)
    #     fontsize="small"
    # else:
    #     sns.set_palette(sns.color_palette("Paired"), n_colors=n_tools)
    #     fontsize="medium"

    for metric, description in metrics.items():

        _df = df.copy()
        _df[metric] = pd.to_numeric(_df[metric])
        _df['tool'] = df.tool + " (mean=" + df["avg_{}".format(metric)].astype(str) + ")"
        
        if metric != "fraction_nan":
            _df = _df.sort_values("avg_{}".format(metric), ascending=False)
        else:
            _df = _df.sort_values("avg_{}".format(metric))

        sns.catplot(x="bin",
            y=metric,
            kind='point',
            order=[i[0] for i in filter_intronic_bins if i[0] not in bins_to_exclude and i[0] in bins_with_scores],
            data=_df,
            linestyles="-", 
            scale=0.75, 
            aspect=1.5,
            palette=a,
            markers='s',
            fontsize=fontsize,
            legend=False,
            dodge=True,
            hue="tool")

        plt.subplots_adjust(right=0.7)
        plt.legend(loc="upper center",
        bbox_to_anchor=(1.15, 1),
        borderaxespad=0.,
        prop=dict(size=7))

        plt.yticks(np.arange(0, 1.05, 0.1))

        plt.xlabel("Intron bin (bp)")
        plt.ylabel(description)
        plt.tight_layout()
        out = fname + '_' + variant_class + "_" + metric + '.pdf'

        plt.savefig(out, bbox_inches="tight")
        plt.close()
    