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
            for container in ax.containers:
                ax.bar_label(container)
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
                                 s=10,
                                 hue="outcome",
                                 palette={"Benign": "skyblue",
                                          "Pathogenic": "chocolate"})
            ax.get_legend().set_title('')

        ax.set(xlabel='Variant intronic offset', ylabel="{} allele frequency".format(af_column))
        plt.tight_layout()
        out = os.path.join(outdir, 'AF_distribution.pdf')
        plt.savefig(out)
        plt.close()


def plot_donor_vs_acceptor(df: pd.DataFrame, 
                           metric: str,
                           fname: str,
                           title: str):

    """
    Plot performance of individual tools on intronic variants
    associated with splicing donor vs splicing acceptors

    :param pd.DataFrame data: Df with stats for each tool
    :pararm str metric: Metric to rank tools
    :param str fname: Output file
    :param str title: Title of the plot
    """
    fig, ax = plt.subplots()
    plt.rcParams.update({'font.size': 10}) 
   
    if len(df['Splice site'].unique().tolist()) == 1:
        df = df.sort_values(metric, ascending=False)
        g = sns.catplot(data=df, 
                        kind="bar",
                        x=metric, 
                        y="Tool",
                        hue="Splice site",
                        palette={'Donor': 'darkred', 'Acceptor': 'skyblue'}, 
                        edgecolor='k',
                        legend_out=True,
                        alpha=.8, 
                        aspect=.8,
                        height=8)

        g.despine(left=True)

    # If both donors and acceptors    
    else:
        sns.stripplot(x=metric, y="Tool", hue='Splice site', linewidth=1, size=10, data=df)
        sns.boxplot(x=metric, y="Tool", data=df, color='white')
   
    plt.legend(loc="upper right",
               bbox_to_anchor=(1.1, 1),
               borderaxespad=0,
               prop=dict(size=8))
    plt.title(title) 
    plt.ylabel('')
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    plt.close()

def plot_metrics_by_bin_split_ss(df: pd.DataFrame, 
                        fname: str, 
                        aggregate_classes: bool = False):
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

    df['fraction_scored'] = 1 - df.fraction_nan
    metrics = {"F1": "F1_Score",
               "weighted_F1": "F1 score (weighted)",
               "fraction_scored": "Fraction scored"}

    avg = df.groupby(['tool', 'variant_class']).agg({'F1':'mean',
                             'weighted_F1': 'mean',
                             'fraction_scored': 'mean'})
    
    avg = avg.round(2)
    avg.columns = ["avg_{}".format(x) for x in avg.columns]

    df = pd.merge(df, avg, on=["tool", "variant_class"], how='left')

    n_tools = df.tool.unique().size 
    if n_tools > 10:
        pal = sns.diverging_palette(250, 30, l=65, n=n_tools, center="dark")
        fontsize="small"
    else:
        pal = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=False, reverse=True, n_colors=n_tools)
        fontsize="medium"
    
    categories=["501-1000", "201-500", "41-200", "11-40"]
    categories.append("1-10") if aggregate_classes else categories.extend(["3-10", "1-2"])
   
    hue_name = 'tool'
    legend_out = False
    for metric, description in metrics.items():
        _df = df.copy()
        d = df[df.variant_class.str.contains('donor')].copy()
        a = df[df.variant_class.str.contains('acceptor')].copy()
        a.bin = pd.Categorical(a.bin, 
                      categories=categories,
                      ordered=True)
        a['tool_with_metric'] = a.tool + " (avg_acc=" + a["avg_{}".format(metric)].astype(str) + ")"
        d['tool_with_metric'] = d.tool + " (avg_don=" + d["avg_{}".format(metric)].astype(str) + ")"
        
        _df = pd.concat([a, d])
        _df = _df.replace({'all_donor_related': 'Donor', 'all_acceptor_related': 'Acceptor'}).sort_values('variant_class')
        _df[metric] = pd.to_numeric(_df[metric])

        # Add avg metric per donor/acceptor
        mean_map = {}
        for _, row in _df.iterrows():
            tool = row.tool
            bin = row.bin
            val = row['avg_{}'.format(metric)]

            if tool in mean_map.keys():
                if bin in mean_map[tool].keys():
                    if ";" in  mean_map[tool][bin]:
                        continue
                    mean_map[tool][bin] = mean_map[tool][bin]  + ";{})".format(val)
                else:
                    mean_map[tool].update({bin: " ({}".format(val)})
            else:
                mean_map[tool] = {bin: " ({}".format(val)}
                         
        mean_df = pd.DataFrame.from_dict(mean_map, orient='index').rename_axis('tool').reset_index()
        mean_df = mean_df.melt(id_vars='tool', value_vars=list(mean_df), var_name='bin', value_name='tool_full')
        _df = pd.merge(_df, mean_df, on=['tool', 'bin'], how='left')
        _df['tool_with_metric'] = _df.tool + _df.tool_full
        hue_name = 'tool_with_metric'
        legend_out = True   
        # Comment code above (from # Add avg ..) to just display tool names in the middle of the two facets
        
        _df = _df.sort_values("avg_{}".format(metric), ascending=False)

        g = sns.catplot(x="bin",
            y=metric,
            kind='point',
            col='variant_class',
            col_order=['Acceptor', 'Donor'],
            order=[i[0] for i in filter_intronic_bins if i[0] not in bins_to_exclude and i[0] in bins_with_scores],
            data=_df,
            sharex=False,
            sharey=False,
            linestyles="-", 
            scale=0.75, 
            aspect=1.2,
            palette=pal,
            markers='s',
            fontsize=fontsize,
            legend= legend_out,
            legend_out = legend_out,
            dodge=False,
            hue=hue_name)

        g.axes[0,0].set_xlabel('Intron bin (bp)')
        g.axes[0,1].set_xlabel('Intron bin (bp)')
        g.axes[0,0].set_ylabel(description)
        g.axes[0,1].set_ylabel(description)

        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
        g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
        
        # Invert acceptor axis
        for i, ax in enumerate(g.axes[0]):
           if i == 0:
               ax.invert_xaxis()
        
        if legend_out:
            # Change label name
            g._legend.set_title("Tool (Mean acceptor; Mean donor)")
        else:
            # Increase space between the plots
            plt.subplots_adjust(hspace=0.4, wspace=0.4)   
            # Place legend in between the plots   
            g.fig.get_axes()[0].legend(loc='upper right', bbox_to_anchor=(1.3, 1)) 
         
        g.set(yticks=np.arange(0, 1.05, 0.1))
        out = fname + '_all_split_ss_' + metric + '.pdf'
        plt.savefig(out, bbox_inches="tight")
        plt.close()
        
def plot_metrics_by_bin(df: pd.DataFrame, 
                        fname: str, 
                        aggregate_classes: bool):
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
    
    df['fraction_scored'] = 1 - df.fraction_nan
    df['weighted_roc_auc'] = df.apply(lambda x: None if x.roc_auc is None else x.fraction_scored * x.roc_auc, axis = 1)
    df['weighted_ap_score'] = df.apply(lambda x: None if x.ap_score is None else x.fraction_scored * x.ap_score, axis = 1)
    
    # Keep only tools with auc and ap scores in all bins
    na_rows = df.loc[df.isnull()['weighted_roc_auc'] | df.isnull()['weighted_ap_score']]
    tools_to_remove = na_rows.tool.unique()
    
    _df = df[~df.tool.isin(tools_to_remove)]
    _avg = _df.groupby('tool').agg({'weighted_roc_auc': 'mean',
                                   'weighted_ap_score': 'mean'})
    _avg = _avg.rename(columns={'weighted_roc_auc': 'avg_weighted_roc_auc_all_bins',
                                'weighted_ap_score': 'avg_weighted_ap_score_all_bins'})
    
    metrics = {"F1": "F1_Score",
               "weighted_F1": "F1 score (weighted)",
               "fraction_scored": "Fraction scored",
               "roc_auc": "auROC",
               "weighted_roc_auc": "auROC (weighted)",
               "weighted_roc_auc_all_bins": "auROC (weighted)",
               "ap_score": "Average precision",
               "weighted_ap_score": "Average precision (weighted)",
               "weighted_ap_score_all_bins": "Average precision (weighted)"}
    
    avg = df.groupby('tool').agg({'F1':'mean',
                             'weighted_F1': 'mean',
                             'fraction_scored': 'mean',
                             'roc_auc': 'mean',
                             'ap_score': 'mean',
                             'weighted_roc_auc': 'mean',
                             'weighted_ap_score': 'mean'})
    
    avg = avg.round(3)
    avg.columns = ["avg_{}".format(x) for x in avg.columns]
    df = pd.merge(df, avg, on="tool", how='left').merge(_avg, on='tool', how='left')

    n_tools = df.tool.unique().size 
    if n_tools > 10:
        pal = sns.diverging_palette(250, 30, l=65, n=n_tools, center="dark")
        fontsize="small"
    else:
        pal = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=False, reverse=True, n_colors=n_tools)
        fontsize="medium"

    n_tools = df.tool.unique().size 

    for _metric, description in metrics.items():
        
        _df = df.copy()
        if "_all_bins" in _metric:
            _df = _df[~df['avg_{}'.format(_metric)].isnull()]   
        metric = _metric.replace("_all_bins", "") 
            
        _df[metric] = pd.to_numeric(_df[metric])
        _df['tool'] = _df.tool + " (mean=" + round(df["avg_{}".format(_metric)], 3).astype(str) + ")"
        
        _df = _df.sort_values("avg_{}".format(_metric), ascending=False)
        g = sns.catplot(x="bin",
            y=metric,
            kind='point',
            order=[i[0] for i in filter_intronic_bins if i[0] not in bins_to_exclude and i[0] in bins_with_scores],
            data=_df,
            linestyles="-", 
            scale=0.75, 
            aspect=1.5,
            palette=pal,
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
        
        if '_all_bins' in _metric:
            out = fname + '_' + variant_class + "_" + metric + '_tools_predicting_all_bins.pdf'    
        else:
            out = fname + '_' + variant_class + "_" + metric + '.pdf'

        plt.savefig(out, bbox_inches="tight")
        plt.close()
    