import os.path
import argparse
import numpy as np
import pandas as pd
from datasets import preprocess
from latex import generate_datasets_table

from filters import *
from thresholds import threshold_list
from tools import apply_tool_predictions

from plots.threshold_analysis import generate_threshold_analysis
from plots.performance_comparison import generate_performance_comparison
from plots.heatmap import generate_heatmap
from plots.predictions import *
from plots.ml import generate_ml_analysis, generate_ml_feature_correlation
from datasets import vep,vep_cleaning

PAPER_LOCATION = "../../Paper_OptimizingPipelines"

def main():
    """ Generates all the tables and plots for the paper """
    parser = argparse.ArgumentParser(description='Script to trigger the full benchmark analysis')
    parser.add_argument("--out_dir", help='Path to store all the output results. Default: "../../Paper_OptimizingPipelines"')
    parser.add_argument("--location", default="HGVSc", choices=("HGVSc","Consequence"),
                        help='VCF field to extract location of the variant')
    referencesetmode = parser.add_argument_group('Tools performance on reference variant sets')
    referencesetmode.add_argument("--limitAnalysis", metavar='dataset', type=lambda s: list(map(str, s.split(","))), default=[],
                        help='Set the analysis to be performed, other than clinvar. Default: Run clinvar. '
                             'Choices: [hcm_beatriz,unifun,swissvar,humsavar,cardiomyopathy]')

    runtoolsmode = parser.add_argument_group('Tools performance on an unseen VCF')
    runtoolsmode.add_argument('-vcf_file', dest='vcf_file', help="VCF to evaluate tools performance")
    runtoolsmode.add_argument('--top_tools', metavar='top_tools', help="If set, perform analysis on the set of tools represented in the file.")
    runtoolsmode.add_argument('--n_top_tools', metavar='n_top_tools', type=int, default=5, help="If set, it refers to the number of top tools to perform the analysis. Default_5.")
    runtoolsmode.add_argument('--plot_tool', metavar='plot_specific_tool', type=lambda s: list(map(str, s.split(","))), default=[],
                              help="Plot score distribution of given tools. Use ',' to separate multiple tools.")
    args = parser.parse_args()

    datasets = ['unifun', 'swissvar', 'humsavar', 'cardiomyopathy', 'hcm_beatriz']
    if args.out_dir :
        PAPER_LOCATION = args.out_dir
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)
            if not args.vcf_file:
                os.mkdir(os.path.join(args.out_dir,"figures"))
                os.mkdir(os.path.join(args.out_dir,"datasets"))
    else:
        PAPER_LOCATION = "../../Paper_OptimizingPipelines"

    if args.top_tools and not args.vcf_file:
        print("--vcf_file is required when --toptools is set.")
        exit(1)

    if args.plot_tool and not set(args.plot_tool).issubset([i[0] for i in threshold_list]):
        print("Please set valid tools in the --plot_tool argument")
        exit(1)

    if args.vcf_file and os.path.isfile(args.vcf_file):
        df = vep.get_df_ready(args.vcf_file,False,True, args.location)
        df = vep_cleaning(df)
        df = apply_tool_predictions(df,threshold_list).set_index('id')
        df_t = pd.concat([df[[tool + "_prediction" for tool, *args in threshold_list]], df["HGVSc"],df["location"].to_frame()], axis=1).copy()
        df_t = df_t.rename(columns={col: col.split('_')[0] for col in df_t.columns})
        df_new = df_t.drop(["location", "HGVSc"], axis = 1).apply(lambda x: x.value_counts(True, dropna=False), axis=1).fillna(0).sort_values([True],ascending=False)
        df_new.rename(columns={False: 'is_benign', True: 'is_pathogenic', np.nan: "unpredictable"}, inplace=True)
        plot_area(df_new, PAPER_LOCATION)
        df_new["unpredictable"] *= 100
        df_t.to_csv(os.path.join(PAPER_LOCATION, "predictions.tsv"), sep="\t")
        plot_heatmap(df_new,PAPER_LOCATION, False)
        plot_heatmap(df_new[df_new.is_pathogenic > 0.75], PAPER_LOCATION,True)

        if args.top_tools and os.path.isfile(args.top_tools):
            tools=pd.read_csv(args.top_tools, sep="\t")['tool'].head(args.n_top_tools).tolist()
            tools.append("location")
            df_top=df_t[tools]
            booldic = {True: 1, False: -1}
            plot_heatmap_toptools(df_top.replace(booldic), filters, PAPER_LOCATION)

        if args.plot_tool:
            for tool in args.plot_tool:
                plot_tool_score_distribution(df, tool, PAPER_LOCATION)


    elif args.limitAnalysis and all(x in datasets for x in args.limitAnalysis):
        for analysis in args.limitAnalysis:
            datasets = preprocess(args.location, dataset=analysis)
            generate_datasets_table(datasets, filters_var_type, PAPER_LOCATION)
            #new_thresholds = generate_threshold_analysis(datasets['3s_l'], filters, threshold_list, '3s_l', PAPER_LOCATION, 100)
            generate_performance_comparison(datasets[analysis], filters_var_type, filters, threshold_list, analysis, PAPER_LOCATION)
            #generate_performance_comparison(datasets[analysis], filters_var_type, filters, threshold_list, analysis, PAPER_LOCATION,
            #                                new_thresholds=new_thresholds[1])

            df = apply_tool_predictions(datasets[analysis], threshold_list)
            generate_heatmap(df, filters_var_type, filters, threshold_list, analysis, PAPER_LOCATION)
            #generate_ml_feature_correlation(df, analysis, threshold_list, PAPER_LOCATION)
            #generate_ml_analysis(df, filters, threshold_list,analysis, PAPER_LOCATION)

    elif not args.limitAnalysis:
        datasets = preprocess(args.location)
        generate_datasets_table(datasets, filters_var_type, PAPER_LOCATION)

        #new_thresholds = generate_threshold_analysis(datasets['3s_l'], filters, threshold_list, '3s_l', PAPER_LOCATION,1000)
        generate_performance_comparison(datasets['3s_l'], filters_var_type, filters, threshold_list, "3s_l", PAPER_LOCATION)
        #generate_performance_comparison(datasets['3s_l'], filters_var_type, filters, threshold_list, "3s_l", PAPER_LOCATION,
        #                                new_thresholds=new_thresholds[1])

        d3s_l = apply_tool_predictions(datasets['3s_l'], threshold_list)
        generate_heatmap(d3s_l, filters_var_type, filters, threshold_list, "3s_l", PAPER_LOCATION)
        #generate_ml_feature_correlation(d3s_l, "3s_l", threshold_list, PAPER_LOCATION)
        #generate_ml_analysis(d3s_l, filters, threshold_list, "3s_l", PAPER_LOCATION)

    else:
        print("Please limit your analysis to one (or more) of the following dataset options:\t{}".format(datasets))
        exit(1)


if __name__ == '__main__':
    main()
