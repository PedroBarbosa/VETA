import os.path
import argparse
import osutils
from latex import generate_datasets_table
from predictions.scores_evaluation_unseen_vcf import inspect_predictions
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
from filters import *
from thresholds import threshold_list
from predictions.apply_tools import apply_tool_predictions
from predictions.performance_comparison import generate_performance_comparison
from plots.plot_heatmap import generate_heatmap
from preprocessing import preprocess
from predictions.threshold_analysis import generate_threshold_analysis
from plots.ml.plot_feature_correlation import generate_ml_feature_correlation
from plots.ml import generate_ml_analysis
import warnings
warnings.filterwarnings("ignore")
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))


def run_standard_analysis(dataset, all_data, out_dir, clinvarStars, thresholdAnalysis, machineLearning, skipHeatmap, new_thresholds_done):
    if dataset == "clinvar":
        out_dir = os.path.join(out_dir, dataset, clinvarStars)
    else:
        out_dir = os.path.join(out_dir, dataset)

    generate_datasets_table(all_data, filters_var_type, out_dir)
    generate_performance_comparison(all_data[dataset], filters_var_type, filters, threshold_list, dataset, out_dir)
    df = apply_tool_predictions(all_data[dataset], threshold_list)
    if skipHeatmap is False:
        generate_heatmap(df, filters_var_type, filters, threshold_list, dataset, out_dir)

    if thresholdAnalysis and new_thresholds_done is False:
        new_thresholds = generate_threshold_analysis(all_data['3s_l'], filters, threshold_list, '3s_l', out_dir, 100)
        generate_performance_comparison(all_data[dataset], filters_var_type, filters, threshold_list, dataset, out_dir,
                                    new_thresholds=new_thresholds[1])

    if machineLearning:
        generate_ml_feature_correlation(df, dataset, threshold_list, out_dir)
        generate_ml_analysis(df, filters, threshold_list, dataset, out_dir)

    return True if thresholdAnalysis else False


def main():

    """ Framework to analyse variant prediction methods """
    parser = argparse.ArgumentParser(description='Script to trigger the full benchmark analysis')
    parser.add_argument('-o', '--out_dir', help='Path to store all the output results. Default: out_VETA')
    parser.add_argument('-l', '--location', default="HGVSc", choices=("HGVSc","Consequence"),
                        help='VCF field to extract location of the variant. Default: "HGVSc')

    referencesetmode = parser.add_argument_group('Tools performance on reference variant sets.')
    referencesetmode.add_argument('-d', '--datasets', metavar='datasets', type=lambda s: list(map(str, s.split(","))),
                                  default=['clinvar'], help='Run the analysis on the given dataset(s). Default: clinvar.'
                                'User may provide the path to the dataset name, where \'benign\' and \'pathogenic\' tags'
                                ' must exist in the filenames. It is possible to analyze one of the following datasets,'
                                ' located within the project \'datasets\' folder:[clinvar, hcm_manually_curated,unifun,'
                                'swissvar,humsavar]')
    referencesetmode.add_argument('-c', '--clinvarStars', metavar='stars', default='3s_l',
                                  help='Level of filtering when running clinvar analysis. By default, a high confidence'
                                       'clinvar subset (3starts with likely annotations) is used for tool performance'
                                       'evaluation and reference threshold analysis. Value: "3s_l". To see all the'
                                       'possible filtering levels, run the software with the argument '
                                       '--listClinvarLevels.')

    referencesetmode.add_argument('--listClinvarLevels', action='store_true', help='List the available clinvar filtering'
                                        'levels to run the analysis')

    referencesetmode.add_argument('--thresholdAnalysis', action="store_true", help="Enable"
                            " threshold analysis using Clinvar database. It does not depend on the --limitAnalysis "
                            "argument.")
    referencesetmode.add_argument('--machineLearning', action="store_true", help="Enable machine learning "
                            "analysis.")
    referencesetmode.add_argument('--phenotypeSpecific', action="store_true", help="Enable phenotype "
                            "specific analysis of Clinvar variants.It does not depend on the --limitAnalysis argument.")
    referencesetmode.add_argument('--skipHeatmap', action="store_true", help="Skip heatmap generation for performance"
                            "analysis. If the dataset is large, this step takes very long time. Argument useful "
                            "in such situations")

    runtoolsmode = parser.add_argument_group('Tools performance on an unseen VCF')
    runtoolsmode.add_argument('-v', '--vcf_file', metavar='vcf', dest='vcf_file', help="VCF file to evaluate tools "
                        "performance")
    runtoolsmode.add_argument('-t', '--top_tools', metavar='top_tools',help="Restrict analysis to the best set of tools"
                        " obtained from a previous run using the reference variant catalogs (e.g. Clinvar)."
                        " It needs the file 'tools_ranking*.csv' that is written when running the aforementioned "
                        "analysis.")
    runtoolsmode.add_argument('-n', '--n_top_tools', metavar='int', type=int, default=5, help="Number of top "
                        "tools selected from the ranking provided in the '--top_tools' argument. Default_5.")
    runtoolsmode.add_argument('--plot_tool', metavar='tool_name', type=lambda s: list(map(str, s.split(","))),
                        default=[], help="Plot score distribution of given tools. Use ',' to separate multiple "
                        "tools.")
    args = parser.parse_args()

    possible_datasets = ['clinvar', 'unifun', 'swissvar', 'humsavar', 'hcm_manually_curated']
    unseen_vcf_analysis=False
    reference_vcf_analysis=True

    ############################
    ##Argparse args processing##
    ############################
    if args.listClinvarLevels:
        # dict = {
        #     'clinvar': df_clinvar_sure,
        #     '1s': filter_clinvar_1_star(df_clinvar_sure),
        #     '2s': filter_clinvar_2_stars(df_clinvar_sure),
        #     '3s': filter_clinvar_3_stars(df_clinvar_sure),
        #     '4s': filter_clinvar_4_stars(df_clinvar_sure),
        #     'clinvar_l': df_clinvar,
        #     '1s_l': filter_clinvar_1_star(df_clinvar),
        #     '2s_l': filterrm _clinvar_2_stars(df_clinvar),
        #     '3s_l': filter_clinvar_3_stars(df_clinvar),
        #     '4s_l': filter_clinvar_4_stars(df_clinvar),
        # }
        print("Possible filtering values:")


    if args.vcf_file:
        osutils.check_file_exists(args.vcf_file)
        unseen_vcf_analysis=True
        reference_vcf_analysis = False

    if args.top_tools:
        osutils.check_file_exists(args.top_tools)

    if args.top_tools and not args.vcf_file:
        logging.error("--vcf_file is required when --top_tools is set.")
        exit(1)

    if args.plot_tool and not set(args.plot_tool).issubset([i[0] for i in threshold_list]):
        logging.error("Please set valid tools names in the --plot_tool argument")
        logging.info("List of available tools:\n{}".format("\n".join([i[0] for i in threshold_list])))
        exit(1)

    ###########################
    ##Output directory setup###
    ###########################
    if args.out_dir:
        OUT_DIR = args.out_dir
    else:
        OUT_DIR = os.path.join(os.getcwd(), "out_VETA")

    if os.path.isdir(OUT_DIR):
        logging.info("Output directory {} exists. Please remove it or select a new one.".format(OUT_DIR))
        exit(1)

    else:
        os.mkdir(OUT_DIR)
        if not args.vcf_file:
            for dataset in args.datasets:
                os.mkdir(os.path.join(OUT_DIR, dataset))
                if dataset == "clinvar":
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars))
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars, "figures"))
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars, "preprocessing"))
                else:
                    os.mkdir(os.path.join(OUT_DIR, dataset, "figures"))
                    os.mkdir(os.path.join(OUT_DIR, dataset, "preprocessing"))

        elif args.datasets and args.vcf_file:
            logging.error("You can't provide --vcf_file and --datasets arguments together. Please choose one of the "
                          "analysis")

    ###########################
    ####Unseen VCF analysis####
    ###########################
    if unseen_vcf_analysis:
        inspect_predictions(args.vcf_file, args.top_tools, args.n_top_tools, args.plot_tool, args.location, args.out_dir)


    ##################################
    #####Reference datasets analysis##
    ##################################
    if reference_vcf_analysis:
        if args.datasets:
            # variable set to run threshold analysis only once, even if multiple datasets are provided
            new_thresholds_done = False
            for dataset in args.datasets:
                if dataset in possible_datasets:
                    all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis, dataset=dataset)
                else:
                    all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,dataset=dataset,
                                          newDataset=True)

                new_thresholds_done = run_standard_analysis(dataset, all_data, OUT_DIR, args.clinvarStars, args.thresholdAnalysis,
                                                            args.machineLearning, args.skipHeatmap, new_thresholds_done)

        else:
            clinvar_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis)
            run_standard_analysis("3s_l", clinvar_data, OUT_DIR, args.thresholdAnalysis,
                                  args.machineLearning, args.skipHeatmap, False)


if __name__ == '__main__':
    main()
