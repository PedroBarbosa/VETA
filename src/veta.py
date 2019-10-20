import os.path
import argparse
import osutils
from latex import generate_datasets_table
from predictions.scores_evaluation_unseen_vcf import inspect_predictions
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from filters import *
from thresholds import subset_toolset_by_scope, threshold_list_complete
from predictions.apply_tools import apply_tool_predictions
from predictions.performance_comparison import *
from plots.plot_heatmap import generate_heatmap
from preprocessing import preprocess, utils
from predictions.threshold_analysis import generate_threshold_analysis
from plots.ml.plot_feature_correlation import generate_ml_feature_correlation
from plots.ml import generate_ml_analysis
import warnings

warnings.filterwarnings("ignore")
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))


def run_standard_analysis(dataset, all_data, out_dir, clinvarStars, intronic_analysis, thresholdAnalysis,
                          scope_to_predict,
                          types_to_analyse,
                          metric_to_evaluate,
                          machineLearning,
                          skipHeatmap,
                          new_thresholds_done,
                          clinvar_from_file=None):
    if clinvar_from_file:
        dataset = "clinvar"

    if dataset == "clinvar" or clinvar_from_file:
        df_to_evaluate = all_data[clinvarStars]
        out_dir = os.path.join(out_dir, dataset, clinvarStars)
    else:
        dataset, clinvar_from_file = utils.check_dataset_arg(dataset)
        df_to_evaluate = all_data[dataset]
        out_dir = os.path.join(out_dir, dataset)

    generate_datasets_table(all_data, filters_var_type, out_dir)
    threshold_list = subset_toolset_by_scope(threshold_list_complete, scope_to_analyse=scope_to_predict)
    variant_types_to_evaluate = subset_variants_by_type(filters_var_type, types_to_analyse)
    generate_performance_comparison(df_to_evaluate, variant_types_to_evaluate, filters, threshold_list,
                                    metric_to_evaluate, dataset, out_dir)

    df = apply_tool_predictions(df_to_evaluate, threshold_list)
    if intronic_analysis:
        threshold_list_no_protein = subset_toolset_by_scope(threshold_list_complete, to_intronic=True)
        perform_intron_analysis(df_to_evaluate, filter_intronic_bins, threshold_list_no_protein, metric_to_evaluate,
                                dataset, out_dir)

    if skipHeatmap is False:
        generate_heatmap(df, variant_types_to_evaluate, filters, threshold_list, dataset, out_dir)

    if thresholdAnalysis and new_thresholds_done is False:
        new_thresholds = generate_threshold_analysis(all_data['3s_l'], filters, threshold_list, out_dir, 100)
        generate_performance_comparison_with_new_thresholds(df_to_evaluate, variant_types_to_evaluate,
                                                            filters, threshold_list, metric_to_evaluate,
                                                            dataset, out_dir, new_thresholds[1])

    if machineLearning:
        logging.info("---------------------------")
        logging.info("Starting ML analysis. ")
        logging.info("---------------------------")
        logging.info("Generating feature correlation matrix from scores.")
        out_dir = os.path.join(out_dir, "machineLearning")
        os.mkdir(out_dir)
        generate_ml_feature_correlation(df, threshold_list, out_dir)
        generate_ml_analysis(df, filters, threshold_list, dataset, out_dir)
        logging.info("All done!")
    return True if thresholdAnalysis else False


def main():
    """ Simple tool to evaluate variant prediction methods """
    parser = argparse.ArgumentParser(description='Script to trigger the full benchmark analysis')
    parser.add_argument('-o', '--out_dir', help='Path to store all the output results. Default: out_VETA')
    parser.add_argument('-l', '--location', default="HGVSc", choices=("HGVSc", "Consequence"),
                        help='VCF field to extract location of the variant. Default: "HGVSc')
    parser.add_argument('-s', metavar='scope_to_predict', help='Restrict analysis to a subset of tools based on'
                                                               ' their scope. Choices: (Conservation, Functional,'
                                                               ' Protein, Splicing)',
                        nargs='+', choices=('Conservation', 'Functional', 'Protein', 'Splicing'))
    parser.add_argument('-t', metavar='types_of_variant', help='Restrict analysis to the given variant types. '
                                                               'By default, performance is measured for all variants '
                                                               'and each subtype.. Choices: (all_types, snps, indels'
                                                               'insertions, deletions, mnps)',
                        nargs='+', choices=('all_types', 'snps', 'indels', 'insertions', 'deletions', 'mnps'))
    parser.add_argument('-i', '--intronic',
                        help='Perform additional analysis of intronic variants extracted from HGVSc '
                             'field.', action='store_true')

    parser.add_argument('--metric', metavar='metric', help="Metric to use to rank tools performance. Default:"
                                                           " weighted_accuracy. Choices: (weighted_accuracy, "
                                                           "accuracy, f1, weighted_f1, coverage)",
                        choices=('weighted_accuracy', 'accuracy', 'f1', 'weighted_f1', 'coverage'),
                        default='weighted_accuracy')

    referencesetmode = parser.add_argument_group('Tools performance on reference variant sets.')
    referencesetmode.add_argument('-d', '--datasets', metavar='datasets', type=lambda s: list(map(str, s.split(","))),
                                  default=['clinvar'],
                                  help='Run the analysis on the given dataset(s). Default: clinvar.'
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

    referencesetmode.add_argument('--listClinvarLevels', action='store_true',
                                  help='List the available clinvar filtering'
                                       'levels to run the analysis')

    referencesetmode.add_argument('--thresholdAnalysis', action="store_true", help="Enable"
                                                                                   " threshold analysis using Clinvar database. It does not depend on the --limitAnalysis "
                                                                                   "argument, and Clinvar 3 stars with likely assignments will be used as the ground truth")
    referencesetmode.add_argument('--machineLearning', action="store_true", help="Enable machine learning "
                                                                                 "analysis.")
    #  referencesetmode.add_argument('--phenotypeSpecific', action="store_true", help="Enable phenotype "
    #                                                                                 "specific analysis of Clinvar variants.It does not depend on the --limitAnalysis argument.")
    referencesetmode.add_argument('--skipHeatmap', action="store_true", help="Skip heatmap generation for performance"
                                                                             " analysis. If the dataset is large, this step takes very long time. Argument useful "
                                                                             "in such situations")

    runtoolsmode = parser.add_argument_group('Tools performance on an unseen VCF')
    runtoolsmode.add_argument('-v', '--vcf_file', metavar='vcf', dest='vcf_file', help="VCF file to evaluate tools "
                                                                                       "performance")
    runtoolsmode.add_argument('-tt', '--top_tools', metavar='top_tools',
                              help="Restrict analysis to the best set of tools"
                                   " obtained from a previous run using the reference variant catalogs (e.g. Clinvar)."
                                   " It needs the file 'tools_ranking*.csv' that is written when running the aforementioned "
                                   "analysis.")
    runtoolsmode.add_argument('-n', '--n_top_tools', metavar='int', type=int, default=5,
                              help="Number of top tools selected from the ranking provided in the '--top_tools'"
                                   " argument. Default: 5.")
    runtoolsmode.add_argument('--plot_tool', metavar='tool_name', type=lambda s: list(map(str, s.split(","))),
                              default=[], help="Plot score distribution of given tools. Use ',' to separate multiple "
                                               "tools.")
    runtoolsmode.add_argument('--labels', metavar='vcf_has_labels', choices=("Benign", "Pathogenic"),
                              help="Input VCF represents a list of labelled variants. Additional metrics will be"
                                   "inferred.")
    args = parser.parse_args()

    possible_datasets = ['clinvar', 'unifun', 'swissvar', 'humsavar', 'hcm_manually_curated']
    clinvar_stars = ['1s', '2s', '3s', '4s', '1s_l', '2s_l', '3s_l', 'clinvar', 'clinvar_l']
    unseen_vcf_analysis = False
    reference_vcf_analysis = True
    is_clinvar_from_file = False

    ############################
    ##Argparse args processing##
    ############################
    if args.listClinvarLevels:
        print("clinvar ---- Whole clinvar with Pathogenic and Benign assignments" + "\n" +
              "clinvar_l ---- Same as clinvar but with likely assignments" + "\n" +
              "1s --- clinvar with 1 star or more" + "\n" +
              "2s --- 2 stars or more" + "\n" +
              "3s --- 3 stars or more" + "\n" +
              "4s --- 4 stars" + "\n" +
              "1s_l ---- 1 star or more with likely assignments" + "\n" +
              "2s_l ---- 2 stars or more with likely assignments" + "\n" +
              "3s_l ---- 3 stars or more with likely assignments" + "\n" +
              "4s_l ---- 4 stars with likely assignments")
        exit(1)

    if args.clinvarStars not in clinvar_stars:
        logging.error(
            "Invalid clinvar stars argument. Please set '--listClinvarLevels' to see available options. If empty, "
            "defaults will be employed (3s_l).")
        exit(1)

    if args.vcf_file:
        if args.datasets and args.datasets[0] != "clinvar":
            logging.error(
                "--datasets and --vcf_file can't be set simultaneously. Please select only one type of analysis.")
            exit(1)
        osutils.check_file_exists(args.vcf_file)
        unseen_vcf_analysis = True
        reference_vcf_analysis = False
        args.datasets = None
    if args.top_tools:
        osutils.check_file_exists(args.top_tools)

    if args.top_tools and not args.vcf_file:
        logging.error("--vcf_file is required when --top_tools is set.")
        exit(1)

    if args.plot_tool and not set(args.plot_tool).issubset([i[0] for i in threshold_list_complete]):
        logging.error("Please set valid tools names in the --plot_tool argument")
        logging.info("List of available tools:\n{}".format("\n".join([i[0] for i in threshold_list_complete])))
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
                dataset, is_clinvar_from_file = utils.check_dataset_arg(dataset)
                os.mkdir(os.path.join(OUT_DIR, dataset))

                if dataset == "clinvar":
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars))
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars, "tools_benchmark"))
                    os.mkdir(os.path.join(OUT_DIR, dataset, args.clinvarStars, "preprocessing"))
                else:
                    os.mkdir(os.path.join(OUT_DIR, dataset, "tools_benchmark"))
                    os.mkdir(os.path.join(OUT_DIR, dataset, "preprocessing"))

        elif args.datasets and args.vcf_file:
            logging.error("You can't provide --vcf_file and --datasets arguments together. Please choose one of the "
                          "analysis")

    ###########################
    ####Unseen VCF analysis####
    ###########################
    if unseen_vcf_analysis:
        variant_types_to_evaluate = subset_variants_by_type(filters_var_type, args.t)
        inspect_predictions(args.vcf_file, args.top_tools, args.n_top_tools, args.plot_tool, args.s,
                            variant_types_to_evaluate, args.labels, args.location, args.intronic, args.metric, OUT_DIR)

    ##################################
    #####Reference datasets analysis##
    ##################################
    if reference_vcf_analysis:
        if args.datasets:
            # variable set to run threshold analysis only once, even if multiple datasets are provided
            new_thresholds_done = False
            for dataset in args.datasets:
                if dataset in possible_datasets:
                    all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
                                          args.intronic,
                                          is_clinvar_from_file,
                                          dataset=dataset)

                else:
                    all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
                                          args.intronic,
                                          is_clinvar_from_file,
                                          dataset=dataset, newDataset=True)

                new_thresholds_done = run_standard_analysis(dataset, all_data, OUT_DIR, args.clinvarStars,
                                                            args.intronic,
                                                            args.thresholdAnalysis,
                                                            args.s,
                                                            args.t,
                                                            args.metric,
                                                            args.machineLearning,
                                                            args.skipHeatmap,
                                                            new_thresholds_done,
                                                            clinvar_from_file=is_clinvar_from_file)

        else:
            clinvar_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
                                      args.intronic, is_clinvar_from_file)
            run_standard_analysis("3s_l", clinvar_data, OUT_DIR, args.clinvarStars, args.intronic,
                                  args.thresholdAnalysis,
                                  args.s,
                                  args.t,
                                  args.metric,
                                  args.machineLearning,
                                  args.skipHeatmap, False)


if __name__ == '__main__':
    main()
