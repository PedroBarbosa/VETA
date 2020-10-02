import os.path
import argparse
import sys
import logging
from src.preprocessing.osutils import print_clinvar_levels
from src.inspect import PredictionsEval
from src.benchmark import BenchmarkTools
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from src.predictions.filters import *

import warnings

warnings.filterwarnings("ignore")
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))


# def run_standard_analysis(dataset, all_data, out_dir, clinvarStars, intronic_analysis, thresholdAnalysis,
#                           scope_to_predict,
#                           types_to_analyse,
#                           metric_to_evaluate,
#                           machineLearning,
#                           skipHeatmap,
#                           new_thresholds_done,
#                           clinvar_from_file=None):
#
#
#     if skipHeatmap is False:
#         generate_heatmap(df, variant_types_to_evaluate, filters, threshold_list, dataset, out_dir)
#
#     if thresholdAnalysis and new_thresholds_done is False:
#         new_thresholds = generate_threshold_analysis(all_data['3s_l'], filters, threshold_list, out_dir, 100)
#         generate_performance_comparison_with_new_thresholds(df_to_evaluate, variant_types_to_evaluate,
#                                                             filters, threshold_list, metric_to_evaluate,
#                                                             dataset, out_dir, new_thresholds[1])
#
#     if machineLearning:
#         logging.info("---------------------------")
#         logging.info("Starting ML analysis. ")
#         logging.info("---------------------------")
#         logging.info("Generating feature correlation matrix from scores.")
#         out_dir = os.path.join(out_dir, "machineLearning")
#         os.mkdir(out_dir)
#         generate_ml_analysis(df, filters, threshold_list, dataset, out_dir)
#         logging.info("All done!")
#     return True if thresholdAnalysis else False


def main():
    parser = argparse.ArgumentParser(prog="veta", description='Simple tool to evaluate variant '
                                                              'prediction methods')

    subparsers = parser.add_subparsers(dest="command")

    # Parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-o', metavar='--out_dir', help='Path to store all the output results. '
                                                               'Default: "out_VETA"')

    parent_parser.add_argument('-s', metavar='--scope_to_predict',
                               help='Restrict analysis to a subset of tools based on their '
                                    'scope. Available options: {%(choices)s}. Default: Tools from '
                                    'all scores are used.',
                               nargs='+', choices=('Conservation', 'Functional', 'Protein', 'Splicing'))

    parent_parser.add_argument('-t', metavar='types_of_variant',
                               help='Restrict analysis to the given variant types. '
                                    'Available options: {%(choices)s}. Default: Performance '
                                    'is measured for all variants and each subtype.',
                               nargs='+', choices=('all_types', 'snps', 'indels', 'insertions', 'deletions', 'mnps'))

    parent_parser.add_argument('-m', metavar='--metric', help='Metric to rank the tools. Available options: '
                                                              '{%(choices)s}. Default: "weighted_accuracy."',
                               choices=('weighted_accuracy', 'accuracy', 'F1', 'weighted_F1', 'coverage'),
                               default='weighted_accuracy')

    parent_parser.add_argument('-l', metavar='--location', default="HGVSc", choices=("HGVSc", "Consequence"),
                               help='VCF field to extract location of the variant. Available options: '
                                    '{%(choices)s}. Default: "HGVSc".')

    parent_parser.add_argument('-g', metavar='--genome', default="hg19", choices=("hg19", "hg38"),
                               help='Genome build of the VCF. Available options: {%(choices)s}. '
                                    'Default: "hg19".')

    parent_parser.add_argument('-i', '--is_intronic', help='Perform additional analysis of '
                                                           'intronic variants extracted from '
                                                           'HGVSc field.', action='store_true')
    parent_parser.add_argument('--config', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                "map_tools2vcf_annotation.txt"),
                               help='Path to the config file that maps tools to the corresponding '
                                    'VCF annotation. Default: \'map_tools2vcf_annotation.txt\' file '
                                    'in the src code directory')

    # Subparsers based on parent
    benchmark_parser = subparsers.add_parser("benchmark", help='Benchmark prediction tools based on labelled data.',
                                             parents=[parent_parser])

    benchmark_parser.add_argument(dest='dataset',
                                  help='Path to the reference dataset. If a single file is '
                                       'given, veta expects the dataset to be from clinvar. If a previous '
                                       'run of VETA was performed, a cached file with the '
                                       'addition of \'tsv\' will exist, thus making processed clinvar data '
                                       'to be loaded much quicker. If a directory is '
                                       'provided instead, veta looks for labelled files within it. The '
                                       '\'*benign*\' (or \'*neutral*\' and \'*pathogenic*\' (or \'*deleterious*\') '
                                       'tags must exist in their names.')

    benchmark_parser.add_argument('-c', metavar='--clinvar_stars', default='3s_l',
                                  help='Level of filtering when dataset refers to the clinvar database. '
                                       'By default, a high confidence clinvar subset (3 stars with likely '
                                       'annotations) is used for performance evaluation and reference '
                                       'threshold analysis (if --do_threshold_analysis is True). '
                                       'Default: "3s_l". All the possible filtering levels are visible '
                                       'with the argument \'--listClinvarLevels\'.')

    benchmark_parser.add_argument('--list_clinvar_levels', action='store_true',
                                  help='List the available clinvar filtering '
                                       'levels to run the analysis.')

    benchmark_parser.add_argument('--do_threshold_analysis', action="store_true",
                                  help="Enable reference thresholds analysis when Clinvar is used. "
                                       "It does not depend on the \'--clinvar_stars\' argument, which "
                                       "means that \'3s_l\' variants will be used as the ground truth "
                                       "to ensure that only good quality variants are used. "
                                       "Default: False")

    benchmark_parser.add_argument('--do_machine_learning', action="store_true",
                                  help="Enable machine learning analysis based on the "
                                       "tools scores to inspect the best predictors and "
                                       "create an ensemble classifier that combines multiple"
                                       " scores. Default: False.")

    benchmark_parser.add_argument('--skip_heatmap', action="store_true",
                                  help="Skip heatmap generation for performance analysis. "
                                       "If the dataset is large, this step takes quite a "
                                       "while. very long time.")

    predict_parser = subparsers.add_parser("inspect", help='Evaluate predictions from unlabelled variants.',
                                           parents=[parent_parser])

    predict_parser.add_argument(dest='vcf', help="VCF file to evaluate tools predictions.")

    predict_parser.add_argument('-b', metavar='--best_tools',
                                help="Restrict analysis to the best set of tools "
                                     "obtained from a previous run using a reference "
                                     "catalog (e.g. Clinvar). It must refer to the file "
                                     "\'tools_ranking*.csv\' that is written when running "
                                     "the aforementioned analysis. Default: Use all tools "
                                     "available in the analysis.")

    predict_parser.add_argument('-n', metavar='--n_best_tools', type=int, default=5,
                                help="Number of best tools selected from the ranking "
                                     "provided in the \'--best_tools\' argument. Default: 5.")

    predict_parser.add_argument('--plot_these_tools', metavar='tool_name', nargs='+',
                                help="Plot scores distribution for the given tools.")

    predict_parser.add_argument('--labels', metavar='vcf_has_labels',
                                choices=("Benign", "Pathogenic"),
                                help="If VCF represents a list of labelled variants, "
                                     "additional metrics will be inferred. This "
                                     "argument asks for the label type of input vcf. "
                                     "Available options: {%(choices)s}. Default: VCF "
                                     "does not have variants with known labels.")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    clinvar_stars = ['1s', '2s', '3s', '4s', '1s_l', '2s_l', '3s_l', 'clinvar', 'clinvar_l']
    is_clinvar_from_file = False

    ##############################
    ## Argparse args processing ##
    ##############################

    if args.command == "inspect":
        PredictionsEval(args.vcf,
                        args.o,
                        args.s,
                        args.t,
                        args.m,
                        args.l,
                        args.g,
                        args.is_intronic,
                        args.b,
                        args.n,
                        args.plot_these_tools,
                        args.labels,
                        args.config
                        )
    elif args.command == "benchmark":
        if args.list_clinvar_levels:
            print_clinvar_levels()

        assert args.c in clinvar_stars, "Set a valid value for the '--clinvar_stars' argument."

        BenchmarkTools(args.dataset,
                       args.o,
                       args.s,
                       args.t,
                       args.m,
                       args.l,
                       args.g,
                       args.is_intronic,
                       args.c,
                       args.do_threshold_analysis,
                       args.do_machine_learning,
                       args.skip_heatmap,
                       args.config)




    # ##################################
    # #####Reference datasets analysis##
    # ##################################
    # if reference_vcf_analysis:
    #     if args.datasets:
    #         # variable set to run threshold analysis only once, even if multiple datasets are provided
    #         new_thresholds_done = False
    #         for dataset in args.datasets:
    #             if dataset in possible_datasets:
    #                 all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
    #                                       args.intronic,
    #                                       is_clinvar_from_file,
    #                                       dataset=dataset)
    #
    #             else:
    #                 all_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
    #                                       args.intronic,
    #                                       is_clinvar_from_file,
    #                                       dataset=dataset, newDataset=True)
    #
    #             new_thresholds_done = run_standard_analysis(dataset, all_data, OUT_DIR, args.clinvarStars,
    #                                                         args.intronic,
    #                                                         args.thresholdAnalysis,
    #                                                         args.s,
    #                                                         args.t,
    #                                                         args.metric,
    #                                                         args.machineLearning,
    #                                                         args.skipHeatmap,
    #                                                         new_thresholds_done,
    #                                                         clinvar_from_file=is_clinvar_from_file)
    #
    #     else:
    #         clinvar_data = preprocess(args.location, ROOT_DIR, args.thresholdAnalysis,
    #                                   args.intronic, is_clinvar_from_file)
    #         run_standard_analysis("3s_l", clinvar_data, OUT_DIR, args.clinvarStars, args.intronic,
    #                               args.thresholdAnalysis,
    #                               args.s,
    #                               args.t,
    #                               args.metric,
    #                               args.machineLearning,
    #                               args.skipHeatmap, False)


if __name__ == '__main__':
    main()
