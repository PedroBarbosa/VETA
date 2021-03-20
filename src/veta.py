import argparse
import logging
import os.path
import sys
import pkgutil, pkg_resources

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
hgvs_logger = logging.getLogger('hgvs')
hgvs_logger.setLevel(logging.CRITICAL)
scikit_logger = logging.getLogger('sklearn')
scikit_logger.setLevel(logging.CRITICAL)

from src.benchmark import BenchmarkTools
from src.inspect import PredictionsEval
from src.preprocessing.osutils import print_clinvar_levels

CONFIG_PATH = pkg_resources.resource_filename('src', 'config/tools_config.txt')


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(prog="veta",
                                     description='Simple tool to evaluate variant '
                                                 'prediction methods')

    subparsers = parser.add_subparsers(dest="command")

    # Parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-o', '--out_dir', metavar='', help='Path to store all the output results. '
                                                                   'Default: "out_VETA"')

    parent_parser.add_argument('-s', '--scopes_to_evaluate', metavar='',
                               help='Restrict analysis to a subset of tools based on their '
                                    'scope. Available options: {%(choices)s}. Default: Tools from '
                                    'all scores are used.',
                               nargs='+', choices=('Conservation', 'Functional', 'Protein', 'Splicing'))

    parent_parser.add_argument('-t', '--types_of_variant', metavar='',
                               help='Restrict analysis to the given variant types. '
                                    'Available options: {%(choices)s}. Default: Performance '
                                    'is measured for all variants and each subtype.',
                               nargs='+', choices=('all_types', 'snps', 'indels', 'insertions', 'deletions', 'mnps'))

    parent_parser.add_argument('-m', '--metric', metavar='',
                               help='Metric to rank the tools. Available options: '
                                    '{%(choices)s}. Default: "weighted_accuracy."',
                               choices=('weighted_accuracy', 'accuracy', 'F1', 'weighted_F1', 'coverage'),
                               default='weighted_accuracy')

    parent_parser.add_argument('-l', '--location', metavar='', default="HGVSc", choices=("HGVSc", "Consequence"),
                               help='VCF field to extract location of the variant. Available options: '
                                    '{%(choices)s}. Default: "HGVSc".')

    parent_parser.add_argument('-g', '--genome', metavar='', default="hg19", choices=("hg19", "hg38"),
                               help='Genome build of the VCF. Available options: {%(choices)s}. '
                                    'Default: "hg19".')

    parent_parser.add_argument('-i', '--do_intronic_analysis', help='Perform additional analysis of '
                                                                    'intronic variants extracted from '
                                                                    'HGVSc field in a bin-based faction. '
                                                                    'Bins are attributed based on the distance '
                                                                    'of the variant to the nearest splice junction',
                               action='store_true')

    parent_parser.add_argument('-a', '--allele_frequency', metavar='', default="gnomADg_AF",
                               help='VCF field (within VEP annotations, or INFO field) '
                                    'that measures frequency of the variant in a '
                                    'population. If it exists in the input data, additional '
                                    'plots will be drawn. If absent, VETA will simply ignore it. '
                                    'Default: "gnomADg_AF". Note: Missing data will be treated as '
                                    'if the variant is absent in population (converted to 0). '
                                    'Fo example, if ExAC frequencies are given, intronic variants '
                                    'will be given 0, but that does not mean that the variant does '
                                    'not exist in gnomAD, for example. In this case, '
                                    'it would be appropriate to just analyze the frequency plots of '
                                    'coding variants.')

    parent_parser.add_argument('--skip_heatmap', action="store_true",
                               help="Skip heatmap generation for performance analysis. "
                                    "If the dataset is large, this step takes quite a "
                                    "while. very long time.")

    parent_parser.add_argument('--config', default=CONFIG_PATH,
                               help='Path to the config file that maps tools to the corresponding '
                                    'VCF annotation. Default: \'tools_config.txt\' file '
                                    'in the src/config github directory')
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
                                       '\'*benign*\' (or \'*neutral*\') and \'*pathogenic*\' (or \'*deleterious*\') '
                                       'tags must exist in their names.')

    benchmark_parser.add_argument('-c', '--clinvar_stars', metavar='', default='3s_l',
                                  help='Level of filtering when dataset refers to the clinvar database. '
                                       'By default, a high confidence clinvar subset (3 stars with likely '
                                       'annotations) is used for performance evaluation and reference '
                                       'threshold analysis (if --do_threshold_analysis is True). '
                                       'Default: "3s_l". All the possible filtering levels are visible '
                                       'with the argument \'--listClinvarLevels\'.')

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

    predict_parser = subparsers.add_parser("inspect", help='Evaluate predictions from unlabelled variants.',
                                           parents=[parent_parser])

    predict_parser.add_argument(dest='vcf', help="VCF file to evaluate tools predictions.")

    predict_parser.add_argument('-b', '--best_tools', metavar='',
                                help="Restrict heatmap analysis to the best set of tools "
                                     "obtained from a previous run using a reference "
                                     "catalog (e.g. Clinvar). It must refer to the file "
                                     "\'tools_ranking*.csv\' that is written when running "
                                     "the aforementioned analysis. Default: Use all tools "
                                     "available in the analysis. It requires '--skip_heatmap' "
                                     "is 'False'")

    predict_parser.add_argument('-n', '--n_best_tools', metavar='', type=int, default=5,
                                help="Number of best tools selected from the ranking "
                                     "provided in the \'--best_tools\' argument. Default: 5.")

    predict_parser.add_argument('--plot_these_tools', metavar='tool_name', nargs='+',
                                help="Plot scores distribution for the given tools.")

    predict_parser.add_argument('--labels', metavar='label_type',
                                choices=("Benign", "Neutral", "Pathogenic",
                                         "Deleterious", "Functional"),
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
    ##############################
    ## Argparse args processing ##
    ##############################
    if args.command == "inspect":
        PredictionsEval(args.vcf,
                        args.out_dir,
                        args.scopes_to_evaluate,
                        args.types_of_variant,
                        args.metric,
                        args.location,
                        args.genome,
                        args.do_intronic_analysis,
                        args.best_tools,
                        args.n_best_tools,
                        args.plot_these_tools,
                        args.labels,
                        args.allele_frequency,
                        args.skip_heatmap,
                        args.config)

    elif args.command == "benchmark":
        if args.clinvar_stars not in clinvar_stars:
            logging.info("Error. Set a valid value for the '--clinvar_stars' argument. Possible options "
                         "(value with a simple description):")
            print_clinvar_levels()

        assert args.clinvar_stars in clinvar_stars, "Set a valid value for the '--clinvar_stars' argument."

        BenchmarkTools(args.dataset,
                       args.out_dir,
                       args.scopes_to_evaluate,
                       args.types_of_variant,
                       args.metric,
                       args.location,
                       args.genome,
                       args.do_intronic_analysis,
                       args.clinvar_stars,
                       args.do_threshold_analysis,
                       args.do_machine_learning,
                       args.allele_frequency,
                       args.skip_heatmap,
                       args.config)


if __name__ == '__main__':
    main()
