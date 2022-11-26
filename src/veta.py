import argparse
import sys
import logging
from importlib import resources
import config
from interrogate import PredictionsEval
from benchmark import BenchmarkTools
from preprocessing.osutils import print_clinvar_levels
import importlib.metadata

__version__ = importlib.metadata.version("veta")
hgvs_logger = logging.getLogger('hgvs')
hgvs_logger.setLevel(logging.CRITICAL)
scikit_logger = logging.getLogger('sklearn')
scikit_logger.setLevel(logging.CRITICAL)
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')

CONFIG_PATH = resources.open_text(config, 'tools_config.txt')

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
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parent_parser.add_argument('-o', '--out_dir', metavar='', help='Path to store all the output results. '
                                                                   'Default: "out_VETA"')

    parent_parser.add_argument('-s', '--scopes_to_evaluate', metavar='',
                               help='Restrict analysis to a subset of tools based on their '
                               'scope. Available options: {%(choices)s}. Default: Tools from '
                               'all scores are used.',
                               nargs='+', choices=('Conservation', 'Whole_genome', 'Protein', 'Splicing'))

    parent_parser.add_argument('-t', '--types_of_variant', metavar='',
                               help='Restrict analysis to the given variant types. '
                               'Available options: {%(choices)s}. Default: Performance '
                               'is measured for all variants and each subtype.',
                               nargs='+', choices=('all_types', 'snps', 'indels', 'insertions', 'deletions', 'mnps'))

    parent_parser.add_argument('-m', '--metric', metavar='',
                               help='Metric to rank the tools. Available options: '
                               '{%(choices)s}. Default: "weighted_accuracy."',
                               choices=('weighted_accuracy', 'accuracy',
                                        'F1', 'weighted_F1', 'coverage',
                                         'norm_mcc', 'weighted_norm_mcc'),

                               default='weighted_accuracy')

    parent_parser.add_argument('-l', '--location', metavar='', default="HGVSc", choices=("HGVSc", "Consequence"),
                               help='VCF field to extract location of the variant. Available options: '
                               '{%(choices)s}. Default: "HGVSc".')

    parent_parser.add_argument('-a', '--aggregate_classes', help='Aggregate specific variant classes into higher level concepts '
                               '(e.g. missense, synonymous variants will be evaluated as coding variants).',
                               action='store_true')
    
    parent_parser.add_argument('-v', '--top_vep_consequence', metavar='', default='gene_body', choices=("first", "gene_body", "smallest_offset"), 
                               help='How to select the top VEP consequence for each variant. Available options: {%(choices)s}. '
                               'Default: "gene_body": First consequence occurring in the body of a gene is selected. '
                               'If "first" is set, first consequence is selected. '
                               'If "smallest_offset" is set, consequence within body genes with the smallest offset is selected (e.g. exonic variants first, then intronic variants closest to a splice site).')
    
    parent_parser.add_argument('-i', '--do_intronic_analysis', help='Perform additional analysis of '
                               'intronic variants extracted from '
                               'HGVSc field in a bin-based faction. '
                               'Bins are attributed based on the distance '
                               'of the variant to the nearest splice junction',
                               action='store_true')

    parent_parser.add_argument('--split_splice_sites', help='When "--do_intronic_analysis" is set, ' 
                               'perform separate analysis for donor and acceptor associated variants.',
                               action='store_true')
        
    parent_parser.add_argument('-af', '--allele_frequency', metavar='', default="gnomADg_AF",
                               help='VCF field (within VEP annotations, or INFO field) '
                               'that measures frequency of the variant in a '
                               'population. If it exists in the input data, additional '
                               'plots will be drawn. If absent, VETA will simply ignore it. '
                               'Default: "gnomADg_AF". Note: Missing data will be treated as '
                               'if the variant is absent in population (converted to 0).')

    parent_parser.add_argument('--skip_heatmap', action="store_true",
                               help="Skip heatmap generation for performance analysis. "
                               "If the dataset is large, this step takes quite a while.")

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

    benchmark_parser.add_argument('-c', '--clinvar_stars', metavar='', default='1s_l',
                                  help='Level of filtering when dataset refers to the clinvar database. '
                                       'By default, a high confidence clinvar subset (3 stars with likely '
                                       'annotations) is used for performance evaluation and reference '
                                       'threshold analysis (if --do_threshold_analysis is True). '
                                       'Default: "1s_l". All the possible filtering levels are visible '
                                       'with the argument \'--listClinvarLevels\'.')

    benchmark_parser.add_argument('--omim_ids', nargs='+', help='When dataset refers to the clinvar database, '
                                  'selects variants belonging to the provided OMIM ids.')
    
    benchmark_parser.add_argument('--medgen_ids', nargs='+', help='When dataset refers to the clinvar database, '
                                  'selects variants belonging to the provided MedGen ids.')
    
    benchmark_parser.add_argument('--mondo_ids', nargs='+', help='When dataset refers to the clinvar database, '
                                  'selects variants belonging to the provided MONDO ids.')
        
    benchmark_parser.add_argument('--do_threshold_analysis', action="store_true",
                                  help="Enable reference thresholds analysis for the input dataset."
                                       "Default: False")

    benchmark_parser.add_argument('--bootstrapping', action="store_true",
                                  help="Enable bootstrapping analysis when '--do_threshold_analysis' is set."
                                       "Default: False")
    
    benchmark_parser.add_argument('--do_machine_learning', action="store_true",
                                  help="Enable machine learning analysis based on the "
                                       "tools scores to inspect the best predictors and "
                                       "create an ensemble classifier that combines multiple"
                                       " scores. Default: False.")

    predict_parser = subparsers.add_parser("interrogate", help='Evaluate predictions from unlabelled variants.',
                                           parents=[parent_parser])

    predict_parser.add_argument(
        dest='vcf', help="VCF file to evaluate tools predictions.")

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

    clinvar_stars = ['1s', '2s', '3s', '4s', '1s_l',
                     '2s_l', '3s_l', '0s', '0s_l']
    ##############################
    ## Argparse args processing ##
    ##############################
    if args.command == "interrogate":

        assert args.metric not in ['F1', 'weighted_F1', 'mcc', 'norm_mcc', 'weighted_norm_mcc'], "Metric provided is not valid for interrogate mode."
     
        PredictionsEval(args.vcf,
                        args.out_dir,
                        args.scopes_to_evaluate,
                        args.types_of_variant,
                        args.metric,
                        args.location,
                        args.aggregate_classes,
                        args.top_vep_consequence,
                        args.do_intronic_analysis,
                        args.plot_these_tools,
                        args.labels,
                        args.allele_frequency,
                        args.skip_heatmap,
                        args.config,
                        interrogate_mode=True
                        )

    elif args.command == "benchmark":
        if args.clinvar_stars not in clinvar_stars:
            logging.info("Error. Set a valid value for the '--clinvar_stars' argument. Possible options "
                         "(value with a simple description):")
            print_clinvar_levels()

        assert args.clinvar_stars in clinvar_stars, "Set a valid value for the '--clinvar_stars' argument."
        phenotype_ids = [args.omim_ids, args.medgen_ids, args.mondo_ids]
        BenchmarkTools(args.dataset,
                       args.out_dir,
                       args.scopes_to_evaluate,
                       args.types_of_variant,
                       args.metric,
                       args.location,
                       args.aggregate_classes,
                       args.top_vep_consequence,
                       args.do_intronic_analysis,
                       args.split_splice_sites,
                       args.clinvar_stars,
                       phenotype_ids,
                       args.do_threshold_analysis,
                       args.bootstrapping,
                       args.do_machine_learning,
                       args.allele_frequency,
                       args.skip_heatmap,
                       args.config)


if __name__ == '__main__':
    main()
