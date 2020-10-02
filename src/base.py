from src.preprocessing.osutils import check_file_exists, setup_output_directory, ensure_folder_exists
from typing import Union, List, Tuple
from src.predictions.filters import update_thresholds, subset_toolset_by_scope,\
    subset_variants_by_type, filters_location
from src.preprocessing.vcf import process_vcf
from src.preprocessing.location import *
import hgvs
from src.preprocessing.utils import *

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
import pandas as pd
from src.preprocessing.clinvar import get_clinvar_cached, remove_clinvar_useless
from collections import defaultdict
import os
TOOLS_CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class Base(object):
    """
    Base class to parse input args and
    create the dataframe for downstream
    analysis
    """

    def __init__(self, vcf: Union[Tuple, str],
                 out_dir: str,
                 scope_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 genome: str = "hg19",
                 is_intronic: bool = False,
                 is_clinvar: bool = False,
                 tools_config: str = TOOLS_CONFIG):

        """
        :param Union[Tuple, str] vcf: Input VCF(s) to analyse.
        :param str out_dir: Path to the output directory.

        :param List scope_to_predict: Restrict analysis
            to a subset of tools based on their scope.
            Available options: ['Conservation', 'Functional',
            'Protein', 'Splicing']. Default: `None` tools
            from all scopes are used.

        :param List types_of_variant: Restrict analysis to the
            given variant types. Available options: ['snps',
            'indels', 'insertions', 'deletions', 'mnps'].
            Default: `None`, Performance is measured for all
            types together and for each subtype separately.

        :param str metric: Metric to rank the tools.
            Available options: ['weighted_accuracy', 'accuracy',
            'F1', 'weighted_F1', 'coverage']. Default:
            `weighted_accuracy`.

        :param str location: VEP annotated VCF field to
            extract location of the variants. Available
            options: ['HGVSc, 'Consequence']. Default: "HGVSc".

        :param str genome: Genome build where variants are mapped
            Available options: ['hg19', 'hg38']. Default: `hg19`.

        :param bool is_intronic: Whether additional analysis of
            intronic variants extracted from HGVSc field will
            be performed

        :param bool is_clinvar: Whether `vcf` is from Clinvar.
            Default: `False`

        :param str tools_config: Path to the tools config
            where available tools in VCF are mapped to the
            corresponding VCF annotation

        :return pd.DataFrame: Processed dataframe
        """
        self.metric = metric
        self.location_from = location
        self.genome = genome
        self.is_intronic = is_intronic
        self.location_filters = filters_location
        self.is_clinvar = is_clinvar
        self.out_dir = setup_output_directory(out_dir)
        ensure_folder_exists(self.out_dir)

        tools_config = self._parse_tools_config(tools_config)
        self.thresholds = update_thresholds(tools_config)
        self.thresholds = subset_toolset_by_scope(self.thresholds, scope_to_predict, self.is_intronic)
        self.available_tools = [t[0] for t in self.thresholds]
        self.tools_config = {k: v for k, v in tools_config.items() if k in self.available_tools}
        self.variant_types = subset_variants_by_type(types_of_variant)

        if self.location_from == "Consequence" and self.is_intronic:
            raise ValueError('If \'--is_intronic\' is set, \'--location\' '
                             'must be \'HGVSc\' because intronic bin analysis '
                             'will be performed based on the HGVS nomenclature.')

        if isinstance(vcf, str):
            check_file_exists(vcf)
            # if it's clinvar, additional
            # processing will be done
            self.df = self.get_df_ready(vcf)

        else:
            _dfs = []
            for f in vcf:
                is_benign = False if "pathogenic" in f or "deleterious" in f else True
                _dfs.append(self.get_df_ready(f, is_benign=is_benign))

            self.df = pd.concat(_dfs)

    def get_df_ready(self, vcf: str, is_benign: bool = False):
        """
        Process input VCF(s) so that df are ready for
        downstream analysis

        :param str vcf: Input VCF file
        :param bool is_benign: Whether `vcf` refers to a set of
            benign variants. Only useful when `vcf` is not from
            Clinvar. Default: `False`
        :return pd.DataFrame: Cleaned dataframe
        """
        logging.info("Processing {} file".format(os.path.basename(vcf)))

        if self.is_clinvar:
            _df = get_clinvar_cached(vcf)
            # If cached file was found, DF is ready
            if _df is not None:
                return _df

        scores = process_vcf(vcf, tools_config=self.tools_config,
                             thresholds=self.thresholds,
                             is_clinvar=self.is_clinvar)
        df = pd.DataFrame.from_dict(scores, orient='index')

        # Fix col names
        df = self._fix_column_names(df)

        if self.is_clinvar:
            df = df.dropna(subset=["CLNSIG", "CLNREVSTAT"])
            df = remove_clinvar_useless(df)
            logging.info("Number of variants after removing "
                         "irrelevant CLNSIG values: {}".format(df.shape[0]))

        # Extract variants location
        if self.location_from == "HGVSc":
            hp = hgvs.parser.Parser()
            df['location'] = df['HGVSc'].apply(get_location, hp=hp)
            if self.is_intronic:
                df[['intron_bin', 'intron_offset']] = df.apply(
                    lambda x: assign_intronic_bins(x['HGVSc'], hp, x['location']), axis=1)

        elif self.location_from == "Consequence":
            df['location'] = df['Consequence'].apply(get_location_from_consequence)

        # Clean predictions scores
        df = self._clean_scores(df.reset_index())

        # Set labels
        if self.is_clinvar:
            df['label'] = (df['CLNSIG'].str.contains('Pathogenic')) | \
                          (df['CLNSIG'].str.contains('Likely_pathogenic'))

            # cache processed df
            df.to_csv(vcf + '.tsv', index=False)

        else:
            df['label'] = False if is_benign else True

        booleanDictionary = {True: 'Pathogenic', False: 'Benign'}
        df["outcome"] = df["label"].map(booleanDictionary)
        return df

    def _fix_column_names(self, df: pd.DataFrame):
        """
        Rename column names
        :param pd.DataFrame df: Input df
        :return pd.DataFrame: Renamed df
        """
        genome_cols = ['hg19.chr', 'hg19.pos'] if self.genome == "hg19" else ['chr', 'pos']
        new_col_names = genome_cols + ['ref', 'alt', 'id', 'type', 'subtype',
                                       'rsID', 'HGVSc', 'Gene', 'Consequence', 'gnomAD_exomes']

        for column in df:

            if isinstance(df[column].iloc[0], (tuple,)):
                if df[column].iloc[0][0] == "gnomADg_AF":
                    new_col_names.append("gnomAD_genomes")
                else:
                    new_col_names.append(df[column].iloc[0][0])
                df[column] = df[column].map(tuple2float)

        rename_dict = {i: j for i, j in zip(list(df), new_col_names)}

        df.rename(columns=rename_dict, inplace=True)
        df['gnomAD_genomes'] = df['gnomAD_genomes'].replace(r'^\s*$', "0", regex=True)
        df['gnomAD_exomes'] = df['gnomAD_exomes'].replace(r'^\s*$', "0", regex=True)
        return df

    def _clean_scores(self, df: pd.DataFrame):
        """
        Transform tools scores so that proper
        evaluations can be performed.

        By default, VETA has some methods
        ready to process the scores of specific
        tools.

        :param pd.DataFrame df: Input dataframe
        with information about each variant.

        :return pd.DataFrame: Processed df.
        """
        # pd.options.display.float_format = '{:,.6f}'.format

        logging.info("Engineering the scores.")
        clean_functions = {
            "GERP": get_top_pred,
            "phyloP": get_top_pred,
            "phastCons": get_top_pred,
            "SiPhy": to_numeric,

            "LRT": to_numeric,
            "Sift": get_top_pred,
            "Polyphen2HVAR": get_top_pred,
            "Polyphen2HDIV": get_top_pred,
            "MutationAssessor": get_top_pred,
            "MutationTaster": get_top_pred,
            "FATHMM": get_top_pred,
            "Provean": get_top_pred,
            "Mutpred": to_numeric,
            "CAROL": process_condel_carol,
            "Condel": process_condel_carol,
            "M-CAP": to_numeric,
            "MetaLR": to_numeric,
            "MetaSVM": to_numeric,
            "REVEL": to_numeric,
            "VEST4": to_numeric,

            "fitCons": to_numeric,
            "LINSIGHT": to_numeric,
            "GWAVA": to_numeric,
            "CADD": to_numeric,
            "Eigen": to_numeric,
            "FATHMM-MKL": to_numeric,
            "FunSeq2": get_top_pred,
            "DANN": to_numeric,
            "ReMM": to_numeric,

            "HAL": process_kipoi_tools,
            "S-CAP": process_scap,
            "MMSplice": process_kipoi_tools,
            "kipoiSplice4": process_kipoi_tools,
            "kipoiSplice4_cons": process_kipoi_tools,
            "TraP": process_trap,
            "SPIDEX": to_numeric,
            "dbscSNV": get_top_pred,
            "MaxEntScan": to_numeric,
            "SpliceAI": process_spliceai
        }

        tools_to_get_min = ['Sift', 'FATHMM', 'Provean']
        tools_to_keep_negatives = ['FATHMM', 'Provean']
        tools_that_require_loc = ["SpliceAI", "TraP"]

        # For each tool belonging to the given scope
        for _tool, info in self.tools_config.items():

            # if it is known tool
            if _tool in clean_functions.keys():

                # get function to apply
                _f = clean_functions[_tool]

                # if the following tools, use custom methods
                # that need variant location
                if _tool in tools_that_require_loc:

                    df[_tool] = df[[_tool] + ['location']].apply(_f, axis=1)

                # apply single function to
                # the target columns
                else:
                    if _tool in tools_to_get_min:
                        absolute = False if _tool in tools_to_keep_negatives else True
                        df[_tool] = df[_tool].apply(_f, is_max=False, absolute=absolute)
                    else:
                        df[_tool] = df[_tool].apply(_f)
        return df

    def _parse_tools_config(self, config: str):
        """
        Parse tools config file to ensure
        correctness for downstream analysis.
        - Empty and commented lines (^#) are ignored.
        - First two fields are required
        - If tool name (1st col) not in the list
        of available tools, it expects to be a new
        custom method so that, 3rd (scope) and 4th
        and 5th columns are required.

        :param str config: Input file
        :return defaultdict: Dict with
        info about each tool
        """
        AVAILABLE_TOOLS_NAMES = ["GERP", "phyloP", "phastCons", "SiPhy", "LRT",
                                 "Sift", "Polyphen2HVAR", "Polyphen2HDIV",
                                 "MutationTaster", "MutationAssessor", "FATHMM",
                                 "Provean", "Mutpred", "CAROL", "Condel",
                                 "M-CAP", "MetaSVM", "MetaLR", "REVEL", "VEST4",
                                 "fitCons", "LINSIGHT", "ReMM", "GWAVA", "FATHMM-MKL",
                                 "Eigen", "FunSeq2", "CADD", "DANN", "HAL", "SPIDEX",
                                 "dbscSNV", "MaxEntScan", "SpliceAI", "S-CAP", "TraP",
                                 "MMSplice", "kipoiSplice4", "kipoiSplice4_cons"]

        AVAILABLE_SCOPES = ["Protein", "Conservation", "Functional", "Splicing"]

        _infile = open(config, 'r')
        data = defaultdict(list)

        for line in _infile:
            line = line.rstrip()
            if line.startswith("#") or not line:
                continue

            _fields = line.split("\t")
            assert len(_fields) >= 2, "First 2 fields are required in the tools config file. " \
                                      "Problematic line: {}".format(line)

            assert _fields[0] not in data.keys(), "Repeated tool name found in the tools " \
                                                  "config file: {}".format(_fields[0])

            data[_fields[0]].append(_fields[1].split(","))
            if _fields[0] not in AVAILABLE_TOOLS_NAMES:
                assert len(_fields) == 5, "\'{}\' tool is not in the list of available " \
                                          "tools. Thus, It is assumed to be a custom method " \
                                          "provided by the user, which means that a reference " \
                                          "threshold (and its directionality) must be provided. " \
                                          "Please set that info (5 columns), or check if the tool name" \
                                          " was just misspelled. List of default tools: " \
                                          "{}".format(_fields[0], '\n' + '\n'.join(AVAILABLE_TOOLS_NAMES))

            if len(_fields) > 2:
                assert _fields[2] in ['>', "<"], "Directionality of the threshold (3th col) in " \
                                                 "the tools config must be set to '>' or '<' at " \
                                                 "line:\n{}".format(line)

                try:
                    _ref_thresh = float(_fields[3])
                except ValueError:
                    raise ValueError("Reference threshold  (4th col) in the tools config "
                                     "file must be a number at line:\n{}".format(line))

                data[_fields[0]].extend([_fields[2], _ref_thresh])

                if len(_fields) > 4:
                    assert _fields[4] in AVAILABLE_SCOPES, "\'{}\' is not in the list of " \
                                                           "available scopes: {} at line:\n{}".format(_fields[4],
                                                                                                      AVAILABLE_SCOPES,
                                                                                                      line)
                    data[_fields[0]].append(_fields[4])

        return data
