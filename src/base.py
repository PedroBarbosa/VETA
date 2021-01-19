import os
from collections import defaultdict
from typing import Union, List, Tuple

import hgvs
import pandas as pd

from src.predictions.filters import update_thresholds, subset_toolset_by_scope, \
    subset_variants_by_type, filters_location
from src.preprocessing.clinvar import get_clinvar_cached, remove_clinvar_useless
from src.preprocessing.location import *
from src.preprocessing.osutils import check_file_exists, setup_output_directory, ensure_folder_exists
from src.preprocessing.utils_tools import *
from src.preprocessing.vcf import process_vcf

TOOLS_CONFIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Base(object):
    """
    Base class to parse input args and
    create the dataframe for downstream
    analysis
    """

    def __init__(self, vcf: Union[Tuple, str],
                 out_dir: str,
                 scopes_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 genome: str = "hg19",
                 is_intronic: bool = False,
                 is_clinvar: bool = False,
                 allele_frequency_col: str = "gnomADg_AF",
                 skip_heatmap: bool = False,
                 tools_config: str = TOOLS_CONFIG):

        """
        :param Union[Tuple, str] vcf: Input VCF(s) to analyse.
        :param str out_dir: Path to the output directory.

        :param List scopes_to_predict: Restrict analysis
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

        :param str allele_frequency_col: VCF field that measures
            population frequency of variants. If not present,
            analysis that depend of the field will be ignored.
            Default: `gnomADg_AF`

        :param bool skip_heatmap: Whether drawing of heatmaps
            should be skipped. Default: `False`

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
        self.allele_frequency_col = allele_frequency_col
        self.skip_heatmap = skip_heatmap

        tools_config = self._parse_tools_config(tools_config)
        self.thresholds = update_thresholds(tools_config)
        self.thresholds = subset_toolset_by_scope(self.thresholds, scopes_to_predict, self.is_intronic)
        self.available_tools = [t[0] for t in self.thresholds]

        # TODO change global thresholds list to remove excess of fields when custom models are provided
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

        # Replace missing AF with 0 for AF plots
        if self.allele_frequency_col in self.df.columns:
            self.df[self.allele_frequency_col] = self.df[self.allele_frequency_col]. \
                apply(lambda x: x[0] if isinstance(x, list) else x).replace(r'^\s*$', "0", regex=True).fillna(0)

        # Replace missing AF if AF is used as a tool itself to evaluate performance
        if self.allele_frequency_col in [v[0][0] for v in self.tools_config.values()]:
            name = [tool for tool, v in self.tools_config.items() if v[0][0] == self.allele_frequency_col][0]
            # If tool name is different than vcf field (it it is the same, it's dealt above)
            if name != self.allele_frequency_col:
                self.df[name] = self.df[name].fillna(0)

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

        scores, fields_missing_in_vcf = process_vcf(vcf, tools_config=self.tools_config,
                                                    thresholds=self.thresholds,
                                                    is_clinvar=self.is_clinvar,
                                                    af_col=self.allele_frequency_col)

        # Remove tool from config if corresponding field is absent from at least one of the VCFs
        self.tools_config = {k: v for k, v in self.tools_config.items() if k not in fields_missing_in_vcf}
        df = pd.DataFrame.from_dict(scores, orient='index')

        # Fix col names
        df = self._fix_column_names(df)

        if self.is_clinvar:
            df = df.dropna(subset=["CLNSIG", "CLNREVSTAT"])
            df = remove_clinvar_useless(df)
            logging.info("Number of variants after removing "
                         "irrelevant CLNSIG values: {}".format(df.shape[0]))

        # Extract variants
        logging.info("Assigning variants location")
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

        else:
            df['label'] = False if is_benign else True

        booleanDictionary = {True: 'Pathogenic', False: 'Benign'}
        df["outcome"] = df["label"].map(booleanDictionary)

        if self.is_clinvar:
            # cache processed df
            print(list(df))
            df.to_csv(vcf + '.tsv', index=False)
        return df

    def _fix_column_names(self, df: pd.DataFrame):
        """
        Rename column names
        :param pd.DataFrame df: Input df
        :return pd.DataFrame: Renamed df
        """
        genome_cols = ['hg19.chr', 'hg19.pos'] if self.genome == "hg19" else ['chr', 'pos']
        new_col_names = genome_cols + ['ref', 'alt', 'id', 'type', 'subtype',
                                       'rsID', 'HGVSc', 'Gene', 'Consequence']

        for column in df:

            if isinstance(df[column].iloc[0], (tuple,)):
                # if df[column].iloc[0][0] == "gnomADg_AF":
                #     new_col_names.append("gnomAD_genomes")
                # else:
                new_col_names.append(df[column].iloc[0][0])
                df[column] = df[column].map(tuple2float)

        rename_dict = {i: j for i, j in zip(list(df), new_col_names)}
        df.rename(columns=rename_dict, inplace=True)
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

        # values: [function_name, get_max, get_absolute_value]
        available_functions = {'to_numeric': [to_numeric, None, None],
                               'top_max': [get_top_pred, True, True],
                               'top_min': [get_top_pred, False, False],
                               'top_min_abs': [get_top_pred, False, True],
                               'spliceai_like': [process_spliceai, None, None],
                               'carol_like': [process_condel_carol, None, None],
                               'kipoi_like': [process_kipoi_tools, None, None],
                               'trap': [process_trap, None, None],
                               'scap': [process_scap, None, None]
                               }

        logging.info("Engineering the scores.")
        clean_functions = {
            "GERP": available_functions['top_max'],
            "phyloP": available_functions['top_max'],
            "phastCons": available_functions['top_max'],
            "SiPhy": available_functions['to_numeric'],

            "LRT": available_functions['to_numeric'],
            "Sift": available_functions['top_min_abs'],
            "Polyphen2HVAR": available_functions['top_max'],
            "Polyphen2HDIV": available_functions['top_max'],
            "MutationAssessor": available_functions['top_max'],
            "MutationTaster": available_functions['top_max'],
            "FATHMM": available_functions['top_min'],
            "Provean": available_functions['top_min'],
            "Mutpred": available_functions['to_numeric'],
            "CAROL": available_functions['carol_like'],
            "Condel": available_functions['carol_like'],
            "M-CAP": available_functions['to_numeric'],
            "MetaLR": available_functions['to_numeric'],
            "MetaSVM": available_functions['to_numeric'],
            "REVEL": available_functions['to_numeric'],
            "VEST4": available_functions['to_numeric'],

            "fitCons": available_functions['to_numeric'],
            "LINSIGHT": available_functions['to_numeric'],
            "GWAVA": available_functions['to_numeric'],
            "CADD": available_functions['to_numeric'],
            "Eigen": available_functions['to_numeric'],
            "FATHMM-MKL": available_functions['to_numeric'],
            "FunSeq2": available_functions['top_max'],
            "DANN": available_functions['to_numeric'],
            "ReMM": available_functions['to_numeric'],

            "HAL": available_functions['kipoi_like'],
            "S-CAP": available_functions['scap'],
            "MMSplice": available_functions['kipoi_like'],
            "kipoiSplice4": available_functions['kipoi_like'],
            "kipoiSplice4_cons": available_functions['kipoi_like'],
            "TraP": available_functions['trap'],
            "SPIDEX": available_functions['to_numeric'],
            "dbscSNV": available_functions['top_max'],
            "MaxEntScan": available_functions['to_numeric'],
            "SpliceAI": available_functions['spliceai_like']
        }

        _functions_that_require_loc = ["process_spliceai", "process_trap"]

        _absent_tools = {_tool: info for _tool, info in self.tools_config.items()
                         if _tool not in clean_functions.keys()}

        # If new custom tool
        if _absent_tools:

            # get function, if provided (5th col in config)
            for _tool, info in _absent_tools.items():
                _function = info[4] if len(info) > 4 else 'to_numeric'
                clean_functions[_tool] = available_functions[_function]

        # For each tool belonging to the given scope
        for _tool, info in self.tools_config.items():

            # if it is known tool
            if _tool in clean_functions.keys():

                # get function to apply
                _f = clean_functions[_tool][0]
                # is_max
                is_max = clean_functions[_tool][1]
                # is_absolute
                is_absolute = clean_functions[_tool][2]

                # if the following tools, use custom methods
                # that need variant location
                if _f.__name__ in _functions_that_require_loc:
                    df[_tool] = df[[_tool] + ['location']].apply(_f, axis=1)

                # apply single function to
                # the target columns
                else:
                    # if max and absolute info are not
                    # required (e.g. to_numeric, and custom
                    # function to specific tools)
                    if all([_v is None for _v in [is_max, is_absolute]]):
                        df[_tool] = df[_tool].apply(_f)
                    else:
                        df[_tool] = df[_tool].apply(_f, is_max=is_max,
                                                    absolute=is_absolute)

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
        AVAILABLE_FUNCTIONS = ['to_numeric', 'top_max', 'top_min', 'top_min_abs',
                               'spliceai_like', 'kipoi_like', 'carol_like', 'trap',
                               'scap']

        _infile = open(config, 'r')
        data = defaultdict(list)


        for line in _infile:
            line = line.rstrip()
            if line.startswith("#") or not line:
                continue

            _fields = line.split("\t")
            assert len(_fields) >= 2, "First 2 fields are required in the tools config file. " \
                                      "Problematic line: {}.\nMake sure fields are tab separated.".format(line)

            assert _fields[0] not in data.keys(), "Repeated tool name found in the tools " \
                                                  "config file: {}".format(_fields[0])

            data[_fields[0]].append(_fields[1].split(","))
            if _fields[0] not in AVAILABLE_TOOLS_NAMES:
                assert len(_fields) >= 5, "\'{}\' tool is not in the list of available " \
                                          "tools. Thus, It is assumed to be a custom method " \
                                          "provided by the user, which means that a reference " \
                                          "threshold, its directionality and scope must be provided. " \
                                          "Please set that info (at least 5 columns), or check if the " \
                                          "tool name was just misspelled. List of default tools: " \
                                          "{}".format(_fields[0], '\n' + '\n'.join(AVAILABLE_TOOLS_NAMES))

            if len(_fields) > 2:
                assert _fields[2] in ['>', "<"], "Directionality of the threshold (3th col) in " \
                                                 "the tools config must be set to '>' or '<' at " \
                                                 "line:\n{}".format(line)

                try:
                    _ref_thresh = float(_fields[3])
                except ValueError:
                    raise ValueError("Reference threshold (4th col) in the tools config "
                                     "file must be a number at line:\n{}".format(line))

                data[_fields[0]].extend([_fields[2], _ref_thresh])

                if len(_fields) > 4:
                    assert _fields[4] in AVAILABLE_SCOPES, "\'{}\' is not in the list of " \
                                                           "available scopes: {} at line:\n{}".format(_fields[4],
                                                                                                      AVAILABLE_SCOPES,
                                                                                                      line)
                    data[_fields[0]].append(_fields[4])

                if len(_fields) > 5:
                    assert _fields[5] in AVAILABLE_FUNCTIONS, "\'{}\' is not in the list of " \
                                                              "available functions: {} at " \
                                                              "line:\n{}".format(_fields[5], AVAILABLE_FUNCTIONS,
                                                                                 line)
                    data[_fields[0]].append(_fields[5])

        if len(data) < 2:
            print(data)
            raise ValueError('Config file requires at least two tools to be compared.')

        return data
