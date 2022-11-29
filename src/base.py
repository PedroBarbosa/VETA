import os
from collections import defaultdict
from posixpath import split
from typing import Union, List, Tuple

import hgvs
import pandas as pd
from importlib import resources
import config

from predictions.filters import update_thresholds, subset_toolset_by_scope, \
    subset_variants_by_type, _extract_possible_filters
from preprocessing.clinvar import get_clinvar_cached, remove_clinvar_useless
from preprocessing.location import *
from preprocessing.osutils import check_file_exists, setup_output_directory, ensure_folder_exists
from preprocessing.utils_tools import *
from preprocessing.vcf import process_vcf

TOOLS_CONFIG = resources.open_text(config, 'tools_config.txt')

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
                 aggregate_classes: bool = False,
                 select_conseqs: str = "gene_body",
                 do_intronic_analysis: bool = False,
                 split_splice_sites: bool = False,
                 is_clinvar: bool = False,
                 allele_frequency_col: str = "gnomADg_AF",
                 skip_heatmap: bool = False,
                 tools_config: str = TOOLS_CONFIG,
                 interrogate_mode: bool = False):

        """
        :param Union[Tuple, str] vcf: Input VCF(s) to analyse.
        :param str out_dir: Path to the output directory.

        :param List scopes_to_predict: Restrict analysis
            to a subset of tools based on their scope.
            Available options: ['Conservation', 'Whole_genome',
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

        :param bool aggregate_classes: Aggregate coding and splice region
            variant classes into more high-level concepts to be analyzed 
            together
        
        :paran str select_conseqs: How to select the top consequence
            block per variant
        
        :param bool do_intronic_analysis: Whether additional analysis of
            intronic variants extracted from HGVSc field will
            be performed
            
        :param bool split_splice_sites: Do separate analysis for donor and 
        acceptor variants when 'do_intronic_analysis' is set
            
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
            
        :param bool interrogate_mode: Whether veta run
            is in interrogate mode
            
        :return pd.DataFrame: Processed dataframe
        """
        self.metric = metric
        self.location_from = location
        self.do_intronic_analysis = do_intronic_analysis
        self.split_splice_sites = split_splice_sites
        self.aggregate_classes = aggregate_classes
        self.select_conseqs = select_conseqs
        self.location_filters = _extract_possible_filters(self.aggregate_classes)

        self.is_clinvar = is_clinvar
        self.is_interrogate_mode = interrogate_mode
        self.out_dir = setup_output_directory(out_dir)
        ensure_folder_exists(self.out_dir)
        self.allele_frequency_col = allele_frequency_col
        self.skip_heatmap = skip_heatmap

        tools_config = self._parse_tools_config(tools_config)
        self.thresholds = update_thresholds(tools_config)
        self.thresholds = subset_toolset_by_scope(self.thresholds, scopes_to_predict)
        self.available_tools = [t[0] for t in self.thresholds]

        # TODO change global thresholds list to remove excess of fields when custom models are provided
        self.tools_config = {k: v for k, v in tools_config.items() if k in self.available_tools}
        self.variant_types = subset_variants_by_type(types_of_variant)

        if self.location_from == "Consequence" and self.do_intronic_analysis:
            raise ValueError('If \'--do_intronic_analysis\' is set, \'--location\' '
                             'must be \'HGVSc\' because intronic bin are assigned '
                             'based on the distance offset extracted from HGVSc expressions.')

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

            self.df = pd.concat(_dfs).reset_index(drop=True)

        # Replace missing AF with 0 for AF plots
        if self.allele_frequency_col in self.df.columns:

            if self.df[self.allele_frequency_col].dtype == object:
                self.df[self.allele_frequency_col] = self.df[self.allele_frequency_col].\
                    apply(lambda x: [0 if v is None else x[0] for v in x][0])

        # Replace missing AF if AF is used as a tool itself to evaluate performance
        if self.allele_frequency_col in [v[0][0] for v in self.tools_config.values()]:
            name = [tool for tool, v in self.tools_config.items() if v[0][0] == self.allele_frequency_col][0]
            # If tool name is different than vcf field (it it is the same, it's dealt above)
            try:
                if name != self.allele_frequency_col:
                    self.df[name] = self.df[name].fillna(0)
            except KeyError:
                raise KeyError('Unknown column {} in the processed dataframe. '
                               'Remove *tsv file and run again Clinvar analysis '
                               'with the added tool.'.format(name))

        if self.is_clinvar:
            # cache processed df
            self.df.to_csv(vcf + '.tsv', index=False)

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

        df, fields_missing_in_vcf = process_vcf(vcf, 
                                                tools_config=self.tools_config,
                                                thresholds=self.thresholds,
                                                select_conseq = self.select_conseqs,
                                                is_clinvar=self.is_clinvar,
                                                af_col=self.allele_frequency_col)

        # Remove tool from config and thresholds if corresponding field is absent from at least one of the VCFs
        self.tools_config = {k: v for k, v in self.tools_config.items() if k not in fields_missing_in_vcf}
        self.thresholds = [x for x in self.thresholds if x[0] not in fields_missing_in_vcf]

        if self.is_clinvar:
            df = df.dropna(subset=["CLNSIG", "CLNREVSTAT"])
            df = remove_clinvar_useless(df)
            logging.info("Number of variants after removing "
                         "irrelevant CLNSIG values: {}".format(df.shape[0]))

        # Extract variants
        logging.info("Assigning variants location/type")
        if self.location_from == "HGVSc":
            hp = hgvs.parser.Parser()

            df['location'] = df.apply(get_location,
                                      hp=hp,
                                      aggregate=self.aggregate_classes,
                                      axis=1,
                                      result_type='expand')
            
            if self.do_intronic_analysis:
                df[['intron_bin', 'intron_offset', 'which_ss']] = df.apply(
                    lambda x: assign_intronic_bins(x['HGVSc'], 
                                                   hp, 
                                                   x['location'], 
                                                   self.aggregate_classes), axis=1)

        elif self.location_from == "Consequence":
            df['location'] = df['Consequence'].apply(get_location_from_consequence, 
                                                     detailed=self.detailed_conseqs)

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
                               'categorical_to_numeric': [categorical_to_numerical, None, None],
                               'carol_like': [process_condel_carol, None, None],
                               'kipoi_like': [process_kipoi_tools, None, None],
                               'spliceai': [process_spliceai, None, None],
                               'trap': [process_trap, None, None],
                               'scap': [process_scap, None, None],
                               'conspliceml': [process_conspliceml, None, None],
                               'conspliceml_like': [process_conspliceml_like, None, None],
                               'pangolin': [process_pangolin, None, None],
                               'spip': [process_spip, None, None],
                               'labranchor': [process_labranchor, None, None],
                               'absplice_dna': [process_absplice, None, None]
                               }

        logging.info("Engineering the scores.")
        clean_functions = {
            "GERP": available_functions['top_max'],
            "phyloP": available_functions['top_max'],
            "phastCons": available_functions['top_max'],
            "SiPhy": available_functions['to_numeric'],
            "CDTS": available_functions['to_numeric'],

            "LRT": available_functions['to_numeric'],
            "Sift": available_functions['top_min_abs'],
            "Polyphen2HVAR": available_functions['top_max'],
            "Polyphen2HDIV": available_functions['top_max'],
            "MutationAssessor": available_functions['top_max'],
            "MutationTaster2": available_functions['top_max'],
            "FATHMM": available_functions['top_min'],
            "Provean": available_functions['top_min'],
            "Mutpred": available_functions['to_numeric'],
            "CAROL": available_functions['carol_like'],
            "Condel": available_functions['carol_like'],
            "M-CAP": available_functions['to_numeric'],
            "MetaLR": available_functions['to_numeric'],
            "MetaSVM": available_functions['to_numeric'],
            "REVEL": available_functions['top_max'],
            "VEST4": available_functions['to_numeric'],
            "CardioBoost": available_functions['to_numeric'],
            "CardioVAI": available_functions['to_numeric'],
            "PrimateAI": available_functions['to_numeric'],
            "VARITY": available_functions['to_numeric'],
            "MutFormer": available_functions['to_numeric'],
            "MutScore": available_functions['to_numeric'],
            "MVP": available_functions['to_numeric'],
            "MTR": available_functions['to_numeric'],
            "MPC": available_functions['to_numeric'],
            "MISTIC": available_functions['top_max'],
            "ClinPred": available_functions['top_max'],
            "cVEP": available_functions['categorical_to_numeric'],
            "EVE_Class20": available_functions['categorical_to_numeric'],
            "EVE_Class50": available_functions['categorical_to_numeric'],
            "EVE_Class90": available_functions['categorical_to_numeric'],
            "EVE": available_functions['to_numeric'],
            
            "fitCons": available_functions['to_numeric'],
            "LINSIGHT": available_functions['to_numeric'],
            "GWAVA": available_functions['to_numeric'],
            "CADD_v1.5": available_functions['to_numeric'],
            "CADD-Splice": available_functions['to_numeric'],
            "Eigen": available_functions['to_numeric'],
            "FATHMM-MKL": available_functions['to_numeric'],
            "FunSeq2": available_functions['top_max'],
            "DANN": available_functions['to_numeric'],
            "ReMM": available_functions['to_numeric'],
            "CAPICE": available_functions['top_max'],

            "HAL": available_functions['kipoi_like'],
            "S-CAP": available_functions['scap'],
            "MMSplice": available_functions['kipoi_like'],
            "kipoiSplice4": available_functions['kipoi_like'],
            "kipoiSplice4_cons": available_functions['kipoi_like'],
            "TraP": available_functions['trap'],
            "SPIDEX": available_functions['to_numeric'],
            "dbscSNV": available_functions['top_max'],
            "MaxEntScan": available_functions['top_max'],
            "SpliceAI": available_functions['spliceai'],
            "SQUIRLS": available_functions['to_numeric'],
            "ConSpliceML": available_functions['conspliceml'],
            "IntSplice2": available_functions['to_numeric'],
            "CI-SpliceAI": available_functions['spliceai'],
            "Pangolin": available_functions['pangolin'],
            "SPiP": available_functions['spip'],
            "MLCsplice": available_functions['top_max'],
            "AbSplice-DNA": available_functions['absplice_dna'],
            "LaBranchoR": available_functions['labranchor']
        }

        _functions_that_require_loc = ["process_trap"]
        _functions_that_require_symbol = ["process_spliceai", "process_conspliceml"]
        _absent_tools = {_tool: info for _tool, info in self.tools_config.items()
                         if _tool not in clean_functions.keys()}

        # If new custom tool
        if _absent_tools:
 
            # get function, if provided (5th col in config)
            for _tool, info in _absent_tools.items():
                if len(info) > 4:    
                    _function = info[4]
                    #self.tools_config[_tool] = info[:-1]
                else:
                    logging.info('Processing function was not provided for {}. Using the default to_numeric function. Be careful.'.format(_tool))
                    _function = 'to_numeric'
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

                elif _f.__name__ in _functions_that_require_symbol:
                    
                    # If CI-SpliceAI
                    if _f.__name__ == "process_spliceai" and _tool != "SpliceAI":
                        df[_tool] = df[[_tool] + ['SYMBOL']].apply(_f, check_gene_name=False, axis=1)
                    else:
                        df[_tool] = df[[_tool] + ['SYMBOL']].apply(_f, axis=1)
                    
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

    def _parse_tools_config(self, config):
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
        AVAILABLE_TOOLS_NAMES = ["GERP", "phyloP", "phastCons", "SiPhy",
                                 "CDTS", "LRT", "Sift", "Polyphen2HVAR",
                                 "Polyphen2HDIV", "MutationTaster2",
                                 "MutationAssessor", "FATHMM", "Provean",
                                 "Mutpred", "CAROL", "Condel", "M-CAP",
                                 "MetaSVM", "MetaLR", "REVEL", "VEST4",
                                 "MVP", "PrimateAI", "cVEP", "VARITY",
                                 "MTR", "MPC", "MutScore", "MutFormer",
                                 "ClinPred", "MISTIC",
                                 "EVE", "EVE_class20", "EVE_class50","EVE_class90",                             
                                 "CardioBoost", "CardioVAI",
                                 "fitCons", "LINSIGHT", "ReMM", "GWAVA", "FATHMM-MKL",
                                 "Eigen", "FunSeq2", "CADD_v1.5", "CADD-Splice", "DANN", "CAPICE",
                                 "HAL", "SPIDEX", "dbscSNV", "MaxEntScan", "SpliceAI", "S-CAP", "ConSpliceML",
                                 "TraP", "MMSplice", "SQUIRLS", "IntSplice2", "kipoiSplice4", "kipoiSplice4_cons",
                                 "CI-SpliceAI", "Pangolin", "SPiP", "LaBranchoR", "MLCsplice", "AbSplice-DNA"]

        AVAILABLE_SCOPES = ["Protein", "Conservation", "Whole_genome", "Splicing"]
        AVAILABLE_FUNCTIONS = ['to_numeric', 'top_max', 'top_min', 'top_min_abs',
                               'categorical_to_numeric', 'kipoi_like', 'carol_like', 
                               'trap', 'scap', 'conspliceml', 'conspliceml_like',
                               'spliceai', 'ci_spliceai',
                               'pangolin', 'spip', 'labranchor']

        data = defaultdict(list)
        if isinstance(config, str):
            _c = open(config, 'r')
            config = _c

        for line in config:
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
    
        if len(data) < 2 and self.is_interrogate_mode is False:
            raise ValueError('Config file requires at least two tools to be compared.')

        return data
