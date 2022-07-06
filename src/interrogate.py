import logging
from collections import defaultdict
from typing import List
import numpy as np

from base import Base
from plots.plots_interrogate_mode import *
from predictions.apply import apply_tool_predictions
from predictions.metrics import generate_statistics
from preprocessing.osutils import check_file_exists
from preprocessing.osutils import ensure_folder_exists


class PredictionsEval(Base):
    """
    Inspect predictions of an unlabelled VCF
    """

    def __init__(self, vcf: str,
                 out_dir: str,
                 scopes_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 aggregate_classes: bool = False,
                 select_conseqs: str = "in_gene_body",
                 do_intronic_analysis: bool = False,
                 best_tools: str = None,
                 n_best_tools: int = None,
                 plot_these_tools: List = None,
                 labels: str = None,
                 allele_frequency_col: str = "gnomADg_AF",
                 skip_heatmap: bool = False,
                 tools_config: str = None):

        """
        ----
        Base args described in Base class
        ----
        :param str best_tools: Restrict analysis to the best set
            of tools obtained from a previous VETA run using a
            reference catalog (e.g. Clinvar). It must refer to the
            file `tools_ranking*.tsv` that is written when running

            the aforementioned analysis. Default: `None`, use all
            tools  available considering the `scope_to_predict`
            argument.

        :param int n_best_tools: Number of best tools selected
            from the ranking provided in the `--best_tools` argument.
            Default: `5`.

        :param List plot_these_tools: Plot scores distribution for the
            given tools. Default: `None`.

        :param str labels: If variants stored in VCF represent a list
            of labelled variants of a given type. If so, additional
            metrics will be inferred. Available types: ['benign',
            'pathogenic']. Default: `None`, variants are unlabelled.
        """

        super().__init__(vcf=vcf,
                         out_dir=out_dir,
                         scopes_to_predict=scopes_to_predict,
                         types_of_variant=types_of_variant,
                         metric=metric,
                         location=location,
                         aggregate_classes = aggregate_classes,
                         select_conseqs = select_conseqs,
                         do_intronic_analysis=do_intronic_analysis,
                         is_clinvar=False,
                         allele_frequency_col=allele_frequency_col,
                         skip_heatmap=skip_heatmap,
                         tools_config=tools_config)

        if best_tools is not None:
            check_file_exists(best_tools)

        if plot_these_tools is not None:
            assert set(plot_these_tools).issubset(self.available_tools), "Invalid tool(s) name(s) in the " \
                                                                         "'--plot_these_tools' argument. " \
                                                                         "List of available tools:\n{}". \
                format(self.available_tools)

        if "F1" in self.metric and self.labels is not None:
            raise ValueError("F1-based metrics are not available in the "
                             "interrogate mode, since all variants refer to a "
                             "single class type (--label set to {})".format(self.labels))

        self.best_tools = check_file_exists(best_tools)
        self.n_best_tools = n_best_tools
        self.plot_these_tools = plot_these_tools
        self.labels = labels

        self.inspect_predictions()

    def inspect_predictions(self):
        """
        Inspect predictions for each
        variant type and produces output
        files and plots

        :return:
        """

        df_with_pred = apply_tool_predictions(self.df, self.thresholds). \
            rename(columns={'index': 'variant_id'})

        for var_type, _func in self.variant_types:

            outdir = os.path.join(self.out_dir, var_type)
            ensure_folder_exists(outdir)
      
            _df_by_type = _func(df_with_pred).copy()

            logging.info("Looking at {} ({} variants)".format(var_type,
                                                              _df_by_type.shape[0]))

            if _df_by_type.shape[0] == 0:
                logging.warning("WARN: There are no {} in the variant set. "
                                "Skipping this analysis.".format(var_type))
                continue
            

            ###################
            ### Predictions ###
            ###################
            # df = apply_tool_predictions(df, threshold_list)
            logging.info("Extracting predictions")

            out_c = ['chr', 'pos', 'ref', 'alt', 'Consequence', 'HGVSc', 'SYMBOL'] + [col for col in _df_by_type.columns if '_prediction' in col]
            out_df = _df_by_type[out_c]
            out_df = out_df.rename(columns={col: col.split('_prediction')[0] for col in out_df.columns}).to_csv(os.path.join(outdir, "individual_predictions.tsv"), index=False, sep="\t")
            
            _df_just_pred = pd.concat([_df_by_type['variant_id'],
                                       _df_by_type[[col for col in _df_by_type.columns if '_prediction' in col]],
                                       _df_by_type["HGVSc"],
                                       _df_by_type["SYMBOL"],
                                       _df_by_type["location"].to_frame()],
                                      axis=1)
            
            _df_just_pred['HGVSc'] = _df_just_pred.apply(lambda x: x.SYMBOL + ":" + x.HGVSc.split(":")[1], axis=1)
            _df_just_pred = _df_just_pred.rename(columns={col: col.split('_prediction')[0]
                                                          for col in _df_just_pred.columns})
       
            ##############
            ### Ratios ###
            ##############
            #try:
            logging.info("Inspecting variants for which a large "
                        "fraction of tools predicts pathogenicity.")

            # TODO calculate ratios more efficiently 
            ratios_df = _df_just_pred.set_index('HGVSc').drop(["location", "variant_id", "SYMBOL"], axis=1).copy()
        
            ratios_df = ratios_df.apply(lambda x: x.value_counts(True, dropna=False),
                                        axis=1).fillna(0).sort_values([True], ascending=False)
            ratios_df = ratios_df.loc[:,~ratios_df.columns.duplicated()]
    
            ratios_df.rename(columns={False: 'Is benign',
                                        True: 'Is pathogenic',
                                        np.nan: "Unpredictable"},
                                inplace=True)
   
            ratios_df = self._fix_col_names(ratios_df)

            _top_predicted_patho = ratios_df[ratios_df['Is pathogenic'] > 0.5]
            _top_predicted_benign = ratios_df[ratios_df['Is benign'] > 0.5]
            if _top_predicted_patho.shape[0] > 0:
                _top_predicted_patho.to_csv(os.path.join(outdir, "top_variant_candidates.tsv"), sep="\t")


            plot_area(ratios_df, outdir)
            ratios_df["Unpredictable"] *= 100

            plot_heatmap(ratios_df, outdir)
            plot_heatmap(_top_predicted_patho, 
                            outdir, 
                            display_annot=True)
                
            plot_heatmap(_top_predicted_patho, 
                        outdir, 
                        display_annot=True,
                        benign_too = _top_predicted_benign)
                
            # except KeyError:
            #     logging.info("No tool has predictions for the given variant type ({}), "
            #                  "analysis is going to be skipped.".format(var_type))
            # if VCF refers to variants of a
            # given label, compute additional
            # metrics
            ###################
            ### Performance ###
            ###################
            if self.labels is not None:
                logging.info("Inspecting tools performance based on the label provided ({})".format(self.labels))
                _df_just_pred['label'] = False if self.labels in ["Benign", "Neutral"] else True
                self.generate_performance_with_label(_df_just_pred,
                                                     outdir=outdir,
                                                     var_type=var_type)


            ###############
            ### Heatmap ###
            ###############
            if self.skip_heatmap is False:
                logging.info("Generating heatmap")
                class_dict = {True: 1, False: 0, np.nan: -1}
                if self.best_tools is not None:
                    tools = pd.read_csv(self.best_tools, sep="\t")['tool'].head(self.n_best_tools).tolist()
                    tools.append("location")

                    if self.labels:
                        tools.append("label")

                    _df = _df_just_pred[tools].replace(class_dict)

                else:
                    _df = _df_just_pred.replace(class_dict)

                plot_heatmap_toptools(_df, filters=self.location_filters, outdir=outdir)

            ##################
            ### Score dist ###
            ##################
            if self.plot_these_tools is not None:
                logging.info("Generating score distribution plots")
                for tool in self.plot_these_tools:
                    # uses original df that contains
                    # scores, not just binary predictions
                    plot_tool_score_distribution(_df_by_type,
                                                 tool=tool,
                                                 thresholds=self.thresholds,
                                                 outdir=outdir)

        if self.do_intronic_analysis:
            logging.info("For now, intronic analysis is not available "
                         "in the inspect mode. Skipping it.")

    def generate_performance_with_label(self,
                                        df_pred: pd.DataFrame,
                                        outdir: str,
                                        var_type:str):

        """
        Evaluate tools performance considering that the
        intput VCF refers to a list of variants with a
        single and known label (e.g. benign or pathogenic)

        :param pd.DataFrame df_pred: Df with predictions
            for each variant
        :param str outdir: Output directory
        :param var_type: Variant types analysed 

        :return:
        """
        assert "F1" not in self.metric, "Can't use F1-based metrics in the inspect mode since it is " \
                                        "necessary to have at least two classes for its calculation." \
                                        "at least "
      
        for _location in self.location_filters:

            outfile = os.path.join(outdir, "statistics_{}_{}.tsv").format(var_type, _location)
            statistics = defaultdict(list)
            
            df = df_pred[df_pred.location == _location].copy()
            if df.empty:
                logging.info("WARN: Input VCF does not have any {} variants. "
                             "Skipping label performance analysis.".format(_location))
                continue

            for tool, *args in self.thresholds:
                try:
                    if np.sum(~df[tool].isnull()) == 0:
                        continue
                    generate_statistics(df, statistics, _location, tool, is_single_label=True)

                except KeyError:
                    continue

                stats_df = pd.DataFrame(statistics).drop(columns=['filter'])
                stats_df.sort_values([self.metric], ascending=False).to_csv(outfile, sep="\t", index=False)
                plot_accuracy(stats_df, self.metric, _location, outdir)

    def _fix_col_names(self, ratios_df: pd.DataFrame):
        """
        Fixes colnames when unpredictable variants were not found

        :param pd.DataFrame ratios_df: Input df
        """

        if "Unpredictable" not in ratios_df.columns:
            for i in range(0, len(ratios_df.columns)):
                if ratios_df.columns[i] not in ['Is benign', 'Is pathogenic']:
                    ratios_df = ratios_df.rename(columns={ratios_df.columns[i]: "Unpredictable"})

        if "Unpredictable" not in ratios_df.columns:
            ratios_df['Unpredictable'] = 0.0
            
        return ratios_df

