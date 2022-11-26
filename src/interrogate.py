import logging
from collections import defaultdict
from typing import List
import numpy as np
import pandas as pd
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
                 plot_these_tools: List = None,
                 labels: str = None,
                 allele_frequency_col: str = "gnomADg_AF",
                 skip_heatmap: bool = False,
                 tools_config: str = None,
                 interrogate_mode: bool = False):

        """
        ----
        Base args described in Base class
        ----
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
                         tools_config=tools_config,
                         interrogate_mode=interrogate_mode)

        if plot_these_tools is not None:
            assert set(plot_these_tools).issubset(self.available_tools), "Invalid tool(s) name(s) in the " \
                                                                         "'--plot_these_tools' argument. " \
                                                                         "List of available tools:\n{}". \
                format(self.available_tools)

        if "F1" in self.metric and self.labels is not None:
            raise ValueError("F1-based metrics are not available in the "
                             "interrogate mode, since all variants refer to a "
                             "single class type (--label set to {})".format(self.labels))

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
        def update_hgvsc(row: pd.Series):
            try:
                return row.SYMBOL + ":" + row.HGVSc.split(":")[1]
            except IndexError:
                if not row.SYMBOL or row.SYMBOL == "":
                    return row.HGVSg
                else:
                    return row.SYMBOL + ":" + row.HGVSg
                    

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

            to_remove_cols = ['variant_id', 'type', 'subtype', 'Existing_variation', 'location', 'label', 'outcome']
            out_c = ['chr', 'pos', 'ref', 'alt', 'Consequence', 'HGVSc', 'HGVSg', 'SYMBOL'] + [col for col in _df_by_type.columns if '_prediction' in col]
            out_c_raw = [col for col in _df_by_type.columns if not '_prediction' in col and col not in to_remove_cols]

            out_df = _df_by_type[out_c]
            out_df_raw = _df_by_type[out_c_raw]
            
            # Add used threshold to the header
            tool_idx = {t[0]: i for i, t in enumerate(self.thresholds)}
            
            for i, _cols in enumerate([out_c, out_c_raw]):
                tool_rename = {}
                for c in _cols:
                    tool = c.split('_')[0]
                    if tool in tool_idx.keys():
                        direction = str(self.thresholds[tool_idx[tool]][1])
                        thresh = str(self.thresholds[tool_idx[tool]][2])
                        tool_rename[c] = tool + ' ({}{})'.format(direction, thresh)
                    else:
                        tool_rename[c] = c
                    
                if i == 0:
                    out_df = out_df.rename(columns=tool_rename) 
                    out_df.to_csv(os.path.join(outdir, "individual_predictions.tsv"), 
                                  index=False, 
                                  sep="\t")
                else:
                    out_df_raw = out_df_raw.rename(columns=tool_rename)
                    out_df_raw.to_csv(os.path.join(outdir, "individual_predictions_raw_scores.tsv"),index=False,
                                      sep="\t")
            
            _df_just_pred = pd.concat([_df_by_type['variant_id'],
                                       _df_by_type[[col for col in _df_by_type.columns if '_prediction' in col]],
                                       _df_by_type["HGVSc"],
                                       _df_by_type["HGVSg"],
                                       _df_by_type["SYMBOL"],
                                       _df_by_type["location"].to_frame()],
                                      axis=1)
                    
            _df_just_pred['HGVSc'] = _df_just_pred.apply(update_hgvsc, axis=1)
            _df_just_pred = _df_just_pred.rename(columns={col: col.split('_prediction')[0]
                                                          for col in _df_just_pred.columns})

            _df_just_pred['tools_preds'] = _df_just_pred.apply(lambda row: row[row == True].index.tolist(), axis=1)
            _df_just_pred['tools_preds'] = _df_just_pred['tools_preds'].apply(lambda x: ';'.join(x))
            any_patho = _df_just_pred[_df_just_pred.tools_preds.apply(lambda x: len(x)) > 0]
            if any_patho.shape[0] > 0:
                any_patho[['variant_id', 'HGVSc', 'HGVSg', 'tools_preds']].to_csv(os.path.join(outdir, 'variants_with_any_pathogenic_pred.tsv'), index=False, sep="\t")
            else:
                logging.info('No variant had any pathogenic prediction by any tool')
            _df_just_pred.drop(columns='tools_preds', inplace=True)
            ##############
            ### Ratios ###
            ##############
            #try:
            logging.info("Inspecting variants for which a large "
                        "fraction of tools predicts pathogenicity.")


            ratios_df = _df_just_pred.set_index('HGVSc').drop(["location", "variant_id", "HGVSg", "SYMBOL"], axis=1).copy()

            n_tools = len(list(ratios_df))
            patho_counts = ratios_df.eq(True).sum(axis=1)
            benign_counts = ratios_df.eq(False).sum(axis=1)
            unpredictable_counts = ratios_df.fillna(-1).eq(-1).sum(axis=1)
            ratios_df = pd.concat([patho_counts, benign_counts, unpredictable_counts], axis=1)
            ratios_df.columns = ['Is pathogenic', 'Is benign', 'Unpredictable']
     
            ratios_df = ratios_df / n_tools
            #ratios_df = self._fix_col_names(ratios_df)

            _top_predicted_patho = ratios_df[ratios_df['Is pathogenic'] > 0.5]
      
            _top_predicted_benign = ratios_df[ratios_df['Is benign'] > 0.5]
            if _top_predicted_patho.shape[0] > 0:
                _top_predicted_patho.to_csv(os.path.join(outdir, "top_variant_candidates.tsv"), sep="\t")

            plot_area(ratios_df, outdir)
        
            ratios_df['Unpredictable'] *= 100
    
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

        os.makedirs(os.path.join(outdir, 'out_performance'), exist_ok=True)
        outdir = os.path.join(outdir, 'out_performance')
        for _location in self.location_filters:

            outfile = os.path.join(outdir, "statistics_{}_{}.tsv").format(var_type, _location)
            statistics = defaultdict(list)

            if _location == "all":
                df = df_pred.copy()
            else:
                df = df_pred[df_pred.location == _location].copy()
                
            if df.empty:
                logging.info("WARN: Input VCF does not have any '{}' variants. "
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

