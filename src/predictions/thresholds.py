import logging
import os
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
from collections import OrderedDict, defaultdict
from tqdm import tqdm
from sklearn.utils import resample
from plots.plots_threshold_analysis import plot_optimal_thresholds, plot_bootstrap_thresholds
from preprocessing.osutils import ensure_folder_exists
from preprocessing.utils_tools import ratio


def fbeta_at_thr(df, tool, beta_values, stats):
    """
    Calculate stats at a single threshold
    """
    tp = np.sum(df['label'].eq(df[tool + '_prediction']) & df['label'])
    tn = np.sum(df['label'].eq(df[tool + '_prediction']) & ~df['label'])
    fp = np.sum(df['label'].ne(df[tool + '_prediction']) & ~df['label'])
    fn = np.sum(df['label'].ne(df[tool + '_prediction']) & df['label'])
    total = df.shape[0]
    correct = tp + tn
    specificity = ratio(tn, tn + fp)
    precision = ratio(tp, tp + fp)
    recall = ratio(tp, tp + fn)
    accuracy = ratio(correct, total)
        
    for _b in beta_values:
        try:
            stats['Fbeta_' + str(_b)].append(round((1.0 + _b ** 2) * (precision * recall) /
                                                        ((precision * _b ** 2) + recall), 3))
        except ZeroDivisionError:
            stats['Fbeta_' + str(_b)].append(0)

    stats['accuracy'].append(accuracy)
    stats['specificity'].append(specificity)
    stats['sensitivity'].append(recall)
    return stats
                    
def classify_at_threshold(df, thr_range, **kwargs):
    """
    Iterate over each putative threshold, classify, and generate metrics
    """
    stats = defaultdict(list)
    for threshold in thr_range:
        if kwargs['direction'] == ">":
            def classification_f(
                x): return x == np.nan and np.nan or x >= threshold
        else:
            def classification_f(
                x): return x == np.nan and np.nan or x <= threshold

        df[kwargs['tool'] + '_prediction'] = df[kwargs['tool']].map(classification_f)
        
        stats = fbeta_at_thr(df, kwargs['tool'], kwargs['beta_values'], stats)
    
    # If evaluating with best thr after bootstrapping
    if len(thr_range) == 1:
        return stats
    
    new_t = {}
    f_betas = defaultdict(list)
    for _beta in kwargs['beta_values']:
        _res_at_each_t = stats['Fbeta_' + str(_beta)]
        f_betas[_beta] = _res_at_each_t
        new_t[_beta] = thr_range[_res_at_each_t.index(max(_res_at_each_t))]
        kwargs['final_thresh'][kwargs['tool']].append((round(float(new_t[_beta]), 3), max(_res_at_each_t)))
        
    return new_t, stats, f_betas

def _calculate_ci(bootst_thr: pd.DataFrame, **kwargs):
    """
    Calculates best threshold by averaging the
    list of thresholds obtained at each bootstrap iteration.
    Additionally, calculates 95% CI and returns Fbeta at best
    threshold
    
    :param pd.DataFrame bootst_thr: Df with list of bootstrap threshold
    for each beta value
    """
    # Info from threshold analysis of the data sample (before bootstrapping)
    ref_thr_info = kwargs['final_thresh'][kwargs['tool']][0]
    adj_thr_info = kwargs['final_thresh'][kwargs['tool']][1:]

    bootst_thr = bootst_thr.iloc[0]
    out = [ref_thr_info]
    for i, b in enumerate([0.5, 1, 1.5]):
  
        # Deriving best thresh by averaging bootstrap statistic values
        # best_thr = round(np.mean(bootst_thr['thr_f{}'.format(b)]), 3)
        #stats = classify_at_threshold(kwargs['df_per_tool'], [best_thr], **kwargs)   
        #fbeta_at_best_thr = stats['Fbeta_{}'.format(b)][0]
        
        dist_thr = bootst_thr['thr_f{}'.format(str(b))]

        # Empirical bootstrap confidence interval
        #ci95_up = np.mean(dist_thr) + 1.96 * np.std(dist_thr) / np.sqrt(len(dist_thr))
        #ci95_down = np.mean(dist_thr) - 1.96 * np.std(dist_thr) / np.sqrt(len(dist_thr))
        
        # Boostrap percentile method
        ci = pd.Series(dist_thr).quantile([0.025, 0.975])
        ci95_down = round(ci[0.025], 3)
        ci95_up = round(ci[0.975], 3)
        
        # Add percentiles to output
        best_thr = adj_thr_info[i][0]
            
        adj_with_percent = ('{} ({},{})'.format(best_thr, str(ci95_down), str(ci95_up)),
                           '{}'.format(adj_thr_info[i][1]))

        out.append(adj_with_percent)
        plot_bootstrap_thresholds(dist_thr, kwargs['tool'], b, kwargs['_loc'], best_thr, kwargs['ref_thresh'], ci95_down, ci95_up, kwargs['outdir'])

    kwargs['final_thresh'][kwargs['tool']] = out

def _do_bootstapping(_iter: int, **kwargs):

    df_per_tool = kwargs['df_per_tool']
    tool = kwargs['tool']
    n_of_points = kwargs['n_of_points']

    max_thr = df_per_tool[tool].max()
    min_thr = df_per_tool[tool].min()

    if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
        logging.info("Something strange in max/min thresholds {} {} "
                     "{}".format(tool, max_thr, min_thr))
        return

    data = resample(df_per_tool, 
                    replace=True, 
                    n_samples=df_per_tool.shape[0], 
                    stratify=df_per_tool.label)
    
    #n_of_points = random.randint(n_of_points, n_of_points + 500)
    #step = (max_thr - min_thr) / float(n_of_points)
    #thr_range = np.arange(min_thr, max_thr, step)

    thr_range = np.sort(np.random.uniform(min_thr, max_thr, size=n_of_points))
    new_t, _, _ = classify_at_threshold(data, thr_range, **kwargs)
    return new_t


def perform_threshold_analysis(dataset: pd.DataFrame,
                               location_filters: list,
                               threshold_list: list,
                               tools_config: dict,
                               outdir: str,
                               f1_at_ref_thresh: dict,
                               n_of_points: int = 100,
                               beta_values: list = [0.5, 1, 1.5],
                               do_bootstrapping: bool = False):
    """
    Do threshold analysis using Clinvar as ground
    truth dataset for threshold evaluation.

    Additionally, several precision/recall trade-offs
    are tested to find the most appropriate threshold
    per tool that maximimizes the number of true positives
    detected.

    :param pd.DataFrame dataset: Input dataset
    :param list location_filters: List of location filters
    :param list threshold_list: List of thresholds
    :param dict tools_config: Dict with config for each tool
    :param str outdir: Output directory
    :param dict f1_at_ref_thresh: Dict with F1 scores
    using the reference thresholds
    :param n_of_points: Default: `100`
    :param list beta_values: List of values to use as beta
        in the F_beta score calculation (precision/recall
        tradeoff). Default: `[0.5, 1, 1.5]`
    :param bool do_bootrstrapping: Set this to True to 
    generate confidence intervals
    :return dict: New thresholds for all the tools at all locations
    """
    logging.info("---------------------------")
    logging.info("Starting threshold analysis")
    logging.info("---------------------------")
    outdir = os.path.join(outdir, "thresholds_analysis")
    ensure_folder_exists(outdir)

    available_for_threshold_analysis = [
        'all', 'coding', 'missense', 'splice_site', 'splice_region', 'intronic', 'deep_intronic']

    # Output dict initialization
    thresholds_to_return = {}
    for beta in beta_values:
        thresholds_to_return[beta] = {}
        for _loc in location_filters:
            if _loc not in available_for_threshold_analysis:
                continue
            thresholds_to_return[beta][_loc] = {}

    for _loc in location_filters:

        if _loc not in available_for_threshold_analysis:
            continue
        
        elif _loc == 'all':
            df_f = dataset.copy()
        
        else:
            df_f = dataset[dataset.location == _loc].copy()
            
        print()
        logging.info("-------------------------")
        logging.info("Looking at {} variants (N={}, {} pos and {} neg)".format(_loc,
                                                                               df_f.shape[0],
                                                                               df_f[df_f.label].shape[0],
                                                                               df_f[~df_f.label].shape[0]))
        logging.info("-------------------------")
        # Do threshold analysis only if dataset has a minimum size
        if df_f.shape[0] < 100:
            logging.info('No minimum number of variants required (N=100)')
            continue

        if not all(i >= 50 for i in df_f.label.value_counts().tolist()):
            logging.warning(
                'No minimum number of variants in the minority class required (N=50)')

            continue

        final_thresholds = OrderedDict()
        for tool, direction, reference_threshold, *args in threshold_list:

            # Threshold analysis is not performed for
            # S-CAP (it uses multiple reference thresholds) and

            # cVEP (predictions are discrete categories that VETA converts to fixed floats)
            if tool not in df_f.columns or tool.lower() in ['scap', 's-cap', 'cvep']:
                continue

            # Do threshold analysis for tools with sufficient prediction power
            df_per_tool = df_f.loc[pd.notnull(df_f[tool]), ][[
                tool, 'label']].copy()

            # Transform TraP scores for proper threshold analysis
            if tool == "TraP":
                df_per_tool['TraP'] = df_per_tool['TraP'].apply(
                    lambda x: x - 1 if x > 1 else x)

                # If all variants, it is not possible to distinguish coding/non-coding thresholds
                if _loc in ["coding", "missense"]:
                    reference_threshold = 0.416
                else:
                    reference_threshold = 0.289

            ratio_predicted = df_per_tool.shape[0] / df_f.shape[0]

            if df_per_tool.shape[0] < 100:
                logging.info(
                    '{} did not predict the minimum number of variants required (N=100)'.format(tool))
                continue
            
            if not all(i >= 50 for i in df_per_tool.label.value_counts().tolist()):
                logging.warning('{} did not predict the minimum number of variants in the minority class required (N=50)'.format(tool))
                continue
            
            if ratio_predicted < 0.5:
                logging.info("{} predicted less than 50% of variants (N={} predictions). "
                             "Threshold analysis will be performed, but be cautious in the "
                             "interpretation.".format(tool, df_per_tool.shape[0]))

            kwargs = {'df_per_tool': df_per_tool,
                      'tool': tool,
                      'n_of_points': n_of_points,
                      'direction': direction,
                      'f1_at_ref_thresh': f1_at_ref_thresh,
                      '_loc': _loc,
                      'beta_values': beta_values,
                      'ref_thresh': reference_threshold,
                      'final_thresh': final_thresholds,
                      'outdir': outdir}
            kwargs['final_thresh'][tool] = [(kwargs['ref_thresh'], kwargs['f1_at_ref_thresh'][_loc][tool])]
            max_thr = df_per_tool[tool].max()
            min_thr = df_per_tool[tool].min()
 
            if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
                logging.info("Something strange in max/min thresholds {} {} "
                    "{}".format(tool, max_thr, min_thr))


            step = (max_thr - min_thr) / float(n_of_points)
            thr_range = np.arange(min_thr, max_thr, step)

            best_t, stats, fbetas = classify_at_threshold(df_per_tool, thr_range, **kwargs)

            for v in [False, True]:
                plot_optimal_thresholds(tool,
                                        _loc,
                                        thr_range,
                                        stats['accuracy'],
                                        stats['sensitivity'],
                                        stats['specificity'],
                                        fbetas,
                                        reference_threshold,
                                        best_t,
                                        outdir,
                                        simple=v)
            
            if do_bootstrapping:
                logging.info("Bootstrapping analysis for {}".format(tool))
                thr_f05, thr_f1, thr_f15 = [], [], []
                with multiprocessing.Pool() as p:
                    new_t = list(tqdm(p.imap(partial(_do_bootstapping, **kwargs), range(1000)), total=1000))    
                    for t in new_t: 
                        thr_f05.append(t[0.5])
                        thr_f1.append(t[1])
                        thr_f15.append(t[1.5])
                    
                cols = ['tool', 'thr_f0.5', 'thr_f1', 'thr_f1.5']
                thr_df = pd.DataFrame([[tool, thr_f05, thr_f1, thr_f15]], columns=cols)
                _calculate_ci(thr_df, **kwargs)

        _outdir_new_thresh = os.path.join(outdir, "new_thresholds")
        _outdir_new_config = os.path.join(outdir, "new_configs")
        ensure_folder_exists(_outdir_new_config)
        ensure_folder_exists(_outdir_new_thresh)


        ###########################
        ### Write new thresholds ##
        ###########################
        final_thresholds = kwargs['final_thresh']
        with open(os.path.join(_outdir_new_thresh, "proposed_thresholds_{}.tsv".format(_loc)), 'w') as out:
            out.write("Tool\tReference_threshold\tF1_score_at_ref_thresh\t{}\n".format("\t".join("Thresh_beta_{}".format(b) +
                                                                                                 "\t" +
                                                                                                 "Fbeta_{}".format(b)
                                                                                                 for b in beta_values)))

            for tool, new_thresholds in final_thresholds.items():
                out.write(
                    "{}\t{}\n".format(tool, '\t'.join(["{}\t{}".format(str(t[0]), str(t[1])) for t in new_thresholds])))

        for i, beta in enumerate(beta_values):

            config = open(os.path.join(_outdir_new_config,
                          "new_config_{}_f{}.txt".format(_loc, beta)), 'w')

            for tool, new_thresh in final_thresholds.items():

                # delete original threshold
                _new_thresh = new_thresh[1:]
                thresholds_to_return[beta][_loc][tool] = _new_thresh

                t = _new_thresh[i][0]
                # if ci are present, just report thrsh
                if isinstance(t, str):
                    t = t.split()[0]
                info = [v for v in threshold_list if v[0] == tool][0]
                vcf_field = ','.join(tools_config[tool][0])
                outline = [tool, vcf_field, info[1], str(t), info[3]]
                config.write('\t'.join(outline) + '\n')

    return thresholds_to_return