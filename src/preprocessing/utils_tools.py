import logging
from typing import Union
import numpy as np
import pandas as pd


def ratio(a, b):
    if b == 0:
        return 0
    return float(a) / float(b)


def tuple2float(x: tuple):
    """
    Returns predictions for a given
    tool.

    :param tupple x: Tuple
    with info about the tool (1st elem)
    and its predictions (2nd elem).

    :return list: List of predictions (2nd elem)
    """

    return x[1] if x[1] != '.' or x[1] != '' or x[1] else np.nan


# def custom_range(series, steps):
#     mi = series.min()
#     ma = series.max()
#     step = (ma - mi) / float(steps)
#     return np.arange(mi, ma, step)


def to_numeric(preds: list, absolute: bool = True):
    """
    Converts a prediction score
    to a numeric value. This is
    for tools that output a single 
    numeric value 
    :param list preds: List with 
    predictions for a single variant.
    List is as long as the number 
    of VCF annotations that compose 
    the tool.

    :param bool absolute: Whether absolute value
    should be returned (e.g. splicing tools that
    predict PSI changes may harbor a negative
    magnitude value. Default: `True`.
    :return float: Single numeric prediction
    """

    # by default, we expect a single
    # VCF annotation for such tools
    # (e.g. CADD)
    if len(preds) > 1:
        raise ValueError('Multiple VCF fields were used in the '
                         'config file in a tool that expects a '
                         'single field with a single numeric value. '
                         'Error: {}'.format(preds))

    _p = preds[0]
    if isinstance(_p, float):
        score = round(_p, 4)

    else:
        try:
            if any([k in _p for k in [',', '|', '&']]):
                raise ValueError('Multiple predictions were found '
                                 'for a tool that expects a single '
                                 'numeric value: {}'.format(_p))

            score = round(float(_p), 4)

        except (TypeError, ValueError):
            return np.nan

    return abs(score) if absolute else score


def categorical_to_numerical(pred: list):
    """
    Parse single prediction that comes as 
    a category and transform it in a numeric
    value. It expects the categorical prediction
    to represent classification labels as seen
    in clinvar. Tools such as cVEP and EVE use
    this taxonomy. 
    
    :param list pred: Single prediction. 
    """

    class_map = {'benign': 0,
                 'likely_benign': 0.25,
                 'uncertain_significance': None,
                 'uncertain': None,
                 'conflicting': None,
                 'likely_pathogenic': 0.75,
                 'pathogenic': 1
                 }
    
    if len(pred) == 1:
        pred = pred[0]

        if pred is None:
            return np.nan

        pred = pred.lower()
        all_p = pred.split(",")
    
        # If multiple predictions
        if len(all_p) > 1:

            is_benign = all(["benign" in x for x in all_p])
            is_patho = all(["pathogenic" in x for x in all_p])
            
            if is_benign and is_patho:
                pred = "conflicting"
            
            else:
                pred = all_p[0]
        
        return class_map[pred]
    else:
        raise ValueError('Unexpected number of different preds found: {}'.format(pred))
        
    
def parse_vep_or_vcfanno(pred: str, is_max: bool,
                         _val: float,
                         absolute: bool,
                         average: bool):
    """
    Parse single prediction that harbors
    multiple values (e.g. VEP uses '&' as
    separator, vcnanno used ',') and returns
    an updated _max

    :param str pred: String with multiple scores
    to extract maximum value from
    :param bool is_max: Whether max value should
    be returned. Default: `True`. If `False`,
    the objective is to return the minimum value

    :param float _val: Top score so far
    :param bool absolute: Whether absolute values
    should be used to calculate max score
    :param bool average: Whether values should
    be average. If `True` average is calculated
    on the real values and only at the end the
    absolute number is calculated, if `absolute`
    is `True`

    :return float: Updated `_val` score
    """
    sep = "," if "," in pred else "&"
    _avg = []

    for s in pred.split(sep):
        if not s or s == 'None' or s == ".":
            continue
        s = abs(float(s)) if absolute else float(s)
        if average:
            _avg.append(float(s))
        else:
            s = abs(float(s)) if absolute else float(s)
            if np.isnan(_val):
                _val = s

            elif (is_max and s > _val) or (not is_max and s < _val):
                _val = s

    if average:
        s = abs(np.mean(_avg)) if absolute else np.mean(_avg)
        if (is_max and s > _val) or (not is_max and s < _val):
            _val = s
    return _val


def get_top_pred(preds: list,
                 is_max: bool = True,
                 absolute: bool = True,
                 average: bool = False):
    """
    Return the top scored predicion based
    on the maximum value found for a
    single tool. It may may be composed by
    multiple scores. Example of such tool
    is 'dbscSNV', which has two scores, 
    the max of them will be returned. 
    CAPICE (that has SNV or Indel scores)
    may also fit for this function.

    :param list preds: List of Predictions
    where each prediction contains a single
    score (or multiple scores separated by
    '&'. ',' or '|')

    :param bool is_max: Whether max value should
    be returned. Default: `True`. If `False`,
    the objective is to return the minimum value

    :param bool absolute: Whether absolute value
    should be returned. Maximum value will consider
    the absolute magnitude over the possible
    multiple predictions (including negative values,
    e.g. dPSI predictors). Default: `True`. If `False`
    maximum positive value will be returned

    :param bool average: Whether values of a single
    prediction (when multiple scores exist) should be
    averaged. If `True`, max between averaged multiple
    predictions will be returned.

    :return float: Top prediction
    """
    if all(not v for v in preds):
        return np.nan

    _top = np.nan
    for _p in preds:
        if _p and np.isnan(_top):
            try:
                _top = abs(float(_p)) if absolute else float(_p)

            # if multiple
            except ValueError:
                _top = parse_vep_or_vcfanno(_p,
                                            is_max=is_max,
                                            _val=_top,
                                            absolute=absolute,
                                            average=average)

        elif _p:
            try:
                if (is_max and float(_p) > _top) or (not is_max and float(_p) < _top):
                    _top = float(_p)

            except ValueError:
                _top = parse_vep_or_vcfanno(_p,
                                            is_max=is_max,
                                            _val=_top,
                                            absolute=absolute,
                                            average=average)

    return round(_top, 4)


def process_condel_carol(preds: list):
    """
    Process Condel and CAROL scores
    so that a numeric prediction
    is returned.

    This method is designed to process
    protein predictions made by some
    plugins in VEP, where an outcome
    category is added to the score.
    (e.g. SIFT scores are also returned
    like this by default, although VETA
    expects the dbNSFP format for that tool)

    :param list preds: List of predictions
    :return float: Final prediction
    """

    if len(preds) > 1:
        raise ValueError('Multiple VCF fields were used in the '
                         'config file in a tool that expects a '
                         'single field with predictions made by '
                         'VEP (e.g. Condel, CAROL, SIFT)'
                         'Error: {}'.format(preds))

    _p = preds[0]

    if not _p:
        return np.nan

    # Function expects prediction to start
    # with the 'Deleterious' or 'Neutral' strings
    elif _p.lower().startswith('deleterious') or _p.lower().startswith('neutral'):
        score = float(_p.split("(")[1].split(")")[0])

    else:
        try:
            score = float(_p)
        except ValueError:
            raise ValueError('Unknown string at the beginning '
                             'when processing Condel-like scores: {}'.format(_p))

    return round(score, 4)


def format_spliceai_fields(record, symbol: str):
    """
    Formats spliceAI scores when multiple INFO
    fields are added.

    :param cyvcf2.Record record: Single cyvcf2 record
    :param str symbol: Gene name extracted from SYMBOL
        field added by VEP
    :return str: SpliceAI output formatted in a single
        field.
    """
    fields = [0] * 10
    gname = record.INFO.get("gene_name")
    if gname and gname == symbol or gname and ',' not in gname:
        fields = [record.ALT[0], gname, record.INFO.get("spliceAI_strand"), record.INFO.get("DS_AG"),
                  record.INFO.get("DS_AL"), record.INFO.get("DS_DG"),
                  record.INFO.get("DS_DL"), record.INFO.get("DP_AG"),
                  record.INFO.get("DP_AL"), record.INFO.get("DP_DG"),
                  record.INFO.get("DP_DL")]

    elif gname and ',' in gname:  # if multiple predictions on different genes
        idx = ""
        for g in gname.split(","):
            if g == symbol:
                idx = gname.split(",").index(g)

        if isinstance(idx, str):
            gene_from_id = record.ID.split(":")[0]
            try:
                idx = gname.split(",").index(gene_from_id)

            except ValueError:
                logging.info("Variant ({} {} {}/{}) has multiple SpliceAI scores, and it was not possible "
                             "to extract the proper gene, both from the SYMBOL or ID fields. First ocurrence"
                             " will be employed. SYMBOL:{}, ID:{}, SpliceAI_genes:{}".format(
                    record.CHROM, record.POS, record.REF, record.ALT[0], symbol, record.ID, gname))
                idx = 0

        fields = [record.ALT[0], gname.split(",")[idx], record.INFO.get("DS_AG")[idx], record.INFO.get("DS_AL")[idx],
                  record.INFO.get("DS_DG")[idx], record.INFO.get("DS_DL")[idx],
                  record.INFO.get("DP_AG")[idx], record.INFO.get("DP_AL")[idx],
                  record.INFO.get("DP_DG")[idx], record.INFO.get("DP_DL")[idx]]

    return '|'.join(map(str, [x for x in fields]))


def process_spliceai(preds: pd.Series, check_gene_name: bool = True):
    """
    Processes a list of SpliceAI/CI-SpliceAI
    predictions (SNV or Indel)
    and returns the top scored field.

    It also sets SpliceAI predictions to null
    in cases where the gene of prediction do not 
    match the gene for the selected  VEP consequence.

    :param pd.Series preds: Series with the
    SpliceAI predictions for the given
    SNV or Indel (One of them will always
    be None) and the location of the variant

    :param bool check_gene_name: Retrieve SpliceAI prediction
    on the same gene assigned in the VEP SYMBOL field. Does not
    work for CI-SpliceAI, that outputs gene IDs.
    
    :return float: Top prediction
    """
    tool = "SpliceAI" if check_gene_name else "CI-SpliceAI"
    if all(v is None or v == "." for v in preds[tool]):
        return np.nan
  
    # Iterate over SNVs and Indels:
    _max = 0
    for _p in preds[tool]:

        if _p is not None and _p != ".":
            
            for ind_pred in _p.split(","):
                
                if check_gene_name:
                    if ind_pred.split("|")[1] == preds.SYMBOL:
                        return max(map(float, ind_pred.split("|")[2:6]))
                else:
                    _new_max = max(map(float, ind_pred.split("|")[2:6]))
                    if _new_max > _max:
                        _max = _new_max
                        
            return _max
    return np.nan

def process_conspliceml_like(preds: pd.Series):
    """
    Process ConSpliceML-like scores (any tool with 'something|pred' format)
    to return a single numeric value for variant
    """
    assert len(preds) == 1, "Multiple VCF fields were wrongly provided for a tool that uses the 'conspliceml_like' function."
    
    if all(v is None or v == "." for v in preds):
        return np.nan
    
    if all("|" not in v for v in preds):
        return np.nan
    
    _preds = preds[0].split(",")
    _max = 0
    for _p in _preds:
       
        new_max = round(float(_p.split("|")[1]), 3)
        if abs(new_max) > _max:
            _max = abs(new_max)

    return _max

def process_conspliceml(preds: pd.Series):
    """
    Process ConSpliceML  scores to return
    a single numeric value for the gene
    where variant occurs
    
    :return float: Prediction value
    """
  
    assert len(preds.ConSpliceML) == 1, "Multiple VCF fields were provided in the config for ConSpliceML."
    
    if all(v is None or v == "." for v in preds.ConSpliceML):
        return np.nan

    _preds = preds.ConSpliceML[0].split(",")
    for _p in _preds:
        
        if _p.split("|")[0] == preds.SYMBOL:
            return round(float(_p.split("|")[1]), 3)

    return np.nan

def process_pangolin(preds: list):
    """
    Process Pangolin scores so that the 
    maximum change in splice site usage 
    is returned (as a positive number)
    
    :param list preds: List with the
    Pangolin predictions for the given variant
    
    :return float: Prediction value to be evaluated
    """
    assert len(preds) == 1, "Multiple VCF fields were provided in the config for Pangolin."
    
    if all(v is None or v == "." for v in preds):
        return np.nan
    
    _preds = preds[0].split(",")
    _max = 0
    for _p in _preds:
        s_increase = abs(float(_p.split("|")[1].split(":")[1]))
        s_decrease = abs(float(_p.split("|")[2].split(":")[1]))
        new_max = abs(max(s_increase, s_decrease))
        if new_max > _max:
            _max = new_max

    return _max

def process_spip(preds: list):
    """
    Process SPiP scores so that the 
    the risk of the variant impacting
    splicing is returned
    
    :param list preds: List with the
    SPiP predictions for the given variant
    
    :return float: Prediction value to be evaluated
    """
    assert len(preds) == 1, "Multiple VCF fields were provided in the config for SPiP."
        
    if all(v is None or v == "." for v in preds):
        return np.nan
    
    _preds = preds[0].split(",")
    _max = 0
    for _p in _preds:
        # T|NM_001142771:g.55626401:G>A|Alter by SPiCE|98.41 % [91.47 % - 99.96 %]|1.000|-|55626401|substitution|G>A|Intron 28| ...
        #ci = float(_p.split("|")[3].split()[0])
        score = float(_p.split("|")[4])
        if score > _max:
            _max = score
    return _max
    
def process_scap(preds: pd.Series):
    """
    Process S-CAP scores so that
    different reference thresholds
    are taken into account
    for pathogenicity prediction

    If prediction is above the reference
    threshold, 1 is summed up to the final
    score. To evaluate predictions, S-CAP scores
    are then compared to 1. Note that ROC
    analysis won't be performed on such
    scores due to this transformation.

    :param pd.Series preds: Series with the
    S-CAP predictions for the given variant
    as well as its location (which won't be used)

    :return float: Prediction value
    """

    is_not_core_region = ["3intronic", "exonic", "5extended", "5intronic"]
    thresholds = {"3intronic": 0.006,
                  "exonic": 0.009,
                  "5extended": 0.005,
                  "5intronic": 0.006,
                  "5core_dominant": 0.034,
                  "5core_recessive": 0.367,
                  "3core_dominant": 0.033,
                  "3core_recessive": 0.264
                  }

    if len(preds) > 1:
        raise ValueError('Multiple VCF fields were used in the '
                         'config file when setting up S-CAP. However '
                         'this tool expects a single field with several '
                         'values to be parsed according the region defined. Error: {}'.format(preds))

    _p = preds[0]

    if not _p or _p == ".":
        return np.nan

    elif "COMMON" in _p:
        return 0

    else:
        try:
            fields = _p.split(":")
            if fields[1] in is_not_core_region:
                # rawscore
                if float(fields[2]) >= thresholds[fields[1]]:
                    return float(fields[3]) + 1
                else:
                    return float(fields[3])

            elif fields[1] == "5core":
                # We do not know whether a variant is
                # dominant and recessive, thus we test
                # for both dom and rec thresholds.
                # If any greater than threshold, prediction is
                # pathogenic
                if float(fields[4]) >= thresholds['5core_dominant']:
                    return float(fields[4]) + 1

                elif float(fields[6]) >= thresholds['5core_recessive']:
                    return float(fields[6]) + 1

                else:
                    return float(fields[4])

            elif fields[1] == "3core":
                if float(fields[4]) >= thresholds['3core_dominant']:
                    return float(fields[4]) + 1

                elif float(fields[6]) >= thresholds['3core_recessive']:
                    return float(fields[6]) + 1

                else:
                    return float(fields[4])

        # if prediction is like '1:exonic:.:.:.:.:.:.'
        except ValueError:
            return np.nan


def process_trap(preds: pd.Series):
    """
    Process TrAP scores so that
    different coding and intronic
    thresholds are taken into account
    for pathogenicity prediction

    If prediction is above the reference
    threshold, 1 is summed up to the final
    score. To evaluate predictions, TraP scores
    are then compared to 1. Note that ROC
    analysis won't be performed on such
    scores due to this transformation.

    :param pd.Series preds: Series with the
    TraP predictions for the given variant
    as well as its location

    :return float: Prediction value
    """
   
    if all(v is None for v in preds.TraP):
        return np.nan

    if len(preds.TraP) > 1:
        raise ValueError('Multiple VCF fields were used in the '
                         'config file when setting up TraP. However '
                         'this tool expects a single field with numeric '
                         'predictions. Error: {}'.format(preds.TraP))

    _p = preds.TraP[0]
    score = ""
    try:
        score = float(_p)
    except ValueError:
        if "," in _p:
            score = max([float(v) for v in _p.split(',') if _p.strip() != '.'])

    if preds.location in ["splice_site", "splice_region", "intronic",
                          "deep_intronic", "5primeUTR", "3primeUTR"]:
        return score + 1 if score > 0.289 else score

    else:
        return score + 1 if score > 0.416 else score


def process_kipoi_tools(x: list, _max: bool = True):
    """
    Process scores given by Kipoi, which usually
    comprise several predictions, one for each
    genomic feature present in the gft used
    that overlaps with the input variant

    For each element (if multiple VCF fields
    define the score) in x, one single prediction
    is returned. Max value from each element is
    retrieved by default (_max).

    :param x: List of lists with the scores
    for a single variant
    :param bool _max: Whether max value should
    be retrieved. Default: `True`. If `False`,
    min is used.

    :return Union[float, List]: Top predictions
    for each element. If tool is defined by one
    single VCF field, returns a single value
    """

    def _get_prediction(elem: Union[str, float]):
        """
        Returns top prediction for a single
        element

        :param str elem: Score(s)
        :return float:
        """
        if isinstance(elem, float):
            return abs(elem) if _max is True else elem

        if not elem or elem == "." or elem == "None":
            return np.nan

        elif "," in elem and _max is True:
            return round(max([abs(float(v)) for v in elem.split(",")]), 4)

        elif "," in elem:
            return round(min([abs(float(v)) for v in elem.split(",")]), 4)

        elif _max is True:
            return abs(float(elem))

        else:
            return float(elem)

    out = []
    for e in x:
        out.append(_get_prediction(e))

    return out[0] if len(out) == 1 else out

def process_labranchor(preds: list):
    """
    Process LaBranchoR scores so that the 
    the risk of the variant impacting
    splicing is returned
    
    :param list preds: List with the
    LaBranchoR predictions for the given variant
    
    :return float: Prediction value to be evaluated
    """
    assert len(preds) == 1, "Multiple VCF fields were provided in the config for LaBranchoR."
    
    if all(v is None or v == "." for v in preds):
        return np.nan
    
    _preds = preds[0].split(",")
    _max = 0
    for _p in _preds:
        positions = map(float, _p.split("|"))
        score = abs(max(positions, key=abs))
        if score > _max:
            _max = score
    return _max


def process_absplice(preds: list):
    """
    Process AbSplice_DNA scores so that the 
    the tissue with the highest score is returned


    :param list preds: List with the per-tissue
    AbSplice_Dna predictions for the given variant
    
    :return float: Prediction value to be evaluated
    """
    assert len(preds) == 1, "Multiple VCF fields were provided in the config for AbSplice_DNA."
    
    if preds[0] is None:
        return np.nan
    
    # If variant overlaps multiple genes, there are multiple predictions
    _preds = preds[0].split(",")
    
    _max = -1
    for _p in _preds:
        
        # index 0 is the gene name
        _p = _p.split("|")[1:]
        
        top_score = max([-1 if x == '' else float(x) for x in _p])
        if top_score > _max:
            _max = top_score
    
    if _max == -1:
        return np.nan 
    return _max