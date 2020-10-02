import numpy as np
import pandas as pd
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


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
        score = round(_p, 3)

    else:
        try:
            if any([k in _p for k in [',', '|', '&']]):
                raise ValueError('Multiple predictions were found '
                                 'for a tool that expects a single '
                                 'numeric value: {}'.format(_p))

            score = round(float(_p), 3)

        except (TypeError, ValueError):
            return np.nan

    return abs(score) if absolute else score


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
    is 'dbscSNV', which has two scores, the
    max of them will be returned

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

    # Before, DANN, ReMM and SPIDEX
    # used to average preds
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
    return round(_top, 3)


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

    :param scores:
    :return:
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

    return round(score, 3)


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


def process_spliceai(preds: pd.Series):
    """
    Processes a list of SpliceAI
    predictions (SNV or Indel)
    and returns the top scored field.

    It also sets SpliceAI predictions to null
    in cases where the variant location
    is not within a gene for the selected
    VEP consequence (which is always the first).
    Ensures fidelity in cases where there may be
    a SpliceAI prediction in a variant located
    within a gene, but another consequence (non-genic)
    was selected.

    :param pd.Series preds: Series with the
    SpliceAI predictions for the given
    SNV or Indel (One of them will always
    be None) and the location of the variant

    :return float: Top prediction
    """

    if all(v is None for v in preds.SpliceAI):
        return np.nan

    # Iterate over SNVs and Indels:
    for _p in preds.SpliceAI:

        if _p is not None:
            if _p.split("|")[1] == "0" and (preds.location == "unknown" or
                                            preds.location == "regulatory_variant" or
                                            preds.location == "mithocondrial"):
                return np.nan

            return max(map(float, _p.split("|")[2:6]))


def process_scap(preds: pd.Series):
    """
    Process S-CAP scores so that
    different reference thresholds
    thresholds are taken into account
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
                # senscore
                if float(fields[3]) > thresholds[fields[1]]:
                    return float(fields[3]) + 1
                else:
                    return float(fields[3])

            elif fields[1] == "5core":
                # We do not know whether a variant is
                # dominant and recessive, thus we test
                # for both dom and rec thresholds.
                # If any greater than threshold, prediction is
                # pathogenic
                if float(fields[5]) > thresholds['5core_dominant']:
                    return float(fields[5]) + 1

                elif float(fields[7]) > thresholds['5core_recessive']:
                    return float(fields[7]) + 1

                else:
                    return float(fields[5])

            elif fields[1] == "3core":
                if float(fields[5]) > thresholds['3core_dominant']:
                    return float(fields[5]) + 1

                elif float(fields[7]) > thresholds['3core_recessive']:
                    return float(fields[7]) + 1

                else:
                    return float(fields[5])

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
    try:
        score = float(_p)
    except ValueError:
        if "," in _p:
            score = max([float(v) for v in _p.split(',') if _p.strip() != '.'])

    if preds.location == "coding":
        return score + 1 if score > 0.416 else score

    else:
        return score + 1 if score > 0.289 else score


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

    def _get_prediction(elem: str):
        """
        Returns top prediction for a single
        element

        :param str elem: Score(s)
        :return:
        """
        if not elem or elem == "." or elem == "None":
            return np.nan

        elif "," in elem and _max is True:
            return round(max([abs(float(v)) for v in elem.split(",")]), 3)

        elif "," in elem:
            return round(min([abs(float(v)) for v in elem.split(",")]), 3)

        else:
            return abs(float(elem))

    out = []
    for e in x:
        out.append(_get_prediction(e))

    return out[0] if len(out) == 1 else out
