from typing import List
import numpy as np
import pandas as pd


def apply_tool_predictions(df: pd.DataFrame, thresholds: List):
    """
    Evaluate tools predictions according to the
    reference thresholds

    :param pd.DataFrame df: Input df
    :param List thresholds: List with the
    reference thresholds

    :return pd.DataFrame: Df with predictions
    for each tool
    """
    for tool, direction, threshold, *args in thresholds:
        if tool not in df.columns:
            # logging.info(tool + " scores were not processed. This tool may not be"
            #                     " in the original VCF header or did not score any"
            #                     " variant.")
            continue

        if direction == ">":
            classification = lambda x: pd.isna(x) and np.nan or x >= threshold
        else:
            classification = lambda x: pd.isna(x) and np.nan or x <= threshold

        df[tool + '_prediction'] = df[tool].apply(classification)
    return df
