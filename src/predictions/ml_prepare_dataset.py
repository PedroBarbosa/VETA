import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn_pandas import DataFrameMapper, CategoricalImputer, gen_features
from imblearn.under_sampling import RandomUnderSampler

MISSING_CONSTANT = -9
forbidden_prefixes = [
    'clinvar',
    'class'
]

label_fields = []

def build_converter(df, threshold_list):
    """ Prepares a mapper between Pandas Dataframe and sklearn matrix """
    fields = []
    for f in df.columns:
        if any([ f.startswith(p) for p in forbidden_prefixes ]):
            continue
        elif f in label_fields:
            fields.extend(gen_features(
                columns = [f],
                classes = [CategoricalImputer, LabelBinarizer]
            ))
        elif df[f].dtype == np.float64:
            fields.append(([f], SimpleImputer(strategy="most_frequent")))
        elif df[f].dtype == bool:
            df[f] = pd.to_numeric(df[f])
            fields.append(([f], None))
        #else:
        #    print(f, df[f].head())
    fields = []
    for t, _, _, _, _ in threshold_list:
        fields.append(([t], SimpleImputer(strategy="median")))

    mapper = DataFrameMapper(fields)
    return mapper

def convert_data(df, mapper):
    X = mapper.fit_transform(df)
    y = df['class'].values
    return X, y

def prepare_dataset(df, threshold_list):
    for tool in threshold_list:
        df[[tool[0]]] = df[[tool[0]]].fillna(value=MISSING_CONSTANT)
    mapper = build_converter(df, threshold_list)
    X, y = convert_data(df, mapper)
    return X, y

def undersample(X, y):
    rus = RandomUnderSampler(random_state=0)
    X_resampled, y_resampled = rus.fit_sample(X, y)
    return X_resampled, y_resampled