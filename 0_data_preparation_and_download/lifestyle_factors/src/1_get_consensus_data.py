import os
import pandas as pd
import numpy as np


def get_consensus(row):
    values = [i for i in row.unique() if not pd.isnull(i)]
    uniq_val = np.nan
    if len(values)>1:
        uniq_val = "inconsistent"
    elif len(values) == 1:
        uniq_val = values[0]
    return uniq_val

def get_mean(row):
    return row.mean()

def process_sample_info(df, categorical_fields, numerical_fields):

    # Make consensus column for categorical fields
    for field in categorical_fields:
        field_columns = df.loc[:, df.columns.str.startswith(field)].columns.values
        df[field] = df.loc[:, field_columns].apply(get_consensus, axis=1)
        df = df.drop(columns=field_columns)

    # Get mean value for numerical fields
    for field in numerical_fields:
        field_columns = df.loc[:, df.columns.str.startswith(field)].columns.values
        df[field] = df.loc[:, field_columns].apply(get_mean, axis=1)
        df = df.drop(columns=field_columns)

    # filter samples with nan values for numerical fields
    df = df.loc[~df.loc[:, numerical_fields].isna().any(axis=1)]

    return df


if __name__ == "__main__":

    binary_fields = ["met"]

    integer_fields = [
        "sleep", "tv", "computer", "cookedvegetable", "salad", "freshfruit", 
        "driedfruit", "oilyfish", "nonoilyfish", "procmeat", "poultry", "beef",
        "mutton", "pork", "bread", "cereal", "tea", "coffee", "water"
        ]

    categorical_fields = ["alcohol", "smokecurr", "smokepast"]

    dfs = []
    for i in range(11):
        df = pd.read_csv(f"/data6/deepro/ukb_bmi/0_data_preparation_and_download/lifestyle_factors/data/lifestyle_raw/lifestyle_block{i}.csv.gz")
        filtered_df = process_sample_info(df, categorical_fields + binary_fields + integer_fields, [])
        dfs.append(filtered_df)

    df = pd.concat(dfs, axis=0)
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/lifestyle_factors/data/lifestyle_processed/filtered_lifestyle_info.csv.gz"
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    df.to_csv(save_file, index=False)
