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

def process_sample_info(df: pd.DataFrame, categorical_fields: dict, numerical_fields: dict):
    # Remove samples with >=10 third-degree relatives
    df = df.loc[df.genetic_kinship_to_other_participants != "Ten or more third-degree relatives identified"]

    # Make consensus column for categorical fields
    for field, field_measures in categorical_fields.items():
        df[field] = df.loc[:, field_measures].apply(get_consensus, axis=1)
        df = df.drop(columns=field_measures)

    # Get mean value for numerical fields
    for field, field_measures in numerical_fields.items():
        df[field] = df.loc[:, field_measures].apply(get_mean, axis=1)
        df = df.drop(columns=field_measures)

    # filter samples with nan values for numerical fields
    df = df.loc[~df.loc[:, numerical_fields.keys()].isna().any(axis=1)]

    return df


if __name__ == "__main__":

    numerical_fields = {
        "bmi": ["bmi0", "bmi1", "bmi2"],
        "age": ["age_assessment0", "age_assessment1", "age_assessment2"],
        }

    categorical_fields = {
        "menopause": ["menopause0", "menopause1", "menopause2"],
        "ethnic_background": ["ethnic_background0", "ethnic_background1", "ethnic_background2"]
        }

    dfs = []
    for i in range(11):
        df = pd.read_csv(f"/data6/deepro/bmi_project/0_data_preparation_and_download/phenotype/data/bmi_raw/bmi_block{i}.csv.gz")
        filtered_df = process_sample_info(df, categorical_fields, numerical_fields)
        dfs.append(filtered_df)

    df = pd.concat(dfs, axis=0)
    save_file = "/data6/deepro/bmi_project/0_data_preparation_and_download/phenotype/data/bmi_processed/filtered_bmi_info.csv.gz"
    df.to_csv(save_file, index=False)
