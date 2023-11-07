import os
import pandas as pd
import numpy as np


def binarize_lifestyle(df, cat_encoding_dict, numerical_fields, combine_field_dict):
    for field, pattern in cat_encoding_dict.items():
        df[field] = df[field].str.fullmatch(rf"{pattern}", na=False).astype(int)

    # encode all non numeric values as numeric ones
    df.loc[:, numerical_fields] = df.loc[:, numerical_fields].replace(["Less than an hour a day", "Less than one", "Do not know", "Prefer not to answer", "inconsistent"], [0, 0, 0, 0, 0]).astype(float)
    for num_field in numerical_fields:
        num_field_quantile_thresh = df[num_field].quantile(q=0.95)
        df[num_field] = (df[num_field]>num_field_quantile_thresh).astype(int)
    
    # combine fields to one numeric value
    for new_field, old_fields in combine_field_dict.items():
        df[new_field] = (df.loc[:, old_fields]==1).any(axis=1).astype(int)
        df = df.drop(columns=old_fields)
    return df

if __name__=="__main__":
    cat_encoding_dict = {
        "alcohol": "Daily or almost daily",
        "smokecurr": "Yes, on most or all days",
        "smokepast": "Yes, on most or all days",
        "met": "No",
        "oilyfish": "5-6 times a week|Once or more daily",
        "nonoilyfish": "5-6 times a week|Once or more daily",
        "procmeat": "5-6 times a week|Once or more daily",
        "poultry": "5-6 times a week|Once or more daily",
        "beef": "5-6 times a week|Once or more daily",
        "mutton": "5-6 times a week|Once or more daily",
        "pork": "5-6 times a week|Once or more daily"
    }

    numerical_fields = [
        'sleep', 'tv', 'computer', 'cookedvegetable', 'salad', 'freshfruit', 
        'driedfruit', 'bread', 'cereal', 'tea', 'coffee', 'water'
        ]
    
    combine_field_dict = {
        "smoke": ["smokecurr", "smokepast"],
        "fish": ["oilyfish", "nonoilyfish"],
        "meat": ["procmeat", "beef", "mutton", "pork"],
        "sedentary": ["tv", "computer"]
    }

    lifestyle_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/lifestyle_factors/data/lifestyle_processed/filtered_lifestyle_info.csv.gz"
    lifestyle_df = pd.read_csv(lifestyle_file)
    binarized_df = binarize_lifestyle(lifestyle_df, cat_encoding_dict, numerical_fields, combine_field_dict)

    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/lifestyle_factors/data/lifestyle_processed/filtered_lifestyle_binarized.csv.gz"
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    binarized_df.rename(columns={"sample_names": "Sample_Name"}).to_csv(save_file, index=False)
