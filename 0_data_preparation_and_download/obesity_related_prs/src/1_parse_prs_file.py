import os
import pandas as pd


def parse_file(prs_dir, save_file):
    dfs = []
    for file in os.scandir(prs_dir):
        if file.name.endswith(".csv.gz"):
            df = pd.read_csv(file.path, index_col=0)
            df = df.dropna(how="all").reset_index()
            dfs.append(df)
    meta_df = pd.concat(dfs)
    meta_df.to_csv(save_file, index=False)
    return 


if __name__ == "__main__":
    prs_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/obesity_related_prs/data/prs_raw"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/obesity_related_prs/data/prs_processed/obesity_related.csv.gz"
    parse_file(prs_dir, save_file)
