import os
import pandas as pd


def main(in_dir, in_sub_dir, infilename, wes_table, store_dir):

    # read wes file
    wes_df = pd.read_csv(wes_table, sep="\t", index_col=0, low_memory=False)
    print("wes loaded")

    # create categorical tables for rarecomb input
    lifestyle_dfs = []
    for isd in in_sub_dir:
        sub_dir = os.path.join(in_dir, isd)
        infilepath = os.path.join(sub_dir, infilename)
        # read encoded extreme lifestyle factors
        lifestyle_df = pd.read_csv(infilepath, index_col=0)
        lifestyle_dfs.append(lifestyle_df)

    # create meta table for rarecomb input
    # concatenate lifestyle factors
    lifestyle_df = pd.concat(lifestyle_dfs, axis=1)
    # merge lifestyle factors with wes data
    df = wes_df.merge(lifestyle_df, left_index=True, right_index=True)
    df.index.name = "Sample_Name"
    # fill na values which is a result of pd concat, rarecomb can't work with nas
    df = df.fillna(0)
    # save rarecomb table
    wes_with_lifestyle = os.path.join(store_dir, "wes_with_lifestyle.tsv")
    df.to_csv(wes_with_lifestyle, sep="\t")
    # save primary entities
    primary_entities_file = os.path.join(store_dir, "primary.csv")
    with open(primary_entities_file, "w") as f:
        for c in lifestyle_df.columns:
            f.write(f"{c}\n")
    # save secondary entities
    secondary_entities_file = os.path.join(store_dir, "secondary.csv")
    with open(secondary_entities_file, "w") as f:
        for c in wes_df.columns:
            f.write(f"{c}\n")           
    return


if __name__ == "__main__":
    wes_table = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables/wes.tsv"
    in_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/grouped_fields"
    in_sub_dirs = ["physical_activity", "diet"] # "mental_health"
    infilename = "cat_extremes.csv"
    store_dir = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables"
    main(in_dir, in_sub_dirs, infilename, wes_table, store_dir)
