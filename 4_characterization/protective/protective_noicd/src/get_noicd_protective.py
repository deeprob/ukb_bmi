import pandas as pd
import numpy as np


sample_to_icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/ukb30075_icd10.csv"
combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/protective/data/parsed_tables/combo_2.csv"
combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/protective/data/parsed_tables/combo_3.csv"
phenotypes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/protective/data/samples.csv"

sample_to_icd_df = pd.read_csv(sample_to_icd_file, low_memory=False)

def get_protective_combos(sample_to_icd_df, combo_file, pheno_file, ncombos):
    combinations_df = pd.read_csv(combo_file, low_memory=False)
    phenotypes_df = pd.read_csv(pheno_file, low_memory=False, usecols=["eid", "bmi"], dtype={"bmi": np.float64})
    # get unique case samples with combo info
    combinations_df["Case_Samples"] = combinations_df.Case_Samples.str.replace('"', "").str.split(",")
    combinations_df = combinations_df.explode("Case_Samples").astype({"Case_Samples": int})
    # add icd info and bmi info to case samples
    df = combinations_df.merge(sample_to_icd_df, left_on="Case_Samples", right_on="eid").merge(phenotypes_df, left_on="Case_Samples", right_on="eid")
    # select low bmi df
    df_low_bmi = df# .loc[(df.Phenotype=="low_bmi")&(df.Group=="white_british")]
    # select samples without any icd10 code
    icd10_columns = [c for c in df_low_bmi.columns if (c.startswith("41204-")or(c.startswith("41202-")))]
    df_low_bmi_noicd10 = df_low_bmi.loc[df_low_bmi.loc[:, icd10_columns].isna().all(axis=1)].drop(columns=icd10_columns)
    index_cols = ["Phenotype"] +[f"Item_{i}_symbol" for i in range(1, ncombos+1)] + [f"Item_{i}" for i in range(1, ncombos+1)]
    df_low_bmi_noicd10_pivot = df_low_bmi_noicd10.pivot_table(
        index=index_cols, 
        values=["Case_Samples", "bmi"],
        aggfunc=lambda x: ",".join(map(str, x))
    ).reset_index()
    # only keep those combinations with atleast 4 samples
    df_low_bmi_noicd10_pivot = df_low_bmi_noicd10_pivot.loc[df_low_bmi_noicd10_pivot.Case_Samples.apply(lambda x: len(x.split(","))>3)]
    return df_low_bmi_noicd10_pivot

df_low_bmi_noicd10_combo2 = get_protective_combos(sample_to_icd_df, combo2_file, phenotypes_file, 2)
df_low_bmi_noicd10_combo3 = get_protective_combos(sample_to_icd_df, combo3_file, phenotypes_file, 3)


df_low_bmi_noicd10_combo2.to_csv("/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_2.csv", index=False)
df_low_bmi_noicd10_combo3.to_csv("/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_3.csv", index=False)


