import pandas as pd


# get combos that contain genes which can cause both high and low bmi
combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"

protective_combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_2.csv"
protective_combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_3.csv"


def parse_combo_df(combo2_file, combo3_file):
    combo2_df = pd.read_csv(combo2_file, usecols=["Phenotype", "Item_1", "Item_2", "Item_1_symbol", "Item_2_symbol", "Case_Samples"])
    combo2_df["Item_3"] = ""
    combo2_df["Item_3_symbol"] = ""
    combo3_df = pd.read_csv(combo3_file, usecols=["Phenotype", "Item_1", "Item_2", "Item_3", "Item_1_symbol", "Item_2_symbol", "Item_3_symbol", "Case_Samples"])
    combo_df = pd.concat((combo2_df, combo3_df))
    return combo_df.loc[:, ["Phenotype", "Item_1", "Item_2", "Item_3", "Item_1_symbol", "Item_2_symbol", "Item_3_symbol", "Case_Samples"]].reset_index(drop=True)

def get_genes_from_df(df):
    geneset = set([g for g in df.loc[:, [f"Item_{i}_symbol" for i in range(1, 4)]].values.flatten() if g])
    return geneset

def get_geneids_from_df(df):
    geneset = set([g for g in df.loc[:, [f"Item_{i}" for i in range(1, 4)]].values.flatten() if g])
    return geneset

high_bmi_df = parse_combo_df(combo2_file, combo3_file)
low_bmi_df = parse_combo_df(protective_combo2_file, protective_combo3_file)
high_bmi_genes = get_genes_from_df(high_bmi_df)
low_bmi_genes = get_genes_from_df(low_bmi_df)
both_pheno_genes = high_bmi_genes.intersection(low_bmi_genes)
high_bmi_geneids = get_geneids_from_df(high_bmi_df)
low_bmi_geneids = get_geneids_from_df(low_bmi_df)
both_pheno_geneids = high_bmi_geneids.intersection(low_bmi_geneids)

save_file_genes = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/variably_expressive/bothphenogenes.txt"
save_file_geneids = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/variably_expressive/bothphenogeneids.txt"

with open(save_file_genes, "w") as f:
    for g in both_pheno_genes:
        f.write(f"{g}\n")

with open(save_file_geneids, "w") as f:
    for g in both_pheno_geneids:
        f.write(f"{g}\n")
        