import pandas as pd


def read_combos_files(combo_files):
    combos_dfs = list(map(lambda x: pd.read_csv(x), combo_files))
    combo_df = pd.concat(combos_dfs)
    return combo_df

def get_genes_from_df(df):
    geneset = set([g for g in df.loc[:, [f"Item_{i}_symbol" for i in range(1, 4)]].values.flatten() if not pd.isnull(g)])
    return geneset

def get_geneids_from_df(df):
    geneset = set([g for g in df.loc[:, [f"Item_{i}" for i in range(1, 4)]].values.flatten() if not pd.isnull(g)])
    return geneset

def get_gene_specific_combo_info(df, gene):
    df = df.loc[(df.Item_1_symbol==gene)|((df.Item_2_symbol==gene))|((df.Item_3_symbol==gene))].reset_index(drop=True)
    return df

def get_case_samples(df):
    high_pheno_cases = ",".join(df.loc[df.Phenotype=="high_bmi"].Case_Samples.str.strip('"').values.flatten()).split(",")
    low_pheno_cases = ",".join(df.loc[df.Phenotype=="low_bmi"].Case_Samples.str.strip('"').values.flatten()).split(",")
    return high_pheno_cases, low_pheno_cases

if __name__ == "__main__":
    risk_combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
    risk_combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    protective_combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_2.csv"
    protective_combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_3.csv"
    variants_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv"
    selected_variants_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/variably_expressive/gene_variants.csv"
    
    risk_df = read_combos_files([risk_combinations2_file, risk_combinations3_file])
    protective_df = read_combos_files([protective_combinations2_file, protective_combinations3_file])
    risk_genes = get_genes_from_df(risk_df)
    protective_genes = get_genes_from_df(protective_df)
    both_pheno_genes = risk_genes.intersection(protective_genes)
    both_pheno_df = pd.concat((risk_df, protective_df))
    both_pheno_df = both_pheno_df.loc[(both_pheno_df.Item_1_symbol.isin(both_pheno_genes))|((both_pheno_df.Item_2_symbol.isin(both_pheno_genes)))|((both_pheno_df.Item_3_symbol.isin(both_pheno_genes)))].reset_index(drop=True)
    variants_df = pd.read_csv(variants_file, sep="\t", low_memory=False, usecols=["Sample", "variant_id", "Gene", "SYMBOL", "Mut_type"], dtype=str)
    risk_case_samples, protection_case_samples = get_case_samples(both_pheno_df)
    risk_variants_df = variants_df.loc[(variants_df.Sample.isin(risk_case_samples)&(variants_df.SYMBOL.isin(both_pheno_genes)))].reset_index(drop=True)
    protection_variants_df = variants_df.loc[(variants_df.Sample.isin(protection_case_samples)&(variants_df.SYMBOL.isin(both_pheno_genes)))].reset_index(drop=True)
    risk_variants_df["obesity_type"] = "risk"
    protection_variants_df["obesity_type"] = "protection"    
    selected_variants_df = pd.concat((risk_variants_df, protection_variants_df))
    selected_variants_df.to_csv(selected_variants_file, index=False)
    