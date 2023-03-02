import numpy as np
import pandas as pd


def create_oligo_table(combo_dfs, ncombos, profile_df, save_file):
    all_enriched_combos = []
    all_cols = []
    all_gene_cols = []
    all_lifestyle_cols = []
    for combo_df, nc in zip(combo_dfs, ncombos):
        # get unique items that came out of the enriched combinations
        enriched_combos = combo_df.loc[:, [f"Item_{c}" for c in range(1, nc+1)]].values
        columns_to_use = sorted(set(enriched_combos.flatten()))
        gene_columns_to_use = [c for c in columns_to_use if "_ENSG" in c]
        lifestyle_columns_to_use = [c for c in columns_to_use if "_ENSG" not in c]
        all_enriched_combos.extend(enriched_combos)
        all_cols.extend(columns_to_use)
        all_gene_cols.extend(gene_columns_to_use)
        all_lifestyle_cols.extend(lifestyle_columns_to_use)
    # Filter samples based on conditions
    # The conditions are:
    # 1. Samples who have any one of the genes but not the lifestyle factors
    # 2. Samples who have any one of the lifestyle factors but not the genes
    # 3. Samples who have these combinations identified by rarecomb
    # 4. Samples who do not have either the gene or the lifestyle factors
    # condition 1 eval string
    cond1 = "(" + " | ".join([f"({gene} == 1)" for gene in gene_columns_to_use]) + ")" + " & " + "(" + " & ".join([f"({lf} == 0)" for lf in lifestyle_columns_to_use]) + ")"
    # condition 2 eval string
    cond2 = "(" + " | ".join([f"({lf} == 1)" for lf in lifestyle_columns_to_use]) + ")" + " & " + "(" +" & ".join([f"({gene} == 0)" for gene in gene_columns_to_use]) + ")"
    # condition 3 eval string
    cond3 = " | ".join(["(" + " & ".join([f"({c} == 1)" for c in ec_i]) + ")" for ec_i in enriched_combos])
    # condition 4 eval string
    cond4 = " & ".join([f"({lf_g} == 0)" for lf_g in columns_to_use])  
    # create plot df
    conds = [cond1, cond2, cond3, cond4]
    cond_cats = ["genes only", "lifestyles only", "combos", "no gene or lifestyle"]    
    plot_df = pd.DataFrame()
    for cond, condcat in zip(conds, cond_cats):
        pltdf = profile_df.loc[profile_df.eval(cond), "bmi"].to_frame()
        pltdf["category"] = condcat
        plot_df = pd.concat((plot_df, pltdf))
    plot_df.to_csv(save_file, index=False)
    return 


if __name__ == "__main__":
    combo_files = [
        "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_3.csv",
        "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_4.csv",
        ]
    ncombos = [3, 4]
    profile_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables/wes_with_lifestyle.tsv"
    phenotypes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/samples_with_residuals.csv"
    combo_dfs = [pd.read_csv(cf) for cf in combo_files]
    profile_df = pd.read_csv(profile_file, sep="\t", low_memory=False, index_col=0)
    profile_df.columns = [c.lstrip('Input_') for c in profile_df.columns]
    print("read profile df")
    phenotypes_df = pd.read_csv(phenotypes_file, low_memory=False, usecols=["eid", "bmi"], index_col=0, dtype={"eid": str, "bmi": np.float64})
    profile_df = profile_df.merge(phenotypes_df, left_index=True, right_index=True)
    print("merged phenotypes to profile")
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/data/oligogenic/oligo_table.csv"
    create_oligo_table(combo_dfs, ncombos, profile_df, save_file)
