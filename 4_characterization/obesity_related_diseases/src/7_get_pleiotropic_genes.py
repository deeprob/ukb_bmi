import pandas as pd

def get_unique_genes(combo_dfs, ncombos_list):
    # get list of unique genes in combos
    unique_combo_genes = set()
    for combo_df, ncombo in zip(combo_dfs, ncombos_list):
        combo_set = set(combo_df.loc[:, [f"Item_{i}_symbol" for i in range(1, ncombo + 1)]].values.flatten())
        unique_combo_genes.update(combo_set)
    return set(g for g in unique_combo_genes if g)


obesity_risk_combo_file = [
    "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv",
    "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
]

icd_risk_combo_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/merged_combos/risk_combos.csv"
icd_protection_combo_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/merged_combos/protective_combos.csv"

obesity_risk_combo_dfs = list(map(lambda x: pd.read_csv(x), obesity_risk_combo_file))
risk_genes = get_unique_genes(obesity_risk_combo_dfs, [2, 3])

icd_risk_df = pd.read_csv(icd_risk_combo_file)
icd_risk_genes = get_unique_genes([icd_risk_df, ], [3, ])

icd_protection_df = pd.read_csv(icd_protection_combo_file)
icd_protection_genes = get_unique_genes([icd_protection_df, ], [3, ])

pleiotropic_genes = risk_genes.intersection(icd_risk_genes).intersection(icd_protection_genes)
save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/pleiotropic/genes.txt"
with open(save_file, "w") as f:
    for g in pleiotropic_genes:
        f.write(f"{g}\n")
