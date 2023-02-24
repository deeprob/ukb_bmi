import pandas as pd
from statsmodels.stats.multitest import multipletests

def multi_test_hpo_table(hpo_file, save_file):
    hpo_df = pd.read_csv(hpo_file)
    hpo_df = hpo_df.rename(columns={
        "HPO_ID": "ID",
        "Definition": "hpo_phenotype",
        "OR":"oddsratio",
        "OR_lower":"conf_int_lower",
        "OR_upper":"conf_int_upper",
        "p":"pvalue",
        "Bonferroni_p":"Bonferroni_p",
        "Group":"Group",
        "Phenotype":"bmi_status"}
        )
    hpo_df = hpo_df.loc[(hpo_df.Group=="white_british")&(hpo_df.bmi_status=="high")]
    hpo_df['FDR'] = multipletests(hpo_df.pvalue, method='fdr_bh')[1]
    hpo_df = hpo_df.loc[:, ["hpo_phenotype", "oddsratio", "pvalue", "conf_int_lower", "conf_int_upper", "FDR"]]
    hpo_df.to_csv(save_file, index=False)
    return 


if __name__ == "__main__":
    hpo_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/hpo/hpo_enrichment_og.csv"
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/hpo/hpo_enrichment.csv"
    multi_test_hpo_table(hpo_file, save_file)
