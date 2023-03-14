import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
import subprocess
import os
from tqdm import tqdm


def run_fishers_enrichment(contingency, table_out, fishers_bash_path):
    cmd = [
        "bash", fishers_bash_path, 
        contingency, table_out
        ]
    subprocess.run(cmd)
    return

def run_fishers_exact_in_r(contingency_table, fishers_bash_path, tmp_dir="/data5/deepro/tmp"):
    # create temporary contingency file
    contingency_file = os.path.join(tmp_dir, "tmp_contingency.csv")
    with open(contingency_file, "w") as f:
        f.write("data_type,condition,control\n")
        f.write(f"overlap,{contingency_table[0][0]},{contingency_table[0][1]}\n")
        f.write(f"nonoverlap,{contingency_table[1][0]},{contingency_table[1][1]}\n")
    results_file = os.path.join(tmp_dir, "tmp_res.csv")
    # run enrichment
    run_fishers_enrichment(contingency_file, results_file, fishers_bash_path)
    df = pd.read_csv(results_file)
    os.remove(contingency_file)
    os.remove(results_file)
    return df.iloc[0].values

def get_samples_with_lifestyle(lifestyle_df, lifestyle):
    df_with_lifestyle_reading = lifestyle_df.loc[~lifestyle_df[lifestyle].isna()]
    samples_with_lifestyle_reading = list(df_with_lifestyle_reading.eid)
    samples_with_lifestyle = df_with_lifestyle_reading.loc[df_with_lifestyle_reading[lifestyle]==1, "eid"].to_list()
    return set(str(s) for s in samples_with_lifestyle_reading), set(str(s) for s in samples_with_lifestyle)


def create_lifestyle_enrichment_table(combo_dfs, ncombos_list, selected_samples_file, lifestyle_df, fishers_bash_path, save_file):
    # get all samples to study
    with open(selected_samples_file, "r") as f:
        all_samples = set([l.strip() for l in f.readlines()])
    print(len(all_samples))
    
    stats=[]
    for lifestyle in tqdm(lifestyle_df.columns):
        samples_with_lifestyle_reading, samples_with_lifestyle = get_samples_with_lifestyle(lifestyle_df, lifestyle)
        all_samples_this_lifestyle = samples_with_lifestyle_reading.intersection(all_samples)
        lifestyle_samples_selected = samples_with_lifestyle.intersection(all_samples_this_lifestyle)
        # only do this if there are greater than 2000 individuals for a diagnosis
        if len(lifestyle_samples_selected)>2000:
            for combo_df, ncombo in  zip(combo_dfs, ncombos_list):
                for combo_row in combo_df.itertuples():
                    combo = ";".join([getattr(combo_row, f"Item_{c}_symbol") for c in range(1, ncombo+1)])
                    # get samples with combinations that had the lifestyle factor reading
                    combo_samples = set(getattr(combo_row, "Case_Samples").split(","))
                    combo_samples_this_lifestyle = combo_samples.intersection(all_samples_this_lifestyle)
                    noncombo_samples = all_samples.difference(combo_samples)
                    noncombo_samples_this_lifestyle = noncombo_samples.intersection(all_samples_this_lifestyle)
                    combo_with_lifestyle_samples = combo_samples_this_lifestyle.intersection(lifestyle_samples_selected)
                    combo_without_lifestyle_samples = combo_samples_this_lifestyle.difference(lifestyle_samples_selected)
                    noncombo_with_lifestyle_samples = noncombo_samples_this_lifestyle.intersection(lifestyle_samples_selected)
                    noncombo_without_lifestyle_samples = noncombo_samples_this_lifestyle.difference(lifestyle_samples_selected)

                    contingency_table = [[len(combo_with_lifestyle_samples), len(noncombo_with_lifestyle_samples)],
                                        [len(combo_without_lifestyle_samples), len(noncombo_without_lifestyle_samples)]]
                    
                    # result = run_fishers_exact_in_r(contingency_table, fishers_bash_path)
                    result = fisher_exact(contingency_table)
                    oddsratio = result[0]
                    pvalue = result[1]
                    ci_low = np.nan #result[2]
                    ci_high = np.nan #result[3]
                    stats.append([combo, lifestyle, oddsratio, pvalue, ci_low, ci_high, len(combo_with_lifestyle_samples), len(combo_without_lifestyle_samples), len(noncombo_with_lifestyle_samples), len(noncombo_without_lifestyle_samples)])
    # multiple testing for stats df
    stats = pd.DataFrame(stats, columns=['combo_name', 'icd_phenotype', 'oddsratio', 'pvalue', 'conf_int_lower', 'conf_int_upper', 'Num_samples_with_combo_and_phenotype', 'Num_samples_with_combo_and_without_phenotype', 'Num_samples_without_combo_and_with_phenotype', 'Num_samples_without_combo_and_without_phenotype'])
    stats['FDR'] = multipletests(stats.pvalue, method='fdr_bh')[1]
    # save file
    stats.sort_values("FDR").to_csv(save_file, index=False)
    return stats

if __name__ == "__main__":
    lifestyle_info_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/binarized_tables/meta/meta_pheno_table.csv"
    lifestyle_info_cols_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/binarized_tables/meta/meta_pheno_table_cols.csv"
    selected_lifestyle_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/lifestyle_v2.xlsx"
    combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    fishers_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/src/scripts/fishers_exact.sh"
    wes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables/wes.tsv"
    selected_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british/data/cases_controls/cases.txt"
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/data/enrichment/lifestyle/enrichment_specific_risk_combos.csv"
    


    selected_lifestyle_df = pd.read_excel(selected_lifestyle_file)
    selected_lifestyle_df = selected_lifestyle_df[selected_lifestyle_df.shortlist == 'X']
    selected_lifestyle_df = selected_lifestyle_df[selected_lifestyle_df.Num_exome_samples_with_phenotype > 190e3]

    lifestyle_cols_df = pd.read_csv(lifestyle_info_cols_file)
    cols_to_keep = [c for c in lifestyle_cols_df.new if c.split("_")[1] in list(map(str, selected_lifestyle_df.Phenotype_ID))]
    lifestyle_df = pd.read_csv(lifestyle_info_file, usecols=["eid"] + cols_to_keep, low_memory=False)

    combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)

    stats = create_lifestyle_enrichment_table([combinations2_df, combinations3_df], [2, 3], selected_samples_file, lifestyle_df, fishers_bash_path, save_file)
