import os
import json
import numpy as np
import pandas as pd


def join_array(xarr):
    xarr = [i for i in xarr if not pd.isnull(i)]
    return "||".join(xarr)

def get_all_items_combos_samples(combos_files, max_combo=4, firstn=5, genes_of_interest=[]):
    combos_df = pd.concat([pd.read_csv(cf) for cf in combos_files])
    if not genes_of_interest:
        combos_df = combos_df.sort_values("Effect_Size").head(firstn)
    else:
        # get any combination that contains the genes of interest
        combos_df = combos_df.loc[combos_df.loc[:, [f"Item_{i}_symbol" for i in range(1, max_combo+1)]].apply(lambda row: row.str.fullmatch("|".join(genes_of_interest), case=False, na=False)).any(axis=1)]
    # Store the gene ids of all the genes from the selected combinations
    unique_items = set([i for i in combos_df.loc[:, [f"Item_{i}" for i in range(1, max_combo+1)]].values.flatten() if not pd.isnull(i)])
    # get the case samples which contains the combos
    case_samples = set(sum(list(map(lambda v: v.split(","), combos_df.Case_Samples.values)), []))
    # assign each combo a unique name 
    unique_combo_names = list(map(join_array, combos_df.loc[:, [f"Item_{i}" for i in range(1, max_combo+1)]].values))
    return unique_items, unique_combo_names, case_samples

def read_wes_file(wes_file, unique_items):
    cols_to_use = [f"Input_{gid}" for gid in unique_items] + ["Sample_Name"]
    wes_df = pd.read_csv(wes_file, sep="\t",  low_memory=False, index_col=0, usecols=cols_to_use)
    return wes_df

def create_profile(wes_df, case_samples, unique_combos):
    case_samples = list(map(int, case_samples))
    # only keep those individuals who have the combo and are cases
    selected_wes_df = wes_df.loc[wes_df.index.isin(case_samples)]
    cond = " | ".join(["(" + " & ".join([f"(Input_{c} == 1)" for c in ec_i]) + ")" for ec_i in list(map(lambda x: x.split("||"), unique_combos))])
    profile_df = selected_wes_df.loc[selected_wes_df.eval(cond)]
    return profile_df


if __name__ == "__main__":
    wes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables/wes_with_lifestyle.tsv"
    combos_files = [
        "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_3.csv",
        "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_4.csv"
    ]
    save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/data/profiles"
    neuro_genes = ["DNAH1", "DNAH2", "DNAH5", "DNAH7", "DNAH9", "DNAH11", "DNAH17", "APC", "PC", "CAPN1", "COX6B2", "ITPR1", "PRPH", "RYR3", "TUBA3D"] 
    dnah_genes = neuro_genes[:7]
    other_neuro = neuro_genes[7:]
    for goi,name in zip([dnah_genes, other_neuro], ["dnah_genes", "other_neuro_genes"]):
        # read unique genes, combos, obese samples
        unique_items, unique_combos, case_samples = get_all_items_combos_samples(combos_files, genes_of_interest=goi)
        # read the wes profile of obese individuals
        wes_df = read_wes_file(wes_file, unique_items)
        # get profile 
        profile_df = create_profile(wes_df, case_samples, unique_combos)
        # save the profile df
        save_file = os.path.join(save_dir, f"{name}.csv")
        profile_df.to_csv(save_file)
