import os
import argparse
import pandas as pd
from functools import reduce
from scipy.stats import ttest_ind

import utils.parsing as utpa
import utils.plotting as utpl

def get_combo_samples_from_files(combo_files):
    combo_samples = [utpa.get_combo_samples_from_file(cf) for cf in combo_files]
    combo_samples = reduce(lambda x,y: x.union(y), combo_samples)
    return combo_samples

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_files", type=str, help="Filepath of the cohort pheno file", nargs="+")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--combo_dir", type=str, help="Filepath of the combo files dir given by Rarecomb")
    parser.add_argument("--combo_filepre", type=str, help="Prefix of the combo files for groups")
    parser.add_argument("--groups", type=str, help="The different groups to compare", nargs="+", default=[])
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    if len(cli_args.groups)==0:
        cohort_df = pd.read_csv(cli_args.cohort_files[0], usecols=["sample_names", "bmi", "bmi_prs"])
        all_combo_samples = get_combo_samples_from_files(cli_args.combo_files)
        cohort_df["combo_carriers"] = cohort_df.sample_names.astype(str).isin(all_combo_samples)
        combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==True]
        non_combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==False]
        ttest_pval = ttest_ind(non_combo_hit_pheno.bmi, combo_hit_pheno.bmi, alternative="less").pvalue
        fig, ax = utpl.plot_box_single_bmi(cohort_df, ttest_pval)
    
    else:
        assert len(cli_args.cohort_files)==len(cli_args.groups)
        cohort_df_list = []
        ttest_pvals = []
        for group, cohort_file in zip(cli_args.groups, cli_args.cohort_files):
            cohort_df = pd.read_csv(cohort_file, usecols=["sample_names", "bmi", "bmi_prs"])
            combo_files = [cf.path for cf in os.scandir(os.path.join(cli_args.combo_dir, group)) if cf.name.startswith(cli_args.combo_filepre)]
            all_combo_samples = get_combo_samples_from_files(combo_files)
            cohort_df["combo_carriers"] = cohort_df.sample_names.astype(str).isin(all_combo_samples)
            cohort_df["group"] = group
            combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==True]
            non_combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==False]
            ttest_pval = ttest_ind(non_combo_hit_pheno.bmi, combo_hit_pheno.bmi, alternative="less").pvalue
            cohort_df_list.append(cohort_df)
            ttest_pvals.append(ttest_pval)
        cohort_df = pd.concat(cohort_df_list)
        fig, ax = utpl.plot_box_oligo(cohort_df, ttest_pvals, xvar="group", huevar="combo_carriers", order=cli_args.groups, hue_order=[False, True], figsize=(12, 6))

    utpl.save_pdf(cli_args.save_file, fig)
