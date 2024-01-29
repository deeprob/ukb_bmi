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
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi", "bmi_prs"])
    all_combo_samples = get_combo_samples_from_files(cli_args.combo_files)

    cohort_df["combo_carriers"] = cohort_df.sample_names.astype(str).isin(all_combo_samples)
    combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==True]
    non_combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==False]
    ttest_pval = ttest_ind(non_combo_hit_pheno.bmi, combo_hit_pheno.bmi, alternative="less").pvalue
    
    fig, ax = utpl.plot_box_single_bmi(cohort_df, ttest_pval)

    utpl.save_pdf(cli_args.save_file, fig)
