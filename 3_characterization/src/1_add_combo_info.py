import numpy as np
import os
import pandas as pd
import argparse
from functools import reduce


def get_combo_samples(combos, genotype_df):
    samples_per_gene = genotype_df.loc[genotype_df.gene.isin(combos)].samples.str.split(",").values
    samples_per_combo = reduce(lambda a,b: set(a).intersection(set(b)), samples_per_gene)
    return "|".join(sorted(samples_per_combo))

def get_sample_bmi_info(samples, phenotype_df):
    samples = samples.split("|")
    sample_phenotypes = phenotype_df.loc[phenotype_df.sample_names.astype(str).isin(samples)]
    sample_bmi_dict = dict(zip(sample_phenotypes.sample_names.astype(str), sample_phenotypes.bmi))
    sample_bmi_prs_dict = dict(zip(sample_phenotypes.sample_names.astype(str), sample_phenotypes.bmi_prs))
    sample_bmi = "|".join([str(round(sample_bmi_dict[s], 4))  if s in sample_bmi_dict.keys() else "NA" for s in samples])
    sample_bmi_prs = "|".join([str(round(sample_bmi_prs_dict[s], 4)) if s in sample_bmi_dict.keys() else "NA" for s in samples ])
    return pd.Series({"combo_samples_bmi": sample_bmi, "combo_samples_bmi_prs": sample_bmi_prs})

def get_mean_values(sample_bmi_info):
    sample_bmi_info = sample_bmi_info.split("|")
    sample_bmi_info = [float(sb) for sb in sample_bmi_info if sb!="NA" ]
    return np.mean(sample_bmi_info)

def get_combos_with_sample_info(combo_df, phenotype_df, genotype_df):
    combo_df["combos"] = combo_df.uniq_items.apply(lambda x: [g.replace("Input_", "", 1) for g in x.split("|")])
    combo_df["combo_samples"] = combo_df.combos.apply(get_combo_samples, args=(genotype_df, ))
    combo_df = pd.concat((combo_df, combo_df.combo_samples.apply(get_sample_bmi_info, args=(phenotype_df, ))), axis=1)
    combo_df["mean_bmi"] = combo_df.combo_samples_bmi.apply(get_mean_values)
    combo_df["mean_bmi_prs"] = combo_df.combo_samples_bmi_prs.apply(get_mean_values)
    return combo_df.loc[:, ["uniq_items", "combo_samples", "combo_samples_bmi", "combo_samples_bmi_prs", "mean_bmi", "mean_bmi_prs"]]


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs='+')
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--phenotype_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    genotype_df = pd.read_csv(cli_args.genotype_file)
    phenotype_df = pd.read_csv(cli_args.phenotype_file, usecols=["sample_names", "bmi", "bmi_prs"])
    combo_dfs = [pd.read_csv(cf, usecols=["uniq_items"]) for cf in cli_args.combo_files]
    combo_df = pd.concat(combo_dfs).reset_index(drop=True)

    combo_df = get_combos_with_sample_info(combo_df, phenotype_df, genotype_df)
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    combo_df.to_csv(cli_args.save_file, index=False)
