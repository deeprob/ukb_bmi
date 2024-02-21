import os
import argparse
from functools import reduce
import numpy as np
import pandas as pd


def get_combo_samples(combos, genotype_df, cohort_samples):
    samples_per_gene = genotype_df.loc[genotype_df.gene.isin(combos)].samples.str.split(",").values
    samples_per_combo = reduce(lambda a,b: set(a).intersection(set(b)), samples_per_gene)
    samples_per_combo = cohort_samples.intersection(samples_per_combo)
    return "|".join(sorted(samples_per_combo))

def get_combos_with_sample_info(combo_df, cohort_samples, genotype_df):
    combo_df["combos"] = combo_df.uniq_items.apply(lambda x: [g.replace("Input_", "", 1) for g in x.split("|")])
    combo_df["combo_samples"] = combo_df.combos.apply(get_combo_samples, args=(genotype_df, cohort_samples))
    return combo_df.loc[:, ["uniq_items", "combo_samples"]]


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Combo info")
    parser.add_argument("--combo_file", type=str, help="Filepath of the combo files given by Rarecomb")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the lifestyle factor matrix", default="")

    cli_args = parser.parse_args()

    combo_df = pd.read_csv(cli_args.combo_file, usecols=["uniq_items"])
    genotype_df = pd.read_csv(cli_args.genotype_file)
    if cli_args.lifestyle_file:
        lifestyle_df = pd.read_csv(cli_args.lifestyle_file)
        cols_to_melt = list(lifestyle_df.columns)
        lifestyle_df = lifestyle_df.melt(id_vars="Sample_Name", value_vars=cols_to_melt, var_name="gene")
        lifestyle_df = lifestyle_df.loc[lifestyle_df.value>0].drop(columns="value")
        lifestyle_df = lifestyle_df.groupby("gene").agg(lambda x: ",".join(map(str, x))).reset_index().rename(columns={"Sample_Name": "samples"})
        genotype_df = pd.concat((genotype_df, lifestyle_df)).reset_index(drop=True)
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi"])
    cohort_samples = set(cohort_df.sample_names.astype(str))

    combo_df = get_combos_with_sample_info(combo_df, cohort_samples, genotype_df)
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    combo_df.to_csv(cli_args.save_file, index=False)
