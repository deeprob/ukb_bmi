import argparse
import pandas as pd
import numpy as np


import utils.parsing as utpa
import utils.plotting as utpl


def get_gene_samples(row, gene_col, genotype_df):
    if type(row[gene_col])!=str:
        if pd.isnull(row[gene_col]):
            return ""
    gene = row[gene_col].replace("Input_", "")
    all_gene_samples = set(genotype_df.loc[genotype_df.gene==gene, "samples"].values[0].split(","))
    all_combo_samples = set(row.combo_samples.split("|"))
    all_single_gene_samples = all_gene_samples.difference(all_combo_samples)
    return "|".join(all_single_gene_samples)

def get_combo_bmi_info(combo_df, cohort_df, ncombo):
    combo_samples = set("|".join(combo_df.combo_samples.values).split("|"))
    gene_samples = [set("|".join(combo_df[f"Gene{i}_samples"].values).split("|")) for i in range(1, ncombo+1)]
    all_samples = set(cohort_df.sample_names.values)
    non_gene_samples = all_samples.difference(gene_samples[0].union(gene_samples[1]).union(combo_samples))
    gene1_only_samples = gene_samples[0].difference(gene_samples[1])
    gene2_only_samples = gene_samples[1].difference(gene_samples[0])
    return non_gene_samples, gene1_only_samples, gene2_only_samples, combo_samples

def get_bmi(phenotype_df, samples):
    return phenotype_df.loc[phenotype_df.sample_names.isin(samples), "bmi"].values

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Additive model")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--combo_file", type=str, help="Filepath of the combo info file")
    parser.add_argument("--save_file", type=str, help="Filepath where additive plot will be stored")

    cli_args = parser.parse_args()
    genotype_df = pd.read_csv(cli_args.genotype_file)
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi"], dtype={"sample_names": str, "bmi": float})
    combo_gene_df = pd.read_csv(cli_args.combo_file).dropna()

    combo_gene_df[["Gene1", "Gene2"]] = combo_gene_df.uniq_items.str.split("|", expand=True)
    combo_gene_df["Gene1_samples"] = combo_gene_df.apply(get_gene_samples, args=("Gene1", genotype_df), axis=1)
    combo_gene_df["Gene2_samples"] = combo_gene_df.apply(get_gene_samples, args=("Gene2", genotype_df), axis=1)


    non_gene_samples, gene1_samples, gene2_samples, combo_samples = get_combo_bmi_info(combo_gene_df, cohort_df, 2)
    all_bmi_info = [get_bmi(cohort_df, s) for s in [non_gene_samples, gene1_samples, gene2_samples, combo_samples]]

    plot_df = pd.DataFrame({
        "Gene A": ["No" for i in range(len(all_bmi_info[0]))] + ["Yes" for i in range(len(all_bmi_info[1]))] + ["No" for i in range(len(all_bmi_info[2]))] + ["Yes" for i in range(len(all_bmi_info[3]))],
        "Gene B": ["No" for i in range(len(all_bmi_info[0]))] + ["No" for i in range(len(all_bmi_info[1]))] + ["Yes" for i in range(len(all_bmi_info[2]))] + ["Yes" for i in range(len(all_bmi_info[3]))],
        "BMI": np.concatenate(all_bmi_info)
    })

    fig, ax = utpl.get_interaction_plot(plot_df)
    utpl.save_pdf(cli_args.save_file, fig)
