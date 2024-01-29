import argparse
import pandas as pd
from scipy.stats import ttest_ind

import utils.parsing as utpa
import utils.plotting as utpl


def prepare_oligo_table(genotype_df, cohort_df, combo_genes, combo_samples):
    # only keep the combo genes
    combo_gene_df = genotype_df.loc[genotype_df.gene.isin(combo_genes)]
    # set samples as str type
    combo_gene_df.loc[:, "samples"] = combo_gene_df.samples.str.split(",")
    # explode by samples
    combo_gene_df = combo_gene_df.explode("samples")
    # create one hot encoded data
    combo_gene_df = pd.crosstab(combo_gene_df.samples, combo_gene_df.gene)
    # get samples that have a single hit in one of the combo genes
    single_hit_samples = set(combo_gene_df.loc[combo_gene_df.sum(axis=1)<2].index)
    # get single hit sample bmi
    single_hit_pheno = cohort_df.loc[cohort_df.sample_names.isin(single_hit_samples)]
    # get combo samples bmi
    combo_pheno = cohort_df.loc[cohort_df.sample_names.isin(combo_samples)]
    # prepare the table
    single_hit_pheno["mutation"] = "Single Hit"
    combo_pheno["mutation"] = "Combo carriers"
    oligo_table = pd.concat((single_hit_pheno, combo_pheno))
    ttest_pval = ttest_ind(single_hit_pheno.bmi, combo_pheno.bmi, alternative="less").pvalue
    return oligo_table, ttest_pval

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    genotype_df = pd.read_csv(cli_args.genotype_file)
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi"], dtype={"sample_names": str, "bmi": float})
    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    oligo_df, ttest_pval = prepare_oligo_table(genotype_df, cohort_df, combo_genes, combo_samples)
    fig, ax = utpl.plot_box_single_oligo(oligo_df, ttest_pval)

    utpl.save_pdf(cli_args.save_file, fig)
