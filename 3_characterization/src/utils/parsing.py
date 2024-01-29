import pandas as pd
from functools import reduce

def get_gene_set_from_file(gene_file):
    with open(gene_file, "r") as f:
        genes = set([g.strip() for g in f.readlines()])
    return genes


def get_samples_with_gene_mutation(genotype_df, genes):
    all_gene_samples = set(",".join(genotype_df.loc[genotype_df.gene.isin(genes)].samples.values).split(","))
    return all_gene_samples


def get_case_cont_samples(case_cont_df):
    case_samples = set(case_cont_df.loc[case_cont_df.Output_BMI==1, "Sample_Name"].astype("str").values)
    control_samples = set(case_cont_df.loc[case_cont_df.Output_BMI==0, "Sample_Name"].astype("str").values)
    return case_samples, control_samples


def get_combo_samples_from_file(combo_file):
    combo_df = pd.read_csv(combo_file)
    all_combo_samples = set("|".join(combo_df.combo_samples.dropna().values).split("|"))
    return all_combo_samples


def get_combo_genes_from_file(combo_file):
    combo_df = pd.read_csv(combo_file)
    combo_df["uniq_items"] = combo_df.uniq_items.apply(lambda x: x.replace("Input_", ""))
    combo_genes = set("|".join(combo_df.uniq_items.values).split("|"))
    return combo_genes


def get_combo_info_from_files(combo_files):
    combo_genes = [get_combo_genes_from_file(cf) for cf in combo_files]
    combo_genes = reduce(lambda x,y: x.union(y), combo_genes)
    combo_samples = [get_combo_samples_from_file(cf) for cf in combo_files]
    combo_samples = reduce(lambda x,y: x.union(y), combo_samples)
    return combo_genes, combo_samples
