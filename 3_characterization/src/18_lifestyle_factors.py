import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import kstest,ttest_ind

import utils.parsing as utpa
import utils.plotting as utpl

def get_genes_from_combo_items(combo_items, all_lifestyles):
    combo_genes = set([i.replace("Input_", "", 1) for i in combo_items if i.replace("Input_", "", 1) not in all_lifestyles])
    return combo_genes

def get_lifestyles_from_combo_items(combo_items, all_lifestyles):
    combo_lifestyles = set([i.replace("Input_", "", 1) for i in combo_items if i.replace("Input_", "", 1) in all_lifestyles])
    return combo_lifestyles

def get_lifestyle_samples(lifestyle_df, all_lifestyles):
    lifestyle_df = lifestyle_df.melt(id_vars="Sample_Name", value_vars=all_lifestyles, var_name="gene")
    lifestyle_df = lifestyle_df.loc[lifestyle_df.value>0].drop(columns="value")
    lifestyle_df = lifestyle_df.groupby("gene").agg(lambda x: ",".join(map(str, x))).reset_index().rename(columns={"Sample_Name": "samples"})
    all_lifestyle_samples = set(",".join(lifestyle_df.loc[lifestyle_df.gene.isin(all_lifestyles)].samples.values).split(","))
    return all_lifestyle_samples

def get_bmi(phenotype_df, samples):
    samples_bmi = phenotype_df.loc[phenotype_df.sample_names.isin(samples), "bmi"].values
    return samples_bmi

def get_oligogenicity_data(cohort_df, combo_gene_samples, combo_lifestyle_samples, combo_samples):
    # get category wise samples
    all_samples = set(cohort_df.sample_names.astype(str).values)
    gene_only_samples = combo_gene_samples.difference(combo_samples).difference(combo_lifestyle_samples)
    lifestyle_only_samples = combo_lifestyle_samples.difference(combo_samples).difference(combo_gene_samples)
    non_carriers = all_samples.difference(combo_gene_samples).difference(combo_lifestyle_samples)
    # get category wise bmi
    non_carriers_bmi = get_bmi(cohort_df, non_carriers)
    gene_only_bmi = get_bmi(cohort_df, gene_only_samples)
    lifestyle_only_bmi = get_bmi(cohort_df, lifestyle_only_samples)
    combo_bmi = get_bmi(cohort_df, combo_samples)
    # create data dict
    data_dict = {
        "bmi": np.concatenate((non_carriers_bmi, gene_only_bmi, lifestyle_only_bmi, combo_bmi)),
        "category": 
        ["Non carriers" for i in range(len(non_carriers_bmi))] + 
        ["Gene carriers" for i in range(len(gene_only_bmi))] + 
        ["Lifestyle carriers" for i in range(len(lifestyle_only_bmi))] +  
        ["Combo carriers" for i in range(len(combo_bmi))]
    }
    oligo_df = pd.DataFrame(data_dict)
    return oligo_df

def save_list_to_file(mylist, save_file):
    with open(save_file, "w") as f:
        for l in mylist:
            f.write(f"{l}\n")
    return

def get_common_and_unique_genes(combo_genes, compare_combo_genes, save_dir):
    common_genes = combo_genes.intersection(compare_combo_genes)
    distinct_genes = combo_genes.difference(compare_combo_genes)
    common_genes_savefile = os.path.join(save_dir, "common.list")
    distinct_genes_savefile = os.path.join(save_dir, "unique.list")
    save_list_to_file(common_genes, common_genes_savefile)
    save_list_to_file(distinct_genes, distinct_genes_savefile)
    return

def get_common_combo_samples_gene_info(cohort_df, common_combos_df, suffixes):
    gene_combo_samples = set("|".join(common_combos_df[f"combo_samples{suffixes[0]}"].values).split("|"))
    lf_combo_samples = set("|".join(common_combos_df[f"combo_samples{suffixes[1]}"].values).split("|"))
    gene_only_combo_samples = gene_combo_samples.difference(lf_combo_samples)
    cohort_filtered_df = cohort_df.loc[cohort_df.sample_names.isin(gene_combo_samples)]
    cohort_filtered_df["gene_only"] = cohort_filtered_df.sample_names.isin(gene_only_combo_samples)
    gene_only_samples_bmi = cohort_filtered_df.loc[cohort_filtered_df.gene_only==True, "bmi"]
    gene_w_lf_samples_bmi = cohort_filtered_df.loc[cohort_filtered_df.gene_only==False, "bmi"]
    res = ttest_ind(gene_only_samples_bmi, gene_w_lf_samples_bmi, alternative="less")
    return cohort_filtered_df, res.pvalue


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--gene_only_combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--compare_combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the genotype file with sample info")
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the lifestyle file with sample info")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    genotype_df = pd.read_csv(cli_args.genotype_file)
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi_prs", "bmi"])
    cohort_df["sample_names"] = cohort_df.sample_names.astype(str)
    combo_items, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    lifestyle_df = pd.read_csv(cli_args.lifestyle_file)
    all_lifestyles = set(lifestyle_df.columns)
    all_lifestyles.remove("Sample_Name")

    #################
    # oligogenicity #
    #################
    combo_genes = get_genes_from_combo_items(combo_items, all_lifestyles)
    combo_gene_samples = set(",".join(genotype_df.loc[genotype_df.gene.isin(combo_genes)].samples.values).split(","))
    combo_lifestyles = get_lifestyles_from_combo_items(combo_items, all_lifestyles)
    combo_lifestyle_samples = get_lifestyle_samples(lifestyle_df, combo_lifestyles)

    oligogenicity_df = get_oligogenicity_data(cohort_df, combo_gene_samples, combo_lifestyle_samples, combo_samples)
    save_file = os.path.join(cli_args.save_dir, "oligo.csv.gz")
    os.makedirs(cli_args.save_dir, exist_ok=True)
    oligogenicity_df.to_csv(save_file, index=False)

    fig = utpl.create_lifestyle_oligo_plot(oligogenicity_df)
    save_file = os.path.join(cli_args.save_dir, "oligo.pdf")
    utpl.save_pdf(save_file, fig)

    #########################
    # Common obesity combos #
    #########################
    lifestyle_combo_df = utpa.concat_multiple_combo_files(cli_args.combo_files)
    obesity_combo_df = utpa.concat_multiple_combo_files(cli_args.compare_combo_files)
    obesity_combo_genes, obesity_combo_samples = utpa.get_combo_info_from_files(cli_args.compare_combo_files)
    get_common_and_unique_genes(combo_genes, obesity_combo_genes, cli_args.save_dir)

    combo_lifestyles_formatted = [f"Input_{clf}" for clf in combo_lifestyles]
    lifestyle_combo_df["uniq_genes"] = lifestyle_combo_df.uniq_items.apply(lambda x: "|".join(sorted(set(x.split("|")).difference(combo_lifestyles_formatted))))
    common_combos = obesity_combo_df.merge(lifestyle_combo_df, left_on="uniq_items", right_on="uniq_genes", suffixes=('_obesity', '_lifestyle'))
    save_file = os.path.join(cli_args.save_dir, "common_combos.csv")
    common_combos.to_csv(save_file, index=False)

    gene_only_w_lf_compare_df, compare_pval = get_common_combo_samples_gene_info(cohort_df, common_combos, ('_obesity', '_lifestyle'))
    fig, ax = utpl.plot_box_single_bmi(
        gene_only_w_lf_compare_df, compare_pval, 
        xvar="gene_only", ylim=(15, 48), 
        xticklabel=["Genes Only", "Genes\nw Lifestyle"],
        order=["True", "False"],
        ttest_hline=45, ttest_text=46)
    save_file = os.path.join(cli_args.save_dir, "common_combos_compare.pdf")
    utpl.save_pdf(save_file, fig)

    ###########################################################
    # All combos gene samples only versus genes and lifestyle #
    ###########################################################
    gene_only_combo_df = utpa.concat_multiple_combo_files(cli_args.gene_only_combo_files)
    common_combos_gene_only = gene_only_combo_df.merge(lifestyle_combo_df, left_on="uniq_items", right_on="uniq_genes", suffixes=('_gene_only', '_lifestyle'))
    gene_only_w_lf_compare_df, compare_pval = get_common_combo_samples_gene_info(cohort_df, common_combos_gene_only, ('_gene_only', '_lifestyle'))
    fig, ax = utpl.plot_box_single_bmi(
        gene_only_w_lf_compare_df, compare_pval, 
        xvar="gene_only", ylim=(15, 48), 
        xticklabel=["Genes Only", "Genes\nw Lifestyle"],
        order=["True", "False"],
        ttest_hline=45, ttest_text=46)
    save_file = os.path.join(cli_args.save_dir, "gene_only_compare.pdf")
    utpl.save_pdf(save_file, fig)
