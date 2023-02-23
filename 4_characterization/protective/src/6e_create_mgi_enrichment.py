import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import subprocess
import os


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

def create_mgi_enrichment_table(combo_dfs, ncombos_list, mgi_df, genes_study, fishers_bash_path, save_file):
    # restrict to genes present both in mgi and all analyzed in our study
    genes_study = list(set(genes_study).intersection(set(mgi_df.Gene)))
    mgi_df = mgi_df.loc[mgi_df.Gene.isin(genes_study)]
    mgi_genes = set(mgi_df.Gene)
    # get list of unique genes in combos
    unique_combo_genes = set()
    for combo_df, ncombo in zip(combo_dfs, ncombos_list):
        combo_set = set(combo_df.loc[:, [f"Item_{i}_symbol" for i in range(1, ncombo + 1)]].values.flatten())
        unique_combo_genes.update(combo_set)
    # keep only those present in mgi database
    unique_combo_genes = unique_combo_genes.intersection(mgi_genes)
    stats = []
    # only keep phenotypes with at least 10 genes
    trait_counts = mgi_df.Mapped_Trait.value_counts()
    trait_keep = list(trait_counts.loc[trait_counts >= 10].index)
    mgi_df = mgi_df[mgi_df.Mapped_Trait.isin(trait_keep)]
    for trait in mgi_df.Mapped_Trait.unique():
        all_genes_with_trait = set(mgi_df.loc[mgi_df.Mapped_Trait==trait, "Gene"])
        all_genes_without_trait = mgi_genes.difference(all_genes_with_trait)
        combo_genes_with_trait = unique_combo_genes.intersection(all_genes_with_trait)
        combo_genes_without_trait = unique_combo_genes.difference(combo_genes_with_trait)
        noncombo_genes_with_trait = all_genes_with_trait.difference(unique_combo_genes)
        noncombo_genes_without_trait = all_genes_without_trait.difference(unique_combo_genes)

        contingency_table = [[len(combo_genes_with_trait), len(combo_genes_without_trait)],
                            [len(noncombo_genes_with_trait), len(noncombo_genes_without_trait)]]
        
        result = run_fishers_exact_in_r(contingency_table, fishers_bash_path)
        oddsratio = result[0]
        pvalue = result[1]
        ci_low = result[2]
        ci_high = result[3]
        stats.append([trait, oddsratio, pvalue, ci_low, ci_high, len(combo_genes_with_trait), len(combo_genes_without_trait), len(noncombo_genes_with_trait), len(noncombo_genes_without_trait)])
    # multiple testing for stats df
    stats = pd.DataFrame(stats, columns=['mgi_phenotype', 'oddsratio', 'pvalue', 'conf_int_lower', 'conf_int_upper', 'Num_combo_genes_with_phenotype', 'Num_combo_genes_without_phenotype', 'Num_noncombo_genes_with_phenotype', 'Num_noncombo_genes_without_phenotype'])
    stats['FDR'] = multipletests(stats.pvalue, method='fdr_bh')[1]
    # save file
    stats.sort_values("FDR").to_csv(save_file, index=False)
    return stats

if __name__ == "__main__":
    mgi_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/MGI/data/mgi_genes.csv"
    combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    gencode_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_gencode_genes/data/gencode.v39.parsed.genes.csv"
    fishers_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/src/scripts/fishers_exact.sh"
    save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/mgi/mgi_enrichment.csv"

    # preprocess mgi data for compatibility with GWAS enrichment pipeline
    mgi_df = pd.read_csv(mgi_file, low_memory=False, usecols=["Gene", "mpi_description"])
    mgi_df.columns = ["Gene", "Mapped_Trait"]
    # split multi traits record
    mgi_df["Mapped_Trait"] = mgi_df.Mapped_Trait.str.split(",")
    mgi_df = mgi_df.explode("Mapped_Trait")
    mgi_df["Mapped_Trait"] = mgi_df.Mapped_Trait.apply(lambda x: x.strip().split("_")[1])
    # drop duplicate entries
    mgi_df = mgi_df.drop_duplicates()

    combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)
    gencode_df = pd.read_csv(gencode_file).drop_duplicates('gene_id_stripped').set_index('gene_id_stripped', drop=False)
    # get genes being studied
    genes = pd.read_csv("/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01_format2.tsv", sep='\t', nrows=0)
    genes = list(genes.columns)[1:]
    genes = [s.split('_')[1] for s in genes]
    genes = list(gencode_df.loc[genes]['gene_name'])

    stats = create_mgi_enrichment_table([combinations2_df, combinations3_df], [2, 3], mgi_df, genes, fishers_bash_path, save_file)
