import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import subprocess
import os


def run_kegg_enrichment(gene_in, table_out, kegg_bash_path):
    cmd = [
        "bash", kegg_bash_path, 
        gene_in, table_out
        ]
    subprocess.run(cmd)
    return

def create_kegg_enrichment_table(combo_dfs, ncombos_list, kegg_bash_path, save_dir):
    # get list of unique genes in combos
    unique_combo_genes = set()
    for combo_df, ncombo in zip(combo_dfs, ncombos_list):
        combo_set = set(combo_df.loc[:, [f"Item_{i}" for i in range(1, ncombo + 1)]].values.flatten())
        unique_combo_genes.update(combo_set)
    # save genes to file
    gene_save_file = os.path.join(save_dir, "genelist.txt")
    with open(gene_save_file, "w") as f:
        for cg in unique_combo_genes:
            f.write(f"{cg}\n")
    # run enrichment
    enrich_table = os.path.join(save_dir, "kegg_enrichment.csv")
    run_kegg_enrichment(gene_save_file, enrich_table, kegg_bash_path)
    return

if __name__ == "__main__":
    kegg_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/kegg_Catalog/data/kegg_genes_to_traits.csv"
    combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_2.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/parsed_tables/combo_3.csv"
    kegg_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/src/scripts/kegg_enrich.sh"
    save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/protective_noicd/data/enrichment/kegg/"
    
    combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)
    create_kegg_enrichment_table([combinations2_df, combinations3_df], [2, 3], kegg_bash_path, save_dir)
