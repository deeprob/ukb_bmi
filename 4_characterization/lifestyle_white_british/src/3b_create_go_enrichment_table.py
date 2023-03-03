import pandas as pd
import subprocess
import os


def run_gsea_enrichment(gene_in, table_out, gsea_bash_path):
    cmd = [
        "bash", gsea_bash_path, 
        gene_in, table_out
        ]
    subprocess.run(cmd)
    return

def create_gsea_enrichment_table(combo_dfs, ncombos_list, gsea_bash_path, save_dir):
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
    enrich_table = os.path.join(save_dir, "go_enrichment.csv")
    run_gsea_enrichment(gene_save_file, enrich_table, gsea_bash_path)
    return

if __name__ == "__main__":
    gsea_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/gsea_Catalog/data/gsea_genes_to_traits.csv"
    combinations4_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_4.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_3.csv"
    gsea_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/src/scripts/gsea_enrich.sh"
    save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/data/enrichment/go/"
    
    combinations4_df = pd.read_csv(combinations4_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)
    create_gsea_enrichment_table([combinations3_df, combinations4_df], [3, 4], gsea_bash_path, save_dir)
