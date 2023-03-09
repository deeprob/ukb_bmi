import pandas as pd
import subprocess
import os


def run_tissue_enrichment(table_in, gene_list, table_out, fishers_tissue_bash_path):
    cmd = [
        "bash", fishers_tissue_bash_path, 
        table_in, gene_list, table_out
        ]
    subprocess.run(cmd)
    return

def create_tissue_enrichment_table(combo_dfs, ncombos_list, tissue_file, fishers_tissue_bash_path, save_dir):
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
    enrich_table = os.path.join(save_dir, "tissue_enrichment.csv")
    run_tissue_enrichment(tissue_file, gene_save_file, enrich_table, fishers_tissue_bash_path)
    return

if __name__ == "__main__":
    tissue_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_tissue_expression_data/data/gtex_tissue_specific_expression.csv"
    combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
    combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    fishers_tissue_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/src/scripts/fishers_tissue.sh"
    save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/tissue/"
    
    combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
    combinations3_df = pd.read_csv(combinations3_file, low_memory=False)
    create_tissue_enrichment_table([combinations2_df, combinations3_df], [2, 3], tissue_file, fishers_tissue_bash_path, save_dir)
