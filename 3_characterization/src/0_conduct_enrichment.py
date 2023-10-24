import os
import argparse
import pandas as pd
import subprocess
import json
import requests
from typing import List
import time

"""
Resources:
https://maayanlab.cloud/Enrichr/help#api
https://maayanlab.cloud/Enrichr/#libraries
"""

def get_population_geneset(gene_burden_df, save_dir):
    popgenefile = os.path.join(save_dir, "genes_pop.list")
    with open(popgenefile, "w") as f:
        for g in gene_burden_df.gene.unique():
            f.write(f"{g}\n")
    popgenset = set(gene_burden_df.gene.unique())
    return list(popgenset)


def get_combo_genes(combo_dfs: List[pd.DataFrame], save_dir):
    geneset = set()
    for combo_df in combo_dfs:
        genes = set("|".join(combo_df.uniq_items.values).split("|"))
        genes = [g.replace("Input_", "", 1) for g in genes]
        geneset = geneset.union(genes)
    genefile = os.path.join(save_dir, "genes.list")
    with open(genefile, "w") as f:
        for g in geneset:
            f.write(f"{g}\n")
    return list(geneset)

def save_enrichr_results(results, save_file):
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    with open(save_file, "w") as f:
        f.write("Term,p_val,adj_pval,odds_ratio,combined_score,genes\n")
        for lines in results:
            if lines[2]<0.05:
                f.write(f"{lines[1].replace(',', '')},{lines[2]},{lines[6]},{lines[3]},{lines[4]},{'|'.join(lines[5])}\n")
    return

def run_enrichment_helper(study_genes, population_genes, enrich_database, save_file):
    # upload geneset
    base_url = "https://maayanlab.cloud/speedrichr"
    description = f"Rarecomb gene set enrichment in {enrich_database}"
    study_res = requests.post(
        f"{base_url}/api/addList",
        files=dict(list=(None, '\n'.join(study_genes)), description=(None, description))
        )
    if study_res.ok:
        userlist_response = study_res.json()
        userListId=userlist_response["userListId"]
    else:
        raise NotImplementedError(f"{enrich_database} enrichment not completed due to study gene posting failure")

    # upload background
    pop_res = requests.post(
        f"{base_url}/api/addbackground",
        data=dict(background='\n'.join(population_genes))
    )
    if pop_res.ok:
        background_response = pop_res.json()
        backgroundid=background_response["backgroundid"]
    else:
        raise NotImplementedError(f"{enrich_database} enrichment not completed due to pop gene posting failure")

    # run enrichment
    enrich_res = requests.post(
            f"{base_url}/api/backgroundenrich",
            data=dict(
            userListId=userListId,
            backgroundid=backgroundid,
            backgroundType=enrich_database,
            )
        )
    if enrich_res.ok:
        try:
            results = json.loads(enrich_res.text)[enrich_database]
            save_enrichr_results(results, save_file)
        except:
            print(f"Results could not be decoded for {enrich_database}")
    else:
        raise NotImplementedError(f"{enrich_database} enrichment not completed due to run enrichment failure")
    return


def run_enrichment(study_genes, population_genes, save_dir, enrich_database=["Allen_Brain_Atlas_10x_scRNA_2021", "dbGaP", "GO_Biological_Process_2023", "GTEx_Tissues_V8_2023", "GWAS_Catalog_2023", "Human_Phenotype_Ontology", "KEGG_2021_Human", "MGI_Mammalian_Phenotype_Level_4_2021", "OMIM_Disease"]):
    for ed in enrich_database:
        print(ed)
        save_file = os.path.join(save_dir, f"{ed}_enrich.csv")
        run_enrichment_helper(study_genes, population_genes, ed, save_file)
        time.sleep(30)
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs='+')
    parser.add_argument("--gene_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--save_dir", type=str, help="Filepath of the save dir for enrichment results")

    cli_args = parser.parse_args()
    combo_genes_dfs = [pd.read_csv(cgf) for cgf in cli_args.combo_files]
    pop_gene_df = pd.read_csv(cli_args.gene_file)

    study_genes = get_combo_genes(combo_genes_dfs, cli_args.save_dir)
    population_genes = get_population_geneset(pop_gene_df, cli_args.save_dir)

    run_enrichment(study_genes, population_genes, cli_args.save_dir)
