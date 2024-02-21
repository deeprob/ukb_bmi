import os
import argparse
import pandas as pd
import json
import requests
import time

import utils.parsing as utpa
import utils.plotting as utpl

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
        enrich_file = os.path.join(save_dir, f"{ed}_enrich.csv")
        run_enrichment_helper(study_genes, population_genes, ed, enrich_file)
        fig, ax = utpl.create_dot_plot(enrich_file, ncat=20)
        ax.set_title(f"{ed}")
        fig_save_dir = os.path.join(save_dir, "figures")
        os.makedirs(fig_save_dir, exist_ok=True)
        fig_save_file = os.path.join(fig_save_dir, f"{ed}.pdf")
        utpl.save_pdf(fig_save_file, fig)
        time.sleep(30)
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--gene_file", type=str, help="Filepath of the gene file to test for enrichment")
    parser.add_argument("--population_gene_file", type=str, help="Filepath of the population gene file to test for enrichment")
    parser.add_argument("--save_dir", type=str, help="Filepath of the save dir for enrichment results")

    cli_args = parser.parse_args()

    study_genes = utpa.get_gene_set_from_file(cli_args.gene_file)
    pop_genes = utpa.get_gene_set_from_file(cli_args.population_gene_file)

    run_enrichment(study_genes, pop_genes, cli_args.save_dir)
