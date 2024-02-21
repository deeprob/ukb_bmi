import os
import argparse
import requests
import pandas as pd

import utils.parsing as utpa

def get_ppi_score(genes):
    request_url = "https://string-db.org/api/tsv-no-header/ppi_enrichment?"
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : 9606, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    response = requests.post(request_url, data=params)
    number_of_nodes, number_of_edges, average_node_degree, local_clustering_coefficient, expected_number_of_edges, p_value = list(map(float, response.text.strip().split("\n")[0].split("\t")))
    return pd.Series({
        "number_of_nodes": number_of_nodes,
        "number_of_edges": number_of_edges,
        "average_node_degree": average_node_degree,
        "local_clustering_coefficient": local_clustering_coefficient,
        "expected_number_of_edges": expected_number_of_edges,
        "p_value": p_value 
    })


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Combo info")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    combo_df = pd.concat([pd.read_csv(cf).loc[:, ["uniq_items"]] for cf in cli_args.combo_files])
    combo_df["genes"] = combo_df.uniq_items.str.replace("Input_", "").str.split("|")
    combo_df = pd.concat((combo_df, combo_df.genes.apply(get_ppi_score)), axis=1)
    combo_df = combo_df.drop(columns=["genes"])
    os.makedirs(cli_args.save_dir, exist_ok=True)
    save_file = os.path.join(cli_args.save_dir, "combos.csv")
    combo_df.to_csv(save_file)

    combo_genes, _ = utpa.get_combo_info_from_files(cli_args.combo_files)

    all_genes_ppi_ser = get_ppi_score(combo_genes)
    save_file = os.path.join(cli_args.save_dir, "all_combo_genes.csv")
    all_genes_ppi_ser.to_csv(save_file)
