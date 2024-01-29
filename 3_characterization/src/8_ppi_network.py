import requests
import pandas as pd


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
    combo_files = [
        "/data6/deepro/ukb_bmi/2_rarecomb/data/british/combo2.csv",
        "/data6/deepro/ukb_bmi/2_rarecomb/data/british/combo3.csv"
    ]

    combo_df = pd.concat([pd.read_csv(cf).loc[:, ["uniq_items"]] for cf in combo_files])
    combo_df["uniq_items"] = combo_df.uniq_items.str.replace("Input_", "").str.split("|")
    combo_df = pd.concat((combo_df, combo_df.uniq_items.apply(get_ppi_score)), axis=1)
    