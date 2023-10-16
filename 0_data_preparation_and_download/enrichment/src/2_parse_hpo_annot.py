import pandas as pd


def parse_hpo_annot_file(hpo_annot_file, save_file):
    df = pd.read_csv(hpo_annot_file, sep="\t", usecols=["gene_symbol", "hpo_id"])
    df.groupby("gene_symbol").aggregate(lambda x: ";".join(x)).iloc[1:].to_csv(save_file, index=True, sep="\t", header=False)
    return


if __name__=="__main__":
    hpo_annot_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/hpo/genes_to_phenotype.txt"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/hpo/genes_to_phenotype.annot.tab"
    parse_hpo_annot_file(hpo_annot_file, save_file)

    