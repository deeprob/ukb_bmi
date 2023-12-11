import pandas as pd
import numpy as np


def save_genes(save_file, genes):
    with open(save_file, "w") as f:
        for g in genes:
            f.write(f"{g}\n")
    return


if __name__ == "__main__":
    gwas_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/gwas/data/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/gwas/data/gwas_genes.txt"

    df = pd.read_csv(gwas_file, usecols=["DISEASE/TRAIT", "REPORTED GENE(S)", "MAPPED_GENE", "MAPPED_TRAIT"], sep="\t", dtype=str).fillna("")
    bmi_terms = ["body mass index", "obesity"]
    bmi_df = df.loc[df.MAPPED_TRAIT.str.lower().str.contains("|".join(bmi_terms))]

    bmi_genes = [gi for g in bmi_df.MAPPED_GENE.values for gi in g.split(" - ")]
    bmi_genes = set([gi.strip() for g in bmi_genes for gi in g.split(",") if g])


    save_genes(save_file, bmi_genes)
    