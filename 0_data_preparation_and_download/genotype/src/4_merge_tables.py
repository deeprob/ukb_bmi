import os
import pandas as pd



def get_proten_coding_genes(gtf_file, store_file):
    gene_df = pd.read_csv(gtf_file, usecols=[2,3,4,8], names=["product_type", "start", "end", "product_info"], header=None, sep="\t", comment="#")
    # select gene products only
    gene_df = gene_df.loc[gene_df.product_type=="gene"]
    gene_df["gene_length"] = gene_df.end-gene_df.start
    # create gene id, type and name columns
    gene_df["gene_id"] = gene_df.product_info.apply(lambda x: x.split(";")[0].strip("gene_id ").strip('"').split(".")[0])
    gene_df["gene_type"] = gene_df.product_info.apply(lambda x: x.split(";")[1].strip("gene_type ").strip('"'))
    gene_df["gene_name"] = gene_df.product_info.apply(lambda x: x.split(";")[2].strip("gene_name ").strip('"'))
    # get only protein coding genes
    gene_df = gene_df.loc[gene_df.gene_type=="protein_coding"]
    # groupby gene names
    gene_df = gene_df.groupby("gene_name").aggregate({"gene_id": lambda x: "|".join(x), "gene_length": lambda x: list(x)[0]})
    # store to disk
    gene_df.to_csv(store_file, index=True)
    return


def get_genes_from_parsed_gtf(parsed_gtf_file):
    df = pd.read_csv(parsed_gtf_file)
    return set(df.gene_name)

def table_merge(burden_dfs, parsed_gtf_file, save_file):
    df = pd.concat(burden_dfs).groupby("gene").aggregate(lambda x: ",".join(x)).reset_index()
    protein_coding_genes = get_genes_from_parsed_gtf(parsed_gtf_file)
    df = df.loc[df.gene.isin(protein_coding_genes)]
    df["samples"] = df.samples.apply(lambda x: ",".join(list(set(x.split(",")))))
    df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    chr_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/burden_tables_0001/"
    burden_dfs = []
    for chrms in os.listdir(chr_dir):
        chrm_dir = os.path.join(chr_dir, chrms)
        for table in os.listdir(chrm_dir):
            if table.endswith(".tsv"):
                filename = os.path.join(chrm_dir, table)
                df = pd.read_csv(filename, usecols=["gene", "samples"], sep="\t")
                burden_dfs.append(df)  
    gtf_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/gene_annotations/gencode.v44.basic.annotation.gtf"
    parsed_gtf_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/gene_annotations/gencode.v44.basic.annotation.parsed.csv"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/processed_burden/all_gene_burden_0001.csv.gz"
    get_proten_coding_genes(gtf_file, parsed_gtf_file)
    table_merge(burden_dfs, parsed_gtf_file, save_file)
