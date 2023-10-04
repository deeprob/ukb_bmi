import os
import pandas as pd


def table_merge(burden_dfs, save_file):
    df = pd.concat(burden_dfs).groupby("gene").aggregate(lambda x: ",".join(x)).reset_index()
    df[["gene", "strand_info"]] = df.gene.str.split("-", n=1, expand=True)
    df = df.groupby("gene").aggregate({"samples": lambda x: ",".join(x), "strand_info": lambda x: ",".join([i for i in x if i])}).reset_index()
    df["samples"] = df.samples.apply(lambda x: ",".join(list(set(x.split(",")))))
    df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    chr_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/burden_tables/"
    burden_dfs = []
    for chrms in os.listdir(chr_dir):
        chrm_dir = os.path.join(chr_dir, chrms)
        for table in os.listdir(chrm_dir):
            if table.endswith(".tsv"):
                filename = os.path.join(chrm_dir, table)
                df = pd.read_csv(filename, usecols=["gene", "samples"], sep="\t")
                burden_dfs.append(df)  

    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/processed_burden/all_gene_burden.csv.gz"
    table_merge(burden_dfs, save_file)
