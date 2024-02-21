import os
import pandas as pd


def create_burden(block_df):
    # select lof and deleterious missense variants
    block_df = block_df.loc[(block_df.lof==True)|(block_df.splice_lof==True)|((block_df.missense==True)&(block_df.del_score>4))]
    block_df = block_df.groupby("gene").agg({"samples": lambda x: ",".join(x)})
    return block_df


if __name__ == "__main__":
    vcfs_per_chrm = {
        "chr1": 97, "chr2": 71, "chr3": 56, "chr4": 39, "chr5": 43, "chr6": 48, 
        "chr7": 47, "chr8": 35, "chr9": 42, "chr10": 40, "chr11": 57, "chr12": 52, 
        "chr13": 18, "chr14": 30, "chr15": 34, "chr16": 47, "chr17": 56, "chr18": 16, 
        "chr19": 65, "chr20": 25, "chr21": 11, "chr22": 23, "chrX": 24, "chrY": 1
    }
    annot_table_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/exome_annot/data/annot_tables_vep109"
    burden_table_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/exome_annot/data/burden_tables/"
    for chr_num in [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]:
        print(chr_num)
        chr_file_num = vcfs_per_chrm[chr_num]
        for filei in range(chr_file_num):
            block_file = os.path.join(annot_table_dir, f"{chr_num}", f"block_{filei}.tsv.gz")
            block_df = pd.read_csv(block_file, sep="\t", index_col=0)
            block_df = create_burden(block_df)
            burden_block_file = os.path.join(burden_table_dir, f"{chr_num}", f"block_{filei}.tsv")
            os.makedirs(os.path.dirname(burden_block_file), exist_ok=True)
            block_df.to_csv(burden_block_file, sep="\t")
