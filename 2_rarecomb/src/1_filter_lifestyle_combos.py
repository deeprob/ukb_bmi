import os
import argparse
import pandas as pd


def filter_single_gene_combos(genotype_df, combo_df, save_file):
    all_genes = set([f"Input_{gene}" for gene in genotype_df.gene.values])
    combo_df["ngenes"] = combo_df.uniq_items.apply(lambda x: len(all_genes.intersection(x.split("|"))))
    combo_df = combo_df.loc[combo_df.ngenes>1]
    combo_df.to_csv(save_file, index=False)
    return combo_df

def add_only_gene_combo_info(genotype_df, combo_df, save_file):
    all_genes = set([f"Input_{gene}" for gene in genotype_df.gene.values])
    combo_df["uniq_items"] = combo_df.uniq_items.apply(lambda x: "|".join(sorted(set(x.split("|")).intersection(all_genes))))
    combo_df.to_csv(save_file)
    return combo_df


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Filter lifestyle with single genes')
    parser.add_argument("--genotype_file", type=str, help="Filepath of the genotype file with gene burden tables")
    parser.add_argument("--combo_file", type=str, help="Filepath of the rarecomb combo output file")
    parser.add_argument("--save_dir_oligo", type=str, help="Filepath of the save directory after filtering for at least two genes")
    parser.add_argument("--save_dir_gene", type=str, help="Filepath of the save directory after filtering for gene info only")

    cli_args = parser.parse_args()

    genotype_df = pd.read_csv(cli_args.genotype_file)
    combo_df = pd.read_csv(cli_args.combo_file)

    os.makedirs(cli_args.save_dir_oligo, exist_ok=True)
    save_file = os.path.join(cli_args.save_dir_oligo, os.path.basename(cli_args.combo_file))
    combo_oligo_df = filter_single_gene_combos(genotype_df, combo_df, save_file)

    os.makedirs(cli_args.save_dir_gene, exist_ok=True)
    save_file = os.path.join(cli_args.save_dir_gene, os.path.basename(cli_args.combo_file))
    combo_oligo_w_gene_df = add_only_gene_combo_info(genotype_df, combo_oligo_df, save_file)
