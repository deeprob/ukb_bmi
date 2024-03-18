import os
import argparse
import pandas as pd

import utils.parsing as utpa

def save_list_to_file(mylist, save_file):
    with open(save_file, "w") as f:
        for l in mylist:
            f.write(f"{l}\n")
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Combo info")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")
    parser.add_argument("--check_lifestyle", action="store_true")

    cli_args = parser.parse_args()

    all_genes = set(pd.read_csv(cli_args.genotype_file).gene.to_list())
    all_gene_file = os.path.join(cli_args.save_dir, "all_genes.list")
    save_list_to_file(all_genes, all_gene_file)

    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    if cli_args.check_lifestyle:
        combo_genes = combo_genes.intersection(all_genes)
    combo_gene_file = os.path.join(cli_args.save_dir, "combo_genes.list")
    save_list_to_file(combo_genes, combo_gene_file)

