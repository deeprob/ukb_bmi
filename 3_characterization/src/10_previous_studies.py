import argparse
import pandas as pd
import itertools as it
import venn

import utils.parsing as utpa
import utils.plotting as utpl


def get_genes_df(g, gl):
    df = pd.DataFrame(g, columns=["gene"])
    df["label"] = gl
    return df

def get_intersecting_num(gene_sets, prod):
    sign = [-1 if p==0 else 1 for p in prod]
    sign_dict = {idx: s for idx,s in zip(range(len(prod)), sign)}
    # get the positive signs first
    sign_dict = dict(sorted(sign_dict.items(), key=lambda item:item[1], reverse=True))
    for i, (k,v) in enumerate(sign_dict.items()):
        if i == 0:
            assert v==1
            final_set = gene_sets[k]
        else:
            if v==-1:
                final_set = final_set.difference(gene_sets[k])
            else:
                assert v==1
                final_set = final_set.intersection(gene_sets[k])
    return len(final_set)

def create_gene_sets_for_venn(genes, labels):
    genes_df = pd.concat([get_genes_df(g,gl) for g,gl in zip(genes, labels)])
    genes_dict = genes_df.groupby("label").agg({"gene": lambda x: set(x)}).to_dict()["gene"]
    gene_labels, gene_sets = list(genes_dict.keys()), list(genes_dict.values())

    subset_dict = dict()
    for prod in list(it.product([0, 1], repeat=4)):
        if all(p == 0 for p in prod):
            continue
        else:
            subset_dict["".join(list(map(str, prod)))] = get_intersecting_num(gene_sets, prod)

    fig, ax = venn.venn4(subset_dict, names=gene_labels)
    return fig


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--all_studies", type=str, help="Filepath of the genelist files from previous studies", nargs="+")
    parser.add_argument("--all_studies_labels", type=str, help="Names of the genelist files from previous studies", nargs="+")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    genes = list(map(utpa.get_gene_set_from_file, cli_args.all_studies))

    fig = create_gene_sets_for_venn(genes, cli_args.all_studies_labels)

    utpl.save_pdf(cli_args.save_file, fig)
