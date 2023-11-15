import os
import pandas as pd
from pyrarecomb.compare_enrichment import compare_enrichment
from pyrarecomb.compare_enrichment_modifiers import compare_enrichment_modifiers
from scipy.sparse import coo_array
import argparse


def get_sparse_df(df):
    df["samples"] = df.samples.str.split(",")
    df = df.explode("samples")
    samples_index = {s:i for i,s in enumerate(df.samples.unique())}
    gene_index = {g:i for i,g in enumerate(df.gene.unique())}
    row_index = df.samples.map(samples_index)
    col_index = df.gene.map(gene_index)
    data = [1 for i in range(len(row_index))]
    # going back to dense because pandas was failing to handle 
    # sparse representation for stack operation
    sparse_arr = coo_array((data, (row_index, col_index)), shape=(len(samples_index), len(gene_index))).todense()
    df = pd.DataFrame(sparse_arr, index=samples_index.keys(), columns=[f"Input_{g}" for g in gene_index.keys()])
    return df.reset_index().rename(columns={"index": "Sample_Name"})


def create_boolean_input_df(pheno_file, geno_file, lifestyle_file, phenotype):
    pheno_df = pd.read_csv(pheno_file, usecols=["Sample_Name", f"Output_{phenotype}"])
    pheno_dict = {str(sn): bm for sn,bm in zip(pheno_df.Sample_Name, pheno_df[f"Output_{phenotype}"])}
    geno_df = pd.read_csv(geno_file, usecols=["gene", "samples"])
    geno_df = get_sparse_df(geno_df)
    samples_w_geno = geno_df.Sample_Name.isin(pheno_df.Sample_Name.astype(str))
    geno_pheno_df = geno_df.loc[samples_w_geno]
    print(len(geno_pheno_df))
    geno_pheno_df[f"Output_{phenotype}"] = geno_pheno_df.Sample_Name.map(pheno_dict)
    if lifestyle_file:
        lifestyle_df = pd.read_csv(lifestyle_file)
        lifestyle_df["Sample_Name"] = lifestyle_df["Sample_Name"].astype(str)
        lifestyle_df.columns = [f"Input_{lf}" if lf!="Sample_Name" else "Sample_Name" for lf in lifestyle_df.columns]
        lifestyles = [c for c in lifestyle_df.columns if c!="Sample_Name"]
        geno_pheno_df = geno_pheno_df.merge(lifestyle_df, left_on = "Sample_Name", right_on="Sample_Name")
        return geno_pheno_df, lifestyles
    return geno_pheno_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rarecomb pipeline.')
    parser.add_argument("pheno_file", type=str, help="Filepath of the phenotype file with binarized pheno")
    parser.add_argument("geno_file", type=str, help="Filepath of the genotype file with gene burden tables")
    parser.add_argument("save_file", type=str, help="Filepath of the combination save file")
    parser.add_argument("--phenotype", type=str, help="The phenotype binarized in phenofile", default="BMI")
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the lifestyle factors file to add as input", default="")
    parser.add_argument("--ncombo", type=int, help="Number of genes forming a combination to mine", default=2)
    parser.add_argument("--min_indv", type=int, help="Minimum number of individuals satisfying a combination", default=5)
    parser.add_argument("--max_freq", type=float, help="Maximum proportion of individuals with the combination", default=0.95)
    parser.add_argument("--adj_pval_type", type=str, help="Multiple testing method", default="BH")
    parser.add_argument("--log_dir", type=str, help="path to rarecomb log directory", default="")

    cli_args = parser.parse_args()
    os.makedirs(cli_args.log_dir, exist_ok=True)
    
    if not cli_args.lifestyle_file:
        boolean_input_df = create_boolean_input_df(cli_args.pheno_file, cli_args.geno_file, cli_args.lifestyle_file, cli_args.phenotype)
        out_df = compare_enrichment(boolean_input_df, cli_args.ncombo, cli_args.min_indv, cli_args.max_freq, adj_pval_type=cli_args.adj_pval_type, logdir=cli_args.log_dir)
    else:
        boolean_input_df, lifestyles = create_boolean_input_df(cli_args.pheno_file, cli_args.geno_file, cli_args.lifestyle_file, cli_args.phenotype)
        out_df = compare_enrichment_modifiers(boolean_input_df, cli_args.ncombo, cli_args.min_indv, cli_args.max_freq, lifestyles, adj_pval_type=cli_args.adj_pval_type, logdir=cli_args.log_dir)
    
    out_df.to_csv(cli_args.save_file, index=False)
