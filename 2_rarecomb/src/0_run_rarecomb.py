import pandas as pd
from pyrarecomb.compare_enrichment import compare_enrichment
from scipy.sparse import coo_array


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


def create_boolean_input_df(pheno_file, geno_file):
    pheno_df = pd.read_csv(pheno_file, usecols=["Sample_Name", "Output_BMI"])
    pheno_dict = {str(sn): bm for sn,bm in zip(pheno_df.Sample_Name, pheno_df.Output_BMI)}
    geno_df = pd.read_csv(geno_file, usecols=["gene", "samples"])
    geno_df = get_sparse_df(geno_df)
    samples_w_geno = geno_df.Sample_Name.isin(pheno_df.Sample_Name.astype(str))
    geno_pheno_df = geno_df.loc[samples_w_geno]
    print(len(geno_pheno_df))
    geno_pheno_df["Output_BMI"] = geno_pheno_df.Sample_Name.map(pheno_dict)
    return geno_pheno_df


if __name__ == "__main__":
    pheno_file = "/data6/deepro/ukb_bmi/1_data_processing/data/british/case_controls.csv"
    geno_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/processed_burden/all_gene_burden.csv.gz"
    save_file = "/data6/deepro/ukb_bmi/2_rarecomb/data/british/combo2.csv"

    # define all other params
    combo_length = 2
    min_indv_threshold = 5
    max_freq_threshold = 0.25

    boolean_input_df = create_boolean_input_df(pheno_file, geno_file)
    out_df = compare_enrichment(boolean_input_df, combo_length, min_indv_threshold, max_freq_threshold)

    out_df.to_csv(save_file)
