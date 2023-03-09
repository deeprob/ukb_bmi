import pandas as pd


combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
gencode_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_gencode_genes/data/gencode.v39.parsed.genes.csv"

combo2_df = pd.read_csv(combo2_file)
combo3_df = pd.read_csv(combo3_file)
gencode_df = pd.read_csv(gencode_file)

def get_gene_info(genename):
    chrom, start, end = gencode_df.loc[gencode_df.gene_name==genename, ["Chrom", "Start", "End"]].values[0]
    return chrom, start, end

def get_unique_genes(combo_dfs, ncombos_list):
    # get list of unique genes in combos
    unique_combo_genes = set()
    for combo_df, ncombo in zip(combo_dfs, ncombos_list):
        combo_set = set(combo_df.loc[:, [f"Item_{i}_symbol" for i in range(1, ncombo + 1)]].values.flatten())
        unique_combo_genes.update(combo_set)
    return set(g for g in unique_combo_genes if g)


# gene list file
gene_list_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_geisenger_data/data/coordinates_for_target_genes.txt"
unique_genes = get_unique_genes([combo2_df, combo3_df], [2, 3])
gene_2_loc_dict = {g: get_gene_info(g) for g in unique_genes}
df = pd.DataFrame(gene_2_loc_dict).T.reset_index()
df.columns = ["Gene", "chrom", "txStart", "txEnd"]
df.sort_values("Gene").to_csv(gene_list_file, sep="\t", index=False)

# gene combo file
combo_list_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_geisenger_data/data/combo_list.txt"
combo_df = pd.concat((combo2_df.loc[:, [f"Item_{i}_symbol" for i in range(1, 3)]], combo3_df.loc[:, [f"Item_{i}_symbol" for i in range(1, 4)]]))
combo_df.to_csv(combo_list_file, index=False, header=False, sep="\t")

