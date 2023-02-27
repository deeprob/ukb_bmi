import pandas as pd

## Turcot et al
giant_table1 = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GIANT/data/Giant_Turcot_Table1.csv"
giant_table2 = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GIANT/data/Giant_Turcot_Table2.csv"
giant1_df = pd.read_csv(giant_table1)
giant2_df = pd.read_csv(giant_table2)
giant_genes = giant1_df['Gene'].to_list() + giant2_df['Gene'].to_list()
giant_genes = list(set(giant_genes))

## MGI
mgi_obesity_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/MGI/data/mgi_bmi_obesity_genes.txt"
mgi_df = pd.read_csv(mgi_obesity_table, header=None)
mgi_genes = mgi_df[0].to_list()

## GWAS 
gwas_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GWAS_Catalog/data/gwas_genes_to_traits.csv"
gwas_df = pd.read_csv(gwas_table)
gwas_bmi_obesity_genes = gwas_df.loc[gwas_df.Mapped_Trait.str.contains("body mass index|obesity", regex=True)].Gene.to_list()

## Marenne et al
marenne_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/Marenne_2020/data/marenne_genes.csv"
mar_df = pd.read_csv(marenne_table)
mar_genes = mar_df['Gene'].to_list()

## Akbari et al
akbari_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/Akbari_2021/data/akbari_genes.csv"
akb_df = pd.read_csv(akbari_table)
akb_genes = akb_df.Gene.to_list()

## Locke et al
locke_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/Locke_2015/data/Locke_genes_parsed.csv"
locke_df = pd.read_csv(locke_table)
locke_genes = locke_df.Gene.to_list()

combo_nums = [2, 3]

gencode_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_gencode_genes/data/gencode.v39.parsed.genes.csv"
gencode = pd.read_csv(gencode_file)
gencode = gencode[gencode.gene_type == 'protein_coding']
gencode = gencode.drop_duplicates('gene_id_stripped')
gencode = gencode.set_index('gene_id_stripped', drop=False)

df = pd.DataFrame()
for combo_num in combo_nums:
						
    statistics_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_female/data/parsed_tables/combo_{combo_num}.csv"
    app = pd.read_csv(statistics_file)

    # app['Comparison'] = app['Phenotype'].apply(lambda s: s.split('_')[0])
    app['Combo_num'] = combo_num

    app['Group'] = "white british"

    if combo_num == 2:
        app = app.melt(id_vars=['Phenotype', 'Group', 'Combo_num'], value_vars=['Item_1', 'Item_2'])
    if combo_num == 3:
        app = app.melt(id_vars=['Phenotype', 'Group', 'Combo_num'], value_vars=['Item_1', 'Item_2', 'Item_3'])

    app = app.drop('variable', axis=1)
    app.columns = ['Phenotype', 'Group', 'Combo_num', 'Gene']	
    app['Group'] = app['Group'] + '_' + app['Phenotype'] + '_' + app['Combo_num'].astype(str)
    app['Gene_symbol'] = app.Gene.apply(lambda s: gencode.loc[s, 'gene_name'])
    app = app[['Gene', 'Gene_symbol', 'Group']]

    df = df.append(app)





df = df.drop_duplicates()


df = df.groupby(['Gene', 'Gene_symbol']).agg(' & '.join)
df = df.reset_index(drop=False)


def compare_genes(gene, gene_trait_dict):
	if gene in gene_trait_dict.keys():
		return gene_trait_dict[gene]
	return ''

df['Giant'] = df['Gene_symbol'].isin(giant_genes)
df['gwas'] = df['Gene_symbol'].isin(gwas_bmi_obesity_genes)
df['MGI_obesity'] = df['Gene_symbol'].isin(mgi_genes)
df['Marenne'] = df['Gene_symbol'].isin(mar_genes)
df['Akbari'] = df['Gene_symbol'].isin(akb_genes)
df['Locke'] = df['Gene_symbol'].isin(locke_genes)


for col in df.columns[3:]:
	if col == 'mgi':
		continue
	print(col, df[col].sum())

df["any_previous_association"] = (df.iloc[:, 3:]==True).any(axis=1)

df.to_csv("/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_female/data/previous_associations/associations.csv", index=False)
