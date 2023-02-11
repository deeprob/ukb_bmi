#!/bin/python3


import pandas as pd
import numpy as np

from pandarallel import pandarallel
pandarallel.initialize(nb_workers=20)


# this file links MGI codes with mouse genes and human homologs
mgi_g_df = pd.read_csv('HMD_HumanPhenotype.rpt', sep='\t', header=None)
mgi_g_df.columns = ['Gene', 1, 'Mouse_gene', 3, 4, 5]
mgi_g_df = mgi_g_df.set_index([3])

# this file links MGI codes with MP codes
mgi_pg_df = pd.read_csv('MGI_PhenoGenoMP.rpt', sep='\t', header=None)
mgi_pg_df.columns = [0, 1, 2, 'Phenotype_Code', 4, 'Gene_Code']


# this file links MP codes with phenotype descriptions
mgi_p_df = pd.read_csv('VOC_MammalianPhenotype.rpt', sep='\t', header=None)
mgi_p_df = mgi_p_df.set_index(0)
pdict = mgi_p_df[1].to_dict()



# split_rows with multiple gene codes
new_df = []
for i, row in mgi_pg_df.iterrows():
	phenotype = row['Phenotype_Code']
	genes = row['Gene_Code']
	for gene in genes.split('|'):
		new_df.append([phenotype, gene])

new_df = pd.DataFrame(new_df, columns=['Phenotype_Code', 'Gene_Code'])


# add phenotype information
new_df['Phenotype'] = new_df.Phenotype_Code.apply(lambda s: pdict[s] if s in pdict else np.nan)
# drop rows without phenotype
new_df = new_df[~new_df.Phenotype.isna()]


def get_mouse_gene(mgi_code):
	if mgi_code not in mgi_g_df.index:
		return np.nan
	subdf = mgi_g_df.loc[mgi_code]
	if type(subdf) == pd.Series:
		return subdf['Mouse_gene']
	mouse_genes = subdf.Mouse_gene.to_list()
	mouse_genes = list(set(mouse_genes))
	if len(mouse_genes) > 1:
		print(mgi_code, subdf)
	return ','.join(mouse_genes)


def get_human_genes(mgi_code):
	if mgi_code not in mgi_g_df.index:
		return np.nan
	subdf = mgi_g_df.loc[mgi_code]
	if type(subdf) == pd.Series:
		return subdf['Gene']
	human_genes = subdf.Gene.to_list()
	human_genes = list(set(human_genes))
	if len(human_genes) > 1:
		print(mgi_code, subdf)
	return '|'.join(human_genes)


# Add Gene information
new_df['Mouse_Gene'] = new_df.Gene_Code.parallel_apply(lambda s:get_mouse_gene(s))
# drop rows without genes
new_df = new_df[~new_df['Mouse_Gene'].isna()]
new_df['Gene'] = new_df.Gene_Code.parallel_apply(lambda s:get_human_genes(s))



# split_rows with multiple human genes
new_new_df = []
for i, row in new_df.iterrows():
	phenotype_code = row['Phenotype_Code']
	gene_code = row['Gene_Code']
	phenotype = row['Phenotype']
	mouse_gene = row['Mouse_Gene']
	genes = row['Gene']
	for gene in genes.split('|'):
		new_new_df.append([phenotype_code, gene_code, phenotype, mouse_gene, gene])

new_new_df = pd.DataFrame(new_new_df, columns=['Phenotype_Code', 'Gene_Code', 'Phenotype', 'Mouse_Gene', 'Human_Gene'])







new_new_df.to_csv('mgi_genotype_phenotype_parsed.csv', index=False)













