#!/bin/python3




import pandas as pd


df = pd.read_csv('gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv', sep='\t')




df = df[['REPORTED GENE(S)', 'DISEASE/TRAIT', 'MAPPED_TRAIT']]
df.columns = ['Gene', 'Disease_Trait', 'Mapped_Trait']
df = df[(~df['Gene'].isna())]
df = df[(~df['Mapped_Trait'].isna())]
df = df[(df['Gene'] != 'NR')]
df = df[(df['Gene'] != '-')]


def split_multi_genes(df, split_char=','):
	result = []
	# split multiple genes
	for i, row in df.iterrows():
		genes = row['Gene']
		trait = row['Disease_Trait']
		Mapped_Trait = row['Mapped_Trait']
		
		genes = genes.split(split_char)
		genes = list(set(genes))
		for gene in genes:
			gene = gene.strip()
			result.append([gene, trait, Mapped_Trait])
	
	result = pd.DataFrame(result, columns=df.columns)
	return result




df = split_multi_genes(df)
df = df.drop_duplicates()

df.to_csv('gwas_genes_phenotypes.csv', index=False)

exit()
