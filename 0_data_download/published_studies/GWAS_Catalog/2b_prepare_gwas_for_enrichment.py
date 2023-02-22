#!/bin/python3




import pandas as pd


df = pd.read_csv('./data/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv', sep='\t')




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



# get gene universe
genes = list(df.Gene.unique())
phenotypes = list(df.Mapped_Trait.unique())




new_df = pd.DataFrame(index=list(set(genes)))
new_df['Gene'] = new_df.index.to_series()
for phenotype in phenotypes:
	print(phenotype)
	genes_with_phenotype = df[df.Mapped_Trait == phenotype].Gene.to_list()
	new_df[phenotype] = new_df.Gene.apply(lambda s: s in genes_with_phenotype)


for phenotype in phenotypes:
	new_df[phenotype] = new_df[phenotype].astype(int)

new_df.to_csv('./data/gwas_genes_phenotypes_binary.csv', index=False)





