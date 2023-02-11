#!/bin/python3


import pandas as pd



gwas_df = pd.read_csv('gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv', sep='\t')



# information about GWAS is here: https://www.ebi.ac.uk/gwas/docs/methods/curation
gwas_df['REPORTED GENE(S)'].value_counts()
# REPORTED GENE: Gene(s) reported by the author; "intergenic" is used to denote a reported intergenic location (or lack of gene if it appeared that gene information was sought); “NR” is used to denote that no gene location information was reported.
gwas_df['MAPPED_GENE'].value_counts()
# MAPPED GENE(S): Gene(s) mapped to the strongest SNP. If the SNP is located within a gene, that gene is listed, with multiple overlapping genes separated by “, ”. If the SNP is intergenic, the upstream and downstream genes are listed, separated by “ - ”.
gwas_df['DISEASE/TRAIT'].value_counts()
# DISEASE/TRAIT: Description of disease/trait examined in the study
gwas_df['MAPPED_TRAIT'].value_counts()
# MAPPED_TRAIT* +: Mapped Experimental Factor Ontology trait for this study




# only keep bmi
gwas_df = gwas_df[gwas_df['MAPPED_TRAIT'] == 'body mass index']


rgwas_df = gwas_df[['REPORTED GENE(S)', 'DISEASE/TRAIT', 'OR or BETA']]
rgwas_df.columns = ['Gene', 'Disease_Trait', 'oddsratio_beta']
rgwas_df = rgwas_df[(~rgwas_df['Gene'].isna())]
rgwas_df = rgwas_df[(rgwas_df['Gene'] != 'NR')]
rgwas_df = rgwas_df[(rgwas_df['Gene'] != '-')]



def split_multi_genes(df, split_char=','):
	result = []
	# split multiple genes
	for i, row in df.iterrows():
		genes = row['Gene']
		trait = row['Disease_Trait']
		or_beta = row['oddsratio_beta']
		
		genes = genes.split(split_char)
		genes = list(set(genes))
		for gene in genes:
			gene = gene.strip()
			result.append([gene, trait, or_beta])
	
	result = pd.DataFrame(result, columns=df.columns)
	return result




rgwas_df = split_multi_genes(rgwas_df)
rgwas_df = rgwas_df.drop_duplicates()


# rgwas_df.to_csv('gwas_reported_genes.csv', index=False)



mgwas_df = gwas_df[['MAPPED_GENE', 'DISEASE/TRAIT', 'OR or BETA']]
mgwas_df.columns = ['Gene', 'Disease_Trait', 'oddsratio_beta']
mgwas_df = mgwas_df[(~mgwas_df['Gene'].isna())]


mgwas_df = split_multi_genes(mgwas_df)
mgwas_df = split_multi_genes(mgwas_df, split_char=' - ')
mgwas_df = mgwas_df[mgwas_df.Gene != 'NA']
mgwas_df = mgwas_df.drop_duplicates()


# mgwas_df.to_csv('gwas_mapped_genes.csv', index=False)


gwas_df = rgwas_df.append(mgwas_df)
gwas_df = gwas_df.drop_duplicates()
gwas_df = gwas_df.set_index('Gene', drop=False)

new = []
for gene in gwas_df.Gene.unique():
	traits = gwas_df.loc[gene]['Disease_Trait']
	if type(traits) != str:
		traits = ', '.join(traits)
	new.append([gene, traits])



new = pd.DataFrame(new, columns=['Gene', 'Disease_Trait'])


new.to_csv('gwas_genes.csv', index=False)














