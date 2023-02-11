#!/bin/python3



# this script gets genes associated with lipidemia and cholesterol in MGI

import pandas as pd


df = pd.read_csv('mgi_genotype_phenotype_parsed.csv')

# From https://www.cell.com/ajhg/fulltext/S0002-9297(22)00265-8:
# "Third, we extracted 1,115 genes associated with “cholesterol” or “lipidemia” phenotypes in mouse knockouts from the Mouse Genome Informatics (MGI) database"
def is_cholesterol_lipidemia_phenotype(s):
	if 'cholesterol' in s:
		return True
	if 'lipidemia' in s:
		return True
	return False



df = df[df.Phenotype.apply(is_cholesterol_lipidemia_phenotype)]


genes = df.Human_Gene.to_list()
genes = list(set(genes))


with open('cholesterol_and_lipidemia_genes.txt', 'w') as f:
	for gene in genes:
		f.write(gene + '\n')















