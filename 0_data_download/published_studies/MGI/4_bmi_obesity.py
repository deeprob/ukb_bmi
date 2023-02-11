#!/bin/python3



# this script gets genes associated with obesity and bmi

import pandas as pd


df = pd.read_csv('mgi_genotype_phenotype_parsed.csv')

def is_cholesterol_lipidemia_phenotype(s):
	if s == 'submission towards male mice':
		# need to be careful of this phenotype because it shows up when searching for bmi
		return False
	if 'bmi' in s.lower():
		return True
	if 'obesity' in s.lower():
		return True
	return False



df = df[df.Phenotype.apply(is_cholesterol_lipidemia_phenotype)]

for item in df.Phenotype.unique():
	print(item)

genes = df.Human_Gene.to_list()
genes = list(set(genes))


with open('mgi_bmi_obesity_genes.txt', 'w') as f:
	for gene in genes:
		f.write(gene + '\n')















