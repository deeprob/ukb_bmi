#!/bin/python3

import pandas as pd
import sys


rarecomb_outfile = sys.argv[1]
save_file = sys.argv[2]
phenotype = sys.argv[3]
combos = sys.argv[4]


filename = '/data5/bx_reference/hg38/annotations/gene_annotations/GENCODE39/gencode.v39.parsed.genes.tsv'
gencode = pd.read_csv(filename, sep='\t')
gencode['gene_id_stripped'] = gencode['gene_id'].apply(lambda s: s.split('.')[0])
gencode = gencode.drop_duplicates('gene_id_stripped')
gencode = gencode.set_index('gene_id_stripped', drop=False)


def convert_item_name(item):
	if item.startswith('Input_ENSG'):
		# this is a gene
		item = item[len('Input_'):]
		item = gencode.at[item, 'gene_name']
	else:
		# this is a lifestyle factor, return as it
		item = item[len('Input_'):]
	return 'Input_' + item


df = pd.read_csv(rarecomb_outfile)
df['Phenotype'] = phenotype

for i in range(1, int(combos)+1):
	df[f'Item_{i}_symbol'] = df[f'Item_{i}'].apply(convert_item_name)


if 'Effect_Size' in df.columns:
	# rearrange columns
	columns = ['Phenotype'] + [s for s in df.columns if s.endswith('_symbol')] + [s for s in df.columns if (not s.endswith('_symbol')) and (not s == 'Phenotype')]
	df = df[columns]
else:
	# if Effect_Size not in columns, this means there were no significant combos
	# drop all
	df.drop(df.index, inplace=True)

# filter
if len(df)>0:
	df = df[df.Control_pvalue_more > 0.05]

	columns_to_change = [s for s in df.columns if s.startswith('Item_')]
	for column in columns_to_change:
		df[column] = df[column].apply(lambda s: s[len('Input_'):])

	df.to_csv(save_file, index=False)
