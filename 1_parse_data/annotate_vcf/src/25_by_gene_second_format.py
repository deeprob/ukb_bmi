#!/bin/python3


# need to get the input into a table with the format
# 	Sample_Name Input_1 Input_2	Output_1	Output2
# Example input files
#	 Dropbox/Vijay/16p12/Input_Files
# Example output is here
#    Dropbox/Vijay/Freq_Itemset/Submission_Files/Genome Res Revision/Revision submission files

import pandas as pd
import sys



root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data'

infile = f'{root_dir}/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv'
outfile = f'{root_dir}/data/variants_by_gene/lof_missense_pred_freq_0.01_format2.tsv'




indf = pd.read_csv(infile, sep='\t')
indf['Sample'] = indf['Sample'].astype(int)


filename = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/samples.csv'
samples_df = pd.read_csv(filename)
samples = samples_df['Sample'].astype(int).to_list()


genes = list(set(indf['Gene']))
genes = sorted(genes)

outdf = pd.DataFrame(index=samples)
outdf['Sample'] = outdf.index.to_series()

i = 0
for gene in genes:
	if i % 1000 == 0:
		print(f'{i} {gene}')
	samples_with_variant_in_gene = indf[indf['Gene'] == gene].Sample.to_list()
	
	outdf[f'Input_{gene}'] = outdf['Sample'].isin(samples_with_variant_in_gene)
	outdf[f'Input_{gene}'] = outdf[f'Input_{gene}'].astype(int)
	
	i = i + 1




outdf.to_csv(outfile, index=False, sep='\t')



