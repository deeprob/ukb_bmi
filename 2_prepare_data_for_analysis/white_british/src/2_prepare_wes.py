#!/bin/python3


import pandas as pd


final_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british/data/samples.csv"

samples_df = pd.read_csv(final_samples_file)
samples_to_use = samples_df['eid'].astype(str).to_list()


var_by_gene_file = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01_format2.tsv'
parsed_wes_file = '/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british/data/tables/wes.tsv'




fin = open(var_by_gene_file, 'r')
fout = open(parsed_wes_file, 'w')

# read and write out header with one change
line = fin.readline()
sline = line.split('\t')
sline[0] = 'Sample_Name'
new_line = '\t'.join(sline)
fout.write(new_line)


# read in every line and write out
# if it is a sample we are using for this analysis
for line in fin:
	sample = line.split('\t')[0]
	if sample in samples_to_use:
		fout.write(line)


fin.close()
fout.close()
