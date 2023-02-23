#!/bin/python3



import pandas as pd



############################################
# get a list of the genes needed for model #
###########################################C

sig_combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_female/data/parsed_tables/combo_2.csv"
sig_combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_female/data/parsed_tables/combo_3.csv"

sig_combos_df2 = pd.read_csv(sig_combo2_file)
sig_combos_df3 = pd.read_csv(sig_combo3_file)


sig_combos_df2['Combo_Name'] = sig_combos_df2['Item_1'] + ';' + sig_combos_df2['Item_2']
sig_combos_df3['Combo_Name'] = sig_combos_df3['Item_1'] + ';' + sig_combos_df3['Item_2'] + ';' + sig_combos_df3['Item_3']


sig_combos_df = sig_combos_df2.append(sig_combos_df3)

all_genes = sig_combos_df['Combo_Name'].to_list()
all_genes = [item for sublist in all_genes for item in sublist.split(';')]


##################################
# load in only the genes we need #
##################################

variant_by_gene_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01_format2.tsv"
columns = ['Sample'] + ['Input_' + s for s in all_genes]
wdf = pd.read_csv(variant_by_gene_file, sep='\t', usecols=columns)


# rename columns
columns = [s if s=='Sample' else s[len('Input_'):] for s in wdf.columns]
wdf.columns = columns

################################################
# Filter for white british samples and add BMI #
################################################

bmi_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british_female/data/samples_with_residuals.csv"
samples_df = pd.read_csv(bmi_file)
samples_df = samples_df.set_index('eid')

# filter for white british samples
wdf = wdf[wdf.Sample.isin(samples_df.index)]

# add bmi
wdf['bmi'] = wdf.Sample.map(samples_df.bmi)


outfile = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_female/data/additive/input_table.csv"
wdf.to_csv(outfile, index=False)
