#!/bin/python3



import pandas as pd



############################################
# get a list of the genes needed for model #
###########################################C

sig_combo4_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_4.csv"
sig_combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/lifestyle_white_british/data/parsed_tables/combo_3.csv"

sig_combos_df4 = pd.read_csv(sig_combo4_file)
sig_combos_df3 = pd.read_csv(sig_combo3_file)

sig_combos_df4['Combo_Name'] = sig_combos_df4['Item_1'] + ';' + sig_combos_df4['Item_2'] + ';' + sig_combos_df4['Item_3'] + ';' + sig_combos_df4['Item_4']
sig_combos_df3['Combo_Name'] = sig_combos_df3['Item_1'] + ';' + sig_combos_df3['Item_2'] + ';' + sig_combos_df3['Item_3']


sig_combos_df = sig_combos_df4.append(sig_combos_df3)

all_genes = sig_combos_df['Combo_Name'].to_list()
all_genes = [item for sublist in all_genes for item in sublist.split(';')]


##################################
# load in only the genes we need #
##################################

variant_by_gene_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/tables/wes_with_lifestyle.tsv"
columns = ['Sample_Name'] + ['Input_' + s for s in all_genes]
wdf = pd.read_csv(variant_by_gene_file, sep='\t', usecols=columns)


# rename columns
columns = [s if s=='Sample_Name' else s[len('Input_'):] for s in wdf.columns]
wdf.columns = columns

################################################
# Filter for white british samples and add BMI #
################################################

bmi_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/lifestyle_white_british/data/samples_with_residuals.csv"
samples_df = pd.read_csv(bmi_file)
samples_df = samples_df.set_index('eid')

# filter for white british samples
wdf = wdf[wdf.Sample_Name.isin(samples_df.index)]

# add bmi
wdf['bmi'] = wdf.Sample_Name.map(samples_df.bmi)


outfile = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/lifestyle_white_british/data/additive/input_table.csv"
wdf.to_csv(outfile, index=False)
