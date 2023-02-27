import numpy as np
import pandas as pd
from scipy.stats import kstest


wes_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv"
phenotype_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/non_white_british/data/samples.csv"


wes_df = pd.read_csv(wes_file, sep='\t')
phenotype_df = pd.read_csv(phenotype_file)

# filter for non-white british samples
# and remove sample-gene duplicates - this is necessary for mapping 
# number of samples with specific combinations 
wes_df = wes_df[wes_df.Sample.isin(phenotype_df.eid)]
wes_df = wes_df.drop_duplicates(['Sample', 'Gene']).copy()

# load combinations of length 2 and 3
combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
sig_combos_df2 = pd.read_csv(combo2_file)
sig_combos_df3 = pd.read_csv(combo3_file)
sig_combos_df2['Combo_Name'] = sig_combos_df2['Item_1'] + '_' + sig_combos_df2['Item_2']
sig_combos_df3['Combo_Name'] = sig_combos_df3['Item_1'] + '_' + sig_combos_df3['Item_2'] + '_' + sig_combos_df3['Item_3']
sig_combos_df2['Combo_Name_symbol'] = sig_combos_df2['Item_1_symbol'] + '_' + sig_combos_df2['Item_2_symbol']
sig_combos_df3['Combo_Name_symbol'] = sig_combos_df3['Item_1_symbol'] + '_' + sig_combos_df3['Item_2_symbol'] + '_' + sig_combos_df3['Item_3_symbol']
sig_combos_df = sig_combos_df2.append(sig_combos_df3)

sig_combos_df = sig_combos_df[['Phenotype', 'Combo_Name', 'Combo_Name_symbol']]
sig_combos_df = sig_combos_df.reset_index(drop=True).copy()

# filter wes_df for only genes that we need
all_genes = sig_combos_df['Combo_Name'].to_list()
all_genes = [item for sublist in all_genes for item in sublist.split('_')]
wes_df = wes_df[wes_df.Gene.isin(all_genes)].copy()

sig_combos_df['Num_individuals_with_combination'] = np.nan
sig_combos_df['Mean_bmi_carriers'] = np.nan
sig_combos_df['Mean_bmi_noncarriers'] = np.nan
sig_combos_df['ks_test_pvalue'] = np.nan
sig_combos_df['directionally_consistent'] = ''
sig_combos_df['directionally_consistent_and_signficant'] = ''

for i, row in sig_combos_df.iterrows():
    combination_genes = row['Combo_Name'].split('_')
    sub_wes_df = wes_df[wes_df.Gene.isin(combination_genes)]
    counts = sub_wes_df.Sample.value_counts()
    # mapping samples to combo is possible because 
    # genes with unique variant info is no longer present - so it's one variant per gene
    samples_with_combo = list(counts[counts >= len(combination_genes)].index)

    if len(samples_with_combo) == 0:
        sig_combos_df.at[i, 'Num_individuals_with_combination'] = len(samples_with_combo)
        sig_combos_df.at[i, 'Mean_bmi_carriers'] = np.nan
        sig_combos_df.at[i, 'Mean_bmi_noncarriers'] = np.nan
        sig_combos_df.at[i, 'ks_test_pvalue'] = np.nan
        sig_combos_df.at[i, 'directionally_consistent'] = np.nan
        sig_combos_df.at[i, 'directionally_consistent_and_signficant'] = np.nan
        continue

    bmi_with_combo = phenotype_df[phenotype_df.eid.isin(samples_with_combo)].bmi.to_list()
    bmi_without_combo = phenotype_df[~phenotype_df.eid.isin(samples_with_combo)].bmi.to_list()

    pvalue = kstest(bmi_with_combo, bmi_without_combo).pvalue

    sig_combos_df.at[i, 'Num_individuals_with_combination'] = len(samples_with_combo)
    sig_combos_df.at[i, 'Mean_bmi_carriers'] = np.mean(bmi_with_combo)
    sig_combos_df.at[i, 'Mean_bmi_noncarriers'] = np.mean(bmi_without_combo)
    sig_combos_df.at[i, 'ks_test_pvalue'] = pvalue
    sig_combos_df.at[i, 'directionally_consistent'] = np.mean(bmi_with_combo) > np.mean(bmi_without_combo)
    sig_combos_df.at[i, 'directionally_consistent_and_signficant'] = (np.mean(bmi_with_combo) > np.mean(bmi_without_combo)) and (pvalue < 0.05)


save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/validation/dir_consistency.csv"
sig_combos_df.to_csv(save_file, index=False)
