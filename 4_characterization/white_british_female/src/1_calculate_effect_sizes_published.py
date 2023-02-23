import pandas as pd
from statsmodels.stats.proportion import proportion_effectsize


phdf = pd.read_csv("/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british_female/data/samples_with_residuals.csv")
wdf = pd.read_csv("/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv", sep='\t')

# add decile info
phdf['bmi_decile'] = pd.qcut(phdf['bmi_residuals'], q=10)
phdf = phdf.sort_values('bmi_residuals')
lowest_group = phdf['bmi_decile'].unique()[[0,1,2]]
highest_group = phdf['bmi_decile'].unique()[[-2, -1]]
phdf['bmi_lowest_group'] = phdf['bmi_decile'].isin(lowest_group)
phdf['bmi_highest_group'] = phdf['bmi_decile'].isin(highest_group)


#########
# Giant #
#########

giant_table1 = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GIANT/data/Giant_Turcot_Table1.csv"
giant_table2 = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GIANT/data/Giant_Turcot_Table2.csv"
giant_effect_size_save = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_female/data/effect_sizes/giant.csv"

giant1_df = pd.read_csv(giant_table1)
giant2_df = pd.read_csv(giant_table2)


high_genes = giant1_df[giant1_df.beta > 0]['Gene'].to_list()

high_bmi_samples = phdf[phdf.bmi_highest_group].eid.to_list()
low_bmi_samples = phdf[phdf.bmi_lowest_group].eid.to_list()


high_hits_df = wdf[wdf.SYMBOL.isin(high_genes)]
# drop duplicated
high_hits_df = high_hits_df.drop_duplicates(['Sample', 'SYMBOL'])


stats = []
for gene in high_genes:
	samples_with_hit = high_hits_df[high_hits_df.SYMBOL == gene].Sample.to_list()
	samples_with_hit = list(set(samples_with_hit))
	proportion_low_bmi = len(list(set(samples_with_hit) & set(low_bmi_samples))) / len(low_bmi_samples)
	proportion_high_bmi = len(list(set(samples_with_hit) & set(high_bmi_samples))) / len(high_bmi_samples)
	
	result = proportion_effectsize(proportion_high_bmi, proportion_low_bmi)
	stats.append([gene, proportion_high_bmi, proportion_low_bmi, result])


stats = pd.DataFrame(stats, columns=['Gene', 'proportion_high_bmi', 'proportion_low_bmi', 'effect_size'])

stats.to_csv(giant_effect_size_save, index=False)


#########
# Akbari #
#########


akbari_table = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/Akbari_2021/data/akbari_genes.csv"
akbari_effect_size_save = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_female/data/effect_sizes/akbari.csv"

akb_df = pd.read_csv(akbari_table)
high_genes = akb_df[akb_df.Beta > 0]['Gene'].to_list()
high_bmi_samples = phdf[phdf.bmi_highest_group].eid.to_list()
low_bmi_samples = phdf[phdf.bmi_lowest_group].eid.to_list()


high_hits_df = wdf[wdf.SYMBOL.isin(high_genes)]
# drop duplicated
high_hits_df = high_hits_df.drop_duplicates(['Sample', 'SYMBOL'])


stats = []
for gene in high_genes:
	samples_with_hit = high_hits_df[high_hits_df.SYMBOL == gene].Sample.to_list()
	samples_with_hit = list(set(samples_with_hit))
	proportion_low_bmi = len(list(set(samples_with_hit) & set(low_bmi_samples))) / len(low_bmi_samples)
	proportion_high_bmi = len(list(set(samples_with_hit) & set(high_bmi_samples))) / len(high_bmi_samples)
	
	result = proportion_effectsize(proportion_high_bmi, proportion_low_bmi)
	stats.append([gene, proportion_high_bmi, proportion_low_bmi, result])


stats = pd.DataFrame(stats, columns=['Gene', 'proportion_high_bmi', 'proportion_low_bmi', 'effect_size'])

stats.to_csv(akbari_effect_size_save, index=False)
