import pandas as pd



phenotypes_data = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_general_pheno/data/ukb_phenotypes_general.csv"
lifestyle_data = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/lifestyle_v2.xlsx'
parsed_lifestyle_data = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/binarized_tables/meta/meta_pheno_table.csv'
final_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british_lifestyle_regressed/data/samples.csv"


df = pd.read_csv(phenotypes_data)
# filter for white british samples
df = df[df.genetic_ethnic_grouping == 1]
# filter for samples with BMI data
df = df[~df.bmi.isna()]



lf_df = pd.read_excel(lifestyle_data)
lf_df = lf_df[lf_df.shortlist == 'X']
lf_df = lf_df[lf_df.Num_exome_samples_with_phenotype > 190e3]


phenotypes_to_keep = lf_df.Phenotype_ID.astype(str).to_list()
def should_keep(s):
	if s == 'eid':
		return True
	s = s.split('_')[1]
	if s in phenotypes_to_keep:
		return True
	return False



lifestyle_df = pd.read_csv(parsed_lifestyle_data)
columns_keep = [s for s in lifestyle_df.columns if should_keep(s)]
lifestyle_df = lifestyle_df[columns_keep]


# remove duplicate categories to get missing counts
tmp_df = lifestyle_df.copy()
cols_keep = []
phenotypes_done = []
for column in tmp_df.columns[1:]:
	if column.split('_')[1] in phenotypes_done:
		continue
	cols_keep.append(column)
	phenotypes_done.append(column.split('_')[1])


tmp_df = tmp_df[cols_keep]
num_missing_per_sample = tmp_df.isna().sum(axis=1)

# remove samples with >60% missing values
cutoff = 0.6 * tmp_df.shape[1]
samples_should_remove = num_missing_per_sample > cutoff
lifestyle_df = lifestyle_df[~samples_should_remove]
# set the rest of the missing values to na
lifestyle_df = lifestyle_df.fillna(0)
# float to integer
for col in lifestyle_df.columns[1:]:
	lifestyle_df[col] = lifestyle_df[col].astype(int)
# intersection
df = df[df.eid.isin(lifestyle_df.eid)]
lifestyle_df = lifestyle_df[lifestyle_df.eid.isin(df.eid)]

df = df.merge(lifestyle_df, on='eid')


df.to_csv(final_samples_file, index=False)
