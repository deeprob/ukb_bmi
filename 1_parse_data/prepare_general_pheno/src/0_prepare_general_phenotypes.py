#!/bin/python3



import pandas as pd
import numpy as np
import math



# Body composition:
# 23099	Body fat percentage
# 23105	Basal metabolic rate
# 48	Waist circumference
# 49	Hip circumference


column_keys = {
	'22006-0.0':'genetic_ethnic_grouping',
	'22001-0.0':'genetic_sex',
	'21001-0.0':'bmi1',
	'21001-1.0':'bmi2',
	'21001-2.0':'bmi3',
	'21002-0.0':'weight1',
	'21002-1.0':'weight2',
	'21002-2.0':'weight3',
	'12144-2.0':'height',
	'21003-0.0':'age_when_attended_assessment_centre1',
	'21003-1.0':'age_when_attended_assessment_centre2',
	'21003-2.0':'age_when_attended_assessment_centre3',
	'21000-0.0':'ethnic_background1',
	'21000-1.0':'ethnic_background2',
	'21000-2.0':'ethnic_background3',
	'22021-0.0':'genetic_kinship_to_other_participants',
	'31-0.0':'sex',
	'34-0.0':'year_of_birth',
	'23099-0.0':'body_fat_percentage1',
	'23099-1.0':'body_fat_percentage2',
	'23105-0.0':'basal_metabolic_rate1',
	'23105-1.0':'basal_metabolic_rate2',
	'48-0.0':'waist_circumfrence1',
	'48-1.0':'waist_circumfrence2',
	'48-2.0':'waist_circumfrence3',
	'49-0.0':'hip_circumfrence1',
	'49-1.0':'hip_circumfrence2',
	'49-2.0':'hip_circumfrence3',
	'21002-0.0':'weight1',
	'21002-1.0':'weight2',
	'21002-2.0':'weight3',
	'51-0.0': 'seated_height1',
	'51-1.0': 'seated_height2',
	'51-2.0': 'seated_height3',
	'2724-0.0': 'had_menopause1',
	'2724-1.0': 'had_menopause2',
	'2724-2.0': 'had_menopause3'
}


columns_wanted = ['eid'] + list(column_keys.keys())

ukb_filename = '/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_field_value/data/ukb30075.csv'
df = pd.read_csv(ukb_filename, encoding='unicode_escape', usecols=columns_wanted)
columns = list(df.columns)


new_columns = [column_keys[s] if s in column_keys.keys() else s for s in columns]
df.columns = new_columns

#=======================
# Remove samples with >=10 third-degree relatives
#=======================


df = df[df.genetic_kinship_to_other_participants != 10]

#=======================
# Add readable ethnicity
#=======================

ethnic_coding_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_general_pheno/data/coding1001.tsv"
coding_df = pd.read_csv(ethnic_coding_file, sep='\t')
coding_df['coding'] = coding_df['coding'].astype(int)
coding_df = coding_df.set_index('coding')


def get_coding_meaning(s):
	if s!=s:
		return ''
	return coding_df.at[int(s), 'meaning']


df['ethnic_background1'] = df['ethnic_background1'].apply(get_coding_meaning)
df['ethnic_background2'] = df['ethnic_background2'].apply(get_coding_meaning)
df['ethnic_background3'] = df['ethnic_background3'].apply(get_coding_meaning)

#=======================
# Add readable menopause
#=======================

menopause_coding_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_general_pheno/data/coding100579.tsv"

coding_df = pd.read_csv(menopause_coding_file, sep='\t')
coding_df['coding'] = coding_df['coding'].astype(int)
coding_df = coding_df.set_index('coding')


def get_coding_meaning(s):
	if s!=s:
		return ''
	return coding_df.at[int(s), 'meaning']


df['had_menopause1'] = df['had_menopause1'].apply(get_coding_meaning)
df['had_menopause2'] = df['had_menopause2'].apply(get_coding_meaning)
df['had_menopause3'] = df['had_menopause3'].apply(get_coding_meaning)


#=======================
# Make consensus column for ethnicity and menopause
#=======================


def get_consensus(row):	
	values = list(row.unique())
	values = [s for s in values if s != '']
			
	if len(values) > 1:
		return 'inconsistent'
	if len(values) == 0:
		return ''
	
	return values[0]



df['ethnic_background'] = df[['ethnic_background1', 'ethnic_background2', 'ethnic_background3']].apply(get_consensus, axis=1)
df['had_menopause'] = df[['had_menopause1', 'had_menopause2', 'had_menopause3']].apply(get_consensus, axis=1)

df = df.drop(['ethnic_background1', 'ethnic_background2', 'ethnic_background3'], axis=1)
df = df.drop(['had_menopause1', 'had_menopause2', 'had_menopause3'], axis=1)



#=======================
# get height
#=======================

# BMI = kg/m2
# field 21001 (weight) is in kgs

def get_height(bmi, weight):
	if bmi!=bmi or weight!=weight:
		return np.nan
	
	height = math.sqrt(weight / bmi)
	return height




for i in range(1, 4):
	df[f'height{i}'] = df[[f'bmi{i}', f'weight{i}']].apply(lambda x: get_height(x[0], x[1]), axis=1)




#=======================
# Get Mean value for BMI and other phenotypes
#=======================

def get_mean_value(row):
	mean = row.mean()
	return mean



df['bmi'] = df[['bmi1', 'bmi2', 'bmi3']].apply(get_mean_value, axis=1)
df['height'] = df[['height1', 'height2', 'height3']].apply(get_mean_value, axis=1)
df['age_when_attended_assessment_centre'] = df[['age_when_attended_assessment_centre1', 'age_when_attended_assessment_centre2', 'age_when_attended_assessment_centre3']].apply(get_mean_value, axis=1)
df['body_fat_percentage'] = df[['body_fat_percentage1', 'body_fat_percentage2']].apply(get_mean_value, axis=1)
df['basal_metabolic_rate'] = df[['basal_metabolic_rate1', 'basal_metabolic_rate2']].apply(get_mean_value, axis=1)
df['waist_circumfrence'] = df[['waist_circumfrence1', 'waist_circumfrence2', 'waist_circumfrence3']].apply(get_mean_value, axis=1)
df['hip_circumfrence'] = df[['hip_circumfrence1', 'hip_circumfrence2', 'hip_circumfrence3']].apply(get_mean_value, axis=1)
df['weight'] = df[['weight1', 'weight2', 'weight3']].apply(get_mean_value, axis=1)
df['seated_height'] = df[['seated_height1', 'seated_height2', 'seated_height3']].apply(get_mean_value, axis=1)

df = df.drop(['bmi1', 'bmi2', 'bmi3', 'height1', 'height2', 'height3', 'age_when_attended_assessment_centre1', 'age_when_attended_assessment_centre2', 'age_when_attended_assessment_centre3', 'body_fat_percentage1', 'body_fat_percentage2', 'basal_metabolic_rate1', 'basal_metabolic_rate2', 'waist_circumfrence1', 'waist_circumfrence2', 'waist_circumfrence3', 'hip_circumfrence1', 'hip_circumfrence2', 'hip_circumfrence3','weight1', 'weight2', 'weight3','seated_height1', 'seated_height2', 'seated_height3'], axis=1)




#=======================
# Add genetic principle components
#=======================


# need to get Genotyping_process_and_sample_QC	22009	Genetic principal components
dwl_df = pd.read_csv(ukb_filename, encoding='unicode_escape', nrows=1)
columns = list(dwl_df.columns)
columns = [s for s in dwl_df if s.startswith('22009-')]
columns_wanted = ['eid'] + columns


dwl_df = pd.read_csv(ukb_filename, encoding='unicode_escape', usecols=columns_wanted)
dwl_df = dwl_df.set_index('eid')


for i in range(1, 41):
	df[f'PC{i}'] = df['eid'].map(dwl_df[f'22009-0.{i}'])



#=======================
# only keep samples with WES
#=======================

wes_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/samples.csv"
wes_df = pd.read_csv(wes_samples_file)

df = df[df.eid.isin(wes_df.Sample)]


#=======================
# Save Table
#=======================



pheno_save_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_general_pheno/data/ukb_phenotypes_general.csv"
df.to_csv(pheno_save_file, index=False)










