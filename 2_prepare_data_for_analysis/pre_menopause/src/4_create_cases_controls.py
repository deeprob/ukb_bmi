#!/bin/python3
import pandas as pd

corrected_file_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/pre_menopause/data/samples_with_residuals.csv"
cases_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/pre_menopause/data/cases_controls/cases.txt"
controls_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/pre_menopause/data/cases_controls/controls.txt"

df = pd.read_csv(corrected_file_samples_file)
df = df.sort_values('bmi_residuals')
df['bmi_decile'] = pd.qcut(df['bmi_residuals'], q=10)

lowest_group = df['bmi_decile'].unique()[[0,1,2]]
highest_group = df['bmi_decile'].unique()[[-2, -1]]

samples_in_lowest_group = df[df['bmi_decile'].isin(lowest_group)]['eid'].astype(str).to_list()
samples_in_highest_group = df[df['bmi_decile'].isin(highest_group)]['eid'].astype(str).to_list()


# write out cases and controls comparing high to middle


# comparing high to low
with open(cases_file, 'w') as f:
	for sample in samples_in_highest_group:
		f.write(sample+'\n')

with open(controls_file, 'w') as f:
	for sample in samples_in_lowest_group:
		f.write(sample+'\n')


