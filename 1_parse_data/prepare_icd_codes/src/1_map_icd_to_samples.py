#!/bin/python3

# This script creates a directory called icd2sample
# Each file in the directory is an ICD10 code (icd2sample/Z899.tsv)
# That has a list of all of the UKBiobank samples
# with a diagnosis for that code

import pandas as pd
import os


sample2icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/ukb30075_icd10.csv"
save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/icd2sample"

df = pd.read_csv(sample2icd_file, low_memory=False)
df = df.set_index('eid')


diagnosis2sample = {}
for column in df.columns:
	print(column)
	codes = df[column]
	codes = codes[~codes.isna()]
	for sample, code in codes.iteritems():
		if code in diagnosis2sample:
			diagnosis2sample[code].update([sample])
		else:
			# updating a set is faster than updating a list
			diagnosis2sample[code] = set([sample])

# set to list
for code in diagnosis2sample:
	diagnosis2sample[code] = list(set(diagnosis2sample[code]))

for code in diagnosis2sample.keys():
	samples_with_diagnosis = diagnosis2sample[code]
	code_save_dir = os.path.join(save_dir, f"{code[:1]}")
	os.makedirs(code_save_dir, exist_ok=True)
	with open(os.path.join(code_save_dir, f'{code}.txt'), 'w') as f:
		for sample in samples_with_diagnosis:
			f.write(f'{sample}\n')
