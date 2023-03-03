#!/bin/python3

# This script only keeps the icd10 diagnosis columns
# in ukb30075.csv

import pandas as pd

ukb_field_value_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_field_value/data/ukb30075.csv"

# get list of columns
df = pd.read_csv(ukb_field_value_file, encoding='unicode_escape')
columns = list(df.columns)



# select columns that I want
# eid 
cols_keep = ['eid']

# I think the combination of the bottom two columns is the Health Episode Statistics (HES) data that all of the other papers are using

# 41202 stands for main diagnosis
cols_keep = cols_keep + [s for s in df.columns if s.startswith('41202-')] 
# 41204 stands for secondary diagnosis
cols_keep = cols_keep + [s for s in df.columns if s.startswith('41204-')] 



df = df[cols_keep]

store_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/ukb30075_icd10.csv"
df.to_csv(store_file, index=False)

