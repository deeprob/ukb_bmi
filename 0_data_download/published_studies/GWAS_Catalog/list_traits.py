#!/bin/python3


import pandas as pd



df = pd.read_csv('gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv', sep='\t')



df = df[['DISEASE/TRAIT', 'MAPPED_TRAIT']]
df.columns = ['Disease_Trait', 'mapped_trait']

df = df.drop_duplicates()
df = df.sort_values('mapped_trait')



df.to_csv('traits.csv', index=False)


