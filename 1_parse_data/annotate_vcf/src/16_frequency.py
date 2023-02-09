#!/bin/python3


import pandas as pd

filename = '/data5/UK_Biobank/annotations/vep/2022_02_27/samples.csv'
samples_df = pd.read_csv(filename)
counts_df = pd.read_csv('../data/counts/exonic_variants_high_impact_moderate_impact_counts.tsv', sep='\t', header=None)
counts_df.columns = ['variant_id', 'cohort_count']


number_of_samples = samples_df.shape[0]





counts_df['cohort_frequency'] = counts_df['cohort_count'] / number_of_samples



counts_df.to_csv('../data/counts/exonic_variants_high_impact_moderate_impact_counts_frequency.tsv', sep='\t', index=False)


















