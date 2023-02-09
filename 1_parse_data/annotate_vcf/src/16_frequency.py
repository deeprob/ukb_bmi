#!/bin/python3


import pandas as pd

filename = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/samples.csv'
samples_df = pd.read_csv(filename)
counts_df = pd.read_csv('/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts.tsv', sep='\t', header=None, low_memory=False)
counts_df.columns = ['variant_id', 'cohort_count']


number_of_samples = samples_df.shape[0]
counts_df['cohort_frequency'] = counts_df['cohort_count'] / number_of_samples
counts_df.to_csv('/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts_frequency.tsv', sep='\t', index=False)
