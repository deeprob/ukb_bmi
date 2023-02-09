#!/bin/python3



import pandas as pd



df = pd.read_csv('/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts_frequency.tsv', sep='\t', low_memory=False)
df = df[df.cohort_frequency <= 0.01]

def chrom2int(s):
	s = s[len('chr'):]
	if s == 'X':
		return 23
	if s == 'Y':
		return 24
	return int(s)

df['chrom'] = df['variant_id'].apply(lambda s: s.split('_')[0])
df['pos'] = df['variant_id'].apply(lambda s: int(s.split('_')[1]))
df['chrom_pos'] = df['chrom'] + '_' + df['pos'].astype(str)
df['chrom_int'] = df['chrom'].apply(chrom2int)
df['dummy'] = '.'

df = df.sort_values(['chrom_int', 'pos'])

df.to_csv('/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts_frequency_sorted.tsv', sep='\t', index=False)

