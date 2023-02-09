#!/bin/python3



import pandas as pd


root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data'


df = pd.read_csv(f'{root_dir}/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense_formatted_pred.tsv', sep='\t')


#==============================
# finalize lof and missense cadd>=25 variants
#==============================

df.to_csv(f'{root_dir}/data/prepared_variants/lof_missense_pred_freq_0.01.tsv', sep='\t', index=False)





