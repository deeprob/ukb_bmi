#!/bin/python3


import pandas as pd
import numpy as np


# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from matplotlib.backends.backend_pdf import PdfPages



root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf'


columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'CADD', 'SIFT_pred', 'LRT_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'MetaSVM_pred', 'MetaLR_pred', 'MetaRNN_pred', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 'BayesDel_noAF_pred', 'ClinPred_pred', 'Sample', 'GT', 'DP', 'AD', 'GQ', 'MIN_DP', 'PL', 'VAF', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'variant_id', 'Mut_type', 'sample_variant_id', 'cohort_count', 'cohort_frequency']


df = pd.read_csv(f'{root_dir}/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense_formatted.tsv', sep='\t', usecols=columns)



#================================
# plot the sum pred for a subset of variants
#================================

# pred columns we will use in prediction
pred_columns = ['SIFT_pred', 'LRT_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'MetaSVM_pred', 'MetaLR_pred', 'MetaRNN_pred', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 'BayesDel_noAF_pred', 'ClinPred_pred']



score2num = {
	'.':0,
	'T':0,
	'D':1,
}

lrt2num = {
	'.':0,
	'N':0,
	'U':0,
	'D':1,
}

ma2num = {
	'.':0,
	'L':0,
	'N':0,
	'M':0.5,
	'H':1
}

provean2num = {
	'.':0,
	'N':0,
	'D':1,
}




def cat2num(s, score_name):
	if score_name == 'LRT_pred':
		return lrt2num[s]
	if score_name == 'MutationAssessor_pred':
		return ma2num[s]
	if score_name == 'PROVEAN_pred':
		return provean2num[s]
	return score2num[s]


# convert scores to integer values
for column in pred_columns:
	df[column + '_int'] = df[column].apply(lambda s: cat2num(s, column))


pred_int_columns = [s+'_int' for s in pred_columns]


# create summed pred score
df['pred_summation'] = 0
for column in pred_int_columns:
	df['pred_summation'] = df['pred_summation'] + df[column]


# filter for variants that are either lof or pred_summation>=5
df = df[(df.pred_summation >= 5) | (df.Mut_type.apply(lambda s: 'lof' in s))]



df.to_csv(f'{root_dir}/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense_formatted_pred.tsv', sep='\t', index=False)





