#!/bin/python3




import pandas as pd


infile = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense.tsv'
outfile = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense_formatted.tsv'

df = pd.read_csv(infile, sep='\t', header=None)


# bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%CSQ\t%CADD_PHRED\t%SIFT_score\t%SIFT_converted_rankscore\t%SIFT_pred\t%SIFT4G_score\t%SIFT4G_converted_rankscore\t%SIFT4G_pred\t%LRT_score\t%LRT_converted_rankscore\t%LRT_pred\t%MutationTaster_score\t%MutationTaster_converted_rankscore\t%MutationTaster_pred\t%MutationAssessor_score\t%MutationAssessor_rankscore\t%MutationAssessor_pred\t%FATHMM_score\t%FATHMM_converted_rankscore\t%FATHMM_pred\t%PROVEAN_score\t%PROVEAN_converted_rankscore\t%PROVEAN_pred\t%MetaSVM_score\t%MetaSVM_rankscore\t%MetaSVM_pred\t%MetaLR_score\t%MetaLR_rankscore\t%MetaLR_pred\t%MetaRNN_score\t%MetaRNN_rankscore\t%MetaRNN_pred\t%MutPred_score\t%MutPred_rankscore\t%MVP_score\t%MVP_rankscore\t%MPC_score\t%MPC_rankscore\t%PrimateAI_score\t%PrimateAI_rankscore\t%PrimateAI_pred\t%DEOGEN2_score\t%DEOGEN2_rankscore\t%DEOGEN2_pred\t%BayesDel_addAF_score\t%BayesDel_addAF_rankscore\t%BayesDel_addAF_pred\t%BayesDel_noAF_score\t%BayesDel_noAF_rankscore\t%BayesDel_noAF_pred\t%ClinPred_score\t%ClinPred_rankscore\t%ClinPred_pred\t%Interpro_domain\t%GTEx_V8_gene\t%GTEx_V8_tissue\t%Aloft_pred\t%Aloft_Confidence\t%DANN_score\t%DANN_rankscore\t%integrated_fitCons_score\t%integrated_fitCons_rankscore\t%integrated_confidence_value\t%phyloP100way_vertebrate\t%phyloP100way_vertebrate_rankscore\t%phyloP30way_mammalian\t%phyloP30way_mammalian_rankscore\t%phastCons100way_vertebrate\t%phastCons100way_vertebrate_rankscore\t%phastCons30way_mammalian\t%phastCons30way_mammalian_rankscore\t%SiPhy_29way_logOdds\t%SiPhy_29way_logOdds_rankscore[\t${sample_name}\t%GT\t%DP\t%AD\t%GQ\t%MIN_DP\t%PL\t%VAF]\n" -o $table_file $out_file

# name columns
columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'CADD', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'SIFT4G_score', 'SIFT4G_converted_rankscore', 'SIFT4G_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'MetaRNN_score', 'MetaRNN_rankscore', 'MetaRNN_pred', 'MutPred_score', 'MutPred_rankscore', 'MVP_score', 'MVP_rankscore', 'MPC_score', 'MPC_rankscore', 'PrimateAI_score', 'PrimateAI_rankscore', 'PrimateAI_pred', 'DEOGEN2_score', 'DEOGEN2_rankscore', 'DEOGEN2_pred', 'BayesDel_addAF_score', 'BayesDel_addAF_rankscore', 'BayesDel_addAF_pred', 'BayesDel_noAF_score', 'BayesDel_noAF_rankscore', 'BayesDel_noAF_pred', 'ClinPred_score', 'ClinPred_rankscore', 'ClinPred_pred', 'Interpro_domain', 'GTEx_V8_gene', 'GTEx_V8_tissue', 'Aloft_pred', 'Aloft_Confidence', 'DANN_score', 'DANN_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_rankscore', 'integrated_confidence_value', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP30way_mammalian', 'phyloP30way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons30way_mammalian', 'phastCons30way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Sample', 'GT', 'DP', 'AD', 'GQ', 'MIN_DP', 'PL', 'VAF', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']
df.columns = columns

df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']


#======================
# do some parsing
#======================





# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
lof_variant_functions = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost']
missense_variant_classes = ['inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant']
def get_one_mutation_type_for_consequence(consequence):
	# sometime there are two consequences sepeareted by & (i.e. missense_variant&splice_region_variant)
	consequence = consequence.split('&')
	
	# lof variants functions are higher priority than missense variant classes
	for variant_class in lof_variant_functions:
		if variant_class in consequence:
			return 'lof'
	for variant_class in missense_variant_classes:
		if variant_class in consequence:
			return 'missense'
	return 'other'


def get_mutation_type(consequences):
	consequences = consequences.split(';')
	mut_types = []
	for consequence in consequences:
		mut_types.append(get_one_mutation_type_for_consequence(consequence))
	
	return ';'.join(mut_types)


df['Mut_type'] = df['Consequence'].apply(lambda s: get_mutation_type(s))


# here, see the uniq mut_types and ensure that they are all missense or lof
# if any are not lof or missense then remove
unique_mut_types = list(df.Mut_type.unique())
for m in unique_mut_types:
	if 'other' in m:
		print(m)
		print('problem')


# see if there are any duplicate variants
df['sample_variant_id'] = df['Sample'].astype(str) + '_' + df['variant_id']
print(df.shape)
print(df['sample_variant_id'].value_counts())
# there is only one duplicate variant
# >>> print(df['sample_variant_id'].value_counts())
# 2401703_chr13_30713841_G_GGT    2
# 2972357_chr2_201680989_C_T      1
# 2192469_chr17_78485738_G_A      1

# drop the duplicate variant
df = df[~df.sample_variant_id.duplicated()]

variant_counts = df['variant_id'].value_counts()
def get_cohort_count(variant_id):
	return variant_counts[variant_id]

df['cohort_count'] = df['variant_id'].apply(get_cohort_count)



# to get the frequency, I need the number of samples in the cohort
# so open a list of samples
samples_df = pd.read_csv('/data5/UK_Biobank/annotations/vep/2022_02_27/samples.csv')
number_of_samples = samples_df.shape[0]


df['cohort_frequency'] = df['cohort_count'] / number_of_samples


print(df.shape)
# filter intracohort count
df = df[df['cohort_frequency'] <= 0.01]
print(df.shape)



df.to_csv(outfile, index=False, sep='\t')







