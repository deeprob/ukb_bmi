#!/bin/python3




import pandas as pd


infile = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare.tsv'
outfile = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare_lof_missense.tsv'




# name columns
# columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'CSQ', 'CADD', 'SIFT_score', 'SIFT_converted_rankscore', 'SIFT_pred', 'SIFT4G_score', 'SIFT4G_converted_rankscore', 'SIFT4G_pred', 'LRT_score', 'LRT_converted_rankscore', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_converted_rankscore', 'MutationTaster_pred', 'MutationAssessor_score', 'MutationAssessor_rankscore', 'MutationAssessor_pred', 'FATHMM_score', 'FATHMM_converted_rankscore', 'FATHMM_pred', 'PROVEAN_score', 'PROVEAN_converted_rankscore', 'PROVEAN_pred', 'MetaSVM_score', 'MetaSVM_rankscore', 'MetaSVM_pred', 'MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred', 'MetaRNN_score', 'MetaRNN_rankscore', 'MetaRNN_pred', 'MutPred_score', 'MutPred_rankscore', 'MVP_score', 'MVP_rankscore', 'MPC_score', 'MPC_rankscore', 'PrimateAI_score', 'PrimateAI_rankscore', 'PrimateAI_pred', 'DEOGEN2_score', 'DEOGEN2_rankscore', 'DEOGEN2_pred', 'BayesDel_addAF_score', 'BayesDel_addAF_rankscore', 'BayesDel_addAF_pred', 'BayesDel_noAF_score', 'BayesDel_noAF_rankscore', 'BayesDel_noAF_pred', 'ClinPred_score', 'ClinPred_rankscore', 'ClinPred_pred', 'Interpro_domain', 'GTEx_V8_gene', 'GTEx_V8_tissue', 'Aloft_pred', 'Aloft_Confidence', 'DANN_score', 'DANN_rankscore', 'integrated_fitCons_score', 'integrated_fitCons_rankscore', 'integrated_confidence_value', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP30way_mammalian', 'phyloP30way_mammalian_rankscore', 'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons30way_mammalian', 'phastCons30way_mammalian_rankscore', 'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'Sample', 'GT', 'DP', 'AD', 'GQ', 'MIN_DP', 'PL', 'VAF']


# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|LoF|LoF_filter|LoF_flags|LoF_info


csq_fields = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','LoF','LoF_filter','LoF_flags','LoF_info']



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


fin = open(infile, 'r')
fout = open(outfile, 'w')


for line in fin:
	sline = line.strip().split('\t')
	csq = sline[6]
	
	
	# a variant may have multiple csq entries sepearated by comma
	# For example, the variant at chr1 46055887 A T
	# has the csq annotation:
	# 	T|missense_variant|MODERATE|P3R3URF-PIK3R3|ENSG00000278139|Transcript|ENST00000540385|protein_coding|7/10||||1004|987|329|N/K|aaT/aaA|||-1||SNV|HGNC|HGNC:54999,T|missense_variant|MODERATE|PIK3R3|ENSG00000117461|Transcript|ENST00000262741|protein_coding|7/10||||1526|849|283|N/K|aaT/aaA|||-1||SNV|HGNC|HGNC:8981
	split_csq = csq.split(',')
	
	# get variant consequence in each csq entry
	consequences = []
	biotypes = []
	for item in split_csq:
		consequence = item.split('|')[1]
		biotype = item.split('|')[7]
		consequences.append(consequence)
		biotypes.append(biotype)
	
	# for each consequence, label it as 'lof', 'missense' or 'other'
	mut_types = []
	for consequence in consequences:
		mut_types.append(get_one_mutation_type_for_consequence(consequence))

	
	
	# only keep csq entries where the variant is lof or missense
	passing_csq_fields = []
	for i in range(len(split_csq)):
		csq_field = split_csq[i]
		mut_type = mut_types[i]
		biotype = biotypes[i]
		# if not protein coding then skip
		if biotype != 'protein_coding':
			continue
		if mut_type == 'lof':
			passing_csq_fields.append(csq_field)
		elif mut_type == 'missense':
			passing_csq_fields.append(csq_field)
	
	# if no csq fields pass the criteria then skip
	if len(passing_csq_fields) == 0:
		continue
	
	passing_csq_fields = [s.split('|') for s in passing_csq_fields]
	
	
	# create a line that only has csq fields that pass
	new_sline = sline
	# new_sline[6] = ','.join(passing_csq_fields)
	del new_sline[6]
	for i in range(len(csq_fields)):
		values = []
		for j in range(len(passing_csq_fields)):
			values.append(passing_csq_fields[j][i])
		values = ';'.join(values)
		new_sline.append(values)
	
	
	new_line = '\t'.join(new_sline) + '\n'
	fout.write(new_line)


fin.close()
fout.close()
	
	
	



