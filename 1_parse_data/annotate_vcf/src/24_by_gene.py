#!/bin/python3




import sys


root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data'


infile = f'{root_dir}/data/prepared_variants/lof_missense_pred_freq_0.01.tsv'
outfile = f'{root_dir}/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv'



# columns:
# Chrom	Pos	Ref	Alt	Qual	Filter	CADD	Sample	GT	DP	AD	GQ	MIN_DP	PL	VAF	variant_id	Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	DISTANCE	STRAND	FLAGS	VARIANT_CLASS	SYMBOL_SOURCE	HGNC_ID	Mut_type	cohort_count	cohort_frequency	sample_variant_id


fin = open(infile, 'r')
fout = open(outfile, 'w')

header = fin.readline().strip().split('\t')

index_Sample = header.index('Sample')
index_Chrom = header.index('Chrom')
index_Pos = header.index('Pos')
index_Ref = header.index('Ref')
index_Alt = header.index('Alt')
index_Gene = header.index('Gene')
index_SYMBOL = header.index('SYMBOL')
index_Mut_type = header.index('Mut_type')
index_variant_id = header.index('variant_id')
index_Qual = header.index('Qual')
index_Filter = header.index('Filter')
index_GT = header.index('GT')
index_DP = header.index('DP')
index_AD = header.index('AD')
index_cohort_count = header.index('cohort_count')
index_cohort_frequency = header.index('cohort_frequency')

# write out new header
new_header = ['Sample','Chrom','Pos','Ref','Alt','Gene','SYMBOL','Mut_type','variant_id','Qual','Filter','GT','DP','AD','cohort_count','cohort_frequency']
new_header = '\t'.join(new_header) + '\n'
fout.write(new_header)



for line in fin:
	sline = line.strip().split('\t')
	
	mut_types = sline[index_Mut_type]
	if ';' not in mut_types:
		outline = [
			sline[index_Sample],
			sline[index_Chrom],
			sline[index_Pos],
			sline[index_Ref],
			sline[index_Alt],
			sline[index_Gene],
			sline[index_SYMBOL],
			sline[index_Mut_type],
			sline[index_variant_id],
			sline[index_Qual],
			sline[index_Filter],
			sline[index_GT],
			sline[index_DP],
			sline[index_AD],
			sline[index_cohort_count],
			sline[index_cohort_frequency],
		]
		outline = '\t'.join(outline) + '\n'
		fout.write(outline)
		continue
	
	mut_types = sline[index_Mut_type].split(';')
	genes = sline[index_Gene].split(';')
	symbols = sline[index_SYMBOL].split(';')
	
	if len(list(set(genes))) == 1 and len(list(set(symbols))) == 1 and len(list(set(mut_types))) == 1:
		# the annotations are all the same gene, symbol, and mutation type
		# this means that the annotations are for different transcripts
		# but we don't need that information right now
		outline = [
			sline[index_Sample],
			sline[index_Chrom],
			sline[index_Pos],
			sline[index_Ref],
			sline[index_Alt],
			genes[0],
			symbols[0],
			mut_types[0],
			sline[index_variant_id],
			sline[index_Qual],
			sline[index_Filter],
			sline[index_GT],
			sline[index_DP],
			sline[index_AD],
			sline[index_cohort_count],
			sline[index_cohort_frequency],
		]
		outline = '\t'.join(outline) + '\n'
		fout.write(outline)
		continue
	
	# for each gene make a map of genes to mut types
	# prioritzing lof mut types
	gene2mut_type = {}
	gene2symbol = {}
	for i in range(len(genes)):
		gene = genes[i]
		symbol = symbols[i]
		mut_type = mut_types[i]
		
		if gene in gene2mut_type:
			if gene2mut_type[gene] != mut_type and mut_type == 'lof':
				gene2mut_type[gene] = mut_type
		else:
			gene2mut_type[gene] = mut_type
			gene2symbol[gene] = symbol 
	
	for gene in gene2mut_type:
		mut_type = gene2mut_type[gene]
		symbol = gene2symbol[gene]
		outline = [
			sline[index_Sample],
			sline[index_Chrom],
			sline[index_Pos],
			sline[index_Ref],
			sline[index_Alt],
			gene,
			symbol,
			mut_type,
			sline[index_variant_id],
			sline[index_Qual],
			sline[index_Filter],
			sline[index_GT],
			sline[index_DP],
			sline[index_AD],
			sline[index_cohort_count],
			sline[index_cohort_frequency],
		]
		outline = '\t'.join(outline) + '\n'
		fout.write(outline)
	
		

fin.close()
fout.close()


