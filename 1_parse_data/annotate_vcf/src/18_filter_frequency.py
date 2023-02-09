#!/bin/python3



import sys


sample = sys.argv[1]


counts_file = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts_frequency_sorted.tsv'
infile = f'/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/annotated_by_sample/{sample[-2:]}/{sample}/{sample}.tsv'
outfile = f'/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/annotated_by_sample/{sample[-2:]}/{sample}/{sample}_filtered.tsv'


def chrom2int(s):
	s = s[len('chr'):]
	if s == 'X':
		return 23
	if s == 'Y':
		return 24
	return int(s)


countsin = open(counts_file, 'r')
fin = open(infile, 'r')
fout = open(outfile, 'w')

# skip first line
countsin.readline()


# two sorted files as input (countsin and fin are both sorted by chromosome and position)
current_chrom = -1
current_pos = -1
current_variant_ids = []
previous_variant_ids = []
for line in fin:
	sline = line.split('\t')
	chrom = sline[0]
	chrom = chrom2int(chrom)
	pos = int(sline[1])
	variant_id = '_'.join(sline[:4])
	
	# ++++++++++
	# if chrom == 17 and pos == 17793787:
	# 	break
	# ++++++++++
	
	
	# already checked for this position
	if variant_id in current_variant_ids:
		fout.write(line)
		continue
	
	# another check that needed to be added in after testing
	if variant_id in previous_variant_ids:
		fout.write(line)
		continue
	
	for counts_line in countsin:
		scounts_line = counts_line.split('\t')
		count_file_chrom = int(scounts_line[6])
		count_file_pos = int(scounts_line[4])
		count_file_variant_id = scounts_line[0]
		
		
		# this file passed chromosome in other file
		if count_file_chrom > chrom:
			break
		
		# this file passed position in the other file
		if count_file_chrom == chrom and count_file_pos > pos:
			break
		
		
		if count_file_chrom == current_chrom and count_file_pos == current_pos:
			current_variant_ids.append(count_file_variant_id)
		else:
			current_variant_ids = [count_file_variant_id]
			current_chrom = count_file_chrom
			current_pos = count_file_pos
		
	if variant_id in current_variant_ids:
		fout.write(line)
	
	if count_file_chrom != current_chrom or count_file_pos != current_pos:
		previous_variant_ids = current_variant_ids
		current_variant_ids = [count_file_variant_id]
		current_chrom = count_file_chrom
		current_pos = count_file_pos






fin.close()
fout.close()




