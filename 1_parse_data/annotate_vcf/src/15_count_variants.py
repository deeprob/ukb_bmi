#!/bin/python3


import pandas as pd

fin = open("/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact.tsv", 'r')


counts = {}



for line in fin:
	sline = line.split('\t')
	variant_id = '_'.join(sline[:4])
	if variant_id in counts:
		counts[variant_id] = counts[variant_id] + 1
	else:
		counts[variant_id] = 1
	




fout = open("/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/counts/exonic_variants_high_impact_moderate_impact_counts.tsv", 'w')
for variant_id in counts:
	outline = f'{variant_id}\t{counts[variant_id]}\n'
	fout.write(outline)

fout.close()
fin.close()


