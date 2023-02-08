#!/bin/bash

invcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep.vcf
outvcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact.vcf

filter_vep -i $invcf \
	-o $outvcf \
	--force_overwrite \
	--only_matched \
	--filter '(BIOTYPE is protein_coding) and (IMPACT is HIGH or IMPACT is MODERATE)'
