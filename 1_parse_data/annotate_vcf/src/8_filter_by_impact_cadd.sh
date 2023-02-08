#!/bin/bash

# next filter for variants that are either HIGH impact (no CADD filter) or 
# MODERATE impact with CADD>=25
# The impact scores line up with the variant classifications we are using 

# see:
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences

# the definitions for lof and missense that I will be using:
# lof_variant_functions = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost']
# missense_variant_classes = ['inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant']

invcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd.vcf.gz
outvcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd_filtered.vcf.gz


/data5/anastasia/sw/bcftools-1.12/bcftools view -i 'CSQ[*]~"HIGH" | CADD_PHRED>=25' $invcf | bgzip > $outvcf
tabix -p vcf $outvcf
