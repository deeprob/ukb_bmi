#!/bin/bash

invcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact.vcf
outvcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd.vcf.gz
cadd_toml_file=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/src/utils/cadd.toml

/data5/software/vcfanno_linux64.1 -p 15 $cadd_toml_file $invcf | bgzip > $outvcf
tabix -p vcf $outvcf
