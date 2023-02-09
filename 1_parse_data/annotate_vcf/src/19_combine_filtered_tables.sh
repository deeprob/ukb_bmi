#!/bin/bash

> /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare.tsv
ls /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/annotated_by_sample | while read i
do
	cat /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/annotated_by_sample/$i/*/*_filtered.tsv >> /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants/exonic_variants_high_impact_moderate_impact_rare.tsv
done
