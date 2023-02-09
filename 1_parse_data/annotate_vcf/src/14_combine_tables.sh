#!/bin/bash

root_dir="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf"
mkdir -p ${root_dir}/data/variants

> ${root_dir}/data/variants/exonic_variants_high_impact_moderate_impact.tsv
ls ${root_dir}/data/vcfs/annotated_by_sample | while read i
do
	cat ${root_dir}/data/vcfs/annotated_by_sample/$i/*/*.tsv >> ${root_dir}/data/variants/exonic_variants_high_impact_moderate_impact.tsv
done
