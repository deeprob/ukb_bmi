#!/bin/bash


bcftools_path="/data5/anastasia/sw/bcftools-1.12/bcftools"
sub_dir="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/tier3"
main_dir="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs"

$bcftools_path merge -m none ${sub_dir}/*.vcf.gz | bgzip > ${main_dir}/all.vcf.gz 

tabix -p vcf ${main_dir}/all.vcf.gz
