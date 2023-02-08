#!/bin/bash

invcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all.vcf.gz 
outvcf=/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep.vcf


# not sure about including --per_gene flag

vep --cache \
	--offline \
	--dir_cache /data5/deepro/tmp/vep_cache \
	--cache_version 105 \
	--assembly GRCh38 \
    --format vcf \
	--gencode_basic \
	--vcf \
	--symbol \
	--biotype \
	--variant_class \
	--no_intergenic \
	--coding_only \
	--force_overwrite \
	--dir /data5/deepro/tmp/.vep \
	-i $invcf \
	-o $outvcf \
	--fork 20
