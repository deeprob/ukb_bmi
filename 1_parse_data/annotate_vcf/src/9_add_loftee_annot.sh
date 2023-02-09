#!/bin/bash


infile="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd.vcf.gz"
outfile="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd_loftee.vcf"

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
	-i $infile \
	-o $outfile \
    --plugin LoF,loftee_path:/data5/deepro/sw/loftee,human_ancestor_fa:/data5/deepro/sw/loftee/loftee_data/human_ancestor.fa.gz,conservation_file:/data5/deepro/sw/loftee/loftee_data/loftee.sql,gerp_bigwig:/data5/deepro/sw/loftee/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
    --dir_plugins /data5/deepro/sw/loftee
