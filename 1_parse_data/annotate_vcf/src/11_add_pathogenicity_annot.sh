#!/bin/bash


meta_vcf="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_vep_protein_coding_and_high_moderate_impact_cadd_loftee.vcf"
annovar_db="/data5/deepro/sw/annovar/humandb"
out_pre="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_annot"

#### Using Annovar to annotate the meta vcf file for pathogenicity ####
perl /data5/deepro/sw/annovar/table_annovar.pl $meta_vcf $annovar_db -vcfinput -protocol dbnsfp42c -operation f -build hg38 -nastring . -out $out_pre

out_vcf="${out_pre}.hg38_multianno.vcf"

bgzip $out_vcf
tabix -p vcf "${out_vcf}.gz"
