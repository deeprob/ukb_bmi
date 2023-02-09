#!/bin/bash

sample_name=$1
sample_in=$2
out_dir=$3

tmp_file="$out_dir"/annotated.vcf.gz
out_file="$out_dir"/annotated_filtered.vcf.gz
table_file="$out_dir"/"$sample_name".tsv

annot_vcf="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/all_annot.hg38_multianno.vcf.gz"


# annotate the variants in a sample
/data5/anastasia/sw/bcftools-1.12/bcftools annotate -a $annot_vcf -c INFO -O z -o $tmp_file $sample_in
# filter the ANNOVAR unannotated samples
/data5/anastasia/sw/bcftools-1.12/bcftools view -i "CSQ!='.'" -O z -o $out_file $tmp_file
# create a table TODO: add multiple vcfs together to see how that works!
/data5/anastasia/sw/bcftools-1.12/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%CSQ\t%CADD_PHRED\t%SIFT_score\t%SIFT_converted_rankscore\t%SIFT_pred\t%SIFT4G_score\t%SIFT4G_converted_rankscore\t%SIFT4G_pred\t%LRT_score\t%LRT_converted_rankscore\t%LRT_pred\t%MutationTaster_score\t%MutationTaster_converted_rankscore\t%MutationTaster_pred\t%MutationAssessor_score\t%MutationAssessor_rankscore\t%MutationAssessor_pred\t%FATHMM_score\t%FATHMM_converted_rankscore\t%FATHMM_pred\t%PROVEAN_score\t%PROVEAN_converted_rankscore\t%PROVEAN_pred\t%MetaSVM_score\t%MetaSVM_rankscore\t%MetaSVM_pred\t%MetaLR_score\t%MetaLR_rankscore\t%MetaLR_pred\t%MetaRNN_score\t%MetaRNN_rankscore\t%MetaRNN_pred\t%MutPred_score\t%MutPred_rankscore\t%MVP_score\t%MVP_rankscore\t%MPC_score\t%MPC_rankscore\t%PrimateAI_score\t%PrimateAI_rankscore\t%PrimateAI_pred\t%DEOGEN2_score\t%DEOGEN2_rankscore\t%DEOGEN2_pred\t%BayesDel_addAF_score\t%BayesDel_addAF_rankscore\t%BayesDel_addAF_pred\t%BayesDel_noAF_score\t%BayesDel_noAF_rankscore\t%BayesDel_noAF_pred\t%ClinPred_score\t%ClinPred_rankscore\t%ClinPred_pred\t%Interpro_domain\t%GTEx_V8_gene\t%GTEx_V8_tissue\t%Aloft_pred\t%Aloft_Confidence\t%DANN_score\t%DANN_rankscore\t%integrated_fitCons_score\t%integrated_fitCons_rankscore\t%integrated_confidence_value\t%phyloP100way_vertebrate\t%phyloP100way_vertebrate_rankscore\t%phyloP30way_mammalian\t%phyloP30way_mammalian_rankscore\t%phastCons100way_vertebrate\t%phastCons100way_vertebrate_rankscore\t%phastCons30way_mammalian\t%phastCons30way_mammalian_rankscore\t%SiPhy_29way_logOdds\t%SiPhy_29way_logOdds_rankscore[\t${sample_name}\t%GT\t%DP\t%AD\t%GQ\t%MIN_DP\t%PL\t%VAF]\n" -o $table_file $out_file

# delete the tmp file
rm $tmp_file
