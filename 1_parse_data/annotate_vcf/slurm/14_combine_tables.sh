#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=meta_tab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=40G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/14_out.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/14_err.log
#SBATCH --nodelist qingyu


echo `date` starting job on $HOSTNAME

bash /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/src/14_combine_tables.sh

echo `date` ending job
