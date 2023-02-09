#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=combine_tables
#SBATCH -o data/logs/10_combine_tables.log
#SBATCH -e data/logs/10_combine_tables.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data5/UK_Biobank/annotations/vep/2022_03_13/

echo `date` starting job on $HOSTNAME

> data/variants/exonic_variants_high_impact_moderate_impact_rare.tsv
ls data/by_sample | while read i
do
	cat data/by_sample/$i/*/*_filtered.tsv >> data/variants/exonic_variants_high_impact_moderate_impact_rare.tsv
done

echo `date` finished




