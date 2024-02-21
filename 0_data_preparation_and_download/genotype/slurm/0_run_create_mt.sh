#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukbmt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=20:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --chdir /data6/ukbiobank/exomes/annot
#SBATCH -o /data6/ukbiobank/exomes/annot/slurm/logs/out_%a.log
#SBATCH -e /data6/ukbiobank/exomes/annot/slurm/logs/err_%a.log
#SBATCH --exclude=durga,ramona,laila
#SBATCH --array 7-8


echo `date` starting job on $HOSTNAME

export HOME="/data6/ukbiobank/exomes/annot"

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data6/ukbiobank/exomes/annot/slurm/files/0_smap.txt)

bash /data6/ukbiobank/exomes/annot/src/0_run_create_mt.sh $LINE

echo `date` ending job on $HOSTNAME
