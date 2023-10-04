#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=pyrarecomb 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=500G
#SBATCH --chdir /data6/deepro/ukb_bmi/2_rarecomb # TODO: set dir to project dir
#SBATCH -o /data6/deepro/ukb_bmi/2_rarecomb/slurm/logs/0_out.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/ukb_bmi/2_rarecomb/slurm/logs/0_err.log # TODO: set slurm input file
#SBATCH --nodelist=sarah # TODO: set nodelist


export HOME="/data6/deepro/ukb_bmi"

source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/dnanexus

echo `date` starting job on $HOSTNAME

python /data6/deepro/ukb_bmi/2_rarecomb/src/0_run_rarecomb.py

echo `date` ending job on $HOSTNAME
