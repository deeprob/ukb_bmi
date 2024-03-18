#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=pyrarecomb 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=1000G
#SBATCH --chdir /data6/deepro/ukb_bmi/2_rarecomb # TODO: set dir to project dir
#SBATCH -o /data6/deepro/ukb_bmi/2_rarecomb/slurm/logs/1_out_%a.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/ukb_bmi/2_rarecomb/slurm/logs/1_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona # TODO: set nodelist
#SBATCH --array 1-2

export HOME="/data6/deepro/ukb_bmi"

source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/pyrarecombtest

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data6/deepro/ukb_bmi/2_rarecomb/slurm/files/1_smap.txt)

python /data6/deepro/ukb_bmi/2_rarecomb/src/1_filter_lifestyle_combos.py $LINE

echo `date` ending job on $HOSTNAME
