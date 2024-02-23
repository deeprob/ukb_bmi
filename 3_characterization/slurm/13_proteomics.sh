#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=goatools 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data6/deepro/ukb_bmi/3_characterization/src # TODO: set dir to project dir
#SBATCH -o /data6/deepro/ukb_bmi/3_characterization/slurm/logs/13_out_%a.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/ukb_bmi/3_characterization/slurm/logs/13_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona # TODO: set nodelist
#SBATCH --array 1-6

export HOME="/data6/deepro/ukb_bmi"

source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/dnanexus

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data6/deepro/ukb_bmi/3_characterization/slurm/files/13_smap.txt)
echo $LINE

python /data6/deepro/ukb_bmi/3_characterization/src/13_proteomics.py $LINE

echo `date` ending job on $HOSTNAME
