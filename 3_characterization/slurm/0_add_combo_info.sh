#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=goatools 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir /data6/deepro/ukb_bmi/3_characterization # TODO: set dir to project dir
#SBATCH -o /data6/deepro/ukb_bmi/3_characterization/slurm/logs/0_out_%a.log # TODO: set slurm output file
#SBATCH -e /data6/deepro/ukb_bmi/3_characterization/slurm/logs/0_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona # TODO: set nodelist
#SBATCH --array 25-28

export HOME="/data6/deepro/ukb_bmi"

source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/dnanexus

echo `date` starting job on $HOSTNAME

 

python /data6/deepro/ukb_bmi/3_characterization/src/0_add_combo_info.py $LINE

echo `date` ending job on $HOSTNAME
