#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rarecomb
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/slurm/logs/0_out_%a.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/slurm/logs/0_err_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/data
#SBATCH --exclude ramona,durga,laila
#SBATCH --array 1-3


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate rarecomb

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/slurm/files/0_smap.txt)
echo $LINE

Rscript /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/src/0_rarecomb.R $LINE

echo `date` ending job
