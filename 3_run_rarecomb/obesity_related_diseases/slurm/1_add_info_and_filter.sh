#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rarecomb
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/slurm/logs/1_out_%a.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/slurm/logs/1_err_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data
#SBATCH --exclude ramona,durga,laila
#SBATCH --array 1-84


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

conda activate ukbiobank

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/slurm/files/1_smap.txt)
echo $LINE

python /data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/src/1_add_info_and_filter.py $LINE

echo `date` ending job
