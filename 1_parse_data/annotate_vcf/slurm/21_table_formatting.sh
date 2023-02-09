#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=meta_tab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=200G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/21_out.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/21_err.log
#SBATCH --nodelist qingyu


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

python /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/src/21_table_formatting.py

echo `date` ending job
