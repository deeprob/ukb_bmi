#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=index
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/2_out.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/2_err.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=6G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data
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

python /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/src/2_index_normalized_vcfs.py


echo `date` finished

