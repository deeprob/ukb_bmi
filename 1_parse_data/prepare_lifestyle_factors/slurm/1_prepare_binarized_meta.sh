#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=pheno_binarize
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/src
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/slurm/logs/out_meta.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/slurm/logs/err_meta.log


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

lifestyle="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/lifestyle_v2.xlsx"
exome="/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_vcf/data/ukb48799.csv"
pheno_info="/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_field_info/data"
pheno_store="/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/data/binarized_tables"

python /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_lifestyle_factors/src/1_prepare_binarized_meta.py $lifestyle $exome $pheno_info $pheno_store

echo `date` ending job
