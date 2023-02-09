#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=meta_tab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data
#SBATCH -o /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/18_out.log
#SBATCH -e /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/logs/18_err.log
#SBATCH --nodelist qingyu
#SBATCH --array 1-1000%100


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

#Set the number of runs that each SLURM task should do
PER_TASK=201

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  #Do your stuff here
  LINE=$(sed -n "$run"p /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/files/13_smap.txt | cut -f1)
  echo $LINE
  python /data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/src/18_filter_frequency.py $LINE
done

echo `date` ending job

