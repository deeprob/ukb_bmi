#!/bin/bash
# set -ue # don't set ue because it interrrupts with conda activate

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

# activate conda environment
conda activate gseanew


contingency_file=$1
out_file=$2

Rscript /data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/src/scripts/fishers_exact.R $contingency_file $out_file
