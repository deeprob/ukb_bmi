#!/bin/bash

chr_num=$1
notebook_path="exome_annot/annot_run/notebooks/chr${chr_num}/Annot_vep109.ipynb"

echo Annotating $chr_num

dx login --token zrKqATrAw7tfABynN0pzqtkBZBZJrRSe 

my_cmd="papermill Annot_vep109.ipynb Annot_vep109_out.ipynb"

dx run dxjupyterlab_spark_cluster \
    -ifeature="HAIL-VEP" \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=960 \
    --destination "exome_annot/annot_run/notebooks/chr${chr_num}" \
    --instance-type "mem2_ssd1_v2_x48" \
    --instance-count 1 -y
