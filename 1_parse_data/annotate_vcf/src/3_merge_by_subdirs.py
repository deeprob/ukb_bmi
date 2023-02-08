import os
import subprocess
import multiprocessing as mp




input_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/by_sample'
output_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/tier3'


def process(i):
    bcftools_path = "/data5/anastasia/sw/bcftools-1.12/bcftools"
    # bcftools merge -m none ${input_dir}/000/*.vcf.gz > vcfs/tier3/000.vcf
    command = f'{bcftools_path} merge -m none {input_dir}/{i:03d}/*.vcf.gz | {bcftools_path} view -G | bgzip > {output_dir}/{i:03d}.vcf.gz ; tabix -p vcf {output_dir}/{i:03d}.vcf.gz'
    subprocess.run(command, shell=True)
    return

pool = mp.Pool(100)
os.makedirs(output_dir, exist_ok=True)
pool.map(process, [s for s in range(0, 1000)])
