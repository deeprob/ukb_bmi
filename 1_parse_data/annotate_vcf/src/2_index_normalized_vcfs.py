import os
import pandas as pd
import subprocess
import multiprocessing as mp

root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/'
input_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/intermediate_files'
output_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/by_sample'


def process(i):
    sample_name = df.at[i, 'Sample']
    invcf = f'{input_dir}/{str(sample_name)[-2:]}/{sample_name}/{sample_name}.1_filter1.vcf'
    outvcf = f'{output_dir}/{str(sample_name)[-3:]}/{sample_name}.vcf.gz'
    command = f'cat {invcf} | bgzip > {outvcf} ; tabix -p vcf {outvcf}'
    subprocess.run(command, shell=True)
    return

for i in range(1000):
	command = f'mkdir -p {output_dir}/{i:03d}'
	subprocess.run(command, shell=True)

df = pd.read_csv(os.path.join(root_dir, "samples.csv"))
num_samples = df.shape[0]

pool = mp.Pool(100)
pool.map(process, [s for s in range(num_samples)])
