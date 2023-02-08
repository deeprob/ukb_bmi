import subprocess
import pandas as pd
import multiprocessing as mp
import os

root_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/'
FASTA='/data5/bx_reference/hg38/Homo_sapiens_assembly38.fasta'
df = pd.read_csv(os.path.join(root_dir, "samples.csv"))
num_samples = df.shape[0]

def filter1(in_vcf_file_path, out_vcf_file_path, log_file_path):
	bcftools_path = "/data5/anastasia/sw/bcftools-1.12/bcftools"
	# this command is first split allele, then apply the filters, then left normalize	
	command = f'{bcftools_path} norm -m-both {in_vcf_file_path} 2> {log_file_path} | {bcftools_path} view -i \'GT="alt" & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25\' 2>> {log_file_path} | {bcftools_path} view -i \'(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9\' 2>> {log_file_path} |  {bcftools_path} norm -f {FASTA} 2>> {log_file_path} > {out_vcf_file_path}'
	subprocess.run(command, shell=True)
	return

def make_intermediate_dir(sample):
	last_two_digits = str(sample)[-2:]
	intermediate_dir = f'{root_dir}/intermediate_files/{last_two_digits}/{sample}'
	command = f'mkdir -p {intermediate_dir}'
	subprocess.run(command, shell=True)
	return intermediate_dir

def make_log_dir(sample):
	last_two_digits = str(sample)[-2:]
	log_dir = f'{root_dir}/logs/intermediate_logs/{last_two_digits}/{sample}'
	command = f'mkdir -p {log_dir}'
	subprocess.run(command, shell=True)
	return log_dir

def process(i):
	sample_name = df.at[i, 'Sample']
	in_vcf_file_path = df.at[i, 'VCF_filename']
	intermediate_dir = make_intermediate_dir(sample_name)
	log_dir = make_log_dir(sample_name)
	filter1(in_vcf_file_path=in_vcf_file_path, out_vcf_file_path=f'{intermediate_dir}/{sample_name}.1_filter1.vcf', log_file_path=f'{log_dir}/{sample_name}.1_filter1.log')
	return



pool = mp.Pool(120)
pool.map(process, [s for s in range(num_samples)])
