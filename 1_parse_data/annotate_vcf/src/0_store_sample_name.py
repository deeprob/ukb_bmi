#!/bin/python3


import pandas as pd
import os


store_dir = '/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data'
bulk_file = '/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_vcf/data/bulk/ukb48799_23151.bulk'
ukb_vcf_dir = '/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/ukb_vcf/data/vcf'

def main():
	with open(bulk_file, 'r') as f:
		bulk_file_lines = f.readlines()
	bulk_file_lines = [s.strip() for s in bulk_file_lines]

	df = pd.DataFrame(index=bulk_file_lines)
	df['Bulk_file_line'] = bulk_file_lines
	df['Sample'] = df['Bulk_file_line'].apply(lambda s: s.split(' ')[0])
	df['Last_two_digits'] = df['Sample'].apply(lambda s: s[-2:])
	df['VCF_filename'] = ukb_vcf_dir + '/' + df['Last_two_digits'] + '/' + df['Sample'] + '_23151_0_0.gz'

	outfile = os.path.join(store_dir, 'samples.csv')
	df.to_csv(outfile, index=False)

if __name__ == "__main__":
	main()
