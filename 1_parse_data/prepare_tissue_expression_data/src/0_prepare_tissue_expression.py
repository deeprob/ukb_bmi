#!/bin/python3



import pandas as pd

gtex_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/tissue_expression/data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
gtex_parsed_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_tissue_expression_data/data/gtex_tissue_specific_expression.csv"
parsed_gencode_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_gencode_genes/data/gencode.v39.parsed.genes.csv"

df = pd.read_csv(gtex_file, skiprows=2, sep='\t')
df['Gene_id'] = df['Name'].apply(lambda s: s.split('.')[0])


all_tissues_in_gtex = ['Adipose - Subcutaneous',
	'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta',
	'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
	'Brain - Anterior cingulate cortex (BA24)',
	'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
	'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
	'Brain - Hippocampus', 'Brain - Hypothalamus',
	'Brain - Nucleus accumbens (basal ganglia)',
	'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
	'Brain - Substantia nigra', 'Breast - Mammary Tissue',
	'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes',
	'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid',
	'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
	'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
	'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex',
	'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland',
	'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary',
	'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
	'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum',
	'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
	'Whole Blood']



tissues_to_use = ['Adipose - Subcutaneous',
	'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta',
	'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
	'Brain - Anterior cingulate cortex (BA24)',
	'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
	'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
	'Brain - Hippocampus', 'Brain - Hypothalamus',
	'Brain - Nucleus accumbens (basal ganglia)',
	'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
	'Brain - Substantia nigra', 'Breast - Mammary Tissue',
	'Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes',
	'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid',
	'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
	'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
	'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex',
	'Kidney - Medulla', 'Liver', 'Lung', 'Minor Salivary Gland',
	'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary',
	'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
	'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum',
	'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina',
	'Whole Blood']


# reorganize
columns = list(df.columns)
columns = ['Gene_id'] + columns[2:-1]
df = df[columns]
df = df.drop_duplicates('Gene_id')
df = df.set_index('Gene_id', drop=True)


def get_z_score(value, mean_value, std_dev):
	if std_dev == 0:
		return 0
	return(value - mean_value) / std_dev


i = 0
num_rows = df.shape[0]
for gene, row in df.iterrows():
	if i % 1000 == 0:
		print("-", end="")
	
	i = i + 1
	mean_value = row.median()
	std_dev = row.std()
	
	zscores = row.apply(lambda s: get_z_score(s, mean_value, std_dev))
	tissue_specific = zscores.apply(lambda s: 1 if s>=2 else 0)
	df.loc[gene] = list(tissue_specific)



for column in df.columns:
	df[column] = df[column].astype(int)


# only keep the protein coding genes from gencode
gencode_df = pd.read_csv(parsed_gencode_file)
gencode_genes = set(gencode_df.loc[gencode_df.gene_type == "protein_coding", "gene_id_stripped"].values)
gencode_genes_in_df = df.index.intersection(gencode_genes)

df.loc[gencode_genes_in_df].to_csv(gtex_parsed_file, index=True)
