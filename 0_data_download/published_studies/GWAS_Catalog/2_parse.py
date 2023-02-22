import pandas as pd


# read GWAS catalogue
gwas_catalog = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/published_studies/GWAS_Catalog/data/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv"
columns_to_keep = ['REPORTED GENE(S)', 'DISEASE/TRAIT', 'MAPPED_TRAIT']
gwas_df = pd.read_csv(gwas_catalog, sep="\t", usecols=columns_to_keep).loc[:, columns_to_keep]
gwas_df.columns = ['Gene', 'Disease_Trait', 'Mapped_Trait']

# filter NA and unknown columns
gwas_df = gwas_df[(~gwas_df['Gene'].isna())]
gwas_df = gwas_df[(~gwas_df['Mapped_Trait'].isna())]
gwas_df = gwas_df[(gwas_df['Gene'] != 'NR')]
gwas_df = gwas_df[(gwas_df['Gene'] != '-')]

# split multi gene records
gwas_df["Gene"] = gwas_df.Gene.str.split(",")
gwas_df = gwas_df.explode("Gene")
gwas_df["Gene"] = gwas_df.Gene.str.strip()

# split multi traits record
gwas_df["Mapped_Trait"] = gwas_df.Mapped_Trait.str.split(",")
gwas_df = gwas_df.explode("Mapped_Trait")
gwas_df["Mapped_Trait"] = gwas_df.Mapped_Trait.str.strip()

# drop duplicate entries
gwas_df = gwas_df.drop_duplicates()
gwas_df.to_csv('./data/gwas_genes_to_traits.csv', index=False)
