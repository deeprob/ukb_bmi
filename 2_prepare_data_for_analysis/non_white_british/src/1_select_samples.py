# for this analysis select non white british samples with whole exome sequencing and bmi data
import pandas as pd


phenotypes_data = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_general_pheno/data/ukb_phenotypes_general.csv"
final_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/non_white_british/data/samples.csv"


phenotype_df = pd.read_csv(phenotypes_data)
# filter for white british samples
phenotype_df = phenotype_df[phenotype_df.genetic_ethnic_grouping != 1]
# filter out samples without sex information (173 samples without genetic sex)
phenotype_df = phenotype_df[~phenotype_df.genetic_sex.isna()]
# filter for samples with BMI data
phenotype_df = phenotype_df[~phenotype_df.bmi.isna()]
phenotype_df.to_csv(final_samples_file, index=False)
