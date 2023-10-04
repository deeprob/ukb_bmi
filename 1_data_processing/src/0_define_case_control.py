import pandas as pd


def create_case_controls(bmi_file):
    df = pd.read_csv(bmi_file, usecols=["sample_names", "bmi_residuals"])
    df = df.rename(columns={"sample_names": "Sample_Name"})
    # can introduce additional filters here based on required cols

    # create deciles
    df["bmi_decile"] = pd.qcut(df.bmi_residuals, q=10, labels=False)
    # select top two and bottom three
    df = df.loc[(df.bmi_decile>7)|(df.bmi_decile<3)] 
    # create output bmi column
    df["Output_BMI"] = df.bmi_decile.apply(lambda x: 1 if x>7 else 0)
    return df.loc[:, ["Sample_Name", "Output_BMI"]]


def create_case_control_df(pheno_file, save_file):
    pheno_df = create_case_controls(pheno_file)
    pheno_df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    bmi_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/phenotype/data/bmi_processed/train_cohort_bmi.csv.gz"
    save_file = "/data6/deepro/ukb_bmi/1_data_processing/data/british/case_controls.csv"

    create_case_control_df(bmi_file, save_file)
