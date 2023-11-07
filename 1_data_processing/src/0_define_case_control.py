import pandas as pd
import os
import argparse

def create_case_controls(bmi_file, case_control_mode):
    df = pd.read_csv(bmi_file, usecols=["sample_names", "bmi_residuals"])
    df = df.rename(columns={"sample_names": "Sample_Name"})
    # can introduce additional filters here based on required cols

    # create deciles
    df["bmi_decile"] = pd.qcut(df.bmi_residuals, q=10, labels=False)
    query_dict = {
        "risk": "(`bmi_decile`>7)|(`bmi_decile`<3)",
        "protective": "(`bmi_decile`>6)|((`bmi_decile`<3)&(`bmi_decile`>0))"
    }
    # select top two and bottom three
    df = df.query(query_dict[case_control_mode])
    # create output bmi column
    case_control_function_dict = {
        "risk": lambda x: 1 if x>7 else 0,
        "protective": lambda x: 0 if x>6 else 1
    }
    df["Output_BMI"] = df.bmi_decile.apply(case_control_function_dict[case_control_mode])
    return df.loc[:, ["Sample_Name", "Output_BMI"]]


def create_case_control_df(pheno_file, save_file, case_control_mode):
    pheno_df = create_case_controls(pheno_file, case_control_mode)
    pheno_df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rarecomb pipeline.')
    parser.add_argument("bmi_file", type=str, help="Filepath of the phenotype file")
    parser.add_argument("save_file", type=str, help="Filepath of the phenotype save file")
    parser.add_argument("--case_control_mode", type=str, help="either risk or protective modes", default="risk")

    cli_args = parser.parse_args()
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    create_case_control_df(cli_args.bmi_file, cli_args.save_file, cli_args.case_control_mode)
