import pandas as pd
import os
import argparse
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression


def filter_data(df, mode):
    filter_mode_dict = {
        "british": ('`ethnic_background` == "British"', '`ethnic_background` != "British"', ["genetic_sex"], ["age"] + [f"genetic_pca{i}" for i in range(1, 40)], ["bmi_prs"]),
        "british_male": ('(`ethnic_background` == "British") & (`genetic_sex` == "Male")', '(`ethnic_background` != "British") & (`genetic_sex` == "Male")', [], ["age"] + [f"genetic_pca{i}" for i in range(1, 40)], ["bmi_prs"]),
        "british_female": ('(`ethnic_background` == "British") & (`genetic_sex` == "Female")', '(`ethnic_background` != "British") & (`genetic_sex` == "Female")', [], ["age"] + [f"genetic_pca{i}" for i in range(1, 40)], ["bmi_prs"])
    }
    train_query, test_query, categorical_cols, numerical_cols, scaled_numerical_cols = filter_mode_dict[mode]
    train_df = df.query(train_query)
    test_df = df.query(test_query)
    return train_df, test_df, categorical_cols, numerical_cols, scaled_numerical_cols

def get_scaled_bmi(df, categorical_cols, numerical_cols, scaled_numerical_cols, save_file):
    # define encoders
    en = LabelEncoder()
    scaler = StandardScaler()
    # select the categorical and numerical columns
    # transform the categorical columns to integer values
    for cat_col in categorical_cols:
        df[cat_col] = en.fit_transform(df[cat_col])
    # scale the numerical columns
    df[numerical_cols] = scaler.fit_transform(df.loc[:, numerical_cols])
    # scale bmi separately
    df["bmi_scaled"] = scaler.fit_transform(df.loc[:, ["bmi"]])
    # Create the target variable (bmi_residuals) using linear regression
    X = df.loc[:, categorical_cols + numerical_cols + scaled_numerical_cols]
    y = df.loc[:, 'bmi_scaled']
    model = LinearRegression()
    model.fit(X, y)
    # save the residuals for bmi
    df['bmi_residuals'] = y - model.predict(X)
    # save file to disk
    df.to_csv(save_file, index=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rarecomb pipeline.')
    parser.add_argument("pheno_file", type=str, help="Filepath of the phenotype file")
    parser.add_argument("filter_mode", type=str, help="Filtration mode for cohort")
    parser.add_argument("save_dir", type=str, help="Filepath of the dir to save cohort data")

    cli_args = parser.parse_args()

    filtered_bmi_file = cli_args.pheno_file
    mode = cli_args.filter_mode
    df = pd.read_csv(filtered_bmi_file, index_col="sample_names")
    # select only the british samples to regress (this is the discovery cohort due to its majority - no internal biases)
    train_df, test_df, categorical_cols, numerical_cols, scaled_numerical_cols = filter_data(df, mode)
    # filter samples that do not have information on the required fields
    required_cols_for_analysis = categorical_cols + numerical_cols + scaled_numerical_cols + ["bmi"]
    train_df = train_df.loc[~train_df.loc[:, required_cols_for_analysis].isna().any(axis=1)]
    # regress and save file to directory as training cohort
    train_save_file = os.path.join(cli_args.save_dir, "train_cohort_bmi.csv.gz") 
    get_scaled_bmi(train_df, categorical_cols, numerical_cols, scaled_numerical_cols, train_save_file)
    # save the non british samples as test cohort
    test_save_file = os.path.join(cli_args.save_dir, "test_cohort_bmi.csv.gz")
    test_df.to_csv(test_save_file, index=True)
