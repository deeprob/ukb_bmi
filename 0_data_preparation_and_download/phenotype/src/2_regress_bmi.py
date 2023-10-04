import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression


def get_scaled_bmi(df, save_file):
    # define encoders
    en = LabelEncoder()
    scaler = StandardScaler()
    # select the categorical and numerical columns
    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    # transform the categorical columns to integer values
    for cat_col in categorical_cols:
        df[cat_col] = en.fit_transform(df[cat_col])
    # scale the numerical columns
    df[numerical_cols] = scaler.fit_transform(df.loc[:, numerical_cols])
    # scale bmi separately
    df["bmi_scaled"] = scaler.fit_transform(df.loc[:, ["bmi"]])
    # Create the target variable (bmi_residuals) using linear regression
    X = df.loc[:, categorical_cols + numerical_cols]
    y = df.loc[:, 'bmi']
    model = LinearRegression()
    model.fit(X, y)
    # save the residuals for bmi
    df['bmi_residuals'] = y - model.predict(X)
    # save file to disk
    df.to_csv(save_file, index=True)
    return


if __name__ == "__main__":
    filtered_bmi_file = "/data6/deepro/bmi_project/0_data_preparation_and_download/phenotype/data/bmi_processed/filtered_bmi_info.csv.gz"
    df = pd.read_csv(filtered_bmi_file, index_col="sample_names")
    # select only the british samples to regress (this is the discovery cohort due to its majority - no internal biases)
    train_df = df.loc[df.ethnic_background=="British"]
    # filter samples that do not have information on the required fields
    required_cols_for_analysis = ["genetic_sex", "age", "bmi"] + [f"genetic_pca{i}" for i in range(1, 40)]
    train_df = train_df.loc[~train_df.loc[:, required_cols_for_analysis].isna().any(axis=1)]
    # regress and save file to directory as training cohort
    train_save_file = "/data6/deepro/bmi_project/0_data_preparation_and_download/phenotype/data/bmi_processed/train_cohort_bmi.csv.gz"
    get_scaled_bmi(train_df, train_save_file)
    # save the non british samples as test cohort
    test_df = df.loc[df.ethnic_background!="British"]
    test_save_file = "/data6/deepro/bmi_project/0_data_preparation_and_download/phenotype/data/bmi_processed/test_cohort_bmi.csv.gz"
    test_df.to_csv(test_save_file, index=True)
