import os
import argparse
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

import utils.plotting as utpl

def divide_into_nocomorbidity(features_df):
    no_comorbidity_df = features_df.loc[features_df.sum(axis=1)==0]
    comorbidity_df = features_df.loc[features_df.sum(axis=1)>0]
    return comorbidity_df, no_comorbidity_df

def cluster_by_comorbidity(train_features_df, test_features_df):
    train_comorbidity_df, train_no_comorbidity_df = divide_into_nocomorbidity(train_features_df)
    test_comorbidity_df, test_no_comorbidity_df = divide_into_nocomorbidity(test_features_df)
    train_no_comorbidity_df["labels"] = 0
    test_no_comorbidity_df["labels"] = 0
    no_comorbidity_df = pd.concat((train_no_comorbidity_df, test_no_comorbidity_df))

    X_train, X_test = train_comorbidity_df.values, test_comorbidity_df.values
    kmeans = KMeans(n_clusters=2, random_state=42, n_init="auto")
    y_train = kmeans.fit_predict(X_train)
    y_test = kmeans.fit_predict(X_test)
    
    train_comorbidity_df["labels"] = y_train
    test_comorbidity_df["labels"] = y_test
    comorbidity_df = pd.concat((train_comorbidity_df, test_comorbidity_df))
    comorbidity_df["labels"] = comorbidity_df.labels.replace(0, 2)
    icd_df = pd.concat((comorbidity_df, no_comorbidity_df))
    return icd_df

def prepare_polynomial_features(data_df, features, y_var, include_bias):
    feature_df = data_df.loc[:, features].dropna()
    X = feature_df.values
    poly = PolynomialFeatures(len(features), interaction_only=True, include_bias=include_bias)
    X_poly = poly.fit_transform(X)
    y = data_df.loc[data_df.index.isin(feature_df.index), y_var].values.reshape(-1,1)
    return X_poly,y, poly.get_feature_names_out(input_features=features)[1:]

def train_model_sm(X, y):
    model = sm.OLS(y, X)
    results = model.fit()
    r2 = results.rsquared
    coefs = results.params[1:]
    conf_ints = results.conf_int()[1:]
    p_vals = results.pvalues[1:]
    return r2, coefs, conf_ints, p_vals


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--icd_features_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--save_dir", type=str, help="Filepath the figures will be stored")

    cli_args = parser.parse_args()

    icd_features_df = pd.read_csv(cli_args.icd_features_file, index_col=0)
    icd_features_df.index = icd_features_df.index.astype(str)

    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi_prs", "bmi", "bmi_residuals", "carrier"], dtype={"sample_names": str})
    cohort_df["bmi_prs_decile"] = pd.qcut(cohort_df.bmi_prs, 10, labels=False)
    cohort_df["bmi_decile"] = pd.qcut(cohort_df.bmi, 10, labels=False)
    cohort_df["bmi_residuals_decile"] = pd.qcut(cohort_df.bmi_residuals, 10, labels=False)

    combo_samples = set(cohort_df.loc[cohort_df.carrier==True, "sample_names"].values)
    noncombo_samples = set(cohort_df.loc[cohort_df.carrier==False, "sample_names"].values)
    combo_features_df = icd_features_df.loc[icd_features_df.index.isin(combo_samples)]
    noncombo_features_df = icd_features_df.loc[icd_features_df.index.isin(noncombo_samples)]

    icd_df = cluster_by_comorbidity(noncombo_features_df, combo_features_df)
    cohort_info_df = cohort_df.merge(icd_df.loc[:, "labels"].to_frame(), left_on="sample_names", right_index=True)
    cohort_info_df["bmi_prs_hue"] = cohort_info_df.bmi_prs_decile.map({
        0: "low", 1: "middle", 
        2: "middle", 3: "middle", 4:"middle", 5: "middle", 6: "middle", 7:"middle",
        8: "middle", 9: "high"})
    fig = utpl.create_variate_interaction_plot(cohort_info_df)
    save_file = os.path.join(cli_args.save_dir, "combo_others_interaction.pdf")
    utpl.save_pdf(save_file, fig)

    X_poly, y, poly_features = prepare_polynomial_features(cohort_info_df, ["carrier", "bmi_prs", "labels"], "bmi", True)
    r2, coefs, conf_ints, p_vals = train_model_sm(X_poly, y)
    df_data = np.concatenate([poly_features.reshape(-1,1), coefs.reshape(-1,1), conf_ints, p_vals.reshape(-1,1)], axis=1)
    coef_df = pd.DataFrame(df_data, columns=["variables", "coefs", "ci_low", "ci_high", "p_value"])
    coef_df["variables"] = coef_df.variables.str.replace(" ", "*")
    fig = utpl.create_variate_interaction_model_plot(coef_df)
    save_file = os.path.join(cli_args.save_dir, "combo_others_interaction_coefs.pdf")
    utpl.save_pdf(save_file, fig)
