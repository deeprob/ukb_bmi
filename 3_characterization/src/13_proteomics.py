import os
import argparse
import statsmodels.api as sm
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import pandas as pd

import utils.parsing as utpa
import utils.plotting as utpl

def prepare_npx_data(df0, df1, df2):
    df = df0.merge(df1, on="sample_names", how="outer", suffixes=("|0", "|1"))
    df2.columns = [f"{c}|2" if c!="sample_names" else "sample_names" for c in df2.columns]
    df = df.merge(df2, on="sample_names", how="outer")
    df = df.set_index("sample_names")
    df.index = df.index.astype(str)
    prefixes = [c+"|" for c in df0.columns if c!="sample_names"]
    grouper = [next(p for p in prefixes if c.startswith(p)) for c in df.columns]
    df_mean = df.groupby(grouper, axis=1).mean()
    df_mean.columns = [c.strip("|") for c in df_mean.columns]
    return df_mean

def prepare_combo_proteomics_data(protein_df, combo_files):
    protein_set = set(protein_df.columns)
    combo_df = pd.concat([pd.read_csv(cf) for cf in combo_files])
    combo_df["genes"] = combo_df.uniq_items.str.replace("Input_", "").str.split("|")
    combo_df["proteomics_data"] = combo_df.genes.apply(lambda g: len(g) == len(set(g).intersection(protein_set)))
    combo_w_proteomics = combo_df.loc[combo_df.proteomics_data==True]
    return combo_w_proteomics

def prepare_polynomial_features(protein_df, genes, y_var, include_bias):
    protein_gene_df = protein_df.loc[:, genes].dropna()
    X = protein_gene_df.values
    poly = PolynomialFeatures(len(genes), interaction_only=True, include_bias=include_bias)
    X_poly = poly.fit_transform(X)
    y = protein_df.loc[protein_df.index.isin(protein_gene_df.index), y_var].values.reshape(-1,1)
    return X_poly,y

def train_model_sk(X, y):
    reg = LinearRegression().fit(X, y)
    r2 = reg.score(X, y)
    coefs = reg.coef_[0]
    int_conf = (np.nan, np.nan)
    int_p_val = np.nan
    return r2, coefs, int_conf, int_p_val

def train_model_sm(X, y):
    model = sm.OLS(y, X)
    results = model.fit()
    r2 = results.rsquared
    coefs = results.params[1:]
    int_conf = results.conf_int()[-1]
    int_p_val = results.pvalues[-1]
    return r2, coefs, int_conf, int_p_val

def train_protein_model(genes, protein_df, y_var, include_bias, model_type):
    model_type_dict = {"sk": train_model_sk, "sm": train_model_sm}
    X_poly, y = prepare_polynomial_features(protein_df, genes, y_var, include_bias)
    r2, coefs, int_conf, int_p_val = model_type_dict[model_type](X_poly, y)
    return r2, coefs, int_conf, int_p_val

def get_protein_coeffs(genes, protein_df, model_type, y_var="bmi_residuals", include_bias=False):
    r2, coefs, int_conf, int_p_val = train_protein_model(genes, protein_df, y_var, include_bias, model_type)
    coefs = dict(zip([f"gene{i}" for i in range(1, len(genes)+1)] + ["interaction"], coefs))
    ser = pd.Series(coefs)
    ser["r2"] = r2
    ser["ci_low"] = int_conf[0]
    ser["ci_high"] = int_conf[1]
    ser["p_val"] = int_p_val
    return ser


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Combo info")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--protein_files", type=str, help="Filepath of the protein expression", nargs="+")
    parser.add_argument("--phenotype_file", type=str, help="Filepath phenotype file for cohort")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    protein_dfs = [pd.read_csv(pf) for pf in cli_args.protein_files]
    phenotype_df = pd.read_csv(cli_args.phenotype_file, usecols=["sample_names", "bmi", "bmi_prs", "bmi_residuals"], dtype={"sample_names": str, "bmi": float, "bmi_prs": float, "bmi_residuals":float})
    
    protein_df = prepare_npx_data(protein_dfs[0], protein_dfs[1], protein_dfs[2])
    combo_w_proteomics = prepare_combo_proteomics_data(protein_df, cli_args.combo_files)
    protein_pheno_df = protein_df.merge(phenotype_df, left_index=True, right_on="sample_names")

    os.makedirs(cli_args.save_dir, exist_ok=True)
    combo_w_proteomics_info = combo_w_proteomics.merge(combo_w_proteomics.genes.apply(get_protein_coeffs, args=(protein_pheno_df, "sm", "bmi_residuals", True)), left_index=True, right_index=True)
    save_file = os.path.join(cli_args.save_dir, "combo_protein_coefs.tsv")
    combo_w_proteomics_info.loc[:, ["uniq_items", "gene1",  "gene2", "interaction", "r2", "ci_low", "ci_high", "p_val"]].to_csv(save_file, index=False)
    
    fig = utpl.create_coefs_plot(combo_w_proteomics_info, figsize=(8, 4))
    save_file = os.path.join(cli_args.save_dir, "combo_protein_coefs.pdf")
    utpl.save_pdf(save_file, fig.figure)
