import os
import argparse
from functools import reduce
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import explained_variance_score

import utils.parsing as utpa
import utils.plotting as utpl

def train_model(df, categorical_cols, numerical_cols, scaled_numerical_cols):
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
    y_pred = model.predict(X)
    return y.values, y_pred

def prepare_genotypes_with_genes(genotype_df, gene_list):
    genotype_df = genotype_df.loc[genotype_df.gene.isin(gene_list)]
    genotype_df.loc[:, "samples"] = genotype_df.samples.str.split(",")
    genotype_df = genotype_df.explode("samples")
    genotype_df = pd.crosstab(genotype_df.samples, genotype_df.gene)
    return genotype_df

def prepare_genotypes_with_combos(genotype_df, combo_genes, combo_samples):
    genotype_df = genotype_df.loc[genotype_df.gene.isin(combo_genes)]
    genotype_df.loc[:, "samples"] = genotype_df.samples.str.split(",")
    genotype_df = genotype_df.explode("samples")
    genotype_df = genotype_df.loc[genotype_df.samples.isin(combo_samples)]
    genotype_df = pd.crosstab(genotype_df.samples, genotype_df.gene)
    return genotype_df

def get_combo_info_from_files(combo_files):
    combo_genes = [utpa.get_combo_genes_from_file(cf) for cf in combo_files]
    combo_genes = reduce(lambda x,y: x.union(y), combo_genes)
    combo_samples = [utpa.get_combo_samples_from_file(cf) for cf in combo_files]
    combo_samples = reduce(lambda x,y: x.union(y), combo_samples)
    return combo_genes, combo_samples


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Variance explained")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--other_study_files", type=str, help="Filepath of the gene files from other studies", nargs="+")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    cohort_df = pd.read_csv(
        cli_args.cohort_file, 
        usecols=["sample_names", "genetic_sex", "age", "ethnic_background"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"], 
        index_col=["sample_names"]).dropna()
    cohort_df.index = cohort_df.index.astype(str)

    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    scaled_numerical_cols = ["bmi_prs"]

    y, y_pred = train_model(cohort_df, categorical_cols, numerical_cols, scaled_numerical_cols)
    prs_exp_variance = explained_variance_score(y, y_pred)
    print(prs_exp_variance)

    genotype_df = pd.read_csv(cli_args.genotype_file)

    other_study_genes = [utpa.get_gene_set_from_file(osf) for osf in cli_args.other_study_files]
    other_study_genes = reduce(lambda x,y: x.union(y), other_study_genes)
    gdf_genes = prepare_genotypes_with_genes(genotype_df, other_study_genes)
    gpdf_genes = cohort_df.merge(gdf_genes, left_index=True, right_index=True, how="left").fillna(0.)
    y, y_pred = train_model(gpdf_genes, categorical_cols, numerical_cols, scaled_numerical_cols + list(gdf_genes.columns))
    prs_rare_exp_variance = explained_variance_score(y, y_pred)
    print(prs_rare_exp_variance)


    combo_genes, combo_samples = get_combo_info_from_files(cli_args.combo_files)
    gdf_combos = prepare_genotypes_with_combos(genotype_df, combo_genes, combo_samples)
    gdf = gdf_combos.merge(gdf_genes, how="outer", left_index=True, right_index=True).fillna(0.)
    gpdf = cohort_df.merge(gdf, left_index=True, right_index=True, how="left").fillna(0.)
    y, y_pred = train_model(gpdf, categorical_cols, numerical_cols, scaled_numerical_cols + list(gdf.columns))
    prs_rare_combo_exp_variance = explained_variance_score(y, y_pred)
    print(prs_rare_combo_exp_variance)

    plot_df = pd.DataFrame(
        {"Model Variates": ["PRS", "PRS + Rare", "PRS + Rare + Combos"], 
         "Variance explained": [prs_exp_variance*100, prs_rare_exp_variance*100, prs_rare_combo_exp_variance*100]}
    )

    fig = utpl.plot_variance_explained(plot_df)
    
    utpl.save_pdf(cli_args.save_file, fig)









    

