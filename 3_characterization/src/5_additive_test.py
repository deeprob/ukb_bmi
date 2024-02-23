import os
import argparse
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder

import utils.parsing as utpa
import utils.plotting as utpl


def prepare_model_data(genotype_df, cohort_df, combo_genes, combo_samples):
    # only keep the combo genes
    combo_gene_df = genotype_df.loc[genotype_df.gene.isin(combo_genes)]
    # set samples as str type
    combo_gene_df.loc[:, "samples"] = combo_gene_df.samples.str.split(",")
    # explode by samples
    combo_gene_df = combo_gene_df.explode("samples")
    # create one hot encoded data
    combo_gene_df = pd.crosstab(combo_gene_df.samples, combo_gene_df.gene)
    # only keep sample that have a mutation in at least one of the combo genes
    combo_gene_df = combo_gene_df.loc[combo_gene_df.sum(axis=1)>0]
    gene_cols = list(combo_gene_df.columns)
    # add sample phenotype info
    combo_gene_df = combo_gene_df.merge(cohort_df, left_index=True, right_index=True)
    # create test data: this will contain samples with the combination
    test_combo_gene_df = combo_gene_df.loc[combo_gene_df.index.isin(combo_samples)]
    # create train data: this will contain sample without any combination
    train_combo_gene_df = combo_gene_df.loc[~combo_gene_df.index.isin(combo_samples)].dropna()
    return train_combo_gene_df, test_combo_gene_df, gene_cols


def prepare_train_test_helper(df, categorical_cols, numerical_cols, scaled_numerical_cols, encoded_categorical_cols):
    X = df.loc[:, categorical_cols + numerical_cols + scaled_numerical_cols + encoded_categorical_cols]
    y = df.loc[:, 'bmi_scaled'].values.reshape(-1,1)    
    return X,y


def prepare_train_test(train_df, test_df, categorical_cols, numerical_cols, scaled_numerical_cols, encoded_categorical_cols):
    
    for cat_col in categorical_cols:
        en = LabelEncoder()
        en.fit(train_df[cat_col])
        train_df[cat_col] = en.transform(train_df[cat_col])
        test_df[cat_col] = en.transform(test_df[cat_col])
    
    if numerical_cols:
        num_col_scaler = StandardScaler()
        num_col_scaler.fit(train_df.loc[:, numerical_cols])
        train_df[numerical_cols] = num_col_scaler.transform(train_df.loc[:, numerical_cols])
        test_df[numerical_cols] = num_col_scaler.transform(test_df.loc[:, numerical_cols])

    bmi_col_scaler = StandardScaler()
    bmi_col_scaler.fit(train_df.loc[:, ["bmi"]])
    train_df["bmi_scaled"] = bmi_col_scaler.transform(train_df.loc[:, ["bmi"]])
    test_df["bmi_scaled"] = bmi_col_scaler.transform(test_df.loc[:, ["bmi"]])
    
    X_train, y_train = prepare_train_test_helper(train_df, categorical_cols, numerical_cols, scaled_numerical_cols, encoded_categorical_cols)
    X_test, y_test = prepare_train_test_helper(test_df, categorical_cols, numerical_cols, scaled_numerical_cols, encoded_categorical_cols)

    return X_train, y_train, X_test, y_test, bmi_col_scaler

def get_model_predictions(genotype_df, cohort_df, combo_genes, combo_samples, group=""):
    train_combo_gene_df, test_combo_gene_df, gene_cols = prepare_model_data(genotype_df, cohort_df, combo_genes, combo_samples)
    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    scaled_numerical_cols = ["bmi_prs"]
    encoded_categorical_cols = gene_cols
    X_train, y_train, X_test, y_test, bmi_col_scaler = prepare_train_test(train_combo_gene_df, test_combo_gene_df, categorical_cols, numerical_cols, scaled_numerical_cols, encoded_categorical_cols)
    model = LinearRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_pred = bmi_col_scaler.inverse_transform(y_pred)
    test_combo_gene_df[f"{pheno_name}_pred"] = y_pred
    ttest_pval = ttest_ind(test_combo_gene_df[f"{pheno_name}_pred"], test_combo_gene_df[f"{pheno_name}"], alternative="less").pvalue
    if group:
        test_combo_gene_df["group"] = group
    return test_combo_gene_df, ttest_pval

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Additive model")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--cohort_files", type=str, help="Filepath of the cohort phenotype file", nargs="+")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--combo_dir", type=str, help="Filepath of the combo files dir given by Rarecomb")
    parser.add_argument("--combo_filepre", type=str, help="Prefix of the combo files for groups")
    parser.add_argument("--groups", type=str, help="The different groups to compare", nargs="+", default=[])
    parser.add_argument("--save_file", type=str, help="Filepath where additive plot will be stored")
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the lifestyle factor matrix", default="")

    cli_args = parser.parse_args()
    pheno_name = "bmi"
    
    genotype_df = pd.read_csv(cli_args.genotype_file)
    if cli_args.lifestyle_file:
        lifestyle_df = pd.read_csv(cli_args.lifestyle_file)
        cols_to_melt = list(lifestyle_df.columns)
        lifestyle_df = lifestyle_df.melt(id_vars="Sample_Name", value_vars=cols_to_melt, var_name="gene")
        lifestyle_df = lifestyle_df.loc[lifestyle_df.value>0].drop(columns="value")
        lifestyle_df = lifestyle_df.groupby("gene").agg(lambda x: ",".join(map(str, x))).reset_index().rename(columns={"Sample_Name": "samples"})
        genotype_df = pd.concat((genotype_df, lifestyle_df)).reset_index(drop=True)

    if len(cli_args.groups)==0:
        cohort_df = pd.read_csv(
            cli_args.cohort_files[0], 
            usecols=["sample_names", "genetic_sex", "age"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"], 
            index_col=["sample_names"]).dropna()
        cohort_df.index = cohort_df.index.astype(str)
        combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
        test_combo_gene_df, ttest_pval = get_model_predictions(genotype_df, cohort_df, combo_genes, combo_samples, group="")
        plot_df = test_combo_gene_df.reset_index().rename(columns={f"{pheno_name}": "Observed", f"{pheno_name}_pred": "Expected"}).melt(id_vars="index")
        fig, ax = utpl.create_additive_plot(plot_df)
        print(ttest_pval)
    else:
        assert len(cli_args.cohort_files)==len(cli_args.groups)
        test_combo_df_list = []
        ttest_pvals = []
        for group, cohort_file in zip(cli_args.groups, cli_args.cohort_files):
            cohort_df = pd.read_csv(
                cohort_file, 
                usecols=["sample_names", "genetic_sex", "age"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"], 
                index_col=["sample_names"]).dropna()
            cohort_df.index = cohort_df.index.astype(str)
            combo_files = [cf.path for cf in os.scandir(os.path.join(cli_args.combo_dir, group)) if cf.name.startswith(cli_args.combo_filepre)]
            combo_genes, combo_samples = utpa.get_combo_info_from_files(combo_files)
            test_combo_gene_df, ttest_pval = get_model_predictions(genotype_df, cohort_df, combo_genes, combo_samples, group=group)
            test_combo_df_list.append(test_combo_gene_df)
            ttest_pvals.append(ttest_pval)
        test_combo_gene_df = pd.concat(test_combo_df_list)
        test_combo_gene_df = test_combo_gene_df.melt(id_vars=["group"], value_vars=[pheno_name, f"{pheno_name}_pred"], var_name="bmi_type", value_name="BMI")
        fig, ax = utpl.plot_box_oligo(test_combo_gene_df, ttest_pvals, xvar="group", yvar="BMI", huevar="bmi_type", hue_order=[f"{pheno_name}_pred", pheno_name], order=cli_args.groups, figsize=(12, 6))

    utpl.save_pdf(cli_args.save_file, fig)
