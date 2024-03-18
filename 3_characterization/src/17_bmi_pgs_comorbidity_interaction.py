import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from scipy.stats.contingency import odds_ratio
from scipy import stats
from functools import reduce

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

from sklearn.cluster import KMeans
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

import utils.parsing as utpa
import utils.plotting as utpl


def get_scaled_bmi(df, categorical_cols, numerical_cols, scaled_numerical_cols):
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
    return df

def create_bmi_res_prs_decile_data(phenotype_df):
    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    scaled_numerical_cols = []#["bmi_prs"]

    phenotype_df = get_scaled_bmi(phenotype_df, categorical_cols, numerical_cols, scaled_numerical_cols)
    phenotype_df["bmi_res_decile"] = pd.qcut(phenotype_df.bmi_residuals, q=10)
    phenotype_df["bmi_res_decile_num"] = pd.qcut(phenotype_df.bmi_residuals, q=10, labels=False)
    phenotype_df["bmi_prs_decile"] = pd.qcut(phenotype_df.bmi_prs, q=10)
    phenotype_df["bmi_prs_decile_num"] = pd.qcut(phenotype_df.bmi_prs, q=10, labels=False)
    phenotype_df["bmi_res_categories"] = phenotype_df.bmi_res_decile_num.map({
        0: "underweight", 
        1:"normal", 2:"normal", 
        3:"overweight", 4:"overweight", 5:"overweight", 6:"overweight", 7:"overweight",
        8:"obese", 9:"severe obesity"
        })
    phenotype_df["bmi_prs_categories"] = phenotype_df.bmi_prs_decile_num.map({
        0: "lowest", 
        1:"middle", 2:"middle", 
        3:"middle", 4:"middle", 5:"middle", 6:"middle", 7:"middle",
        8:"middle", 9:"highest"
        })
    mean_bmi_dict = phenotype_df.groupby("bmi_res_categories")["bmi"].mean().to_dict()
    return phenotype_df, mean_bmi_dict

def create_icd_tree(icd_raw_dir, icd_codes_file):
    icd_samples_df = utpa.create_icd_samples_file(icd_raw_dir)
    icd_codes_df = pd.read_csv(icd_codes_file, usecols=["coding", "meaning", "node_id", "parent_id"], sep="\t")
    icd_codes_df["coding"] = icd_codes_df.coding.str.replace(" ", "")
    pheno_tree, root_pheno, c2nodeid_dict = utpa.create_tree(icd_codes_df, icd_samples_df)
    return pheno_tree, root_pheno, c2nodeid_dict, icd_codes_df

def get_table_icd(samples_of_interest, nonsamples_of_interest, comorbid_samples, field, shortlist):
    table = [
        [len(samples_of_interest.intersection(comorbid_samples)), len(samples_of_interest.difference(comorbid_samples))],
        [len(nonsamples_of_interest.intersection(comorbid_samples)), len(nonsamples_of_interest.difference(comorbid_samples))]
    ]
    df = pd.DataFrame(table, columns=[f"{field}", f"No {field}"], index=shortlist)
    return df

def get_icd_enrich(samples_of_interest, nonsamples_of_interest, all_cohort_samples, pheno_tree, icd_codes_df, c2nodeid_dict, shortlist, mode=""):
    icd_data = []
    for icdc in icd_codes_df.coding.values:
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        comorbid_samples = all_cohort_samples.intersection(comorbid_samples)
        df = get_table_icd(samples_of_interest, nonsamples_of_interest, comorbid_samples, icdc_node.meaning, shortlist)
        res = fisher_exact(df)
        if mode=="lazy":
            or_study = res.statistic
            cil, cih = np.nan, np.nan        
        else:
            or_study = odds_ratio(df)
            cil, cih = or_study.confidence_interval(confidence_level=0.95)
        icdc_node_data = (icdc, icdc_node.meaning, df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
        icd_data.append(icdc_node_data)
    icd_df = pd.DataFrame(icd_data, columns=["coding", "meaning", f"{shortlist[0]}_comorbid", f"{shortlist[0]}_noncomorbid", f"{shortlist[1]}_comorbid", f"{shortlist[1]}_noncomorbid", "odds_ratio", "p_value", "ci_low", "ci_high"])
    icd_df["FDR"] = stats.false_discovery_control(icd_df.p_value)
    return icd_df

def get_level(icd_node, level):
    if icd_node.parent:
        return get_level(icd_node.parent, level+1)
    else:
        return level

def add_more_info(code, pheno_tree, c2nodeid_dict):
    icd_node = pheno_tree.node_dict[c2nodeid_dict[code]]
    level = get_level(icd_node, 0)
    parent = icd_node.parent.meaning
    return pd.Series({"level": level, "parent": parent})

def add_info_and_filter_icd_enrich(icd_df, pheno_tree, c2nodeid_dict, level_filter=2, fdr_filter=0.05):
    icd_df[["level", "parent"]] = icd_df.coding.apply(add_more_info, args=(pheno_tree, c2nodeid_dict))
    # keep only filter level phenotypes
    icd_df = icd_df.loc[icd_df.level==level_filter]
    # keep only fdr<0.05 significant phenotypes
    icd_df = icd_df.loc[icd_df.FDR<fdr_filter]
    # keep only enriched comorbidities
    icd_df = icd_df.loc[icd_df.odds_ratio>1]
    return icd_df

def create_icd_comorbidity_features(icd_codes_of_interest, samples_of_interest, pheno_tree, c2nodeid_dict):
    icd_samples_data_dict = {"icd": [], "samples": []}
    for icd_code in icd_codes_of_interest:
        icd_node = pheno_tree.node_dict[c2nodeid_dict[icd_code]]
        comorbid_cohort_samples = samples_of_interest.intersection(icd_node.samples)
        icd_samples_data_dict["icd"].append(icd_node.code)
        icd_samples_data_dict["samples"].append(comorbid_cohort_samples)

    icd_samples_df = pd.DataFrame(icd_samples_data_dict)
    icd_samples_df = icd_samples_df.explode("samples")
    icd_samples_df = pd.crosstab(icd_samples_df.samples, icd_samples_df.icd)
    # ensure all samples are present
    existing_samples = set(icd_samples_df.index)
    no_comorbidity_samples = samples_of_interest.difference(existing_samples)
    if len(no_comorbidity_samples)>0:
        df_data = np.zeros((len(no_comorbidity_samples), icd_samples_df.shape[1]))
        no_comorbid_df = pd.DataFrame(df_data, index=list(no_comorbidity_samples), columns=icd_samples_df.columns)
        no_comorbid_df.index.name = icd_samples_df.index.name
        icd_samples_df = pd.concat((icd_samples_df, no_comorbid_df))
    return icd_samples_df

def divide_into_nocomorbidity(features_df):
    no_comorbidity_df = features_df.loc[features_df.sum(axis=1)==0]
    comorbidity_df = features_df.loc[features_df.sum(axis=1)>0]
    return comorbidity_df, no_comorbidity_df

def cluster_by_comorbidity(features_df):
    comorbidity_df, no_comorbidity_df = divide_into_nocomorbidity(features_df)
    no_comorbidity_df["labels"] = 0

    X = comorbidity_df.values
    nclust = 2 #find_the_best_cluster(X)
    kmeans = KMeans(n_clusters=nclust, random_state=42, n_init="auto")
    y = kmeans.fit_predict(X)
    comorbidity_df["labels"] = y
    comorbidity_df["labels"] = comorbidity_df.labels.replace(0, nclust)
    icd_df = pd.concat((comorbidity_df, no_comorbidity_df))
    return icd_df.labels

def get_average_comorbidity_per_label(icd_features_df, comorbidity_labels):
    ncomorbidity_df = icd_features_df.sum(axis=1).to_frame()
    ncomorbidity_df.columns = ["ncomorbidity"]
    ncomorbidity_df = ncomorbidity_df.merge(comorbidity_labels, left_index=True, right_index=True).groupby("labels")["ncomorbidity"].mean()
    avg_comorbidities_per_label = ncomorbidity_df.to_dict()
    return avg_comorbidities_per_label
    
def get_rvas_samples(genotype_file, rvas_genes):
    genotype_df = pd.read_csv(genotype_file)
    rvas_geno = genotype_df.loc[genotype_df.gene.isin(rvas_genes)]
    all_rvas_samples = set(",".join(rvas_geno.samples).split(","))
    return all_rvas_samples

def add_carrier_info(phenotype_df, all_factor_samples):
    phenotype_df["carrier"] = phenotype_df.sample_names.isin(all_factor_samples)
    return phenotype_df

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

def create_phenotype_info_file(phenotype_samples_df, all_factor_samples, comorbidity_labels):
    phenotype_samples_df = add_carrier_info(phenotype_samples_df, all_factor_samples)
    phenotype_samples_info_df = phenotype_samples_df.loc[:, ["sample_names", "bmi", "bmi_prs", "bmi_residuals", "carrier", "bmi_res_categories", "bmi_prs_categories"]]
    phenotype_samples_info_df = phenotype_samples_info_df.merge(comorbidity_labels, left_on="sample_names", right_index=True)
    return phenotype_samples_info_df

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--genotype_file", type=str, help="Filepath of genes with rare deleterious mutation samples")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--icd_raw_dir", type=str, help="Filepath of the icd codes to samples files storage dir")
    parser.add_argument("--icd_codes_file", type=str, help="Filepath of the icd codes to parent node file")
    parser.add_argument("--hes_info_file", type=str, help="Filepath of the sample names with icd info")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode")

    cli_args = parser.parse_args()

    # read cohort data and annotate it with categorical info
    cohort_df = pd.read_csv(
        cli_args.cohort_file, 
        usecols=["sample_names", "genetic_sex", "age"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"])
    cohort_df["sample_names"] = cohort_df.sample_names.astype(str)
    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    rvas_samples = get_rvas_samples(cli_args.genotype_file, ["MC4R"])
    phenotype_samples_df, mean_bmi_dict = create_bmi_res_prs_decile_data(cohort_df)
    #####################
    # icd tree creation #
    #####################
    # create icd tree and only keep samples with icd info in cohort df
    pheno_tree, root_pheno, c2nodeid_dict, icd_codes_df = create_icd_tree(cli_args.icd_raw_dir, cli_args.icd_codes_file)
    hes_info_df = pd.read_csv(cli_args.hes_info_file, dtype={"sample_names": str, "hes_info": float})
    all_icd_samples = set(hes_info_df.loc[hes_info_df.hes_info>0, "sample_names"].values)
    print(f"Samples in cohort {len(phenotype_samples_df)}")
    phenotype_samples_df = phenotype_samples_df.loc[phenotype_samples_df.sample_names.isin(all_icd_samples)]
    print(f"Samples in cohort with icd {len(phenotype_samples_df)}")
    print(f"Samples with combo {len(combo_samples)}")
    combo_samples = combo_samples.intersection(all_icd_samples)
    print(f"Samples with combo and icd {len(combo_samples)}")
    rvas_samples = rvas_samples.intersection(all_icd_samples)
    #####################################
    # get obesity related comorbidities #
    #####################################
    # get enriched comorbidities for obese people: these are obesity related comorbidities
    obese_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories=="severe obesity")|(phenotype_samples_df.bmi_res_categories=="obese"), "sample_names"].values)
    notobese_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories!="severe obesity")&(phenotype_samples_df.bmi_res_categories!="obese"), "sample_names"].values)
    all_samples = set(phenotype_samples_df.sample_names.values)
    save_file = os.path.join(cli_args.save_dir, "obese_samples_icd_enrich.csv")
    if (cli_args.debug) and (os.path.exists(save_file)):
        obese_samples_enrich_df = pd.read_csv(save_file)
    else:
        obese_samples_enrich_df = get_icd_enrich(obese_samples, notobese_samples, all_samples, pheno_tree, icd_codes_df, c2nodeid_dict, ["obese", "not obese"])
        obese_samples_enrich_df.to_csv(save_file, index=False)
    obesity_related_comorbidity_df = add_info_and_filter_icd_enrich(obese_samples_enrich_df, pheno_tree, c2nodeid_dict, level_filter=2, fdr_filter=0.05)
    # get rid of obesity terms
    obesity_related_comorbidity_df = obesity_related_comorbidity_df.loc[~obesity_related_comorbidity_df.meaning.str.lower().str.contains("obesity")] 
    # eliminate external factors from ICD file: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=40001
    obesity_related_comorbidity_df = obesity_related_comorbidity_df.loc[~obesity_related_comorbidity_df.parent.str.startswith(("Chapter XV", "Chapter XX", "Chapter XXI"))]
    save_file = os.path.join(cli_args.save_dir, "obesity_related_comorbidities.csv")
    obesity_related_comorbidity_df.to_csv(save_file, index=False)
    print(len(obesity_related_comorbidity_df))
    ###########################################
    # create comorbidity features for samples #
    ###########################################
    # creating a feature file for selected icd comorbidities
    icd_codes_of_interest = set(obesity_related_comorbidity_df.coding)
    samples_of_interest = set(phenotype_samples_df.sample_names.values)
    save_file = os.path.join(cli_args.save_dir, "obesity_related_comorbid_feature_file.csv.gz")
    icd_features_df = create_icd_comorbidity_features(icd_codes_of_interest, samples_of_interest, pheno_tree, c2nodeid_dict)
    icd_features_df.to_csv(save_file)
    #######################################################
    # cluster by comorbidity and assign labels to samples #
    #######################################################
    comorbidity_labels = cluster_by_comorbidity(icd_features_df)
    avg_comorbidity_per_group = get_average_comorbidity_per_label(icd_features_df, comorbidity_labels)
    print(avg_comorbidity_per_group)
    #############################
    # risk factor effect on bmi #
    #############################
    phenotype_samples_combo_info_df = create_phenotype_info_file(phenotype_samples_df, combo_samples, comorbidity_labels)
    save_file = os.path.join(cli_args.save_dir, "all_factors_bmi_info.csv")
    phenotype_samples_combo_info_df.to_csv(save_file, index=False)
    save_file = os.path.join(cli_args.save_dir, "all_factors_bmi_plot.pdf")
    fig = utpl.create_variate_interaction_plot(phenotype_samples_combo_info_df)
    utpl.save_pdf(save_file, fig)

    phenotype_samples_rvas_info_df = create_phenotype_info_file(phenotype_samples_df, rvas_samples, comorbidity_labels)
    save_file = os.path.join(cli_args.save_dir, "all_factors_rvas_bmi_info.csv")
    phenotype_samples_rvas_info_df.to_csv(save_file, index=False)
    save_file = os.path.join(cli_args.save_dir, "all_factors_rvas_bmi_plot.pdf")
    fig = utpl.create_variate_interaction_plot(phenotype_samples_rvas_info_df)
    utpl.save_pdf(save_file, fig)
    ######################
    # risk factors model #
    ######################
    X_poly, y, poly_features = prepare_polynomial_features(phenotype_samples_combo_info_df, ["carrier", "bmi_prs", "labels"], "bmi", True)
    r2, coefs, conf_ints, p_vals = train_model_sm(X_poly, y)
    df_data = np.concatenate([poly_features.reshape(-1,1), coefs.reshape(-1,1), conf_ints, p_vals.reshape(-1,1)], axis=1)
    coef_df = pd.DataFrame(df_data, columns=["variables", "coefs", "ci_low", "ci_high", "p_value"])
    coef_df["variables"] = coef_df.variables.str.replace(" ", "*")
    save_file = os.path.join(cli_args.save_dir, "all_factors_bmi_interaction_plot.pdf")
    fig = utpl.create_variate_interaction_model_plot(coef_df)
    utpl.save_pdf(save_file, fig)

    X_poly, y, poly_features = prepare_polynomial_features(phenotype_samples_rvas_info_df, ["carrier", "bmi_prs", "labels"], "bmi", True)
    r2, coefs, conf_ints, p_vals = train_model_sm(X_poly, y)
    df_data = np.concatenate([poly_features.reshape(-1,1), coefs.reshape(-1,1), conf_ints, p_vals.reshape(-1,1)], axis=1)
    coef_df = pd.DataFrame(df_data, columns=["variables", "coefs", "ci_low", "ci_high", "p_value"])
    coef_df["variables"] = coef_df.variables.str.replace(" ", "*")
    save_file = os.path.join(cli_args.save_dir, "all_factors_rvas_bmi_interaction_plot.pdf")
    fig = utpl.create_variate_interaction_model_plot(coef_df)
    utpl.save_pdf(save_file, fig)   
