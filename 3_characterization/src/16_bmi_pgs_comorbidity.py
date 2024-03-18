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

def get_prs_percent_by_prs(df):
    df1 = df.groupby("bmi_prs_categories")["bmi_res_categories"].value_counts(normalize=True)
    df1 = df1.mul(100)
    df1 = df1.rename('percent').reset_index()
    return df1

def create_bmi_res_cat_prs_cat_df(phenotype_samples_df):
    noncombo = get_prs_percent_by_prs(phenotype_samples_df.loc[phenotype_samples_df.carrier==False])
    noncombo["carrier"] = False
    combo = get_prs_percent_by_prs(phenotype_samples_df.loc[phenotype_samples_df.carrier==True])
    combo["carrier"] = True
    plot_df = pd.concat((combo, noncombo))
    return plot_df

def create_bmi_res_prs_decile_data(phenotype_df, all_combo_samples):
    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    scaled_numerical_cols = []#["bmi_prs"]

    phenotype_df = get_scaled_bmi(phenotype_df, categorical_cols, numerical_cols, scaled_numerical_cols)
    phenotype_df["bmi_res_decile"] = pd.qcut(phenotype_df.bmi_residuals, q=10)
    phenotype_df["bmi_res_decile_num"] = pd.qcut(phenotype_df.bmi_residuals, q=10, labels=False)
    phenotype_df["bmi_prs_decile"] = pd.qcut(phenotype_df.bmi_prs, q=10)
    phenotype_df["bmi_prs_decile_num"] = pd.qcut(phenotype_df.bmi_prs, q=10, labels=False)
    phenotype_df["carrier"] = phenotype_df.sample_names.isin(all_combo_samples)
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

def get_combos_icd_enrichment(phenotype_samples_df, bmi_res_category, pheno_tree, c2nodeid_dict, comorbid_df):
    cat_combo_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories==bmi_res_category)&(phenotype_samples_df.carrier==True), "sample_names"].values)
    cat_noncombo_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories==bmi_res_category)&(phenotype_samples_df.carrier==False), "sample_names"].values)
    cat_samples = set(phenotype_samples_df.loc[phenotype_samples_df.bmi_res_categories==bmi_res_category, "sample_names"].values)
    # individual icd enrichment
    combo_samples_enrich_df = get_icd_enrich(cat_combo_samples, cat_noncombo_samples, cat_samples, pheno_tree, comorbid_df, c2nodeid_dict, ["combo", "non combo"])
    # all comorbidity enrichment
    orc_samples = set()
    for icdc in comorbid_df.sort_values("FDR").coding.values:
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        comorbid_samples = cat_samples.intersection(comorbid_samples)
        orc_samples.update(comorbid_samples)
    df = get_table_icd(cat_combo_samples, cat_noncombo_samples, orc_samples, "obesity related comorbidity", ["Combo", "Non Combo"])
    res = fisher_exact(df)
    or_study = odds_ratio(df)
    cil, cih = or_study.confidence_interval(confidence_level=0.95)
    all_comorbid_data = ("All", "Obesity Related Comorbidities", df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
    all_comorbid_df = pd.Series(all_comorbid_data, index=["coding", "meaning", f"combo_comorbid", f"combo_noncomorbid", f"noncombo_comorbid", f"noncombo_noncomorbid", "odds_ratio", "p_value", "ci_low", "ci_high"])
    return combo_samples_enrich_df, all_comorbid_df

def get_pgs_icd_enrichment(phenotype_samples_df, bmi_res_category, pheno_tree, c2nodeid_dict, comorbid_df):
    cat_pgs_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories==bmi_res_category)&(phenotype_samples_df.bmi_prs_categories=="highest"), "sample_names"].values)
    cat_nonpgs_samples = set(phenotype_samples_df.loc[(phenotype_samples_df.bmi_res_categories==bmi_res_category)&(phenotype_samples_df.bmi_prs_categories!="highest"), "sample_names"].values)
    cat_samples = set(phenotype_samples_df.loc[phenotype_samples_df.bmi_res_categories==bmi_res_category, "sample_names"].values)
    # individual icd enrichment
    pgs_samples_enrich_df = get_icd_enrich(cat_pgs_samples, cat_nonpgs_samples, cat_samples, pheno_tree, comorbid_df, c2nodeid_dict, ["pgs", "non pgs"])
    # all comorbidity enrichment
    orc_samples = set()
    for icdc in comorbid_df.sort_values("FDR").coding.values:
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        comorbid_samples = cat_samples.intersection(comorbid_samples)
        orc_samples.update(comorbid_samples)
    df = get_table_icd(cat_pgs_samples, cat_nonpgs_samples, orc_samples, "obesity related comorbidity", ["pgs", "non pgs"])
    res = fisher_exact(df)
    or_study = odds_ratio(df)
    cil, cih = or_study.confidence_interval(confidence_level=0.95)
    all_comorbid_data = ("All", "Obesity Related Comorbidities", df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
    all_comorbid_df = pd.Series(all_comorbid_data, index=["coding", "meaning", f"pgs_comorbid", f"pgs_noncomorbid", f"nonpgs_comorbid", f"nonpgs_noncomorbid", "odds_ratio", "p_value", "ci_low", "ci_high"])
    return pgs_samples_enrich_df, all_comorbid_df


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--icd_raw_dir", type=str, help="Filepath of the icd codes to samples files storage dir")
    parser.add_argument("--icd_codes_file", type=str, help="Filepath of the icd codes to parent node file")
    parser.add_argument("--hes_info_file", type=str, help="Filepath of the sample names with icd info")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode")

    cli_args = parser.parse_args()
    os.makedirs(cli_args.save_dir, exist_ok=True)

    # read cohort data and annotate it with categorical info
    cohort_df = pd.read_csv(
        cli_args.cohort_file, 
        usecols=["sample_names", "genetic_sex", "age"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"])
    cohort_df["sample_names"] = cohort_df.sample_names.astype(str)
    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    phenotype_samples_df, mean_bmi_dict = create_bmi_res_prs_decile_data(cohort_df, combo_samples)
    save_file = os.path.join(cli_args.save_dir, "phenotype_info.csv.gz")
    cohort_df.loc[:, ["sample_names", "bmi", "bmi_prs", "bmi_residuals", "bmi_res_categories", "bmi_prs_categories", "carrier"]].to_csv(save_file, index=False)
    exit()
    ################################
    # bmi res cat and prs cat plot #
    ################################
    bmi_res_cat_prs_cat_df = create_bmi_res_cat_prs_cat_df(phenotype_samples_df)
    save_file = os.path.join(cli_args.save_dir, "res_cat_prs_cat.pdf")
    fig = utpl.create_prs_cat_bmi_res_cat_plot(bmi_res_cat_prs_cat_df, mean_bmi_dict)
    utpl.save_pdf(save_file, fig)
    #######################################################
    # category wise residual versus prs distribution plot #
    #######################################################
    save_file = os.path.join(cli_args.save_dir, "res_cat_prs_dist.pdf")
    fig, ax = utpl.create_res_cat_prs_dist_plot(phenotype_samples_df, mean_bmi_dict)
    utpl.save_pdf(save_file, fig)
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
    ############################################
    # get severe obesity related comorbidities #
    ############################################
    # get enriched comorbidities for severely obese people: these are severe obesity related comorbidities
    severely_obese_samples = set(phenotype_samples_df.loc[phenotype_samples_df.bmi_res_categories=="severe obesity", "sample_names"].values)
    notseverely_obese_samples = set(phenotype_samples_df.loc[phenotype_samples_df.bmi_res_categories!="severe obesity", "sample_names"].values)
    all_samples = set(phenotype_samples_df.sample_names.values)
    save_file = os.path.join(cli_args.save_dir, "severely_obese_samples_icd_enrich.csv")
    if (cli_args.debug) and (os.path.exists(save_file)):
        severely_obese_samples_enrich_df = pd.read_csv(save_file)
    else:
        severely_obese_samples_enrich_df = get_icd_enrich(severely_obese_samples, notseverely_obese_samples, all_samples, pheno_tree, icd_codes_df, c2nodeid_dict, ["severely obese", "not severely obese"])
        severely_obese_samples_enrich_df.to_csv(save_file, index=False)
    severe_obesity_related_comorbidity_df = add_info_and_filter_icd_enrich(severely_obese_samples_enrich_df, pheno_tree, c2nodeid_dict, level_filter=2, fdr_filter=0.05)
    # get rid of obesity terms
    severe_obesity_related_comorbidity_df = severe_obesity_related_comorbidity_df.loc[~severe_obesity_related_comorbidity_df.meaning.str.lower().str.contains("obesity")] 
    # eliminate external factors from ICD file: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=40001
    severe_obesity_related_comorbidity_df = severe_obesity_related_comorbidity_df.loc[~severe_obesity_related_comorbidity_df.parent.str.startswith(("Chapter XV", "Chapter XX", "Chapter XXI"))]
    save_file = os.path.join(cli_args.save_dir, "severe_obesity_related_comorbidities.csv")
    severe_obesity_related_comorbidity_df.to_csv(save_file, index=False)
    ##################################################################################
    # check severe obesity related comorbidities enrichment in severely obese combos #
    ##################################################################################
    # get individual and all severely obese comorbidity enrichment for combos
    severely_obese_combo_ind_enrich_df, severely_obese_combo_enrich_df = get_combos_icd_enrichment(phenotype_samples_df, "severe obesity", pheno_tree, c2nodeid_dict, severe_obesity_related_comorbidity_df)
    save_file = os.path.join(cli_args.save_dir, "severe_obesity_related_comorbidities_severely_obese_combos.csv")
    severely_obese_combo_ind_enrich_df.to_csv(save_file, index=False)
    #########################################################################
    # check severe obesity related comorbidities enrichment in obese combos #
    #########################################################################
    # get individual and all severely obese comorbidity enrichment for combos in obese group
    obese_combo_ind_enrich_df, obese_combo_enrich_df = get_combos_icd_enrichment(phenotype_samples_df, "obese", pheno_tree, c2nodeid_dict, severe_obesity_related_comorbidity_df)
    save_file = os.path.join(cli_args.save_dir, "severe_obesity_related_comorbidities_obese_combos.csv")
    obese_combo_ind_enrich_df.to_csv(save_file, index=False)
    ########################################################################################
    # check severe obesity related comorbidities enrichment in severely obese pgs carriers #
    ########################################################################################
    # get individual and all severely obese comorbidity enrichment for highest pgs individuals
    severely_obese_pgs_ind_enrich_df, severely_obese_pgs_enrich_df = get_pgs_icd_enrichment(phenotype_samples_df, "severe obesity", pheno_tree, c2nodeid_dict, severe_obesity_related_comorbidity_df)
    save_file = os.path.join(cli_args.save_dir, "severe_obesity_related_comorbidities_severely_obese_pgs.csv")
    severely_obese_pgs_ind_enrich_df.to_csv(save_file, index=False)
    # group all severely obese comorbidity analysis
    all_severely_obese_comorbidity_enrich_df = pd.concat((severely_obese_combo_enrich_df, obese_combo_enrich_df, severely_obese_pgs_enrich_df), axis=1)
    all_severely_obese_comorbidity_enrich_df.columns = ["Severely obese Combo", "Obese Combo", "Severely obese high PGS"]
    save_file = os.path.join(cli_args.save_dir, "all_severe_obesity_related_comorbidities_enrichment.csv")
    all_severely_obese_comorbidity_enrich_df.to_csv(save_file)
