import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
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

def create_bmi_res_prs_decile_data(phenotype_df, all_combo_samples):
    categorical_cols = ["genetic_sex"]
    numerical_cols = ["age"] + [f"genetic_pca{i}" for i in range(1, 40)]
    scaled_numerical_cols = []#["bmi_prs"]

    phenotype_df = get_scaled_bmi(phenotype_df, categorical_cols, numerical_cols, scaled_numerical_cols)
    phenotype_df["bmi_res_decile"] = pd.qcut(phenotype_df.bmi_residuals, q=10)
    phenotype_df["bmi_res_decile_num"] = pd.qcut(phenotype_df.bmi_residuals, q=10, labels=False)
    phenotype_df["bmi_prs_decile"] = pd.qcut(phenotype_df.bmi_prs, q=10)
    phenotype_df["bmi_prs_decile_num"] = pd.qcut(phenotype_df.bmi_prs, q=10, labels=False)
    phenotype_combo_samples_df = phenotype_df.loc[phenotype_df.sample_names.astype(str).isin(list(map(str, all_combo_samples)))]
    phenotype_other_samples_df = phenotype_df.loc[~phenotype_df.sample_names.astype(str).isin(list(map(str, all_combo_samples)))]

    phenotype_other_samples_df["description"] = "Non Combo"
    phenotype_combo_samples_df["description"] = "Combo"
    phenotype_samples_df = pd.concat((phenotype_combo_samples_df, phenotype_other_samples_df))
    return phenotype_samples_df

def get_table_lf(int_decile_combo_lf, int_decile_noncombo_lf, field):
    table = [
        [len(int_decile_combo_lf.loc[int_decile_combo_lf[field]==1]), len(int_decile_combo_lf.loc[int_decile_combo_lf[field]==0])],
        [len(int_decile_noncombo_lf.loc[int_decile_noncombo_lf[field]==1]), len(int_decile_noncombo_lf.loc[int_decile_noncombo_lf[field]==0])]
    ]
    df = pd.DataFrame(table, columns=[f"{field}", f"No {field}"], index=["Combo", "Non Combo"])
    return df

def check_lifestyle_factor_enrichment(phenotype_samples_df, lf_df, decile):
    int_decile = phenotype_samples_df.loc[phenotype_samples_df.bmi_res_decile_num==decile]
    int_decile_combo = int_decile.loc[int_decile.description=="Combo"]
    int_decile_noncombo = int_decile.loc[int_decile.description=="Non Combo"]
    int_decile_combo_lf = int_decile_combo.merge(lf_df, left_on="sample_names", right_on="Sample_Name")
    int_decile_noncombo_lf = int_decile_noncombo.merge(lf_df, left_on="sample_names", right_on="Sample_Name")
    lf_fields_check = ["alcohol", "met", "sleep", "smoke", "meat", "sedentary"]
    all_lf_data = []
    for field in lf_fields_check:
        df = get_table_lf(int_decile_combo_lf, int_decile_noncombo_lf, field)
        res = fisher_exact(df, alternative='less')
        or_study = odds_ratio(df)
        cil, cih = or_study.confidence_interval(confidence_level=0.95)
        lf_data = (field, df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
        all_lf_data.append(lf_data)
    
    lf_enrich_df = pd.DataFrame(all_lf_data, columns=["lifestyle", "combo_lf", "combo_nolf", "noncombo_lf", "noncombo_nolf", "odds_ratio", "p_value", "ci_low", "ci_high"])
    lf_enrich_df["FDR"] = stats.false_discovery_control(lf_enrich_df.p_value)
    return lf_enrich_df

def get_rvas_gene(gene_file):
    with open(gene_file, "r") as f:
        genes = set([l.strip() for l in f.readlines()])
    return genes

def check_rvas_gene_enrichment(phenotype_samples_df, decile, genotype_df, rvas_studies, combo_genes, combo_samples):
    int_decile = phenotype_samples_df.loc[phenotype_samples_df.bmi_res_decile_num==decile]
    int_decile_combo = int_decile.loc[int_decile.description=="Combo"]
    int_decile_noncombo = int_decile.loc[int_decile.description=="Non Combo"]
    rvas_genes = reduce(lambda x,y: x.union(y), [get_rvas_gene(gf) for gf in rvas_studies])
    rvas_other_genes = rvas_genes.difference(combo_genes)
    rvas_geno = genotype_df.loc[genotype_df.gene.isin(rvas_other_genes)]
    all_rvas_samples = set(",".join(rvas_geno.samples).split(","))
    combo_samples = set(int_decile_combo.sample_names.astype(str).values)
    noncombo_samples = set(int_decile_noncombo.sample_names.astype(str).values)
    table = [
        [len(combo_samples.intersection(all_rvas_samples)), len(combo_samples.difference(all_rvas_samples)), ],
        [len(noncombo_samples.intersection(all_rvas_samples)), len(noncombo_samples.difference(all_rvas_samples)), ]
    ]
    df = pd.DataFrame(table, columns=["RVAS", "No RVAS"], index=["Combo", "Non Combo"])
    res = fisher_exact(df, alternative='less')
    or_study = odds_ratio(df)
    cil, cih = or_study.confidence_interval(confidence_level=0.95)
    rvas_data = ("RVAS", df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
    rvas_enrich_df = pd.Series(rvas_data, index=["rvas", "combo_rvas", "combo_norvas", "noncombo_rvas", "noncombo_norvas", "odds_ratio", "p_value", "ci_low", "ci_high"])
    return rvas_enrich_df

def get_table_icd(combo_samples, noncombo_samples, comorbid_samples, field):
    table = [
        [len(combo_samples.intersection(comorbid_samples)), len(combo_samples.difference(comorbid_samples))],
        [len(noncombo_samples.intersection(comorbid_samples)), len(noncombo_samples.difference(comorbid_samples))]
    ]
    df = pd.DataFrame(table, columns=[f"{field}", f"No {field}"], index=["Combo", "Non Combo"])
    return df

def get_icd_enrich_per_decile(phenotype_samples_df, decile, pheno_tree, icd_codes_df, c2nodeid_dict):
    decile_df = phenotype_samples_df.loc[phenotype_samples_df.bmi_res_decile_num==decile]
    combo_samples = set(decile_df.loc[decile_df.description=="Combo", "sample_names"].astype(str).values)
    noncombo_samples = set(decile_df.loc[decile_df.description=="Non Combo", "sample_names"].astype(str).values)

    icd_data = []

    for icdc in icd_codes_df.coding.values:
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        df = get_table_icd(combo_samples, noncombo_samples, comorbid_samples, icdc_node.meaning)
        res = fisher_exact(df)
        or_study = odds_ratio(df)
        cil, cih = or_study.confidence_interval(confidence_level=0.95)
        icdc_node_data = (icdc, icdc_node.meaning, df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
        icd_data.append(icdc_node_data)
    
    icd_df = pd.DataFrame(icd_data, columns=["icd_code", "icd_meaning", "combo_comorbid", "combo_noncomorbid", "noncombo_comorbid", "noncombo_noncomorbid", "odds_ratio", "p_value", "ci_low", "ci_high"])
    icd_df["FDR"] = stats.false_discovery_control(icd_df.p_value)
    return icd_df

def get_prs_percent_by_weight(df):
    df1 = df.groupby("bmi_res_categories")["bmi_prs_categories"].value_counts(normalize=True)
    df1 = df1.mul(100)
    df1 = df1.rename('percent').reset_index()
    mean_bmi_dict = df.groupby("bmi_res_categories")["bmi"].mean().to_dict()
    return df1, mean_bmi_dict

def get_prs_bmi_res_decile_data_weight_percent(phenotype_samples_df):
    phenotype_samples_df["bmi_res_categories"] = phenotype_samples_df.bmi_res_decile_num.map({
        0: "underweight", 
        1:"normal", 2:"normal", 
        3:"overweight", 4:"overweight", 5:"overweight", 6:"overweight", 7:"overweight",
        8:"obese", 9:"severe obesity"
        })
    phenotype_samples_df["bmi_prs_categories"] = phenotype_samples_df.bmi_prs_decile_num.map({
        0: "lowest", 
        1:"middle", 2:"middle", 
        3:"middle", 4:"middle", 5:"middle", 6:"middle", 7:"middle",
        8:"middle", 9:"highest"
        })
    phenotype_combo_samples_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Combo"]
    phenotype_other_samples_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Non Combo"]
    phenotype_combo_samples_plot_df, cbmi_dict = get_prs_percent_by_weight(phenotype_combo_samples_df)
    phenotype_other_samples_plot_df, ncbmi_dict = get_prs_percent_by_weight(phenotype_other_samples_df)
    return phenotype_combo_samples_plot_df, cbmi_dict, phenotype_other_samples_plot_df, ncbmi_dict

def get_prs_percent_by_prs(df):
    df1 = df.groupby("bmi_prs_categories")["bmi_res_categories"].value_counts(normalize=True)
    df1 = df1.mul(100)
    df1 = df1.rename('percent').reset_index()
    mean_bmi_dict = df.groupby("bmi_res_categories")["bmi"].mean().to_dict()
    return df1, mean_bmi_dict

def get_prs_bmi_res_decile_data_prs_percent(phenotype_samples_df):
    phenotype_samples_df["bmi_res_categories"] = phenotype_samples_df.bmi_res_decile_num.map({
        0: "underweight", 
        1:"normal", 2:"normal", 
        3:"overweight", 4:"overweight", 5:"overweight", 6:"overweight", 7:"overweight",
        8:"obese", 9:"severe obesity"
        })
    phenotype_samples_df["bmi_prs_categories"] = phenotype_samples_df.bmi_prs_decile_num.map({
        0: "lowest", 
        1:"middle", 2:"middle", 
        3:"middle", 4:"middle", 5:"middle", 6:"middle", 7:"middle",
        8:"middle", 9:"highest"
        })
    phenotype_combo_samples_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Combo"]
    phenotype_other_samples_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Non Combo"]
    phenotype_combo_samples_plot_df, cbmi_dict = get_prs_percent_by_prs(phenotype_combo_samples_df)
    phenotype_other_samples_plot_df, ncbmi_dict = get_prs_percent_by_prs(phenotype_other_samples_df)
    return phenotype_combo_samples_plot_df, cbmi_dict, phenotype_other_samples_plot_df, ncbmi_dict

def get_table_prs_decile_helper(combo_samples, other_combo_samples, decile_samples):
    combo_in_decile = len(combo_samples.intersection(decile_samples))
    combo_not_in_decile = len(combo_samples.difference(decile_samples))
    other_combo_in_decile = len(other_combo_samples.intersection(decile_samples))
    other_combo_not_in_decile = len(other_combo_samples.difference(decile_samples))
    table = [
        [combo_in_decile, combo_not_in_decile],
        [other_combo_in_decile, other_combo_not_in_decile]
    ]
    res = fisher_exact(table, alternative="greater")
    or_res = odds_ratio(table)
    ci_low, ci_high = or_res.confidence_interval(confidence_level=0.95)
    return combo_in_decile, combo_not_in_decile, other_combo_in_decile, other_combo_not_in_decile, or_res.statistic, ci_low, ci_high, res.pvalue

def get_table_prs_decile(ser, all_combo_samples, decile_samples_bmi_prs):
    if not pd.isnull(ser.combo_samples):
        combo_samples = set(ser.combo_samples.split("|"))
    else:
        combo_samples = set()
    other_combo_samples = all_combo_samples.difference(combo_samples)
    combo_in_decile, combo_not_in_decile, other_combo_in_decile, other_combo_not_in_decile, odds_ratio_bmi_prs, ci_low_bmi_prs, ci_high_bmi_prs, pvalue_bmi_prs = get_table_prs_decile_helper(combo_samples, other_combo_samples, decile_samples_bmi_prs)
    
    return pd.Series({
        "combo_in_decile_bmi_prs": combo_in_decile, "combo_not_in_decile_bmi_prs": combo_not_in_decile, 
        "other_combo_in_decile_bmi_prs": other_combo_in_decile, "other_combo_not_in_decile_bmi_prs": other_combo_not_in_decile,
        "odds_ratio_bmi_prs": odds_ratio_bmi_prs, "ci_low_bmi_prs": ci_low_bmi_prs, "ci_high_bmi_prs": ci_high_bmi_prs, "p_value_bmi_prs": pvalue_bmi_prs
        })

def get_combos_enriched_in_prs_decile(phenotype_samples_df, decile, combo_df, all_combo_samples):
    decile_samples_bmi_prs = set(phenotype_samples_df.loc[phenotype_samples_df.bmi_prs_decile_num==decile, "sample_names"].astype(str).values)
    combo_info_df = combo_df.merge(combo_df.apply(get_table_prs_decile, args=(all_combo_samples, decile_samples_bmi_prs), axis=1), left_index=True, right_index=True)
    return combo_info_df


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the binarized lifestyle factor file")
    parser.add_argument("--rvas_gene_files", type=str, help="Filepath of the binarized lifestyle factor file", nargs="+")
    parser.add_argument("--icd_raw_dir", type=str, help="Filepath of the icd codes to samples files storage dir")
    parser.add_argument("--icd_codes_file", type=str, help="Filepath of the icd codes to parent node file")
    parser.add_argument("--icd_terms_plot", type=str, help="icd codes to plot in comorbidity enrichment", nargs="+")
    parser.add_argument("--save_dir", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    cohort_df = pd.read_csv(
        cli_args.cohort_file, 
        usecols=["sample_names", "genetic_sex", "age"] + [f"genetic_pca{i}" for i in range(1, 40)] + ["bmi_prs", "bmi"])
    cohort_df["sample_names"] = cohort_df.sample_names.astype(str)
    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)

    phenotype_samples_df = create_bmi_res_prs_decile_data(cohort_df, combo_samples)
    fig, ax = utpl.create_decile_rank_plot(phenotype_samples_df)
    save_file = os.path.join(cli_args.save_dir, "prs_by_bmi_res_decile.pdf")
    utpl.save_pdf(save_file, fig)

    fig, ax = utpl.create_regplot(phenotype_samples_df)
    save_file = os.path.join(cli_args.save_dir, "prs_by_bmi_res_regplot.pdf")
    utpl.save_pdf(save_file, fig)

    lf_df = pd.read_csv(cli_args.lifestyle_file)
    lf_df["Sample_Name"] = lf_df.Sample_Name.astype(str)
    lf_enrich_df_tenth = check_lifestyle_factor_enrichment(phenotype_samples_df, lf_df, 9)
    save_file = os.path.join(cli_args.save_dir, "lifestyle_enrich_tenth_decile.csv")
    lf_enrich_df_tenth.to_csv(save_file, index=False)

    genotype_df = pd.read_csv(cli_args.genotype_file)
    rvas_df_tenth = check_rvas_gene_enrichment(phenotype_samples_df, 9, genotype_df, cli_args.rvas_gene_files, combo_genes, combo_samples)
    save_file = os.path.join(cli_args.save_dir, "rvas_enrich_tenth_decile.csv")
    rvas_df_tenth.to_frame().T.to_csv(save_file, index=False)

    icd_samples_df = utpa.create_icd_samples_file(cli_args.icd_raw_dir)
    icd_codes_df = pd.read_csv(cli_args.icd_codes_file, usecols=["coding", "meaning", "node_id", "parent_id"], sep="\t")
    icd_codes_df["coding"] = icd_codes_df.coding.str.replace(" ", "")
    pheno_tree, root_pheno, c2nodeid_dict = utpa.create_tree(icd_codes_df, icd_samples_df)
    icd_df_tenth = get_icd_enrich_per_decile(phenotype_samples_df, 9, pheno_tree, icd_codes_df, c2nodeid_dict)
    save_file = os.path.join(cli_args.save_dir, "icd_enrich_tenth_decile.csv")
    icd_df_tenth.to_csv(save_file)
    icd_codes = cli_args.icd_terms_plot #["I272", "E059", "E06", "Block K90-K93", "I27", "M19", "I259", "M190", "K92", "K804", "M1907"]
    plot_df = icd_df_tenth.loc[icd_df_tenth.icd_code.isin(icd_codes)].sort_values("p_value")
    plot_df = plot_df.loc[~plot_df.ci_high.isin([np.inf, -np.inf])]
    fig = utpl.get_odds_ratio_plot_bmi_extreme(plot_df, xticks=[0, 0.6, 1, 5])
    save_file = os.path.join(cli_args.save_dir, "icd_enrich_tenth_decile.pdf")
    utpl.save_pdf(save_file, fig)


    icd_df_ninth = get_icd_enrich_per_decile(phenotype_samples_df, 8, pheno_tree, icd_codes_df, c2nodeid_dict)
    save_file = os.path.join(cli_args.save_dir, "icd_enrich_ninth_decile.csv")
    icd_df_ninth.to_csv(save_file)
    plot_df = icd_df_ninth.loc[icd_df_ninth.icd_code.isin(icd_codes)].sort_values("p_value")
    plot_df = plot_df.loc[~plot_df.ci_high.isin([np.inf, -np.inf])]
    fig = utpl.get_odds_ratio_plot_bmi_extreme(plot_df, xticks=[0, 0.6, 1, 5])
    save_file = os.path.join(cli_args.save_dir, "icd_enrich_ninth_decile.pdf")
    utpl.save_pdf(save_file, fig)

    icd_df_first = get_icd_enrich_per_decile(phenotype_samples_df, 0, pheno_tree, icd_codes_df, c2nodeid_dict)
    save_file = os.path.join(cli_args.save_dir, "icd_enrich_first_decile.csv")
    icd_df_first.to_csv(save_file)

    phenotype_combo_samples_plot_df, cbmi_dict, phenotype_other_samples_plot_df, ncbmi_dict = get_prs_bmi_res_decile_data_weight_percent(phenotype_samples_df)
    fig = utpl.get_prs_bmi_res_decile_plot(phenotype_combo_samples_plot_df, cbmi_dict)
    save_file = os.path.join(cli_args.save_dir, "bmi_res_prs_decile_compare_combos_weight.pdf")
    utpl.save_pdf(save_file, fig)

    fig = utpl.get_prs_bmi_res_decile_plot(phenotype_other_samples_plot_df, ncbmi_dict)
    save_file = os.path.join(cli_args.save_dir, "bmi_res_prs_decile_compare_others_weight.pdf")
    utpl.save_pdf(save_file, fig)

    phenotype_combo_samples_plot_df, cbmi_dict, phenotype_other_samples_plot_df, ncbmi_dict = get_prs_bmi_res_decile_data_prs_percent(phenotype_samples_df)
    fig = utpl.get_prs_bmi_res_decile_plot(phenotype_combo_samples_plot_df, cbmi_dict)
    save_file = os.path.join(cli_args.save_dir, "bmi_res_prs_decile_compare_combos_prs.pdf")
    utpl.save_pdf(save_file, fig)

    fig = utpl.get_prs_bmi_res_decile_plot(phenotype_other_samples_plot_df, ncbmi_dict)
    save_file = os.path.join(cli_args.save_dir, "bmi_res_prs_decile_compare_others_prs.pdf")
    utpl.save_pdf(save_file, fig)

    combo_df = pd.concat([pd.read_csv(cf, dtype={"uniq_items": str, "combo_samples": str}) for cf in cli_args.combo_files]).reset_index(drop=True)
    combo_info_df_bottom = get_combos_enriched_in_prs_decile(phenotype_samples_df, 0, combo_df, combo_samples)
    save_file = os.path.join(cli_args.save_dir, "combos_in_bottom_prs.csv")
    combo_info_df_bottom.to_csv(save_file, index=False)

    combo_info_df_top = get_combos_enriched_in_prs_decile(phenotype_samples_df, 9, combo_df, combo_samples)
    save_file = os.path.join(cli_args.save_dir, "combos_in_top_prs.csv")
    combo_info_df_top.to_csv(save_file, index=False)
    