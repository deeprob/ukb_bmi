import os
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio
from scipy import stats
from functools import reduce

import utils.parsing as utpa
import utils.plotting as utpl
import plotly.express as px


def get_table_icd(combo_samples, noncombo_samples, comorbid_samples, field):
    table = [
        [len(combo_samples.intersection(comorbid_samples)), len(combo_samples.difference(comorbid_samples))],
        [len(noncombo_samples.intersection(comorbid_samples)), len(noncombo_samples.difference(comorbid_samples))]
    ]
    df = pd.DataFrame(table, columns=[f"{field}", f"No {field}"], index=["Combo", "Non Combo"])
    return df

def get_icd_enrich(phenotype_samples_df, pheno_tree, icd_codes_df, c2nodeid_dict):
    combo_samples = set(phenotype_samples_df.loc[phenotype_samples_df.carrier==True, "sample_names"].astype(str).values)
    noncombo_samples = set(phenotype_samples_df.loc[phenotype_samples_df.carrier==False, "sample_names"].astype(str).values)
    all_cohort_samples = set(phenotype_samples_df.sample_names)
    icd_data = []

    for icdc in icd_codes_df.coding.values:
        icdc_node = pheno_tree.node_dict[c2nodeid_dict[icdc]]
        comorbid_samples = icdc_node.get_samples()
        comorbid_samples = all_cohort_samples.intersection(comorbid_samples)
        df = get_table_icd(combo_samples, noncombo_samples, comorbid_samples, icdc_node.meaning)
        res = fisher_exact(df)
        or_study = odds_ratio(df)
        cil, cih = or_study.confidence_interval(confidence_level=0.95)
        icdc_node_data = (icdc, icdc_node.meaning, df.iloc[0,0], df.iloc[0,1], df.iloc[1,0], df.iloc[1,1], or_study.statistic, res.pvalue, cil, cih)
        icd_data.append(icdc_node_data)
    
    icd_df = pd.DataFrame(icd_data, columns=["icd_code", "icd_meaning", "combo_comorbid", "combo_noncomorbid", "noncombo_comorbid", "noncombo_noncomorbid", "odds_ratio", "p_value", "ci_low", "ci_high"])
    icd_df["FDR"] = stats.false_discovery_control(icd_df.p_value)
    return icd_df

def create_phenotype_samples_df(phenotype_df, combo_samples):
    phenotype_df["carrier"] = phenotype_df.sample_names.isin(combo_samples)
    return phenotype_df

def get_level(icd_node, level):
    if icd_node.parent:
        return get_level(icd_node.parent, level+1)
    else:
        return level

def add_more_info(pheno_tree, code, odds, fdr, c2nodeid_dict):
    icd_node = pheno_tree.node_dict[c2nodeid_dict[code]]
    icd_node.fdr = fdr
    icd_node.odds = odds
    icd_node.level = get_level(icd_node, 0)
    return


def get_icd_enrich_plot(icd_df, pheno_tree, root_pheno, c2nodeid_dict, save_dir):
    for num, row in icd_df.iterrows():
        code, odds, fdr = row.icd_code, row.odds_ratio, row.FDR
        add_more_info(pheno_tree, code, odds, fdr, c2nodeid_dict)
    root_pheno.fdr=1
    root_pheno.level=0
    root_pheno.odds=1

    df_columns = ["level", "node_id", "num_samples", "parent", "meaning", "coding", "fdr", "odds_ratio"]
    df_data_dict = {c:[] for c in df_columns}

    def save_data(curr_node):
        df_data_dict["level"].append(curr_node.level)
        df_data_dict["node_id"].append(curr_node.node)
        df_data_dict["num_samples"].append(len(curr_node.samples))
        if curr_node.parent:
            df_data_dict["parent"].append(curr_node.parent.meaning)
        else:
            df_data_dict["parent"].append("")
        df_data_dict["meaning"].append(curr_node.meaning)
        df_data_dict["coding"].append(curr_node.code)
        df_data_dict["fdr"].append(curr_node.fdr)
        df_data_dict["odds_ratio"].append(curr_node.odds)
        return

    def add_all_parents(curr_node, saved_nodes):
        if curr_node.parent:
            if curr_node.parent not in saved_nodes:
                save_data(curr_node.parent)
                saved_nodes.add(curr_node.parent)
                return add_all_parents(curr_node.parent, saved_nodes)
        return saved_nodes
    
    saved_nodes = set()
    for i,row in icd_df.loc[icd_df.FDR<0.05].iterrows():
        code = row.icd_code
        icd_node = pheno_tree.node_dict[c2nodeid_dict[code]]
        if icd_node not in saved_nodes:
            save_data(icd_node)
            saved_nodes.add(icd_node)
            saved_nodes.update(add_all_parents(icd_node, saved_nodes))

    # filter the phenotypes
    max_level=2
    df_data = pd.DataFrame(df_data_dict)
    df_data = df_data.loc[(df_data.level>0)&(df_data.level<max_level+1)].replace("Root Phenotype", "")
    df_data = df_data.loc[~((df_data.level==max_level)&(df_data.fdr>0.01))]
    df_data["fdr_rank"] = df_data.fdr.rank()
    df_data["fdr_rank_inv"] = df_data.fdr_rank.apply(lambda x: 1/x)
    save_file = os.path.join(save_dir, "icd_enrichment_filtered.csv")
    df_data.to_csv(save_file, index=False)

    fig2 = px.sunburst(
        df_data, names="meaning", parents="parent", values='fdr_rank_inv', color='meaning', 
        width=1000, height=850, )
        # color_continuous_scale='RdBu', range_color=[0., 0.01], color_continuous_midpoint=1e-9)

    save_file = os.path.join(save_dir, "icd_enrichment_filtered.pdf")
    fig2.update_layout(title_text="Two-level Sunburst Diagram", font_size=10)
    fig2.write_image(save_file)
    return df_data

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

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--icd_raw_dir", type=str, help="Filepath of the icd codes to samples files storage dir")
    parser.add_argument("--icd_codes_file", type=str, help="Filepath of the icd codes to parent node file")
    parser.add_argument("--hes_info_file", type=str, help="Filepath of the in patient info file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names",  "bmi", "bmi_prs", "bmi_residuals"], dtype={"sample_names": str, "bmi": float, "bmi_prs": float, "bmi_residuals":float})
    cohort_df["sample_names"] = cohort_df.sample_names.astype(str)
    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    phenotype_samples_df = create_phenotype_samples_df(cohort_df, combo_samples)
    icd_samples_df = utpa.create_icd_samples_file(cli_args.icd_raw_dir)
    icd_codes_df = pd.read_csv(cli_args.icd_codes_file, usecols=["coding", "meaning", "node_id", "parent_id"], sep="\t")
    icd_codes_df["coding"] = icd_codes_df.coding.str.replace(" ", "")
    pheno_tree, root_pheno, c2nodeid_dict = utpa.create_tree(icd_codes_df, icd_samples_df)
    hes_info_df = pd.read_csv(cli_args.hes_info_file, dtype={"sample_names": str, "hes_info": float})
    all_icd_samples = set(hes_info_df.loc[hes_info_df.hes_info>0, "sample_names"].values)
    print(f"Samples in cohort {len(phenotype_samples_df)}")
    phenotype_samples_df = phenotype_samples_df.loc[phenotype_samples_df.sample_names.isin(all_icd_samples)]
    print(f"Samples with icd {len(phenotype_samples_df)}")
    print(f"Samples with combo {len(combo_samples)}")
    combo_samples = combo_samples.intersection(all_icd_samples)
    print(f"Samples with combo and icd {len(combo_samples)}")

    # enrichment of icd codes within combo samples
    icd_df = get_icd_enrich(phenotype_samples_df, pheno_tree, icd_codes_df, c2nodeid_dict)
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    icd_df.to_csv(cli_args.save_file)
    save_dir = os.path.dirname(cli_args.save_file)
    icd_enrich_df = get_icd_enrich_plot(icd_df, pheno_tree, root_pheno, c2nodeid_dict, save_dir)

    # creating a feature file for selected icd comorbidities
    icd_enrich_df = icd_enrich_df.loc[icd_enrich_df.level==2].sort_values("fdr_rank")
    icd_codes_of_interest = set(icd_enrich_df.coding)
    samples_of_interest = set(phenotype_samples_df.sample_names.values)
    icd_samples_df = create_icd_comorbidity_features(icd_codes_of_interest, samples_of_interest, pheno_tree, c2nodeid_dict)
    save_file = os.path.join(save_dir, "icd_features.csv.gz")
    icd_samples_df.to_csv(save_file)
    save_file = os.path.join(save_dir, "cohort_pheno.csv.gz")
    phenotype_samples_df.to_csv(save_file)
