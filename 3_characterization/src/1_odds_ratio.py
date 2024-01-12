import os
import pandas as pd
import argparse
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

import utils.parsing as utpa
import utils.plotting as utpl

def create_contingency_table(study_samples, case_samples, control_samples):
    all_study_case_samples = case_samples.intersection(study_samples)
    all_nonstudy_case_samples = case_samples.difference(study_samples)
    all_study_cont_samples = control_samples.intersection(study_samples)
    all_nonstudy_cont_samples = control_samples.difference(study_samples)
    contingency_table = np.array([[len(all_study_case_samples), len(all_nonstudy_case_samples)], [len(all_study_cont_samples), len(all_nonstudy_cont_samples)]])
    return contingency_table


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Odds ratio")
    parser.add_argument("--case_controls_file", type=str, help="Filepath of the case control file")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--other_study_files", type=str, help="Filepath of the gene files from other studies", nargs="+")
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")


    cli_args = parser.parse_args()

    case_cont_df = pd.read_csv(cli_args.case_controls_file)
    combo_dfs = [pd.read_csv(cf) for cf in cli_args.combo_files]
    genotype_df = pd.read_csv(cli_args.genotype_file)

    case_samples, control_samples = utpa.get_case_cont_samples(case_cont_df)
    other_genes = [utpa.get_gene_set_from_file(gf) for gf in cli_args.other_study_files]
    other_samples = [utpa.get_samples_with_gene_mutation(genotype_df, gs) for gs in other_genes]
    combo_samples = [utpa.get_combo_samples_from_file(cf) for cf in cli_args.combo_files]

    cont_tables = [create_contingency_table(sg, case_samples, control_samples) for sg in other_samples + combo_samples]
    odds_ratios = [odds_ratio(ct) for ct in cont_tables]
    fishers = [fisher_exact(ct) for ct in cont_tables]

    data_dict = {"study_type": [], "odds_ratio": [], "ci_low": [], "ci_high": [], "pvalue": []}
    for or_study, fe_study, st in zip(odds_ratios, fishers, ["Akbari et. al.", "Turcot et. al. ", "Digenic Combinations", "Trigenic Combinations"]):
        data_dict["study_type"].append(st)
        data_dict["odds_ratio"].append(or_study.statistic)
        cil, cih = or_study.confidence_interval(confidence_level=0.95)
        data_dict["ci_low"].append(cil)
        data_dict["ci_high"].append(cih)
        data_dict["pvalue"].append(fe_study.pvalue)

    df = pd.DataFrame(data_dict)
    df["study_type_group"] = ["RVAS", "RVAS", "RareComb", "RareComb"]

    fig = utpl.create_odds_ratio_plot(df)
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    utpl.save_pdf(cli_args.save_file, fig.figure)
