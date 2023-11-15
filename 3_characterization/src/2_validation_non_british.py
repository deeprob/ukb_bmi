import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import kstest,ttest_ind, ks_2samp
import seaborn as sns
from matplotlib.ticker import MultipleLocator
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 14, 'axes.linewidth': 2, 'xtick.major.width': 1.5, 'xtick.major.size': 7, 'ytick.major.width': 1.5, 'ytick.major.size': 7})
from matplotlib.backends.backend_pdf import PdfPages


def save_pdf(save_file, fig):
    pdf = PdfPages(save_file)
    pdf.savefig(fig, bbox_inches='tight', dpi=300)
    pdf.close()
    return


def plot_box(boxdf):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x="combo_carriers",
        y="bmi",
        order=[False, True],
        hue_order=[False,True],
        gap=0.5,
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[7]],
        dodge=False, width=0.75, linewidth=2.25, fliersize=0, capprops={'color':'none'}, 
        boxprops={ 'edgecolor':'k'},  # 'facecolor':'none',
        whiskerprops={'color':'k'}, medianprops={'color':'k'}) # 

    # Adjust Axis
    # ax.set_yticks([-0.02, 0, 0.02, 0.04])
    ax.set_xlim((-1, 2))
    ax.set_ylim((12, 45))
    # ax.set_ylabel('Percentage')
    ax.set_xticklabels(["Non\ncarrier", "Carrier"], rotation=45, ha="center", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("BMI")
    ax.hlines(43, 0, 1, color="k")
    ttest_pval = ttest_ind(non_combo_hit_pheno.bmi, combo_hit_pheno.bmi, alternative="less").pvalue
    ax.text(0.5, 44, f"P = {ttest_pval:.2E}", ha="center", va="bottom", fontsize=14)

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--combo_samples_file", type=str, help="Filepath of the cohort samples file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi", "bmi_prs"])
    combo_samples_df = pd.read_csv(cli_args.combo_samples_file)

    all_combo_samples = set("|".join(combo_samples_df.combo_samples.values).split("|"))
    cohort_df["combo_carriers"] = cohort_df.sample_names.astype(str).isin(all_combo_samples)
    combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==True]
    non_combo_hit_pheno = cohort_df.loc[cohort_df.combo_carriers==False]

    fig, ax = plot_box(cohort_df)

    save_pdf(cli_args.save_file, fig)
