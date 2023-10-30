import os
import argparse
from functools import reduce
import numpy as np
import pandas as pd
from scipy.stats import kstest,ttest_ind
import itertools as it
from collections import Counter
import seaborn as sns
from matplotlib.ticker import MultipleLocator
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 14, 'axes.linewidth': 2, 'xtick.major.width': 1.5, 'xtick.major.size': 7, 'ytick.major.width': 1.5, 'ytick.major.size': 7})
from matplotlib.backends.backend_pdf import PdfPages
from functools import reduce
from scipy.stats import kstest,ttest_ind


def save_pdf(save_file, fig):
    pdf = PdfPages(save_file)
    pdf.savefig(fig, bbox_inches='tight', dpi=300)
    pdf.close()
    return

def get_combo_samples(combos, genotype_df):
    samples_per_gene = genotype_df.loc[genotype_df.gene.isin(combos)].samples.str.split(",").values
    samples_per_combo = reduce(lambda a,b: set(a).intersection(set(b)), samples_per_gene)
    return "|".join(sorted(samples_per_combo))

def get_sample_bmi_info(samples, phenotype_df):
    samples = samples.split("|")
    sample_phenotypes = phenotype_df.loc[phenotype_df.sample_names.astype(str).isin(samples)]
    sample_bmi_dict = dict(zip(sample_phenotypes.sample_names.astype(str), sample_phenotypes.bmi))
    sample_bmi_prs_dict = dict(zip(sample_phenotypes.sample_names.astype(str), sample_phenotypes.bmi_prs))
    sample_bmi = "|".join([str(round(sample_bmi_dict[s], 4))  if s in sample_bmi_dict.keys() else "NA" for s in samples])
    sample_bmi_prs = "|".join([str(round(sample_bmi_prs_dict[s], 4)) if s in sample_bmi_dict.keys() else "NA" for s in samples ])
    return pd.Series({"combo_samples_bmi": sample_bmi, "combo_samples_bmi_prs": sample_bmi_prs})

def get_mean_values(sample_bmi_info):
    sample_bmi_info = sample_bmi_info.split("|")
    sample_bmi_info = [float(sb) for sb in sample_bmi_info if sb!="NA" ]
    return np.mean(sample_bmi_info)

def get_combos_with_sample_info(combo_df, phenotype_df, genotype_df):
    combo_df["combos"] = combo_df.uniq_items.apply(lambda x: [g.replace("Input_", "", 1) for g in x.split("|")])
    combo_df["combo_samples"] = combo_df.combos.apply(get_combo_samples, args=(genotype_df, ))
    combo_df = pd.concat((combo_df, combo_df.combo_samples.apply(get_sample_bmi_info, args=(phenotype_df, ))), axis=1)
    combo_df["mean_bmi"] = combo_df.combo_samples_bmi.apply(get_mean_values)
    combo_df["mean_bmi_prs"] = combo_df.combo_samples_bmi_prs.apply(get_mean_values)
    return combo_df.loc[:, ["uniq_items", "combo_samples", "combo_samples_bmi", "combo_samples_bmi_prs", "mean_bmi", "mean_bmi_prs"]]

def create_decile_rank_plot(phenotype_samples_df):
    fig,ax = plt.subplots(2, 1, figsize=(10, 12), sharex=True, height_ratios=(4, 4))
    ## BMI prs box plots
    sns_ax = sns.boxplot(
        phenotype_samples_df, x="bmi_decile", y="bmi_prs", hue="description", hue_order=["Non Combo", "Combo"],
        width=0.65, linewidth=1.25, fliersize=0, capprops={'color':'none'}, boxprops={'edgecolor':'k'},
        medianprops={'color':'k'},
        linecolor='k',
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[7]],
        legend=True, gap=0.25, ax=ax[0])
    ax[0].set_ylim(-3.5, 3.5)
    ax[0].spines[['right', 'top']].set_visible(False)
    ax[0].hlines([3 for i in range(10)], [(i-0.25+i*0.001953125) for i in range(0, 10)], [(i+0.25+i*0.001953125) for i in range(0,10)], color="k")
    for i, (psdi, psd) in enumerate(phenotype_samples_df.groupby("bmi_decile", observed=False)):
        ttest_res = ttest_ind(psd.loc[psd.description=="Combo", "bmi_prs"], psd.loc[psd.description=="Non Combo", "bmi_prs"], alternative="less")
        ttest_pval = ttest_res.pvalue
        if ttest_pval<0.05:
            pval_text = "*"
            if ttest_pval<0.001:
                pval_text = "**"
                # pval_text = f"P = {ttest_pval:.1E}"
            ax[0].text(0.+i, 3.05, pval_text, ha="center", va="bottom", fontsize=14)
    h,l = ax[0].get_legend_handles_labels()
    ax[0].legend_.remove()
    fig.legend(h,l, ncol=2, loc=(0.35, 0.5025), frameon=False)
    ax[0].set_ylabel("BMI PRS")
    
    phenotype_combo_samples_decile_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Combo"].groupby("bmi_decile", observed=True).agg({
        "sample_names": "count",
        "bmi_prs": "median",
        "bmi": "mean"}
        ).reset_index()
    phenotype_combo_samples_decile_df["decile_rank"] = range(10)
    decile_labels =  [" - ".join([str(round(v.left, 2)),str(round(v.right, 2))]) for v in  phenotype_combo_samples_decile_df.bmi_decile.values]

    ## Combo per decile barplot
    ax[1].bar(phenotype_combo_samples_decile_df.decile_rank, phenotype_combo_samples_decile_df.sample_names, width=0.75, color=sns.color_palette("Reds", 15).as_hex()[:10], edgecolor="k")
    xticklabels = decile_labels
    ax[1].set_xticks(range(10), xticklabels, rotation=45, ha="center", fontsize=11)
    ax[1].set_xlabel("BMI")
    ax[1].set_ylabel("Number of individuals")
    # g.bar_label([g.containers[i] for i in range])
    rects = ax[1].patches
    # Make some labels.
    for rect in rects:
        height = rect.get_height()
        ax[1].text(
            rect.get_x() + rect.get_width() / 2, height + 5, f"{height}", ha="center", va="bottom"
        )
    ax[1].set_ylim(-100, 3500)
    ax[1].spines[['right', 'top']].set_visible(False)

    return fig,ax

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs='+')
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--phenotype_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()

    genotype_df = pd.read_csv(cli_args.genotype_file)
    phenotype_df = pd.read_csv(cli_args.phenotype_file, usecols=["sample_names", "bmi", "bmi_prs"])
    combo_dfs = [pd.read_csv(cf, usecols=["uniq_items"]) for cf in cli_args.combo_files]
    combo_df = pd.concat(combo_dfs).reset_index(drop=True)

    combo_df = get_combos_with_sample_info(combo_df, phenotype_df, genotype_df)
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    combo_df.to_csv(cli_args.save_file, index=False)

    phenotype_df["bmi_decile"] = pd.qcut(phenotype_df.bmi, q=10)

    cohort_df = combo_df
    all_combo_samples = set("|".join(cohort_df.combo_samples.values).split("|"))

    phenotype_combo_samples_df = phenotype_df.loc[phenotype_df.sample_names.astype(str).isin(list(map(str, all_combo_samples)))]
    phenotype_other_samples_df = phenotype_df.loc[~phenotype_df.sample_names.astype(str).isin(list(map(str, all_combo_samples)))]

    phenotype_other_samples_df["description"] = "Non Combo"
    phenotype_combo_samples_df["description"] = "Combo"
    phenotype_samples_df = pd.concat((phenotype_combo_samples_df, phenotype_other_samples_df))

    fig, ax = create_decile_rank_plot(
        phenotype_samples_df
        )

    save_name = os.path.splitext(os.path.basename(cli_args.save_file))[0]
    save_file = os.path.join(os.path.dirname(cli_args.save_file), "figures", f"{save_name}_combo_per_decile.pdf")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    save_pdf(save_file, fig)
