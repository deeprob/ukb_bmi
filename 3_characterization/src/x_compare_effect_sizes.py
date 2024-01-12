import numpy as np
import os
import argparse
import pandas as pd
from statsmodels.stats.proportion import proportion_effectsize
from scipy.stats import kstest,ttest_ind
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
    pdf.savefig(fig, bbox_inches='tight')
    pdf.close()
    return


def read_genes(gene_file):
    with open(gene_file, "r") as f:
        genes = [g.strip() for g in f.readlines()]
    return set(genes)


def get_overlapping_samples(sample_set, overlap_set):
    sample_set = set(sample_set.split(","))
    return len(sample_set.intersection(overlap_set))/len(overlap_set)


def get_effect_size(gene_list, case_control_df, genotype_df, combo_samples, save_file):
    genotype_df = genotype_df.loc[genotype_df.gene.isin(gene_list)]
    case_samples = set(case_control_df.loc[case_control_df.Output_BMI==1].Sample_Name.astype(str).values)
    control_samples = set(case_control_df.loc[case_control_df.Output_BMI==0].Sample_Name.astype(str).values)
    # from case and control, eliminate samples who have the combinations
    case_samples = case_samples.difference(combo_samples)
    control_samples = control_samples.difference(combo_samples)
    genotype_df["prop_case_samples"] = genotype_df.samples.apply(get_overlapping_samples, args=(case_samples, ))
    genotype_df["prop_control_samples"] = genotype_df.samples.apply(get_overlapping_samples, args=(control_samples, ))
    genotype_df["Effect_Size"] = genotype_df.apply(lambda x: proportion_effectsize(x.prop_case_samples, x.prop_control_samples), axis=1)
    genotype_df.loc[:, ["gene", "Effect_Size"]].to_csv(save_file, index=False)
    return

def create_efs_table(files, filenames):
    dfs = [pd.read_csv(fl).loc[:, "Effect_Size"].to_frame() for fl in files]
    for ifn, fn in enumerate(filenames):
        dfs[ifn]["Description"] = fn
    df = pd.concat(dfs)
    return df

def plot_efs(boxdf, order):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    pnts = np.linspace(0, np.pi * 2, 24)
    circ = np.c_[np.sin(pnts) / 2, -np.cos(pnts) / 2]
    vert = np.r_[circ, circ[::-1] * .7]

    open_circle = matplotlib.path.Path(vert)

    # Box Plot
    sns_strip = sns.stripplot(
        data=boxdf,
        color="lightgrey",
        palette=["#008176", "#08509b", "#ff7f0e", "#a30f15", ],  # '#D1245D', '#00ADEE', '#D1245D', '#00ADEE' "#0000a7"   
        x="Description",
        y="Effect_Size",
        order=order,
        orient="v",
        s=4,
        marker=open_circle,
        alpha=0.75, linewidth=0.1, facecolor=(0,0,0,0), dodge=False, ax=ax, jitter=0.15,  # ec='none', 
        )

    sns_box = sns.boxplot(
        data=boxdf,
        palette=["#008176", "#00ADEE",  "#ff7f0e", '#D1245D', "#FF3688", "#0000a7" "#eecc16"], 
        x="Description",
        y="Effect_Size",
        order=order,
        hue_order=order,
        dodge=False, width=0.5, linewidth=2.25, fliersize=0, capprops={'color':'none'}, boxprops={'facecolor':'none', 'edgecolor':'k'}) # 


    # Adjust Axis
    ax.set_yticks([-0.02, 0, 0.02, 0.04])
    ax.set_ylim((-0.025, 0.05))
    # ax.set_ylabel('Percentage')
    ax.set_xticklabels(order, rotation=45, )
    ax.set_xlabel("")

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--gene_files", type=str, help="Filepath of rare variant genes in other studies", nargs='+')
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs='+')
    parser.add_argument("--genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("--case_control_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("--save_files", type=str, help="Filepath where combo info will be stored", nargs='+')
    parser.add_argument("--lifestyle_file", type=str, help="Filepath of the lifestyle factor matrix", default="")

    cli_args = parser.parse_args()

    case_control_df = pd.read_csv(cli_args.case_control_file)
    genotype_df = pd.read_csv(cli_args.genotype_file)
    if cli_args.lifestyle_file:
        lifestyle_df = pd.read_csv(cli_args.lifestyle_file)
        cols_to_melt = list(lifestyle_df.columns)
        lifestyle_df = lifestyle_df.melt(id_vars="Sample_Name", value_vars=cols_to_melt, var_name="gene")
        lifestyle_df = lifestyle_df.loc[lifestyle_df.value>0].drop(columns="value")
        lifestyle_df = lifestyle_df.groupby("gene").agg(lambda x: ",".join(map(str, x))).reset_index().rename(columns={"Sample_Name": "samples"})
        genotype_df = pd.concat((genotype_df, lifestyle_df)).reset_index(drop=True)

    combo_dfs = [pd.read_csv(cf).loc[:, ["Case_Samples", "Control_Samples"]] for cf in cli_args.combo_files]
    combo_df = pd.concat(combo_dfs)
    combo_samples = set("|".join(combo_df.values.flatten().astype(str)).split("|"))

    for glf,sf in zip(cli_args.gene_files, cli_args.save_files):
        gene_list = read_genes(glf)
        os.makedirs(os.path.dirname(sf), exist_ok=True)
        get_effect_size(gene_list, case_control_df, genotype_df, combo_samples, sf)

    files = cli_args.save_files + cli_args.combo_files
    filenames = [os.path.splitext(os.path.basename(f))[0] for f in files]
    efs_df = create_efs_table(files, filenames).reset_index(drop=True)
    fig,ax =  plot_efs(efs_df, filenames)

    figsave_file = os.path.join(os.path.dirname(cli_args.save_files[0]), "figures", "efs.pdf")
    os.makedirs(os.path.dirname(figsave_file), exist_ok=True)
    save_pdf(figsave_file, fig)
