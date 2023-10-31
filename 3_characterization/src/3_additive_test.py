import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
import matplotlib.patches as mpatches
matplotlib.rcParams['font.sans-serif'] = "Arial" # missing fonts:: https://alexanderlabwhoi.github.io/post/2021-03-missingfont/
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 14, 'axes.linewidth': 2, 'xtick.major.width': 1.5, 'xtick.major.size': 7, 'ytick.major.width': 1.5, 'ytick.major.size': 7})
from matplotlib.backends.backend_pdf import PdfPages


def save_pdf(save_file, fig):
    pdf = PdfPages(save_file)
    pdf.savefig(fig, bbox_inches='tight')
    pdf.close()
    return


def get_uniq_items(ser):
    return set("|".join(ser.values).split("|"))

def prepare_model_data(genotype_df, phenotype_df, combo_genes, combo_samples):
    # only keep the combo genes
    combo_gene_df = genotype_df.loc[genotype_df.gene.isin([g.replace("Input_", "") for g in combo_genes])]
    # set samples as str type
    combo_gene_df.loc[:, "samples"] = combo_gene_df.samples.str.split(",")
    # explode by samples
    combo_gene_df = combo_gene_df.explode("samples")
    # create one hot encoded data
    combo_gene_df = pd.crosstab(combo_gene_df.samples, combo_gene_df.gene)
    # only keep sample that have a mutation in at least one of the combo genes
    combo_gene_df = combo_gene_df.loc[combo_gene_df.sum(axis=1)>0]
    # add sample phenotype info
    combo_gene_df = combo_gene_df.merge(phenotype_df, left_on="samples", right_on="sample_names").set_index("sample_names")
    # create test data: this will contain samples with the combination
    test_combo_gene_df = combo_gene_df.loc[combo_gene_df.index.isin(combo_samples)]
    # create train data: this will contain sample without any combination
    train_combo_gene_df = combo_gene_df.loc[~combo_gene_df.index.isin(combo_samples)].dropna()
    return train_combo_gene_df, test_combo_gene_df


def get_features_labels(df, pheno_name):
    y = df[pheno_name].values.reshape(-1,1)
    X = df.loc[:, [c for c in df.columns if c!=pheno_name]].values
    return X,y

def create_additive_plot(plot_df):
    fig,ax = plt.subplots(1,1,figsize=(6,4))

    g = sns.histplot(
        data=plot_df, x='value', 
        hue='variable', hue_order=["Observed", "Expected"],
        edgecolor="k",
        linewidth=0.45,
        palette=["whitesmoke", "whitesmoke"], kde=True, # sns.color_palette("Reds", 15).as_hex()[7]
        element="bars", fill=True, bins=150, line_kws={"linewidth":2, "linestyle":"solid"}, legend=False,
        ax=ax)
    g.axes.lines[1].set_color("#a30f15")
    g.axes.lines[0].set_color("#08509b")
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.set_xlim(20, 45)
    g.set_xlabel("BMI")
    ax.spines[['right', 'top']].set_visible(False)

    # legend
    red_patch = mpatches.Patch(color='#a30f15', label='Observed')
    blue_patch = mpatches.Patch(color='#08509b', label='Expected')
    ax.legend(handles=[blue_patch, red_patch], handlelength=2, handleheight=0.001, frameon=False)
    return fig, ax

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("genotype_file", type=str, help="Filepath of the gene burden file")
    parser.add_argument("phenotype_file", type=str, help="Filepath of the cohort phenotype file")
    parser.add_argument("combo_gene_file", type=str, help="Filepath of the combo info file")
    parser.add_argument("save_file", type=str, help="Filepath where combo info will be stored")

    cli_args = parser.parse_args()
    pheno_name = "bmi"

    genotype_df = pd.read_csv(cli_args.genotype_file)
    phenotype_df = pd.read_csv(cli_args.phenotype_file, usecols=["sample_names", pheno_name])
    phenotype_df.loc[:, "sample_names"] = phenotype_df.sample_names.astype(str)
    combo_gene_df = pd.read_csv(cli_args.combo_gene_file)

    combo_genes = get_uniq_items(combo_gene_df.uniq_items)
    combo_samples = get_uniq_items(combo_gene_df.combo_samples)

    train_combo_gene_df, test_combo_gene_df = prepare_model_data(genotype_df, phenotype_df, combo_genes, combo_samples)
    X_train, y_train = get_features_labels(train_combo_gene_df, pheno_name)
    X_test, y_test = get_features_labels(test_combo_gene_df, pheno_name)
    scaler = StandardScaler()
    scaler.fit(y_train)
    y_train_scaled = scaler.transform(y_train)

    model = LinearRegression()
    model.fit(X_train, y_train_scaled)
    y_pred = model.predict(X_test)
    y_pred = scaler.inverse_transform(y_pred)

    test_combo_gene_df[f"{pheno_name}_pred"] = y_pred
    os.makedirs(os.path.dirname(cli_args.save_file), exist_ok=True)
    test_combo_gene_df.loc[:, [f"{pheno_name}", f"{pheno_name}_pred"]].to_csv(cli_args.save_file, index=True)

    plot_df = test_combo_gene_df.reset_index().rename(columns={f"{pheno_name}": "Observed", f"{pheno_name}_pred": "Expected"}).melt(id_vars="sample_names")
    fig, ax = create_additive_plot(plot_df)

    save_name = os.path.splitext(os.path.basename(cli_args.save_file))[0]
    save_file = os.path.join(os.path.dirname(cli_args.save_file), "figures", f"{save_name}_expected_observed.pdf")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    save_pdf(save_file, fig)
