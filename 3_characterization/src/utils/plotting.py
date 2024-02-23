import os
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial" # missing fonts:: https://alexanderlabwhoi.github.io/post/2021-03-missingfont/
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams.update({'font.size': 14, 'axes.linewidth': 2, 'xtick.major.width': 1.5, 'xtick.major.size': 7, 'ytick.major.width': 1.5, 'ytick.major.size': 7})
from matplotlib.backends.backend_pdf import PdfPages
import forestplot as fp
import upsetplot

from scipy.stats import ttest_ind


##################
# plotting utils #
##################
def save_pdf(save_file, fig):
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    pdf = PdfPages(save_file)
    pdf.savefig(fig, bbox_inches='tight',dpi=300)
    pdf.close()
    return


###################
# odds ratio plot #
###################
def create_odds_ratio_plot(df, figsize=(6,3)):
    fig = fp.forestplot(
        df,  # the dataframe with results data
        estimate="odds_ratio",  # col containing estimated effect size 
        ll="ci_low", hl="ci_high",
        pval="pvalue",
        decimal_precision=3,
        varlabel="study_type",  # column containing variable label
        groupvar="study_type_group",  # Add variable groupings
        color_alt_rows=True,
        table=True,
        annote=["est_ci"],
        annoteheaders=["Est. (95% Conf. Int.)"],
        # group ordering
        group_order=["RVAS", "RareComb"],
        sort=True, # sort in ascending order (sorts within group if group is specified)
        logscale=True, 
        ylabel="Est. (95% Conf. Int.)",  # y-label title
        xlabel="Odds Ratio",  # x-label title
        xticks=[0.1, 1, 10, 100],
        # xticks=[-5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50],
        figsize=figsize,
        **{"marker": "D",  # set maker symbol as diamond
            "markersize": 75,  # adjust marker size
            "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
            "xlinecolor": "k",  # gray color for x-reference line
            "xline": 1,
            "xlinewidth": 1,
            "xtick_size": 12,  # adjust x-ticker fontsize
            "lw": 3,
            }       
        )
    return fig


###########################
# variance explained plot #
###########################
def plot_variance_explained(plot_df):
    fig, axes = plt.subplots(1, 1, figsize=(6, 2))
    axes.barh(
        plot_df["Model Variates"],
        plot_df["Variance explained"],
        height=0.5,
        edgecolor="k", linewidth=2,
        color="white"
        )
    axes.set_xlim(0, 100)
    axes.set_ylim(-0.5, 2.5)
    axes.set_xlabel("BMI Variance Explained (%)")
    axes.set_ylabel("")
    axes.spines[['top', 'right']].set_visible(False)

    rects =axes.patches
    # Make some labels.
    for rect in rects:
        width = rect.get_width()
        axes.text(
            width + 1, rect.get_y() + rect.get_height() / 2,  f"{round(width, 3)}%", ha="left", va="center"
        )
    return fig


###########################
# combo carriers bmi plot #
###########################
def plot_box_single_bmi(
        boxdf, ttest_pval, 
        xvar="combo_carriers", yvar="bmi", ylim=(12, 45), 
        xticklabel=["Non\ncarrier", "Carrier"], ttest_hline=43, ttest_text=44):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
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
    ax.set_ylim(ylim)
    # ax.set_ylabel('Percentage')
    ax.set_xticklabels(xticklabel, rotation=45, ha="center", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("BMI")
    ax.hlines(ttest_hline, 0, 1, color="k")
    ax.text(0.5, ttest_text, f"P = {ttest_pval:.2E}", ha="center", va="bottom", fontsize=14)

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax


###################
# oligogenic plot #
###################
def plot_box_oligo(
        boxdf, ttest_pvals,
        xvar="mutation", yvar="bmi", huevar=None, 
        order=["Single Hit", "Combo carriers"], hue_order=["Single Hit", "Combo carriers"], 
        figsize=(4, 6),  ylim=(10, 50),
        ttest_hline=45, ttest_text=46):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=figsize)
    dodge = False if len(ttest_pvals)==1 else True
    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
        hue=huevar,
        order=order,
        hue_order=hue_order,
        gap=0.5,
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[7]],
        dodge=dodge, width=0.75, linewidth=2.25, fliersize=0, capprops={'color':'none'}, 
        boxprops={ 'edgecolor':'k'},  # 'facecolor':'none',
        whiskerprops={'color':'k'}, medianprops={'color':'k'}) # 

    # Adjust Axis
    # ax.set_yticks([-0.02, 0, 0.02, 0.04])
    ax.set_xlim((-1, len(ttest_pvals)+1))
    ax.set_ylim(ylim)
    # ax.set_ylabel('Percentage')
    # ax.set_xticklabels(xticklabel, rotation=45, ha="center", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("BMI")
    shift = 0.25 if len(ttest_pvals)>1 else 1
    textshift = 0 if len(ttest_pvals)>1 else 0.5
    for i, ttest_pval in enumerate(ttest_pvals):
        ax.hlines(ttest_hline, i-shift, i+shift, color="k")
        ax.text(i+textshift, ttest_text, f"P = {ttest_pval:.2E}", ha="center", va="bottom", fontsize=14)

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax


#######################
# additive model plot #
#######################
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


####################
# interaction plot #
####################
def get_interaction_plot(plot_df):
    fig, axes = plt.subplots(figsize=(10,8))

    sns.pointplot(
        data=plot_df,
        x="Gene A", y="BMI", hue="Gene B",
        markers=["o", "s"], linestyles=[":", "-"], errorbar=("se", 10),
        palette=["navy", "red"], markersize=15,
        dodge=0.075,
        ax=axes
    )

    axes.set_ylim(26, 32, auto=False)

    axes.set_xticklabels(axes.get_xticklabels(), fontsize=18)
    axes.set_yticklabels(axes.get_yticklabels(), fontsize=18)

    axes.set_xlabel("Gene A", fontsize=20)
    axes.set_ylabel("BMI", fontsize=20)

    # Remove Spines
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)


    plt.setp(axes.get_legend().get_texts(), fontsize='14') # for legend text
    plt.setp(axes.get_legend().get_title(), fontsize='18') # for legend title
    return fig, axes


########################
# allelic effects plot #
########################
def plot_box_allelic(
        boxdf, ttest_pval, ttest_pval_cont,
        xvar="mutation", yvar="bmi", ylim=(10, 50), 
        xticklabel=["Controls", "Different\nvariants", "Same\nvariants"], ttest_hline=45, ttest_text=46):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
        order=["Controls", "Different variants", "Same variants"],
        hue_order=["Controls", "Different variants", "Same variants"],
        gap=0.5,
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[3], sns.color_palette("Reds", 15).as_hex()[9]],
        dodge=False, width=0.75, linewidth=2.25, fliersize=0, capprops={'color':'none'}, 
        boxprops={ 'edgecolor':'k'},  # 'facecolor':'none',
        whiskerprops={'color':'k'}, medianprops={'color':'k'}) # 

    # Adjust Axis
    # ax.set_yticks([-0.02, 0, 0.02, 0.04])
    ax.set_xlim((-1, 3))
    ax.set_ylim(ylim)
    # ax.set_ylabel('Percentage')
    ax.set_xticklabels(xticklabel, rotation=45, ha="center", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("BMI")
    ax.hlines(ttest_hline + 5, 0, 2, color="k")
    ax.text(0.5, ttest_text + 5, f"P = {ttest_pval_cont:.2E}", ha="center", va="bottom", fontsize=14)
    ax.hlines(ttest_hline, 1, 2, color="k")
    ax.text(0.5, ttest_text, f"P = {ttest_pval:.2E}", ha="center", va="bottom", fontsize=14)

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax

def plot_box_variant_types(
        boxdf,
        xvar="vtype_parsed", yvar="bmi", ylim=(10, 50), 
        xticklabel=["missense only", "both lof and missense", "lof only"]):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
        order=xticklabel,
        hue_order=xticklabel,
        gap=0.5,
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[3], sns.color_palette("Reds", 15).as_hex()[9]],
        dodge=False, width=0.75, linewidth=2.25, fliersize=0, capprops={'color':'none'}, 
        boxprops={ 'edgecolor':'k'},  # 'facecolor':'none',
        whiskerprops={'color':'k'}, medianprops={'color':'k'}) # 

    # Adjust Axis
    # ax.set_yticks([-0.02, 0, 0.02, 0.04])
    ax.set_xlim((-1, 3))
    ax.set_ylim(ylim)
    # ax.set_ylabel('Percentage')
    ax.set_xticklabels(xticklabel, rotation=45, ha="center", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("BMI")

    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False);
    return fig, ax


############################
# bmi residual decile plot #
############################
def create_decile_rank_plot(phenotype_samples_df):
    fig,ax = plt.subplots(2, 1, figsize=(10, 12), sharex=True, height_ratios=(4, 4))
    ## BMI prs box plots
    sns_ax = sns.boxplot(
        phenotype_samples_df, x="bmi_res_decile", y="bmi_prs", hue="description", hue_order=["Non Combo", "Combo"],
        width=0.65, linewidth=1.25, fliersize=0, capprops={'color':'none'}, boxprops={'edgecolor':'k'},
        medianprops={'color':'k'},
        linecolor='k',
        palette= ["whitesmoke", sns.color_palette("Reds", 15).as_hex()[7]],
        legend=True, gap=0.25, ax=ax[0])
    ax[0].set_ylim(-3.5, 3.5)
    ax[0].spines[['right', 'top']].set_visible(False)
    ax[0].hlines([3 for i in range(10)], [(i-0.25+i*0.001953125) for i in range(0, 10)], [(i+0.25+i*0.001953125) for i in range(0,10)], color="k")
    for i, (psdi, psd) in enumerate(phenotype_samples_df.groupby("bmi_res_decile", observed=False)):
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
    
    phenotype_combo_samples_decile_df = phenotype_samples_df.loc[phenotype_samples_df.description=="Combo"].groupby("bmi_res_decile", observed=True).agg({
        "sample_names": "count",
        "bmi_prs": "median",
        "bmi_residuals": "mean"}
        ).reset_index()
    phenotype_combo_samples_decile_df["decile_rank"] = range(10)
    decile_labels =  [" - ".join([str(round(v.left, 2)),str(round(v.right, 2))]) for v in  phenotype_combo_samples_decile_df.bmi_res_decile.values]

    ## Combo per decile barplot
    ax[1].bar(phenotype_combo_samples_decile_df.decile_rank, phenotype_combo_samples_decile_df.sample_names, width=0.75, color=sns.color_palette("Reds", 15).as_hex()[:10], edgecolor="k")
    xticklabels = decile_labels
    ax[1].set_xticks(range(10), xticklabels, rotation=45, ha="center", fontsize=11)
    ax[1].set_xlabel("BMI Residuals")
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

def create_regplot(phenotype_samples_df):
    fig, axes = plt.subplots(figsize=(6,4))

    sns.regplot(
        phenotype_samples_df.loc[phenotype_samples_df.description=="Non Combo"], x="bmi_residuals", y="bmi_prs", scatter=False, 
        line_kws={"ls": "--", "color":"navy"}, ax=axes)
    sns.regplot(
        phenotype_samples_df.loc[phenotype_samples_df.description=="Combo"], x="bmi_residuals", y="bmi_prs", scatter=False, 
        line_kws={"ls": "-", "color":"darkred"}, ax=axes, ci=99)

    axes.set_xlabel("BMI Residuals", fontsize=18)
    axes.set_ylabel("BMI PRS", fontsize=18)

    # axes.set_xlim(10, 80, auto=False)
    axes.set_ylim(-1.5, 3, auto=False)

    axes.vlines(5, -1.5, 3, linestyles=["dashdot"], color="k")
    axes.axvspan(5, 10, alpha=0.25, facecolor="red")
    axes.text(5, 2.75, f"Extreme Obesity", ha="left", va="top", fontsize=14)

    axes.set_xticklabels(axes.get_xticklabels(), fontsize=14)
    axes.set_yticklabels(axes.get_yticklabels(), fontsize=14)

    # Remove Spines
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    return fig, axes

def get_odds_ratio_plot_bmi_extreme(
        plot_df, varlabel="icd_meaning", 
        right_annotate_extra=["combo_comorbid", "combo_noncomorbid"],
        right_annoteheaders_extra=["Comorbid Carriers", "Non-Comorbid Carriers"],
        xticks=[-5, 1, 5, 10, 20, 30, 40, 50]
    ):
    fig = fp.forestplot(plot_df,  # the dataframe with results data
                estimate="odds_ratio",  # col containing estimated effect size 
                ll="ci_low", hl="ci_high",
                pval="p_value",
                decimal_precision=3,
                varlabel=varlabel,  # column containing variable label
                color_alt_rows=True,
                table=True,
                annote=["est_ci"],
                annoteheaders=["Est. (95% Conf. Int.)"],
                rightannote=["formatted_pval"]+right_annotate_extra,
                right_annoteheaders=["P-value"]+right_annoteheaders_extra,
                # group ordering
                sort=False, # sort in ascending order (sorts within group if group is specified)
                logscale=False, 
                ylabel="Est. (95% Conf. Int.)",  # y-label title
                xlabel="Odds Ratio",  # x-label title
                xticks=xticks,
                # xticks=[-5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50],
                figsize=(6, 9 ),
                **{"marker": "D",  # set maker symbol as diamond
                    "markersize": 75,  # adjust marker size
                    "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
                    "xlinecolor": "k",  # gray color for x-reference line
                    "xline": 1,
                    "xlinewidth": 1,
                    "xtick_size": 12,  # adjust x-ticker fontsize
                    "lw": 3,
                    }       
                )

    return fig.figure

def get_prs_bmi_res_decile_plot(plot_df, bmi_dict):
    g = sns.catplot(data=plot_df, 
        x="bmi_res_categories",
        y="percent",
        hue="bmi_prs_categories", 
        kind="bar",
        height=4, aspect=2,
        order=["underweight", "normal", "overweight", "obese", "severe obesity"],
        hue_order=["lowest", "middle", "highest"],
        palette=["silver", "gray", "black"]
        )

    g.ax.set_ylim(0,100)
    g.ax.set_xticklabels([f"{l.get_text()}\n{round(bmi_dict[l.get_text()], 2)}" for l in g.ax.get_xticklabels()])

    rects = g.ax.patches
    for rect in rects:
        height = rect.get_height()
        g.ax.text(
            rect.get_x() + rect.get_width() / 2, height + 1, f"{round(height, 2)}%", ha="center", va="bottom"
        )
    return g.figure


###################
# enrichment plot #
###################
def create_dot_plot(go_file, terms_col="Term", odds_ratio_col="odds_ratio", gene_col="genes", gene_pval_col="adj_pval", ncat=20, figsize=(5,7)):
    # read and parse go file
    go_df = pd.read_csv(go_file)
    go_df[terms_col] = go_df[terms_col].apply(lambda x: x[0].upper() + x[1:])
    go_df["gene_counts"] = go_df[gene_col].apply(lambda x: len(x.split("|")))
    go_df_adj = go_df.loc[go_df[gene_pval_col]<0.05]
    plot_df = go_df_adj.sort_values(gene_pval_col).head(ncat)
    
    plot_df = plot_df.sort_values(odds_ratio_col, ascending=False)
    fig, axes = plt.subplots(figsize=figsize)
    norm = plt.Normalize(0, 0.05) # plot_df['qvalue'].min(), plot_df['qvalue'].max()
    sns_ax = sns.scatterplot(
        data=plot_df, x=odds_ratio_col, y=terms_col, 
        size="gene_counts", hue=gene_pval_col, ax=axes, sizes=(100, 300), palette='RdBu', hue_norm=norm, linewidth=0.5, edgecolor="k"
        )
    sns_ax.legend(loc='center left', bbox_to_anchor=(1.5, 0.5), ncol=1)
    sns_ax.set_xlabel("Gene ratio")
    sns_ax.set_ylabel("")
    sns_ax.set_axisbelow(True)
    sns_ax.yaxis.grid(ls="--", which="major")
    # Add a colorbar
    sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
    sm.set_array([])
    sns_ax.figure.colorbar(sm, ax=axes, shrink=0.25, aspect=5, ticks=[0, 0.01, 0.05])
    # axes.margins(x=0.1, y=0.1)
    return fig, axes


####################
# proteomics plots #
####################
def bmires_bin_plot_for_proteins(protein_gene_pheno_df, gene1, gene2):
    g = sns.catplot(
        data=protein_gene_pheno_df,
        x="bmi_res_bins",
        row=f"level_{gene1}",
        col=f"level_{gene2}",
        row_order=["above", "median", "below"],
        col_order=["below", "median", "above"],
        kind="count",
        stat="proportion",
        sharex=True,
        sharey=True,
        )
    return g.figure


def create_coefs_plot(df, figsize=(6,3)):
    fig = fp.forestplot(
        df,  # the dataframe with results data
        estimate="interaction",  # col containing estimated effect size 
        ll="ci_low", hl="ci_high",
        pval="p_val",
        decimal_precision=3,
        varlabel="uniq_items",  # column containing variable label
        color_alt_rows=True,
        table=True,
        annote=["est_ci"],
        annoteheaders=["Est. (95% Conf. Int.)"],
        # group ordering
        sort=True,
        ylabel="Est. (95% Conf. Int.)",  # y-label title
        xlabel="Interaction coefficients",  # x-label title
        xticks=[-0.2, 0, 0.2, 0.4, 0.6],
        figsize=figsize,
        **{"marker": "D",  # set maker symbol as diamond
            "markersize": 75,  # adjust marker size
            "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
            "xlinecolor": "k",  # gray color for x-reference line
            "xline": 0,
            "xlinewidth": 1,
            "xtick_size": 12,  # adjust x-ticker fontsize
            "lw": 3,
            }       
        )
    return fig

def plot_box_protein_levels(
        boxdf, genes, protein_levels, sample_id, bmi, bmi_prs, bmi_residuals,
        xvar="genes", yvar="npx"
        ):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
        order=genes,
        gap=0.5,
        color= "whitesmoke",
        dodge=False, width=0.75, linewidth=2.25, fliersize=0, capprops={'color':'none'}, 
        boxprops={ 'edgecolor':'k'},  # 'facecolor':'none',
        whiskerprops={'color':'k'}, medianprops={'color':'k'}) # 

    # Adjust Axis
    ax.set_xlabel("")
    ax.set_ylabel("NPX")
    for i,pl in enumerate(protein_levels):
        ax.plot(i, pl, 'go', markersize=10)

    ax.set_ylim(-4, 4)
    # Remove Spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    print(bmi, bmi_prs, bmi_residuals)
    bmi, bmi_prs, bmi_residuals = list(map(lambda x: round(x, 3), [bmi, bmi_prs, bmi_residuals]))
    ax.set_title(f"{sample_id}\nBMI:{bmi}|PRS:{bmi_prs}|Res:{bmi_residuals}");
    return fig, ax

##############
# upset plot #
##############

def get_upset_plot(parsed_upset_df, figsize=(14, 6)):
    fig,ax = plt.subplots(1,1, figsize=figsize)
    upsetplot.plot(parsed_upset_df.counts, show_counts=True, fig=fig, element_size=None)
    ax.axis("off")
    return fig