import os
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
def plot_box_single_oligo(
        boxdf, ttest_pval, 
        xvar="mutation", yvar="bmi", ylim=(10, 50), 
        xticklabel=["Single\nHit", "Combo\ncarriers"], ttest_hline=45, ttest_text=46):
    # Define Canvas
    fig,ax = plt.subplots(1, 1, figsize=(4, 6))

    # Box Plot
    sns_box = sns.boxplot(
        data=boxdf,
        x=xvar,
        y=yvar,
        order=["Single Hit", "Combo carriers"],
        hue_order=["Single Hit", "Combo carriers"],
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

