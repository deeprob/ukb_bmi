import os
import matplotlib
import matplotlib.pyplot as plt
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
