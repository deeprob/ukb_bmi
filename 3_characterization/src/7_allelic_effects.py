import os
import argparse
import pandas as pd
from scipy.stats import ttest_ind
from functools import reduce

import utils.parsing as utpa
import utils.plotting as utpl


def get_case_cont_samples_from_df(df):
    combo_case_samples = set("|".join(df.Case_Samples.values).split("|"))
    combo_control_samples = set("|".join(df.Control_Samples.dropna().values).split("|"))
    return combo_case_samples, combo_control_samples

def get_vtype_df(block_df, combo_genes, combo_samples, combo_case_samples, combo_control_samples):
    # select lof and deleterious missense variants
    block_df = block_df.loc[(block_df.lof==True)|(block_df.splice_lof==True)|((block_df.missense==True)&(block_df.del_score>4))]
    block_df["samples"] = block_df.samples.str.split(",")
    block_df = block_df.explode("samples")
    block_df = block_df.loc[(block_df.samples.isin(combo_samples.union(combo_case_samples).union(combo_control_samples))) & (block_df.gene.isin(combo_genes))]
    block_df["variant"] = block_df.gene + "_" + block_df.locus + "_" + block_df.alleles
    return block_df

def create_variant_df(annot_table_dir, combo_genes, combo_samples, combo_case_samples, combo_control_samples):
    vcfs_per_chrm = {
        "chr1": 97, "chr2": 71, "chr3": 56, "chr4": 39, "chr5": 43, "chr6": 48, 
        "chr7": 47, "chr8": 35, "chr9": 42, "chr10": 40, "chr11": 57, "chr12": 52, 
        "chr13": 18, "chr14": 30, "chr15": 34, "chr16": 47, "chr17": 56, "chr18": 16, 
        "chr19": 65, "chr20": 25, "chr21": 11, "chr22": 23, "chrX": 24, "chrY": 1
    }
    variant_df = pd.DataFrame()
    for chr_num in [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]:
        chr_file_num = vcfs_per_chrm[chr_num]
        for filei in range(chr_file_num):
            block_file = os.path.join(annot_table_dir, f"{chr_num}", f"block_{filei}.tsv.gz")
            block_df = pd.read_csv(block_file, sep="\t", index_col=0)
            block_df = get_vtype_df(block_df, combo_genes, combo_samples, combo_case_samples, combo_control_samples)
            variant_df = pd.concat((variant_df, block_df))
    return variant_df 

def get_combo_case_loci(variant_case_df):
    combo_case_loci = set()
    for gn, gdf in variant_case_df.groupby("samples"):
        gdf = gdf.sort_values(["variant"])
        combo_case_locus = "|".join((gdf.variant).values)
        combo_case_loci.update([combo_case_locus])
    return combo_case_loci

def get_same_variant_samples(variant_subset_df, combo_case_variant_loci):
    same_variant_samples = set()
    diff_variant_samples = set()
    for gn, gdf in variant_subset_df.groupby("samples"):
        gdf = gdf.sort_values(["variant"])
        combo_cont_locus = "|".join((gdf.variant).values)
        if combo_cont_locus in combo_case_variant_loci:
            same_variant_samples.update([gn])
        else:
            diff_variant_samples.update([gn])
    return same_variant_samples, diff_variant_samples

def get_samples_with_same_variants(ser, variant_df, add_control_sample):
    genes = ser.uniq_items.replace("Input_", "").split("|")
    case_samples = set(ser.Case_Samples.split("|"))
    all_samples = set(ser.combo_samples.split("|"))
    variant_case_df = variant_df.loc[(variant_df.gene.isin(genes))&(variant_df.samples.isin(case_samples))]
    combo_case_variant_loci = get_combo_case_loci(variant_case_df)

    if add_control_sample:
        control_samples = set(ser.Control_Samples.split("|"))
        other_samples = all_samples.difference(case_samples.union(control_samples))
        variant_control_df = variant_df.loc[(variant_df.gene.isin(genes))&(variant_df.samples.isin(control_samples))]
        control_samples_same_variants, control_samples_other_variants = get_same_variant_samples(variant_control_df, combo_case_variant_loci)
        ser = pd.Series({
            "control_samples_same_variants": "|".join(control_samples_same_variants), 
            "control_samples_other_variants": "|".join(control_samples_other_variants)
            })
    else:
        other_samples = all_samples.difference(case_samples)

    variant_other_df = variant_df.loc[(variant_df.gene.isin(genes))&(variant_df.samples.isin(other_samples))]
    other_samples_same_variants, other_samples_other_variants = get_same_variant_samples(variant_other_df, combo_case_variant_loci)
    ser["other_samples_same_variants"] = "|".join(other_samples_same_variants)
    ser["other_samples_other_variants"] = "|".join(other_samples_other_variants)
    return ser 

def compare_allelic_variants(variant_df, combo_info_df, cohort_df, all_case_samples, control_sample_flag):
    combo_info_df = combo_info_df.merge(combo_info_df.apply(get_samples_with_same_variants, args=(variant_df, control_sample_flag), axis=1), left_index=True, right_index=True)
    sv_samples = set("|".join(combo_info_df.loc[combo_info_df.other_samples_same_variants!="", "other_samples_same_variants"].values).split("|"))
    dv_samples = set("|".join(combo_info_df.loc[combo_info_df.other_samples_other_variants!="", "other_samples_other_variants"].values).split("|"))
    pheno_sv_df = cohort_df.loc[cohort_df.sample_names.isin(sv_samples)]
    pheno_dv_df = cohort_df.loc[cohort_df.sample_names.isin(dv_samples)]
    pheno_cont_df = cohort_df.loc[~cohort_df.sample_names.isin(sv_samples.union(dv_samples).union(all_case_samples))]
    pheno_sv_df["mutation"] = "Same variants"
    pheno_dv_df["mutation"] = "Different variants"
    pheno_cont_df["mutation"] = "Controls"
    ttest_pval = ttest_ind(pheno_dv_df.bmi, pheno_sv_df.bmi, alternative="less").pvalue
    ttest_pval_cont = ttest_ind(pheno_cont_df.bmi, pheno_sv_df.bmi, alternative="less").pvalue
    allelic_df = pd.concat([pheno_sv_df, pheno_dv_df, pheno_cont_df])
    return allelic_df, ttest_pval, ttest_pval_cont

def collapse_variant_info_by_type(vtypes):
    vtypes = set(vtypes)
    if "lof" in vtypes:
        vtype = "lof"
    else:
        vtype = "missense"
    return vtype

def collapse_variant_info_by_gene(df):
    df["variant"] = df.gene + "_" + df.locus + "_" + df.alleles
    df["vtype"] = df.apply(lambda x: "lof" if (x.lof or x.splice_lof) else "missense", axis=1)
    return pd.Series({"variant": ",".join(df.variant.values), "vtype": collapse_variant_info_by_type(df.vtype.values)})

def collapse_variant_info_by_samples(ser, cohort_df):
    cvg_info = [collapse_variant_info_by_gene(cohort_df.loc[(cohort_df.samples==ser.combo_samples)&(cohort_df.gene==g)]) for g in ser.combo_genes]
    return reduce(lambda x,y: x + "|" + y, cvg_info)

def set_vtype(val):
    vtypes = set(val.split("|"))
    if len(vtypes)==1:
        if "lof" in vtypes:
            vtype = "lof only"
        else:
            vtype = "missense only"
    else:
        vtype = "both lof and missense"
    return vtype

def group_vtype_by_samples(vtypes):
    vtypes = set(vtypes)
    if "lof only" in vtypes:
        vtype = "lof only"
    elif "both lof and missense" in vtypes:
        vtype = "both lof and missense"
    else:
        vtype = "missense only"
    return vtype

def compare_variant_types(combo_df, variant_df, cohort_df):
    combo_df["combo_genes"] = combo_df.uniq_items.str.replace("Input_", "").str.split("|") 
    combo_df["combo_samples"] = combo_df.combo_samples.str.split("|")
    combo_df = combo_df.explode("combo_samples")
    combo_vinfo_df = pd.concat([combo_df, combo_df.apply(collapse_variant_info_by_samples, args=(variant_df, ), axis=1)], axis=1)
    combo_vinfo_df["vtype_parsed"] = combo_vinfo_df.vtype.apply(set_vtype)
    combo_vtypeinfo_df = combo_vinfo_df.groupby("combo_samples").agg({"vtype_parsed": group_vtype_by_samples}).reset_index()
    plot_df = combo_vtypeinfo_df.merge(cohort_df, left_on="combo_samples", right_on="sample_names")
    return plot_df


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Enrichment analysis")
    parser.add_argument("--annot_table_dir", type=str, help="Filepath of the dir where annotation tables are stored")
    parser.add_argument("--cohort_file", type=str, help="Filepath of the cohort pheno file")
    parser.add_argument("--case_control_file", type=str, help="Filepath of the case control pheno file")
    parser.add_argument("--rarecomb_files", type=str, help="Filepath of the combo files given by Rarecomb", nargs="+")
    parser.add_argument("--combo_files", type=str, help="Filepath of the combo info files", nargs="+")
    parser.add_argument("--save_dir", type=str, help="Filepath where allelic effects info will be stored")
    parser.add_argument("--add_control_sample", action="store_true", help="Flag to add control sample info")

    cli_args = parser.parse_args()

    combo_info_df = pd.concat([pd.read_csv(cif) for cif in cli_args.combo_files])
    combo_df = pd.concat([pd.read_csv(cf, usecols=["uniq_items", "Case_Samples", "Control_Samples"]) for cf in cli_args.rarecomb_files])
    combo_info_df = combo_info_df.merge(combo_df, on="uniq_items")
    combo_info_df = combo_info_df.fillna("")

    combo_genes, combo_samples = utpa.get_combo_info_from_files(cli_args.combo_files)
    combo_case_samples, combo_control_samples = get_case_cont_samples_from_df(combo_df)

    variant_df = create_variant_df(cli_args.annot_table_dir, combo_genes, combo_samples, combo_case_samples, combo_control_samples)
    cohort_df = pd.read_csv(cli_args.cohort_file, usecols=["sample_names", "bmi"], dtype={"sample_names": str, "bmi": float})

    case_cont_df = pd.read_csv(cli_args.case_control_file, dtype={"Sample_Name": str, "Output_BMI": float})
    all_case_samples = set(case_cont_df.loc[case_cont_df.Output_BMI==1, "Sample_Name"].values)
    boxdf, ttest_pval, ttest_pval_cont = compare_allelic_variants(variant_df, combo_info_df, cohort_df, all_case_samples, cli_args.add_control_sample)
    fig, ax = utpl.plot_box_allelic(boxdf, ttest_pval, ttest_pval_cont)
    save_file = os.path.join(cli_args.save_dir, "allelic_effect.pdf")
    utpl.save_pdf(save_file, fig)

    plot_df = compare_variant_types(combo_info_df, variant_df, cohort_df)
    fig, ax = utpl.plot_box_variant_types(plot_df)
    save_file = os.path.join(cli_args.save_dir, "variant_type.pdf")
    utpl.save_pdf(save_file, fig)
