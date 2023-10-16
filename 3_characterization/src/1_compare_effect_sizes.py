import pandas as pd
from statsmodels.stats.proportion import proportion_effectsize


def read_genes(gene_file):
    with open(gene_file, "r") as f:
        genes = [g.strip() for g in f.readlines()]
    return set(genes)


def get_overlapping_samples(sample_set, overlap_set):
    sample_set = set(sample_set.split(","))
    return len(sample_set.intersection(overlap_set))/len(overlap_set)


def get_effect_size(gene_list, case_control_df, burden_df, combo_samples, save_file):
    burden_df = burden_df.loc[burden_df.gene.isin(gene_list)]
    case_samples = set(case_control_df.loc[case_control_df.Output_BMI==1].Sample_Name.astype(str).values)
    control_samples = set(case_control_df.loc[case_control_df.Output_BMI==0].Sample_Name.astype(str).values)
    # from case and control, eliminate samples who have the combinations
    case_samples = case_samples.difference(combo_samples)
    control_samples = control_samples.difference(combo_samples)
    burden_df["prop_case_samples"] = burden_df.samples.apply(get_overlapping_samples, args=(case_samples, ))
    burden_df["prop_control_samples"] = burden_df.samples.apply(get_overlapping_samples, args=(control_samples, ))
    burden_df["Effect_Size"] = burden_df.apply(lambda x: proportion_effectsize(x.prop_case_samples, x.prop_control_samples), axis=1)
    burden_df.loc[:, ["gene", "Effect_Size"]].to_csv(save_file, index=False)
    return


if __name__=="__main__":
    gene_list_files = [
        "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/akbari_2021/data/akbari_genes.list",
        "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/turcot_2018/data/turcot_genes.list",
        "/data6/deepro/ukb_bmi/3_characterization/data/enrichment/british/genes.list"
    ]
    save_files = [
        "/data6/deepro/ukb_bmi/3_characterization/data/effect_sizes/akbari_genes.csv",
        "/data6/deepro/ukb_bmi/3_characterization/data/effect_sizes/turcot_genes.csv",
        "/data6/deepro/ukb_bmi/3_characterization/data/effect_sizes/study_genes.csv"
    ]

    combo_files = [
        "/data6/deepro/ukb_bmi/2_rarecomb/data/british/combo2.csv",
        "/data6/deepro/ukb_bmi/2_rarecomb/data/british/combo3.csv"
    ]
    case_control_file = "/data6/deepro/ukb_bmi/1_data_processing/data/british/case_controls.csv"
    burden_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/genotype/data/processed_burden/all_gene_burden.csv.gz"

    case_control_df = pd.read_csv(case_control_file)
    burden_df = pd.read_csv(burden_file)
    combo_dfs = [pd.read_csv(cf).loc[:, ["Case_Samples", "Control_Samples"]] for cf in combo_files]
    combo_df = pd.concat(combo_dfs)
    combo_samples = set("|".join(combo_df.values.flatten().astype(str)).split("|"))

    for glf,sf in zip(gene_list_files, save_files):
        gene_list = read_genes(glf)
        get_effect_size(gene_list, case_control_df, burden_df, combo_samples, sf)
