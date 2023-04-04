import numpy as np
import pandas as pd
import itertools as it
from collections import Counter
from scipy.stats import chi2_contingency


variants_file = "/data5/UK_Biobank/annotations/vep/2022_03_13/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv"
phenotypes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british/data/samples_with_residuals.csv"
# icd_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_icd_codes/data/ukb30075_icd10.csv"
combinations_files = [
    "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv",
    "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
    ]
lowbmi_combo_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/icd/combosamples_bmidecile1to3.txt"
all_sample_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/icd/allsamples_bmidecile1to3.txt"
combo_samples_per_decile_save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/combo_samples_per_decile/combo_samples_per_decile.csv"
decile_cut_save_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/combo_samples_per_decile/decile_rank_to_bmi.csv"


def get_merged_combo_file(combo_files):
    combo_dfs = [pd.read_csv(cf) for cf in combo_files]
    combo_df = pd.concat(combo_dfs)
    return combo_df


combo_df = get_merged_combo_file(combinations_files)
variants_df = pd.read_csv(variants_file, sep="\t", low_memory=False, usecols=["Sample", "variant_id", "Gene", "SYMBOL", "Mut_type"], dtype=str)
phenotypes_df = pd.read_csv(phenotypes_file, low_memory=False, usecols=["eid", "bmi_residuals"], dtype={"eid": str, "bmi_residuals": np.float64})
# icd_df = pd.read_csv(icd_file, low_memory=False, dtype={"eid": str})


variant_type_dict = {"eid":[], "vtype":[]}
for uci in combo_df.unique_combo_id:
    genes = uci.split("_")
    svdf = variants_df.loc[(variants_df.Gene.isin(genes))]
    for group in svdf.groupby("Sample"):
        # if the set of genes is present in the sample variant df, that means the sample contains the combo
        if set(genes).issubset(set(group[1].Gene)):
            variant_type_dict["eid"].append(group[0])
            variant_type_dict["vtype"].append("-".join(sorted(group[1].Mut_type.unique())))

# adding Quantile_rank column to the DataFrame
decile_cut = pd.qcut(phenotypes_df['bmi_residuals'], 10, labels = False, retbins=True)
decile_bins = ["-".join([str(int(decile_cut[1][i])), str(int(decile_cut[1][i+1]))]) for i in range(10)]
decile_map_dict = dict(zip(range(10), decile_bins))
phenotypes_df['Decile_rank'] = decile_cut[0]
phenotypes_df["Decile_label"] = phenotypes_df.Decile_rank.map(decile_map_dict)



vtype_bmi_df = pd.DataFrame(variant_type_dict).merge(phenotypes_df, left_on="eid", right_on="eid")
vtype_bmi_df.to_csv(combo_samples_per_decile_save_file, index=False)

with open(lowbmi_combo_samples_file, "w") as f:
    for sample in vtype_bmi_df.loc[vtype_bmi_df.Decile_rank<3, "eid"]:
        f.write(f"{sample}\n")

with open(all_sample_file, "w") as f:
    for sample in phenotypes_df.loc[phenotypes_df.Decile_rank<3, "eid"]:
        f.write(f"{sample}\n")

# profile_df = vtype_bmi_df.merge(icd_df, left_on="eid", right_on="eid")
# icd_columns = list(icd_df.columns)[1:]
# vtype_bmi_df["Diagnosis"] = ~profile_df.loc[:, icd_columns].isna().all(axis=1)
# all_icd_list = sum([icds.split(",") for icds in profile_df.loc[profile_df.Decile_rank<=2, icd_columns].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1).values.flatten() if icds], [])
# all_icd_list = sum([icds.split(",") for icds in profile_df.loc[profile_df.Decile_rank<=2, icd_columns].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1).values.flatten() if icds], [])
# icd_bmi_df = vtype_bmi_df.groupby("Decile_rank").agg({"Diagnosis": sum, "eid": len})
# icd_bmi_df["No_diagnosis"] = icd_bmi_df.eid - icd_bmi_df.Diagnosis
# chi2_contingency(icd_bmi_df.loc[:, ["No_diagnosis", "Diagnosis"]].T.values)

