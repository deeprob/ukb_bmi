import pandas as pd
import numpy as np


variants_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv"
combinations2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/protective/data/parsed_tables/combo_2.csv"
combinations3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/protective/data/parsed_tables/combo_3.csv"
phenotypes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/non_white_british/data/samples.csv"

variants_df = pd.read_csv(variants_file, sep="\t", low_memory=False, usecols=["Sample", "variant_id", "Gene", "SYMBOL", "Mut_type"], dtype=str)
combinations2_df = pd.read_csv(combinations2_file, low_memory=False)
combinations3_df = pd.read_csv(combinations3_file, low_memory=False)
phenotypes_df = pd.read_csv(phenotypes_file, low_memory=False, usecols=["eid", "bmi"], dtype={"eid": str, "bmi": np.float64})

# get all single genes
def get_unique_genes_from_combo_df(combo_df, ncombos):
    gene_set = set()
    for i in range(1, ncombos+1):
        gene_set.update(list(combo_df[f"Item_{i}_symbol"]))
    return gene_set


def get_samples_with_hits(variants_df, gene_set):
    return set(variants_df.loc[variants_df.SYMBOL.isin(gene_set)].Sample)

def check_combo_func(gene_list, combo_list):
    for combo in combo_list:
        if all(gene in gene_list for gene in combo):
            return True
    return False

def get_single_combo_phenos(variants_df, phenotypes_df, combinations_dfs, ncombos, ):
    # get unique genes in combos
    gene_set = set()
    all_combos_sorted = []
    for cdf,nc in zip(combinations_dfs, ncombos):
        gs = get_unique_genes_from_combo_df(cdf, nc)
        gene_set.update(gs)
        # sort combos by their gene names for consistency of pipeline
        acs = sorted(map(lambda x: tuple(sorted(x)), zip(*[cdf[f"Item_{c}_symbol"].to_list() for c in range(1, nc+1)])))
        all_combos_sorted.append(acs)
    all_combos_sorted = sum(all_combos_sorted, [])
    # store gene to their combo mappings
    gene_to_combo_map ={g:[] for g in gene_set}
    for g in gene_set:
        for c in all_combos_sorted:
            if g in c:
                gene_to_combo_map[g].append(c)
    # get samples with a single hit of the gene and those with combinations
    all_accepted_samples = set()
    all_rejected_samples = set()
    for g in gene_set:
        # select samples with variants in those genes; these samples may have single or multiple hits
        samples_list = variants_df.loc[variants_df.SYMBOL.isin([g])].Sample
        samples_df = variants_df.loc[variants_df.Sample.isin(samples_list)].groupby("Sample").aggregate(lambda x: list(x))
        # exclude samples which have the combinations for these genes
        rejected_series = samples_df.SYMBOL.apply(check_combo_func, args=(gene_to_combo_map[g], )) 
        rejected_samples = samples_df.loc[rejected_series].index.to_list()
        accepted_samples = samples_df.loc[~rejected_series].index.to_list()
        all_accepted_samples.update(accepted_samples)
        all_rejected_samples.update(rejected_samples)
    # get rid of samples with other combo hits from single gene hits
    all_accepted_samples =  all_accepted_samples - all_rejected_samples
    single_hit_pheno = phenotypes_df.loc[phenotypes_df.eid.isin(all_accepted_samples)]
    combo_hit_pheno = phenotypes_df.loc[phenotypes_df.eid.isin(all_rejected_samples)]
    non_combo_hit_pheno = phenotypes_df.loc[~phenotypes_df.eid.isin(all_rejected_samples)]
    return single_hit_pheno, combo_hit_pheno, non_combo_hit_pheno

single_hit_pheno, combo_hit_pheno, non_combo_hit_pheno = get_single_combo_phenos(variants_df, phenotypes_df, (combinations2_df, combinations3_df), (2,3))

combo_hit_pheno["vtype"] = "Combinations"
non_combo_hit_pheno["vtype"] = "No combinations"

valid_df = pd.concat((non_combo_hit_pheno, combo_hit_pheno), axis=0)

save_table = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/protective/data/validation/valid_mean.csv"
valid_df.to_csv(save_table, index=False)
