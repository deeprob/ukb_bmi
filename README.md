# Non-additive interactions of rare variants and lifestyle factors contribute to obesity
Codebase for manuscript titled "Non-additive interactions of rare variants and lifestyle factors contribute to obesity" by ;

# Repository guidelines
This repo contains step-by-step guide to the analysis carried out for "Non-additive interactions of rare variants and lifestyle factors contribute to obesity". A description of each sub dir is as follows:

# Links to Figures and Supplemetary tables
## Figures
1. Figure 1A: "" # TODO: make a study design figure
2. Figure 1B: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/1_effect_sizes_comparison.ipynb"
3. Figure 1C: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/2_single_hit_vs_combos.ipynb"
4. Figure 1D: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/3_additivity_model.ipynb"
5. Figure 1E: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/4_validation_mean.ipynb"
6. Figure 2A: "//data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/6f_all_pheno_enrichment.ipynb"
7. Figure 2B: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/6b_go_enrichment.ipynb"
8. Figure 2C: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/6b_go_enrichment.ipynb"
9. Figure 2D: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/6b_go_enrichment.ipynb"
10. Supp. Figure 1A: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/notebooks/6c_kegg_enrichment.ipynb"


## Tables
1. Supp. Table 1: "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_2.csv"
2. Supp. Table 2: "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british/data/parsed_tables/combo_3.csv"
3. Supp. Table 3: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/previous_associations/associations.csv"
4. Supp. Table 4: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/effect_sizes/{akbari|giant}.csv"
5. Supp. Table 5: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/validation/dir_consistency.csv"
6. Supp. Table 6: "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_lifestyle_regressed/data/parsed_tables/combo_2.csv"
7. Supp. Table 7: "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_lifestyle_regressed/data/parsed_tables/combo_3.csv"
8. Supp. Table 8: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/gwas/gwas_enrichment.csv"
9. Supp. Table 9: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/go/go_enrichment.csv"
10. Supp. Table 10: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/kegg/kegg_enrichment.csv"
11. Supp. Table 11: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/hpo/hpo_enrichment.csv"
12. Supp. Table 12: "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british/data/enrichment/mgi/mgi_enrichment.csv"

# TODO list of figures
## Result 1
1. Study design - Figure 1A :: **Major update**
    - make from scratch
2. Mean effect size - Figure 1B
    - add line with p-value
    - add proper font and labels
3. Oligogenic inheritance - Figure 1C
    - add individual datapoints as a heatmap
    - add line with p-value
    - draft add t-test p-value
    - add proper font and labels
4. Non-additive figure - Figure 1D
    - add individual datapoints as a heatmap
    - add line with p-value
    - add proper font and labels
5. Non-white British validation
    - add datapoints
    - add significance line
    - check with bmi residuals
    - add proper font and labels
6. Lifestyle regressed out
    - rerun using new metatable

## Result 2
1. Gwas, hpo, mgi figure - Figure 2A
    - add proper font and labels
2. GO enrichment figures - Figure 2B; Supp Figure 1A; Supp Figure 1B; 
    - add proper font and labels
3. KEGG pathway enrichment - Figure 2C
    - add proper font and labels
4. Tissue specific enrichment :: **Major update**
    - conduct analysis - COMPLETED on 03/08
    - make heatmap figure

## Result 3
1. Overlapping combos upset plot - Figure 3A
    - add proper font and labels

2. Combinations in opposite sex
    - combine the two figures
    - add proper font and labels

3. Effect sizes, oligogenicity and non-additivity
    - add proper font and labels
    - add significance values 
    - show datapoints

4. Non-white british validation
    - add proper font and labels
    - add datapoints
    - add significance values

5. Enrichment map figure
    - add all nodes to some cluster
    - get rid of node overlap
    - only include molecular functions that are not the same as biological processes
    - add proper font and labels

6. Tissue specific enrichment
    - analysis for individual cohorts as well

## Result 4
1. Protective combos oligogenicity and non additivity - Supp. Figures 3A and 3B
    - add proper fonts and labels

2. Non white british validation 
    - try with bmi residuals
    - remove this analysis

3. Enrichment figure - Figure 4a :: **Major update**
    - create an enrichment map figure 

4. Risk and protection venn diagram - Figure 4b:
    - do this using matplotlib venn, more flexibility
    - add proper fonts and labels

5. Allelic heterogeneity figure - Figure 4c :: **Major update**
    - create the figure outline
    - check other papers
    - think of a table to represent this

6. Variant classes analysis :: **Major update**
    - can variant class explain the variably expressive genes
    - conduct analysis protective
    - conduct analysis low bmi associated
    - this can also be a figure

## Result 5
1.  Change organization of the written portion :: **Major update** - COMPLETED on 03/09
    - separate risk and protection
    - start with risk 
    - add examples of combinations that added more risk towards relevant disease
    - add protection combinations mining
    - examples of protective combinations
    - go to overlaps between all, then specifics
    - then go to pleiotropy

2. Overlap between all venn - Figure 5A
    - do this using matplotlib venn, more flexibility
    - add proper fonts and labels

3. Upset plot in ICD risk - Supp. Figure 5A
    - add proper fonts and labels

4. Upset plot in ICD protection - Figure 5B
    - add proper fonts and labels

5. Individual profile figure - Figure 5C
    - add proper fonts and labels

6. Obesity related in non-obese individuals
    - this is not required

## Result 6
1. Lifestyle factor enrichment :: **Major update**
    - new analysis

2. Oligogenic inheritance pattern
    - add significance line
    - add proper fonts and labels

3. Non additive figure
    - add datapoints
    - add significance line
    - add proper fonts and labels

4. Enrichment plot - Figure 6C :: **Major update**
    - create an enrichment map figure

5. Add example about the lifestyle factors
    - take examples which are not very general

## Send data to geisenger for validation
1. Add unique genes with their start and end coordinates from
    - obesity risk - COMPLETED on 03/08
    - obesity protection
    - low bmi association

2. Add combo file in format for 
    - obesity risk - COMPLETED on 03/08
    - obesity protection
    - low bmi association

3. Annotate the variants already provided with a flag sign - COMPLETED on 03/08

4. Send email

## Abstract

## Introduction

## Discussion

## Methods

## Repository code needs to be more streamlined
    - reuse code chunks
    - remove redundancy
    - make the interface simpler
    - start with run rarecomb dir
