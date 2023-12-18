# Characterization

## Analysis scripts
1. Combo information
    Given rarecomb output, gene burden file and cohort file, for each combination, adds the uniq items, all samples with combo in a cohort, original samples discovered in cases, original samples discovered in controls
2. Odds Ratio
    Given monogenic or oligogenic factors, gene burden file, and case control file, calculates odds of getting the phenotype and its significance
3. Variance explained
    Given PGS scores, monogenic and oligogenic factors along with covariates and gene burden file, calculates explained variance of different factors together.
4. Oligo pattern
    Given gene burden file, phenotype file, combo info file and an optional lifestyle factor file, calculates the effect of single hit versus combo hit
5. Additive model expected values
    Given genotype file, phenotype file and combo info file and an optional lifestyle factor file, calculated expected value from an additive model
6. Interaction table
    Given genotype file, phenotype file and combo info file, creates digenic interaction tables
7. Previous known association
    Given known gene lists, creates a table overlapping genes in each study 
8. Enrichment analysis
    Given combo info file, performs enrichment analysis of the genes in combo
9. PPI
    Given combo info file, checks PPI scores from string database
10. Variants per gene combo with sample info
    Given combo info file and variant annotation tables, gets the variant info of combos in cases, controls and others
11. 

## Enrichment Analysis

Resources
1. https://autobencoder.com/2021-10-03-gene-conversion/
2. https://pypi.org/project/biomart/
3. https://github.com/openvax/pyensembl
4. https://github.com/lukebfunk/pyclusterprofiler
5. https://github.com/tanghaibao/goatools

