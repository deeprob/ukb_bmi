#!/bin/R

library(glue)
library(RareComb)



args = commandArgs(trailingOnly=T)
wes_file = args[1]
cases_file = args[2]
controls_file = args[3]
combos = as.numeric(args[4])
output_file = args[5]
primary_entities = args[6]
secondary_entities = args[7]

gencode_file = "/data5/bx_reference/hg38/annotations/gene_annotations/GENCODE39/gencode.v39.parsed.genes.tsv"
gencode = read.table(gencode_file, sep='\t', header=T)
gencode[,'gene_id_stripped'] = unlist(lapply(gencode[,'gene_id'], function(x) strsplit(x, '.', fixed=T)[[1]][1]))
gencode = gencode[!duplicated(gencode$gene_id_stripped),]
gencode = gencode[,c('gene_id_stripped', 'Chrom', 'Start', 'End')]
colnames(gencode) = c('Gene', 'Chrom', 'Start', 'End')


df = read.table(wes_file, sep='\t', header=T)
# primary and secondary entities
pent = readLines(primary_entities)
sent = readLines(secondary_entities)


# get list of cases and controls
cases = scan(cases_file)
controls = scan(controls_file)
all_samples = c(cases, controls)



# only keep samples that are cases or controls
df = df[df$Sample_Name %in% all_samples,]

# add column for output
df[, 'Output_phenotype'] = as.integer(df$Sample_Name %in% cases)


combo_length=combos
max_freq_threshold=0.25
pval_filter_threshold=0.05
adj_pval_type='bonferroni'
min_power_threshold=0.7
sample_names_ind='Y'
quiet=FALSE

if (combos >2){
	min_indv_threshold = 7
} else {
	min_indv_threshold = 5
}

# print(min_indv_threshold)

result = compare_enrichment_modifiers(
    df, 
    gene_coordinates_df=gencode, 
    primary_input_entities=pent, 
    secondary_input_entities=sent, 
    combo_length=combo_length, 
    min_indv_threshold=min_indv_threshold, 
    max_freq_threshold=max_freq_threshold, 
    input_format='Input_', 
    output_format='Output_', 
    pval_filter_threshold=pval_filter_threshold, 
    adj_pval_type=adj_pval_type, 
    min_power_threshold=min_power_threshold, 
    sample_names_ind=sample_names_ind, 
    ld_block_size=0, 
    quiet=quiet
    )

write.csv(result, output_file, row.names=F)

rm(df)
rm(pent)
rm(sent)
rm(result)
