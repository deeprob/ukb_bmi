#!/bin/python3


args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=3) {
  stop("expression table, current genelist and output filename must be given", call.=FALSE)
}

table_file = args[1]
gene_file = args[2]
out_file = args[3]


# load in tissue expression data and define gene universe
pref_df = read.csv(table_file, check.names=F)
rownames(pref_df) = pref_df$Gene_id
cell_types = colnames(pref_df)
cell_types = cell_types[2:length(cell_types)]
gene_universe = unique(rownames(pref_df))

print(length(gene_universe))

# load the gene list
genes = read.table(gene_file, header=FALSE)

all_stats = data.frame()

pref_df = pref_df[gene_universe,]
# sig_genes_df = sig_genes_df[sig_genes_df$Gene_id %in% gene_universe,]


# for each variant class do a fishers exact test for each cell type
columns=c('cell_type', 'oddsratio', 'pvalue', 'conf_int_lower', 'conf_int_upper', 'sig_and_pref', 'sig_and_not_pref', 'not_sig_and_pref', 'not_sig_and_pref')
stats = data.frame(matrix(ncol=length(columns), nrow=0))
colnames(stats) = columns
for (cell_type in cell_types) {
	genes_with_second_hits = genes$V1
	genes_for_cell_type = rownames(pref_df[pref_df[,cell_type] == 1,])

	second_hit_and_pref = unique(intersect(genes_with_second_hits, genes_for_cell_type))
	second_hit_and_not_pref = unique(genes_with_second_hits[!(genes_with_second_hits %in% genes_for_cell_type)])
	not_second_hit_and_pref = unique(genes_for_cell_type[!(genes_for_cell_type %in% genes_with_second_hits)])
	not_second_hit_and_not_pref_tmp = unique(gene_universe[!(gene_universe %in% genes_with_second_hits)])
	not_second_hit_and_not_pref = unique(not_second_hit_and_not_pref_tmp[!(not_second_hit_and_not_pref_tmp %in% genes_for_cell_type)])
	
	contingency_table = matrix(nrow=2, ncol=2)
	contingency_table[1,1] = length(second_hit_and_pref)
	contingency_table[1,2] = length(second_hit_and_not_pref)
	contingency_table[2,1] = length(not_second_hit_and_pref)
	contingency_table[2,2] = length(not_second_hit_and_not_pref)
	
	result = fisher.test(contingency_table)
	oddsratio = result$estimate
	pvalue = result$p.value
	conf_int_lower = result$conf.int[1]
	conf_int_upper = result$conf.int[2]
	
	stats[nrow(stats) + 1,] = c(cell_type, oddsratio, pvalue, conf_int_lower, conf_int_upper, length(second_hit_and_pref), length(second_hit_and_not_pref), length(not_second_hit_and_pref), length(not_second_hit_and_not_pref))

	}
	

# FDR correction for each variant class (BH correction)
stats[,'FDR'] = p.adjust(stats$pvalue, method='BH')

stats = stats[order(stats$FDR),]

all_stats = rbind(all_stats, stats)
		

all_stats <- all_stats[order(all_stats$FDR, all_stats$pvalue),]

write.csv(all_stats, out_file, row.names=F)





