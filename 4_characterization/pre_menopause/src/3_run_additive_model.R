#!/bin/R




library(glue)

sig_combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/pre_menopause/data/parsed_tables/combo_2.csv"
# sig_combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/pre_menopause/data/parsed_tables/combo_3.csv"
variant_by_gene_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/pre_menopause/data/additive/input_table.csv"

sig_combos_df2 = read.csv(sig_combo2_file)
# sig_combos_df3 = read.csv(sig_combo3_file)

wdf = read.csv(variant_by_gene_file)



sig_combos_df2[, 'Observed_value'] = NA
sig_combos_df2[, 'Observed_sd'] = NA
sig_combos_df2[, 'Expected_value'] = NA
sig_combos_df2[, 'diff_observed_expected'] = NA
sig_combos_df2[, 'z_score'] = NA
sig_combos_df2[, 'pvalue'] = NA

for (i in 1:nrow(sig_combos_df2)) {
	genes = sig_combos_df2[i, c('Item_1', 'Item_2')]
	genes = as.character(genes)
	
	# to calculate expected additive values, we cannot use the
	# individuals with both hits when fitting the model
	# because it will skew the coefficients
	#
	# in other words, we want to calculate the coefficients based
	# on individuals who carry only one of the genes that are
	# part of the combination
	input_df = wdf[!((wdf[genes[1]] == 1) & (wdf[genes[2]] == 1)),]
	formula = glue('bmi ~ {glue_collapse(genes, sep = " + ")}')
	mod = lm(formula, data=input_df)
	
	data = data.frame(1,1)
	colnames(data) = genes
	estimated = predict(mod, newdata=data, se.fit=T)
	estimated_value = estimated$fit
	# print(estimated)
	
	observed = wdf[(wdf[genes[1]] == 1) & (wdf[genes[2]] == 1),]
	observed_value = mean(observed$bmi)
	observed_std = sd(observed$bmi)
	
	z_score = (observed_value - estimated_value) / observed_std
	# two tailed z-test
	pvalue = 2*pnorm(-abs(z_score))
	
	# add statistics to table
	sig_combos_df2[i, 'Observed_value'] = observed_value
	sig_combos_df2[i, 'Observed_sd'] = observed_std
	sig_combos_df2[i, 'Expected_value'] = estimated_value
	sig_combos_df2[i, 'diff_observed_expected'] = observed_value-estimated_value
	sig_combos_df2[i, 'z_score'] = z_score
	sig_combos_df2[i, 'pvalue'] = pvalue
}

combo2_outfile = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/pre_menopause/data/additive/combo2.csv"
write.csv(sig_combos_df2, combo2_outfile, row.names=F)


# #############################################
# # then same thing for combinations of three #
# #############################################



# sig_combos_df3[, 'Observed_value'] = NA
# sig_combos_df3[, 'Observed_sd'] = NA
# sig_combos_df3[, 'Expected_value'] = NA
# sig_combos_df3[, 'diff_observed_expected'] = NA
# sig_combos_df3[, 'z_score'] = NA
# sig_combos_df3[, 'pvalue'] = NA


# # do the regression
# for (i in 1:nrow(sig_combos_df3)) {
# 	genes = sig_combos_df3[i, c('Item_1', 'Item_2', 'Item_3')]
# 	genes = as.character(genes)
	
# 	input_df = wdf[!((wdf[genes[1]] == 1) & (wdf[genes[2]] == 1) & (wdf[genes[3]] == 1)),]
# 	formula = glue('bmi ~ {glue_collapse(genes, sep = " + ")}')
# 	mod = lm(formula, data=input_df)
	
# 	data = data.frame(1,1,1)
# 	colnames(data) = genes
# 	estimated = predict(mod, newdata=data, se.fit=T)
# 	estimated_value = estimated$fit
# 	# print(estimated)
	
# 	observed = wdf[(wdf[genes[1]] == 1) & (wdf[genes[2]] == 1) & (wdf[genes[3]] == 1),]
# 	observed_value = mean(observed$bmi)
# 	observed_std = sd(observed$bmi)
	
# 	z_score = (observed_value - estimated_value) / observed_std
# 	# two tailed z-test
# 	pvalue = 2*pnorm(-abs(z_score))
	
# 	# add statistics to table
# 	sig_combos_df3[i, 'Observed_value'] = observed_value
# 	sig_combos_df3[i, 'Observed_sd'] = observed_std
# 	sig_combos_df3[i, 'Expected_value'] = estimated_value
# 	sig_combos_df3[i, 'diff_observed_expected'] = observed_value-estimated_value
# 	sig_combos_df3[i, 'z_score'] = z_score
# 	sig_combos_df3[i, 'pvalue'] = pvalue

# }

# combo3_outfile = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/pre_menopause/data/additive/combo3.csv"
# write.csv(sig_combos_df3, combo3_outfile, row.names=F)
