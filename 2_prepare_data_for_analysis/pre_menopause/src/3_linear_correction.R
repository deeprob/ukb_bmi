#!/bin/R
library(glue)

final_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/pre_menopause/data/samples.csv"
corrected_file_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/pre_menopause/data/samples_with_residuals.csv"


df = read.csv(final_samples_file)
binary_input_columns = c('genetic_sex')
continuous_discrete_input_columns = c('age_when_attended_assessment_centre')
pcs = as.character(1:40)
pcs = glue('PC{pcs}')
continuous_discrete_input_columns = c(continuous_discrete_input_columns, pcs)


# scale the continuous and discrete inputs
for (col in continuous_discrete_input_columns){
	df[, glue('{col}_scaled')] = scale(df[, col])
}

# scale bmi
df[, 'bmi_scaled'] = scale(df[, 'bmi'])
formula = glue('bmi_scaled ~ {glue_collapse(binary_input_columns, sep = " + ")} + {glue_collapse(continuous_discrete_input_columns, sep = "_scaled + ")}_scaled')
mod = lm(formula, data=df, x=T)
residuals = resid(mod)
residuals = as.numeric(residuals)
df$bmi_residuals = residuals

write.csv(df, corrected_file_samples_file, row.names=F)
