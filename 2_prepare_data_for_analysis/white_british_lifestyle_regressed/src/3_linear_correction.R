#!/bin/R
library(glue)

final_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british_lifestyle_regressed/data/samples.csv"
corrected_file_samples_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/white_british_lifestyle_regressed/data/samples_with_residuals.csv"

		
df = read.csv(final_samples_file)
binary_input_columns = c('genetic_sex')
continuous_discrete_input_columns = c('age_when_attended_assessment_centre')
pcs = as.character(1:40)
pcs = glue('PC{pcs}')
lifestyle_columns = c('Input_1289_low', 'Input_1289_high', 'Input_1299_low', 'Input_1299_high', 'Input_1309_low', 'Input_1309_high', 'Input_1319_low', 'Input_1319_high', 'Input_1329_low', 'Input_1329_high', 'Input_1339_low', 'Input_1339_high', 'Input_1349_low', 'Input_1349_high', 'Input_1359_low', 'Input_1359_high', 'Input_1369_low', 'Input_1369_high', 'Input_1379_low', 'Input_1379_high', 'Input_1389_low', 'Input_1389_high', 'Input_1408_low', 'Input_1408_high', 'Input_1418_1', 'Input_1418_2', 'Input_1418_3', 'Input_1418_4', 'Input_1418_5', 'Input_1418_6', 'Input_1428_0', 'Input_1428_1', 'Input_1428_2', 'Input_1428_3', 'Input_1438_low', 'Input_1438_high', 'Input_1448_1', 'Input_1448_2', 'Input_1448_3', 'Input_1448_4', 'Input_1458_low', 'Input_1458_high', 'Input_1478_low', 'Input_1478_high', 'Input_1488_low', 'Input_1488_high', 'Input_1498_low', 'Input_1498_high', 'Input_1518_low', 'Input_1518_high', 'Input_1528_low', 'Input_1528_high', 'Input_1538_0', 'Input_1538_1', 'Input_1538_2', 'Input_1548_low', 'Input_1548_high', 'Input_6144_1', 'Input_6144_2', 'Input_6144_3', 'Input_6144_4', 'Input_6144_5', 'Input_1110_low', 'Input_1110_high', 'Input_2237_low', 'Input_2237_high', 'Input_2492_low', 'Input_2492_high', 'Input_6154_1', 'Input_6154_2', 'Input_6154_3', 'Input_6154_4', 'Input_6154_5', 'Input_6154_6', 'Input_6155_1', 'Input_6155_2', 'Input_6155_3', 'Input_6155_4', 'Input_6155_5', 'Input_6155_6', 'Input_6155_7', 'Input_6179_1', 'Input_6179_2', 'Input_6179_3', 'Input_6179_4', 'Input_6179_5', 'Input_6179_6', 'Input_137_low', 'Input_137_high', 'Input_1920_low', 'Input_1920_high', 'Input_1930_low', 'Input_1930_high', 'Input_1940_low', 'Input_1940_high', 'Input_1950_low', 'Input_1950_high', 'Input_1960_low', 'Input_1960_high', 'Input_1970_low', 'Input_1970_high', 'Input_1980_low', 'Input_1980_high', 'Input_1990_low', 'Input_1990_high', 'Input_2000_low', 'Input_2000_high', 'Input_2010_low', 'Input_2010_high', 'Input_2020_low', 'Input_2020_high', 'Input_2030_low', 'Input_2030_high', 'Input_2040_low', 'Input_2040_high', 'Input_2050_low', 'Input_2050_high', 'Input_2060_low', 'Input_2060_high', 'Input_2070_low', 'Input_2070_high', 'Input_2080_low', 'Input_2080_high', 'Input_2090_low', 'Input_2090_high', 'Input_2100_low', 'Input_2100_high', 'Input_6145_1', 'Input_6145_2', 'Input_6145_3', 'Input_6145_4', 'Input_6145_5', 'Input_6145_6', 'Input_864_low', 'Input_864_high', 'Input_874_low', 'Input_874_high', 'Input_884_low', 'Input_884_high', 'Input_904_low', 'Input_904_high', 'Input_924_low', 'Input_924_high', 'Input_943_low', 'Input_943_high', 'Input_1070_low', 'Input_1070_high', 'Input_1080_low', 'Input_1080_high', 'Input_1090_low', 'Input_1090_high', 'Input_1100_low', 'Input_1100_high', 'Input_1100_5', 'Input_6162_1', 'Input_6162_2', 'Input_6162_3', 'Input_6162_4', 'Input_6164_1', 'Input_6164_2', 'Input_6164_3', 'Input_6164_4', 'Input_6164_5', 'Input_1160_low', 'Input_1160_high', 'Input_1170_low', 'Input_1170_high', 'Input_1180_low', 'Input_1180_high', 'Input_1190_low', 'Input_1190_high', 'Input_1200_low', 'Input_1200_high', 'Input_1210_low', 'Input_1210_high', 'Input_1220_low', 'Input_1220_high', 'Input_1031_low', 'Input_1031_high', 'Input_2110_low', 'Input_2110_high', 'Input_6160_1', 'Input_6160_2', 'Input_6160_3', 'Input_6160_4', 'Input_6160_5')
continuous_discrete_input_columns = c(continuous_discrete_input_columns, pcs)
binary_input_columns = c(binary_input_columns, lifestyle_columns)


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
