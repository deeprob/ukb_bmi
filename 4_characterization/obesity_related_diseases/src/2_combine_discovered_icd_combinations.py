import os
import pandas as pd

icd_codes_of_interest = ['I10', 'E780', 'R074', 'I251', 'I259', 'E039', 'E11', 'Block M15-M19', 'K80', 'K81', 'K82', 'F32', 'F33', 'G30']
combos = [2,3,4]

risk_dfs = []
protective_dfs = []

# risk towards the disease
for icd_code in icd_codes_of_interest:
    for combo in combos:
        icd_code = icd_code.replace(" ", "")
        parsed_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/parsed_tables/{icd_code}_combo_{combo}.csv"
        if os.path.exists(parsed_out_file):
            df = pd.read_csv(parsed_out_file)
            risk_dfs.append(df)

risk_df = pd.concat(risk_dfs).reset_index(drop=True)
risk_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/merged_combos/risk_combos.csv"
risk_df.to_csv(risk_file, index=False)

# protective from the disease
for icd_code in icd_codes_of_interest:
    for combo in combos:
        icd_code = icd_code.replace(" ", "")
        parsed_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/parsed_tables/{icd_code}_protective_combo_{combo}.csv"
        if os.path.exists(parsed_out_file):
            df = pd.read_csv(parsed_out_file)
            protective_dfs.append(df)

protective_df = pd.concat(protective_dfs).reset_index(drop=True)
protective_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/merged_combos/protective_combos.csv"
protective_df.to_csv(protective_file, index=False)
