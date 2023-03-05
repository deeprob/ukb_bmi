smap_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/slurm/files/1_smap.txt"

icd_codes_of_interest = ['I10', 'E780', 'R074', 'I251', 'I259', 'E039', 'E11', 'Block M15-M19', 'K80', 'K81', 'K82', 'F32', 'F33', 'G30']

combos = [2,3,4]
with open(smap_file, "w") as sfile:
    # risk towards the disease
    for icd_code in icd_codes_of_interest:
        for combo in combos:
            icd_code = icd_code.replace(" ", "")
            rarecomb_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/statistics/{icd_code}_combo_{combo}.csv"
            parsed_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/parsed_tables/{icd_code}_combo_{combo}.csv"
            sfile.write(" ".join([rarecomb_out_file, parsed_out_file, icd_code, str(combo)]))
            sfile.write("\n")
    # protective from the disease
    for icd_code in icd_codes_of_interest:
        for combo in combos:
            icd_code = icd_code.replace(" ", "")
            rarecomb_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/statistics/{icd_code}_protective_combo_{combo}.csv"
            parsed_out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/parsed_tables/{icd_code}_protective_combo_{combo}.csv"
            sfile.write(" ".join([rarecomb_out_file, parsed_out_file, icd_code, str(combo)]))
            sfile.write("\n")
