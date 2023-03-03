smap_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/slurm/files/0_smap.txt"
wes_file = "/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/tables/wes.tsv"

icd_codes_of_interest = ['I10', 'E780', 'R074', 'I251', 'I259', 'E039', 'E11', 'Block M15-M19', 'K80', 'K81', 'K82', 'F32', 'F33', 'G30']

combos = [2,3,4]
with open(smap_file, "w") as sfile:
    for icd_code in icd_codes_of_interest:
        for combo in combos:
            icd_code = icd_code.replace(" ", "")
            cases_file = f"/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/cases_controls/cases_{icd_code}.txt"
            controls_file = f"/data5/deepro/ukbiobank/papers/bmi_project/2_prepare_data_for_analysis/obesity_related_diseases/data/cases_controls/controls_{icd_code}.txt"
            out_file = f"/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/obesity_related_diseases/data/statistics/{icd_code.strip()}_combo_{combo}.csv"
            sfile.write(" ".join([wes_file, cases_file, controls_file, str(combo), out_file]))
            sfile.write("\n")

