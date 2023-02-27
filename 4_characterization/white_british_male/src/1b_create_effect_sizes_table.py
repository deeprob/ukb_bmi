import pandas as pd

combo2_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/data/parsed_tables/combo_2.csv"
combo3_file = "/data5/deepro/ukbiobank/papers/bmi_project/3_run_rarecomb/white_british_male/data/parsed_tables/combo_3.csv"
akbari_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_male/data/effect_sizes/akbari.csv"
turcot_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_male/data/effect_sizes/giant.csv"

combo2_df = pd.read_csv(combo2_file)
combo3_df = pd.read_csv(combo3_file)
akbari_df = pd.read_csv(akbari_file)
turcot_df = pd.read_csv(turcot_file)

combo2_efs = combo2_df.Effect_Size.to_frame(name="Effect Size")
combo2_efs["Description"] = "Digenic combinations"
combo3_efs = combo3_df.Effect_Size.to_frame(name="Effect Size")
combo3_efs["Description"] = "Trigenic combinations"
akbari_efs = akbari_df.effect_size.to_frame(name="Effect Size")
akbari_efs["Description"] = "Akbari et. al."
turcot_efs = turcot_df.effect_size.to_frame(name="Effect Size")
turcot_efs["Description"] = "Turcot et. al."

efs_df = pd.concat([combo2_efs, combo3_efs, akbari_efs, turcot_efs], axis=0)

save_table = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/white_british_male/data/effect_sizes/all.csv"

efs_df.to_csv(save_table, index=False)
