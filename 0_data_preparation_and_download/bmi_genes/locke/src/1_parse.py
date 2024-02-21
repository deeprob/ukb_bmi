import pandas as pd

def save_genes(save_file, genes):
    with open(save_file, "w") as f:
        for g in genes:
            f.write(f"{g}\n")
    return


if __name__=="__main__":
	locke_files = [
		"/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/locke/data/Table_1.xlsx",
		"/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/locke/data/Table_2.xlsx",
		"/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/locke/data/Extended_Data_Table_2.xlsx"
	]
    
	save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/locke/data/locke_genes.txt"
	
	locke_df = pd.concat([pd.read_excel(lf)  for lf in locke_files], axis=0)
	genes = [gi.split("(")[0] for g in locke_df["Notable gene(s)"].str.replace(" ", "").str.split(";") for gi in g if gi]

	save_genes(save_file, genes)
