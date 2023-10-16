import pandas as pd


if __name__=="__main__":
    turcot_files = [
        "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/turcot_2018/data/Giant_Turcot_Table1.csv",
        "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/turcot_2018/data/Giant_Turcot_Table2.csv"
    ]
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/turcot_2018/data/turcot_genes.list"
    df = pd.concat([pd.read_csv(tf) for tf in turcot_files])
    df.Gene.drop_duplicates().to_csv(save_file, index=False, header=False)
