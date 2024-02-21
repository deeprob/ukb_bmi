import pandas as pd


if __name__=="__main__":
    akbari_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/akbari_2021/data/akbari_genes.csv"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/akbari_2021/data/akbari_genes.list"
    df = pd.read_csv(akbari_file)
    df.Gene.drop_duplicates().to_csv(save_file, index=False, header=False)
