import os
import dxpy


bmi_project_id = "project-GQpgZf8JX0KKbFBGK9yff4Zg"
lifestyle_dir = "/phenotype_processing/lifestyle_info/"
store_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/lifestyle_factors/data/lifestyle_raw/"


def download_folder(project_id, local_dir, dna_nexus_dir):
    os.makedirs(local_dir, exist_ok=True)
    dxpy.download_folder(project_id, local_dir, dna_nexus_dir, overwrite=True)
    return


if __name__ == "__main__":
    download_folder(bmi_project_id, store_dir, lifestyle_dir)
