import os
import dxpy


oligo_project_id = "project-Gg0J78jJfzq8XV4v1G08XByZ"
prs_dir = "/fields/data/prs_obesity_related_info/"
store_dir = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/obesity_related_prs/data/prs_raw"


def download_folder(project_id, local_dir, dna_nexus_dir):
    os.makedirs(local_dir, exist_ok=True)
    dxpy.download_folder(project_id, local_dir, dna_nexus_dir, overwrite=True)
    return


if __name__ == "__main__":
    download_folder(oligo_project_id, store_dir, prs_dir)
