import os
import dxpy


bmi_project_id = "project-GQpgZf8JX0KKbFBGK9yff4Zg"
bmi_dir = "/exome_annot/annot_run/notebooks"
store_dir = "../data/annot_tables_vep109"

def download_folder(project_id, local_dir, dna_nexus_dir):
    os.makedirs(local_dir, exist_ok=True)
    dxpy.download_folder(project_id, local_dir, dna_nexus_dir, overwrite=True)
    return


if __name__ == "__main__":
    for chrm in [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]:
        print(chrm)
        dnanexus_burden_dir = os.path.join(bmi_dir, chrm, "annot_tables_vep109")
        local_burden_dir = os.path.join(store_dir, chrm)
        os.makedirs(local_burden_dir, exist_ok=True)
        download_folder(bmi_project_id, local_burden_dir, dnanexus_burden_dir)
