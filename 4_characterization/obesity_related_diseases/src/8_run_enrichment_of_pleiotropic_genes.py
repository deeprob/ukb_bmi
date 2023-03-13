import pandas as pd
import subprocess
import os


def run_gsea_enrichment(gene_in, table_out, gsea_bash_path):
    cmd = [
        "bash", gsea_bash_path, 
        gene_in, table_out
        ]
    subprocess.run(cmd)
    return

def create_gsea_enrichment_table(gene_file, gsea_bash_path, save_dir):
    # run enrichment
    enrich_table = os.path.join(save_dir, "go_enrichment.csv")
    run_gsea_enrichment(gene_file, enrich_table, gsea_bash_path)
    return

if __name__ == "__main__":
    gene_file = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/pleiotropic/geneids.txt"
    gsea_bash_path = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/src/scripts/gsea_enrich.sh"
    save_dir = "/data5/deepro/ukbiobank/papers/bmi_project/4_characterization/obesity_related_diseases/data/pleiotropic/go/"
    
    create_gsea_enrichment_table(gene_file, gsea_bash_path, save_dir)
