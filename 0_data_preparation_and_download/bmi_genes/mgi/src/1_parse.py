import pandas as pd

def save_genes(save_file, genes):
    with open(save_file, "w") as f:
        for g in genes:
            f.write(f"{g}\n")
    return


if __name__=="__main__":
    human_mouse_map_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/mgi/data/HMD_HumanPhenotype.rpt"
    mgi_pheno_desc_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/mgi/data/VOC_MammalianPhenotype.rpt"
    gene_pheno_map_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/mgi/data/MGI_PhenoGenoMP.rpt"
    save_file = "/data6/deepro/ukb_bmi/0_data_preparation_and_download/bmi_genes/mgi/data/mgi_genes.txt"

    human_to_mouse_gene_df = pd.read_csv(
        human_mouse_map_file, sep="\t", header=None, 
        names=["gene", "entrez", "mouse_gene", "mgi_id", "mgi_phenotypes"], 
        usecols=[0,1,2,3,4]
        )
    
    mouse_to_human_gene = dict(zip(human_to_mouse_gene_df.mgi_id, human_to_mouse_gene_df.gene))

    pheno_desc_df = pd.read_csv(mgi_pheno_desc_file, sep="\t", header=None, names=["mgi_phenotype", "description", "long_description"]).dropna()
    bmi_terms = ["body mass index", "obesity"]
    mgi_codes = pheno_desc_df.loc[pheno_desc_df.description.str.lower().str.contains("|".join(bmi_terms)), "mgi_phenotype"].values

    mouse_gene_to_pheno = pd.read_csv(gene_pheno_map_file, header=None, sep="\t", usecols=[3,5], names=["mgi_pheno", "mgi_gene"])
    mouse_genes = set("|".join(mouse_gene_to_pheno.loc[mouse_gene_to_pheno.mgi_pheno.isin(mgi_codes), "mgi_gene"]).split("|"))
    human_genes = set([mouse_to_human_gene[mg] for mg in mouse_genes if mg in mouse_to_human_gene])

    save_genes(save_file, human_genes)
