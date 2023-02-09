import os
import pandas as pd

df = pd.read_csv('/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/samples.csv')
slurm_file = os.path.join("/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/slurm/files", "13_smap.txt")


if __name__=="__main__":
    sf = open(slurm_file, "w")
    for row in df.itertuples():
        sample_name = str(row.Sample)
        in_vcf = f'/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/by_sample/{sample_name[-3:]}/{sample_name}.vcf.gz'
        out_vcf_pre = f'/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/vcfs/annotated_by_sample/{sample_name[-2:]}/{sample_name}'
        os.makedirs(out_vcf_pre, exist_ok=True)
        sf.write("\t".join([sample_name, in_vcf, out_vcf_pre]))
        sf.write("\n")
    sf.close()
