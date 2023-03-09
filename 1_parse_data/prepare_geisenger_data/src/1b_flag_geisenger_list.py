import pandas as pd


svariant_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_geisenger_data/data/variants_list.tsv"
gvariant_file = "/data5/deepro/ukbiobank/papers/bmi_project/0_data_download/geisenger_data/data/ghs-175k-rare-vep-mod-and-high-in-bmi-genes.txt"
gflagged_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_geisenger_data/data/ghs-175k-rare-vep-mod-and-high-in-bmi-genes_flagged.txt"

svariant_df = pd.read_csv(svariant_file, sep="\t", low_memory=False, dtype=str)
gvariant_df = pd.read_csv(gvariant_file, sep="\t", low_memory=False, header=None, names=["Chrom", "Pos", "Ref", "Alt", "gene", "transcript"], dtype=str)

svariant_df = svariant_df.set_index(["Chrom", "Pos", "Ref", "Alt"])
gvariant_df = gvariant_df.set_index(["Chrom", "Pos", "Ref", "Alt"])

gvariant_df["flag"] = False
gvariant_df.loc[gvariant_df.index.intersection(svariant_df.index), "flag"] = True

gvariant_df.to_csv(gflagged_file, header=False, index=True, sep="\t")
