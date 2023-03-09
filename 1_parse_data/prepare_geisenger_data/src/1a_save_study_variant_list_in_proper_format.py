import pandas as pd
import numpy as np


our_variant_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/annotate_vcf/data/variants_by_gene/lof_missense_pred_freq_0.01.tsv"
save_unique_variants_file = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_geisenger_data/data/variants_list.tsv"


variant_df = pd.read_csv(our_variant_file, sep="\t", low_memory=False)
# pivot table and include only unique variants
variant_pivot = variant_df.pivot_table(index=["Chrom", "Pos", "Ref", "Alt"], values=["Gene", "Sample", "SYMBOL"], aggfunc=lambda x: ",".join(list(map(str, x))))
variant_pivot = variant_pivot.reset_index()
variant_pivot["Chrom"] = variant_pivot.Chrom.str.replace("chr", "")
# save the variant list
variant_pivot.to_csv(save_unique_variants_file, index=False, sep="\t")
