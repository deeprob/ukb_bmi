import pandas as pd



filename = '/data5/bx_reference/hg38/annotations/gene_annotations/GENCODE39/gencode.v39.parsed.genes.tsv'
gencode = pd.read_csv(filename, sep='\t')
gencode['gene_id_stripped'] = gencode['gene_id'].apply(lambda s: s.split('.')[0])
# gencode = gencode.drop_duplicates('gene_id_stripped')
# gencode = gencode.set_index('gene_id_stripped', drop=False)

outfile = "/data5/deepro/ukbiobank/papers/bmi_project/1_parse_data/prepare_gencode_genes/data/gencode.v39.parsed.genes.csv"
gencode.to_csv(outfile, index=False)

