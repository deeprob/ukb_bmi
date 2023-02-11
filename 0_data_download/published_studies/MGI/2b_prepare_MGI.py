#!/bin/python3


import pandas as pd




mgi_df = pd.read_csv('HMD_HumanPhenotype.rpt', sep='\t', header=None)
mgi_df.columns = ['Gene', 1, 'Mouse_gene', 3, 'MPI', 5]


mgi_phenotypes_df = pd.read_csv('VOC_MammalianPhenotype.rpt', sep='\t', header=None)
mgi_phenotypes_df = mgi_phenotypes_df.set_index(0)
mgi_phenotypes_df['mpi_description'] = mgi_phenotypes_df.index.to_series() + '_' + mgi_phenotypes_df[1]
pdict = mgi_phenotypes_df['mpi_description'].to_dict()


# drop genes without phenotypes
mgi_df = mgi_df[~mgi_df.MPI.isna()]
mgi_df = mgi_df[['Gene', 'Mouse_gene', 'MPI']]







def get_descriptions(mpis):
	mpis = mpis.split(', ')
	result = []
	for mpi in mpis:
		result.append(pdict[mpi])
	return ', '.join(result)






mgi_df['mpi_description'] = mgi_df.MPI.apply(get_descriptions)







mgi_df.to_csv('mgi_genes.csv', index=False)













