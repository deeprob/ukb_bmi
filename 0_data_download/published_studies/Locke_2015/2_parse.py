#!/bin/python3




import pandas as pd



df1 = pd.read_excel('Table_1.xlsx')
df2 = pd.read_excel('Table_2.xlsx')
dfe2 = pd.read_excel('Extended_Data_Table_2.xlsx')


print(df1)
print(df2)
print(dfe2)

# df1['P value'] = df1['P value'].astype(float)
# df2['P value'] = df2['P value'].astype(float)
# dfe2['P value'] = dfe2['P value'].astype(float)

df = df1.append(df2)
df = df.append(dfe2)

print(df)

# split genes (one gene per line)
new_df = []
for i, row in df.iterrows():
	genes = row['Notable gene(s)']
	beta = row['Beta']
	pvalue = row['P value']
	analysis = row['Analysis']
	for gene in genes.split(';'):
		gene = gene.split('(')[0]
		new_df.append([gene, beta, analysis])
	

new_df = pd.DataFrame(new_df, columns=['Gene', 'Beta', 'Analysis'])


print(new_df)

new_df.to_csv('Locke_genes_parsed.csv', index=False)







