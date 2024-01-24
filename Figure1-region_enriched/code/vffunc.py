import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load vf data
vf_name = pd.read_csv('../data/PATRIC_sp_gene.txt', sep='\t')
vf_name['PATRIC ID'] = vf_name['PATRIC ID'].str.replace('|','_')
vf = pd.read_csv('../data/patricMap.csv', index_col=0)
vf = vf.set_index('igc2_id')
vf.vf_id = vf.vf_id.str.extract(r'(.*peg.\d+)', expand=True)

# Load gct and msp
genebag = pd.read_csv('../../../data/FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)

# Load westernisation
sortmat = pd.read_csv('../results/westernspecies.csv', index_col=0)

nwvals = sortmat.sort_values('NW').tail(100).index
wvals = sortmat.sort_values('W').tail(100).index

genebag['nwCount'] = 0
for i in nwvals:
    genebag.loc[i, 'nwCount'] = genebag.loc[i, 'nwCount'] + 1

pseudogct = genebag.loc[genebag.nwCount > 0].set_index('gene_name')
vfjoin = pseudogct.join(vf, how='inner')
nwsumvf = vfjoin.groupby('vf_id').sum()

genebag['wCount'] = 0
for i in wvals.index:
    genebag.loc[i, 'wCount'] = genebag.loc[i, 'wCount'] + 1

pseudogct = genebag.loc[genebag.wCount > 0].set_index('gene_name')
vfjoin = pseudogct.join(vf, how='inner')
wsumvf = vfjoin.groupby('vf_id').sum()

merged = pd.merge(nwsumvf['nwCount'], wsumvf['wCount'], on='vf_id', how='outer')
merged.fillna(0, inplace=True)

cutoff = merged[merged.sum(axis=1) > 500]
percentile = cutoff.T.div(cutoff.sum(axis=1)).T
extremes = percentile.loc[abs(percentile['nwCount'].sub(percentile['wCount'])).sort_values().tail(18).index]
sort = extremes.sort_values('nwCount')

sort.join(vf_name.set_index('PATRIC ID')['Gene'])
vf_name = vf_name.set_index('PATRIC ID')['Gene']

# Plot
plt.bar(x = sort.index, height=sort.iloc[:,0])
plt.bar(x = sort.index, height=sort.iloc[:,1], bottom = sort.iloc[:,0])

plt.xticks(rotation=90)
plt.ylim(0,1)
left=0.07
right=0.5
bottom=0.17
top=0.31
plt.subplots_adjust(left=left, right=right, top=top+0.2, bottom=bottom+0.2)

plt.savefig('../results/Figure_1f3.svg')
plt.show()
