import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

vf_name = pd.read_csv('../data/PATRIC_sp_gene.txt', sep='\t')
vf_name['PATRIC ID'] = vf_name['PATRIC ID'].str.replace('|','_')
vf = pd.read_csv('../data/patricMap.csv', index_col=0)
vf = vf.set_index('igc2_id')
vf.vf_id = vf.vf_id.str.extract(r'(.*peg.\d+)', expand=True)

genebag = pd.read_csv('../../../data/FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)
country_codes = pd.read_csv('../data/countrycodes.tsv', sep='\t',  index_col=1)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0).T

meta = pd.read_csv('../data/unique_metadata.csv').set_index('country')
meta = meta.join(country_codes).set_index('secondary_sample_accession')

regionmetamsp= msp.join(meta[['westernised','Country']])
metamsp= regionmetamsp.groupby('Country').mean()

countrymap = meta[['Country', 'westernised']].groupby('Country').first()

from scipy.stats import zscore
zscores = metamsp.apply(zscore)

zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.join(countrymap)
dft = df.T
var = 'westernised'
sortmat = df.groupby(var).mean().T
nwvals = dft.loc[sortmat.sort_values('NW').tail(100).index].infer_objects()
wvals = dft.loc[sortmat.sort_values('W').tail(100).index].infer_objects()

genebag['nwCount'] = 0
for i in nwvals.index:
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
newvf_name = pd.read_csv('vftable.csv',sep='\t',index_col=0)
sort = sort.join(newvf_name).set_index('gene_name')

plt.bar(x = sort.index, height=sort.iloc[:,0])
plt.bar(x = sort.index, height=sort.iloc[:,1], bottom = sort.iloc[:,0])

plt.xticks(rotation=90)
plt.ylim(0,1)
left=0.07
right=0.5
bottom=0.17
top=0.31
plt.subplots_adjust(left=left, right=right, top=top+0.2, bottom=bottom+0.2)

plt.savefig('../results/Figure_1f3.pdf')
plt.show()
