import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cazy_dorin = pd.read_csv('../../../Dorines_supplementary.table.3_CAZymes_subsConversion.tsv', sep='\t', index_col=0) 
cazy_name = pd.read_csv('../../../ATLAS/data/CAZyDB.07292021.fam-activities.txt', sep='\t', index_col=0)
genebag = pd.read_csv('../../../../data/FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)
cazyme = pd.read_csv('../../../../data/FMT/gutdataverse_files/IGC2_vs_cazy.table', sep='\t', index_col=0)
country_codes = pd.read_csv('../../../ATLAS/data/countrycodes.tsv', sep='\t',  index_col=1)
msp = pd.read_csv('../../../ATLAS/data/vect_atlas.csv', index_col=0).T
meta = pd.read_csv('../../../ATLAS/data/unique_metadata.csv').set_index('country')

cazyme.columns = ['cazyme']
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
cazy = pseudogct.join(cazyme, how='inner')
nwsumcazy = cazy.groupby('cazyme').sum()

genebag['wCount'] = 0
for i in wvals.index:
    genebag.loc[i, 'wCount'] = genebag.loc[i, 'wCount'] + 1
pseudogct = genebag.loc[genebag.wCount > 0].set_index('gene_name')
cazy = pseudogct.join(cazyme, how='inner')
wsumcazy = cazy.groupby('cazyme').sum()

merged = pd.merge(nwsumcazy['nwCount'], wsumcazy['wCount'], on='cazyme', how='outer')
merged.fillna(0, inplace=True)

cutoff = merged[merged.sum(axis=1) > 90]
percentile = cutoff.T.div(cutoff.sum(axis=1)).T
extremes = percentile.loc[abs(percentile['nwCount'].sub(percentile['wCount'])).sort_values().tail(18).index]
sort = extremes.sort_values('nwCount')
sort = sort.join(cazy_dorin, how='inner').drop_duplicates()

plt.bar(x = sort.index, height=sort.nwCount)
plt.bar(x = sort.index, height=sort.wCount, bottom = sort.nwCount)
plt.xticks(rotation=90)

plt.ylim(0,1)
left=0.07
right=0.5
bottom=0.17
top=0.31
adj = 0.1
plt.subplots_adjust(left=left, right=right, top=top+adj, bottom=bottom+adj)

#plt.tight_layout()
plt.savefig('../results/Figure_1f1.pdf')
plt.show()
