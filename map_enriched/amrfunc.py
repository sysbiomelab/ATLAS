import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

card = pd.read_csv('../data/hs_10_4_igc2.CARD.tsv', sep='\t')
card_name = pd.read_csv('../data/aro.tsv', sep='\t', index_col=0)
card_name.index = card_name.index.str.replace('ARO:','')
card["gene_name"] = card.ORF_ID.str.split(expand=True)[0]
card["gene_name"] = card["gene_name"].str[0:-2]
card = card.set_index('gene_name')

genebag = pd.read_csv('../../../FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)
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
nwvals = dft.loc[sortmat.sort_values('NW').tail(200).index].infer_objects()
wvals = dft.loc[sortmat.sort_values('W').tail(200).index].infer_objects()
#cardvar='ARO'
cardvar='Best_Hit_ARO'

genebag['nwCount'] = 0
for i in nwvals.index:
    genebag.loc[i, 'nwCount'] = genebag.loc[i, 'nwCount'] + 1
pseudogct = genebag.loc[genebag.nwCount > 0].set_index('gene_name')
cardjoin = pseudogct.join(card, how='inner')
nwsumcard = cardjoin.groupby(cardvar).sum()

genebag['wCount'] = 0
for i in wvals.index:
    genebag.loc[i, 'wCount'] = genebag.loc[i, 'wCount'] + 1
pseudogct = genebag.loc[genebag.wCount > 0].set_index('gene_name')
cardjoin = pseudogct.join(card, how='inner')
wsumcard = cardjoin.groupby(cardvar).sum()

merged = pd.merge(nwsumcard['nwCount'], wsumcard['wCount'], on=cardvar, how='outer')
merged.fillna(0, inplace=True)

cutoff = merged[merged.sum(axis=1) > 2]
percentile = cutoff.T.div(cutoff.sum(axis=1)).T
extremes = percentile.loc[abs(percentile['nwCount'].sub(percentile['wCount'])).sort_values().tail(18).index]
sort = extremes.sort_values('nwCount')

sort.index = sort.index.astype(str)

plt.bar(x = sort.index, height=sort.iloc[:,0])
plt.bar(x = sort.index, height=sort.iloc[:,1], bottom = sort.iloc[:,0])
plt.xticks(rotation=90)

plt.ylim(0,1)
left=0.07
right=0.5
bottom=0.17
top=0.31
plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
#plt.tight_layout()
plt.savefig('../results/Figure_1f2.pdf')
plt.show()
