import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load card data
card = pd.read_csv('../data/hs_10_4_igc2.CARD.tsv', sep='\t')
card_name = pd.read_csv('../data/aro.tsv', sep='\t', index_col=0)
card_name.index = card_name.index.str.replace('ARO:','')
card["gene_name"] = card.ORF_ID.str.split(expand=True)[0]
card["gene_name"] = card["gene_name"].str[0:-2]
card = card.set_index('gene_name')

# Load gct and msp
genebag = pd.read_csv('../../../data/FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)

# Load westernisation
sortmat = pd.read_csv('../results/westernspecies.csv', index_col=0)

nwvals = sortmat.sort_values('NW').tail(100).index
wvals = sortmat.sort_values('W').tail(100).index

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

# Plot
plt.bar(x = sort.index, height=sort.iloc[:,0])
plt.bar(x = sort.index, height=sort.iloc[:,1], bottom = sort.iloc[:,0])
plt.xticks(rotation=90)

plt.ylim(0,1)
left=0.07
right=0.5
bottom=0.17
top=0.31
plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
plt.savefig('../results/Figure_1f2.svg')
plt.show()
