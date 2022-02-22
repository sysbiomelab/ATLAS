import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as PCA
from sklearn.preprocessing import StandardScaler

country_codes = pd.read_csv('../data/countrycodes.tsv', sep='\t',  index_col=1)
taxo = pd.read_csv('../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0)

msptaxo = msp.join(taxo['species']).groupby('species').sum().T

meta = pd.read_csv('../data/unique_metadata.csv').set_index('country')
meta = meta.join(country_codes).set_index('secondary_sample_accession')

#regionmetamsptaxo = msptaxo.join(meta[['westernised','Country']])
regionmetamsptaxo = msp.T.join(meta[['westernised','Country']])
metamsptaxo = regionmetamsptaxo.groupby('Country').mean()

countrymap = meta[['Country', 'westernised']].groupby('Country').first()

from scipy.stats import zscore

'''
zscores = pd.DataFrame(index=metamsptaxo.index)
for j in metamsptaxo.columns:
    data = metamsptaxo.loc[:,j]
    newrow = []
    for i in data:
        newrow.append((i - data.mean())/data.std())
    zscores[j] = newrow
'''
zscores = metamsptaxo.apply(zscore)

zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.join(countrymap)
dft = df.T
var = 'westernised'
sortmat = df.groupby(var).mean().T
sortmat = sortmat.loc[~sortmat.index.str.contains('unclassified')]
nwvals = dft.loc[sortmat.sort_values('NW').tail(20).index]
wvals = dft.loc[sortmat.sort_values('W').tail(20).index]

plotdf = pd.concat([nwvals, wvals]).infer_objects()
plotdf = plotdf.astype('float').infer_objects()
'''
sortmat['difference'] = sortmat.W - sortmat.NW
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
s = scaler.fit_transform(sortmat)
scaled = pd.DataFrame(index=sortmat.index,columns=sortmat.columns,data=s)
scaled['1mindiff'] = 1 - scaled.difference
'''

#sns.diverging_palette(220, 20, as_cmap=True)
sns.clustermap(data=plotdf, cmap='vlag',center = 0 , yticklabels=True)

plt.savefig('../results/Figure_1d.pdf')
plt.show()
