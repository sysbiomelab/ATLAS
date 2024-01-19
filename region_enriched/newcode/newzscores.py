#!/usr/bin/env python
import pandas as pd
import functions
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import zscore
import skbio

taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../../data/vect_atlas.csv', index_col=0)
countries = pd.read_csv('../../data/countries.csv', index_col=2)
meta = pd.read_csv('../../data/basicMetadata.Theo.SL.2022.08.18.csv').set_index('Geography')

msptaxo = msp.join(taxo['species']).groupby('species').sum().T
meta = meta.loc[meta.type == 'Control']
westernmeta = meta.join(countries['westernised']).reset_index().set_index('sample.ID')

nmsptaxo = msptaxo.T.loc[~msptaxo.columns.str.contains('unclassified')]

regionmetamsptaxo = nmsptaxo.T.join(westernmeta[['westernised','Geography']], how='inner')

### western non significance
df, grouping = regionmetamsptaxo.drop('Geography',axis=1), 'westernised'
b = functions.MANNWHIT(df, grouping)
#b.apply(np.log).sort_values().dropna().plot.barh()

metamsptaxo = regionmetamsptaxo.groupby('Geography').mean()

zscores = metamsptaxo.apply(zscore)

zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.copy()

dft = df.T
var = 'westernised'
sortmat = df.groupby(var).mean().T
sortmat = sortmat.loc[~sortmat.index.str.contains('unclassified')]
'''
enriched = zscores.T.loc[(zscores > 2.5).any()] 
new = zscores.T.join(taxo.set_index('species')['gp']).reset_index().set_index(['gp', 'index']).sort_index()
new.to_csv('countryzscores.csv')
enriched = new > 2.5
enriched.to_csv('enrichedcountryzscores.csv')

'''
nwvals = dft.loc[sortmat.sort_values('NW').tail(20).index]
wvals = dft.loc[sortmat.sort_values('W').tail(20).index]

plotdf = pd.concat([nwvals, wvals]).infer_objects()
plotdf = plotdf.astype('float').infer_objects()

'''
### Country beta difference
sortmat['difference'] = sortmat.W - sortmat.NW
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
s = scaler.fit_transform(sortmat)
scaled = pd.DataFrame(index=sortmat.index,columns=sortmat.columns,data=s)
scaled['1mindiff'] = 1 - scaled.difference

Ar_dist = distance.squareform(distance.pdist(metamsptaxo, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
import plotly.express as px

fig = px.scatter(pd.DataFrame(data=Ar_dist,columns=metamsptaxo.index, index=metamsptaxo.index))
sns.clustermap(pd.DataFrame(data=Ar_dist,columns=metamsptaxo.index, index=metamsptaxo.index))
'''

#sns.diverging_palette(220, 20, as_cmap=True)
sns.clustermap(data=plotdf, cmap='vlag',center = 0 , yticklabels=True)
#sns.heatmap(zscores.T.loc[(sortmat.T > 0.65).any()],xticklabels=True, yticklabels=True, cmap='coolwarm')
sns.heatmap(zscores.loc[(zscores > 3).any()],xticklabels=True, yticklabels=True, cmap='coolwarm')

plt.savefig('../results/Figure_1d.pdf')
plt.show()
