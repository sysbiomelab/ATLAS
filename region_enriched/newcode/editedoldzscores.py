import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as PCA
from sklearn.preprocessing import StandardScaler
import skbio

country_codes = pd.read_csv('../data/countrycodes.tsv', sep='\t',  index_col=1)
taxo = pd.read_csv('../../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../../data/vect_atlas.csv', index_col=0)

msptaxo = msp.join(taxo['species']).groupby('species').sum().T

meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv').set_index('country')
meta = meta.join(country_codes).set_index('secondary_sample_accession')
meta = meta.loc[meta.health_status == 'H']

nmsptaxo = msptaxo.T.loc[~msptaxo.columns.str.contains('unclassified')]
regionmetamsptaxo = nmsptaxo.T.join(meta[['westernised','Country']], how='inner')
#regionmetamsptaxo = msp.T.join(meta[['westernised','Country']])
metamsptaxo = regionmetamsptaxo.groupby('Country').mean()

countrymap = meta[['Country', 'westernised']].groupby('Country').first()

from scipy.stats import zscore

zscores = metamsptaxo.apply(zscore)

zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.join(countrymap)

dft = df.T
var = 'westernised'
sortmat = df.groupby(var).mean().T
sortmat = sortmat.loc[~sortmat.index.str.contains('unclassified')]
'''
#enriched = zscores.T.loc[(zscores > 2.5).any()] 
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
sortmat['difference'] = sortmat.W - sortmat.NW
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
s = scaler.fit_transform(sortmat)
scaled = pd.DataFrame(index=sortmat.index,columns=sortmat.columns,data=s)
scaled['1mindiff'] = 1 - scaled.difference

Ar_dist = distance.squareform(distance.pdist(metamsptaxo, metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
import plotly.express as px

fig = px.scatter(pd.DataFrame(data=Ar_dist,columns=metamsptaxo.columns, index=metamsptaxo.columns))
sns.clustermap(pd.DataFrame(data=Ar_dist,columns=metamsptaxo.index, index=metamsptaxo.index))
'''

#sns.diverging_palette(220, 20, as_cmap=True)
sns.clustermap(data=plotdf, cmap='vlag',center = 0 , yticklabels=True, method='spearman')
functions.clustermap(plotdf)
#sns.heatmap(zscores.T.loc[(sortmat.T > 0.65).any()],xticklabels=True, yticklabels=True, cmap='coolwarm')

plt.savefig('../results/Figure_1d.pdf')
plt.show()
