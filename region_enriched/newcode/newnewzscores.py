#!/usr/bin/env python
from scipy.stats import zscore
import functions
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skbio

### Load data
taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../../data/vect_atlas.csv', index_col=0)
countries = pd.read_csv('../../data/countries.csv', index_col=2)
meta = pd.read_csv('../../data/basicMetadata.Theo.SL.2022.08.18.csv').set_index('Geography')
#df = functions.enterotypes(msp, taxo)

### Merge taxo information
msptaxo = msp.join(taxo['species']).groupby('species').sum().T

### Select only classified healthy smpales
meta = meta.loc[meta.type == 'Control']
westernmeta = meta.join(countries['westernised']).reset_index().set_index('sample.ID')
nmsptaxo = msptaxo.T.loc[~msptaxo.columns.str.contains('unclassified')]
regionmetamsptaxo = nmsptaxo.T.join(westernmeta[['westernised','Geography']], how='inner')
#regionmetamsptaxo.to_csv('../results/healthymsp.csv')
#regionmetamsptaxo.join(meta, how='inner')['BioProject'].reset_index().rename(columns={'index':'sample.ID','BioProject':'project.ID'}).to_csv('../results/monoclemeta.csv')

### Calculate difference between western/non-western
df, grouping = regionmetamsptaxo.drop('Geography',axis=1), 'westernised'
qvalues = functions.MANNWHIT(df, grouping)
group1 = df.groupby(grouping).get_group(df[grouping].unique()[0]).drop(grouping, axis=1)
group2 = df.groupby(grouping).get_group(df[grouping].unique()[1]).drop(grouping, axis=1)
dat1 = group1.agg(['mean','std']).T
dat2 = group2.agg(['mean','std']).T
fcdf = pd.DataFrame()
fcdf['Proportion_Western'] = group1.nunique() / group1.shape[0]
fcdf['Mean_Western'] = dat1['mean']
fcdf['Std_Western'] = dat1['std']
fcdf['Proportion_Nonwestern'] = group2.nunique() / group2.shape[0]
fcdf['Mean_Nonwestern'] = dat2['mean']
fcdf['Std_Nonwestern'] = dat2['std']
fcdf['FC(Western_vs_Nonwestern)'] = fcdf['Mean_Western'].div(fcdf['Mean_Nonwestern'])
outdf = fcdf.join(qvalues.to_frame())
outdf.loc[outdf[0] < 0.05 ,'Significant'] = True
outdf.loc[outdf[0] > 0.05 ,'Significant'] = False
outdf.rename(columns={0:'M.W.W. qvalue'}, inplace=True)
outdf.to_csv('../results/westernnonfc.csv')

### Compute regional specific enrichment
metamsptaxo = regionmetamsptaxo.groupby('Geography').mean()
zscores = metamsptaxo.apply(zscore)
zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.copy()
dft = df.T
var = 'westernised'
sortmat = df.groupby(var).mean().T
sortmat = df.groupby(level=0).mean().T
sortmat = sortmat.loc[~sortmat.index.str.contains('unclassified')]

### Compute and plot
plt.rcParams["figure.figsize"] = (6,5)
model = functions.PCOA(regionmetamsptaxo.drop(['Geography','westernised'],axis=1))
fig, ax = plt.subplots()
functions.PCplot(model.join(regionmetamsptaxo), 'westernised', ax=ax)
ax.set_ylabel('PCoA2')
ax.set_xlabel('PCoA1')
plt.tight_layout()
plt.savefig('../results/pcoa.pdf')
plt.show()

### Compute and plot
plt.rcParams["figure.figsize"] = (6,5)
model = functions.UMAP(regionmetamsptaxo.drop(['Geography','westernised'],axis=1))
fig, ax = plt.subplots()
functions.PCplot(model.join(regionmetamsptaxo), 'westernised', ax=ax)
ax.set_ylabel('UMAP2')
ax.set_xlabel('UMAP1')
plt.tight_layout()
plt.savefig('../results/umap.pdf')
plt.show()

### Compute and plot
plt.rcParams["figure.figsize"] = (6,5)
model = functions.TSNE(regionmetamsptaxo.drop(['Geography','westernised'],axis=1))
fig, ax = plt.subplots()
functions.PCplot(model.join(regionmetamsptaxo), 'westernised', ax=ax)
ax.set_ylabel('TSNE2')
ax.set_xlabel('TSNE1')
plt.tight_layout()
plt.savefig('../results/tsne.pdf')
plt.show()

### Compute and plot
plt.rcParams["figure.figsize"] = (6,5)
model = functions.NMDS(regionmetamsptaxo.drop(['Geography','westernised'],axis=1))
fig, ax = plt.subplots()
functions.PCplot(model.join(regionmetamsptaxo), 'westernised', ax=ax)
ax.set_ylabel('NMDS2')
ax.set_xlabel('NMDS1')
plt.tight_layout()
plt.savefig('../results/nmds.pdf')
plt.show()

### Compute permanova
permanovaPval = PERMANOVA(regionmetamsptaxo.drop(['Geography','westernised'],axis=1), regionmetamsptaxo['westernised'])
#permanovaPval = functions.PERMANOVA(regionmetamsptaxo.drop(['Geography','westernised'],axis=1), regionmetamsptaxo['westernised'])

### Select top 20 enriched in western and not
nwvals = dft.loc[sortmat.sort_values('NW').tail(20).index]
wvals = dft.loc[sortmat.sort_values('W').tail(20).index]
plotdf = pd.concat([nwvals, wvals]).infer_objects()
plotdf = plotdf.astype('float').infer_objects()

### Plot and save clustermap
sns.clustermap(data=plotdf, cmap='vlag',center = 0 , yticklabels=True)
functions.clustermap(plotdf.T)
#sns.heatmap(zscores.T.loc[(sortmat.T > 0.65).any()],xticklabels=True, yticklabels=True, cmap='coolwarm')
sns.heatmap(zscores.loc[(zscores > 3).any()],xticklabels=True, yticklabels=True, cmap='coolwarm')
plt.savefig('../results/enrichedHeatmap.pdf')
plt.show()
