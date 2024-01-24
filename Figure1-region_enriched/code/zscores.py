from scipy.spatial import distance
from scipy.stats import zscore
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skbio

# Load data
taxo = pd.read_csv('../../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../../data/msp.csv', index_col=0)
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
meta = pd.read_csv('../../data/sampleID.csv', index_col=0)
wes = pd.read_csv('../../data/westernized.tsv', sep='\t', index_col=0).westernised.to_frame()

# Filter metadata and species
meta = meta.loc[meta.Disease == 'Healthy']
meta = meta.reset_index().set_index('Geography')
meta = meta.join(wes['westernised'])
meta = meta.reset_index().set_index('sample.ID')
meta = meta.rename(columns={'index':'Geography'})
nmsptaxo = msptaxo.T.loc[~msptaxo.columns.str.contains('unclassified')]
regionmetamsptaxo = nmsptaxo.T.join(meta[['westernised','Geography']], how='inner')
metamsptaxo = regionmetamsptaxo.groupby('Geography').mean()

# Calcualte z-scores
zscores = metamsptaxo.apply(zscore)
zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.join(countrymap)
df = zscores.join(wes)
dft = df.T
sortmat = df.groupby('westernised').mean().T
sortmat.to_csv('../results/westernspecies.csv')

# Calculate enrichment
nwvals = dft.loc[sortmat.sort_values('NW').tail(20).index]
wvals = dft.loc[sortmat.sort_values('W').tail(20).index]
plotdf = pd.concat([nwvals, wvals]).infer_objects()
plotdf = plotdf.astype('float').infer_objects()

# Plot
sns.clustermap(data=plotdf, cmap='vlag', center=0, yticklabels=True)
plt.savefig('../results/Figure_1d.svg')
plt.show()
