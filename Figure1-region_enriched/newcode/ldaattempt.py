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

regionmetamsptaxo = msptaxo.join(meta[['westernised','Country']])
metamsptaxo = regionmetamsptaxo.groupby('Country').mean()

countrymap = meta[['Country', 'westernised']].groupby('Country').first()

zscores = pd.DataFrame(index=metamsptaxo.index)
for j in metamsptaxo.columns:
    data = metamsptaxo.loc[:,j]
    newrow = []
    for i in data:
        newrow.append((i - data.mean())/data.std())
    zscores[j] = newrow


zscores.replace([np.inf, -np.inf], np.nan, inplace=True)
zscores.dropna(axis=1, inplace=True)
df = zscores.join(countrymap).T

var = 'westernised'

'''
scaledDf = StandardScaler().fit_transform(df.drop(var, axis=0))
labscaledDf = pd.DataFrame(scaledDf, index=df.drop(var, axis=0).index,columns=df.drop(var, axis=0).columns)
'''
lda = PCA().fit(scaledDf, df.loc['westernised'])
results = lda.transform(scaledDf)
#coef = pd.DataFrame(lda.coef_)
df['LDA_score'] = results[:,0]

sns.clustermap(zscores)

plt.show()
