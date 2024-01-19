#!/usr/bin/env python
import functions
import numpy as np
import pandas as pd
#from skbio.stats.composition import multiplicative_replacement as mult

### Load data
t = pd.read_csv('../data/gutTaxo.csv', index_col=0)
m = pd.read_csv('../../data/vect_atlas.csv', index_col=0).T

meta = pd.read_csv('../../data/basicMetadata.Theo.SL.2022.08.18.csv').set_index('Geography')

### Merge taxo information
mT = m.T.join(t['genus']).groupby('genus').sum()
unclass = mT[mT.index.str.contains("unclassified")].sum()
mT = mT[~mT.index.str.contains("unclassified")].T
mT[-1] = unclass
df = mT.div(mT.sum(axis=1), axis=0)
#df = pd.DataFrame(mult(df), index=df.index, columns=df.columns)
df = df.loc[df.sum(axis=1) != 0, df.sum(axis=0) !=0]
df.to_csv('../data/genusrelabund.csv')
