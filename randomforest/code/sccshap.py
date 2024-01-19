#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import functions
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap

### Load data
guttaxo = pd.read_csv('../../data/gutTaxo.csv', index_col=0)
guttaxo['species'] = guttaxo.species.str.replace(r'/.*','',regex=True).str.replace(r'&.*','', regex=True)
msp = pd.read_csv('../../data/vect_atlas.csv', index_col=0).T
allmeta = pd.read_csv('../../data/basicMetadata.Theo.SL.2022.08.18.csv', index_col=0)
hgmi = pd.read_csv('../../data/hgmiedges.csv', index_col=1)
meta = hgmi.join(allmeta).reset_index().set_index('hgmi_id')
meta.loc[meta['Geography'] == 'United States of America (the)', 'Geography'] = 'USA'

meta['Project details'] = meta['Disease'] +' ' +meta['Geography'] + ' ' + meta['BioProject']
var = 'Disease'

### Shap calculation for each disease
shaps = pd.DataFrame()
hgmaid = meta.index.unique()[13]
for hgmaid in meta.index.unique():
    df = msp.join(meta.loc[hgmaid].set_index('sample.ID')[var],how='inner').set_index(var)
    df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns).reset_index()
    df['Disease'] = df['Disease'].str.replace('Healthy','0')
    final = functions.RFC(df.drop(var, axis=1), df[var])
    shapvals = functions.SHAP_bin(df.drop(var, axis=1), final)
    shaps[hgmaid] = shapvals

plotdf = shaps.T.join(meta.loc[meta.Disease != 'Healthy'].groupby(level=0).first()['Project details']).set_index('Project details').T.join(guttaxo['species']).reset_index().set_index(['index','species']).T

plotdf.T.to_csv('../results/shaps.csv')
plotdf = plotdf.droplevel(0, axis=1)

plotdf = plotdf.loc[:, ~plotdf.columns.str.contains('unclassified')]

val = 0.008
#fplotdf = plotdf.loc[:,plotdf.abs().gt(val).any(axis=0)].T

functions.setupplot()
#sig = pd.DataFrame(index=plotdf.index, columns=plotdf.columns).fillna(False)
#i = plotdf.index[0]
#for i in plotdf.index:
#    sig.loc[i, plotdf.loc[i].abs().sort_values().tail(5).index] = True
#sig = sig.loc[fplotdf.columns, fplotdf.index,].T
#fplotdf = plotdf.loc[:, sig.any()].T
fplotdf = plotdf.loc[:, plotdf.abs().gt(val).any()].T
#fsig = sig.loc[:, sig.any()].T

functions.clustermap(
        fplotdf,
        fplotdf.abs().gt(val),
        figsize=(5.2,8),
        )
plt.savefig('../results/shapannot.svg')
plt.show()
