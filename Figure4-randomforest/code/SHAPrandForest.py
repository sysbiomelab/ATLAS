#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from sklearn.preprocessing import StandardScaler
import functions as f
import matplotlib.pyplot as plt
import pandas as pd

# Load data
taxo = pd.read_csv('../../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../../data/msp.csv', index_col=0)
meta = pd.read_csv('../../data/sampleID.csv', index_col=0)

# add species and standard scale
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
msptaxo = pd.DataFrame(StandardScaler().fit_transform(msptaxo), index=msptaxo.index, columns=msptaxo.columns)

# filter metadata
conditioncount = meta.reset_index().groupby(['BioProject', 'Disease']).nunique()['sample.ID']

# needs to have a disease
prjs = conditioncount.groupby(level=0).count().gt(1)
prjs = prjs.loc[prjs]

# minimum have to have more than 20 samples
nprjs = conditioncount.groupby(level=0).min().gt(10)
nprjs = nprjs.loc[nprjs]

# needs to have a healthy arm
nnprjs = conditioncount.xs('Healthy', level=1).gt(0)
nnprjs = nnprjs.loc[nnprjs]

# final filtering
finalindex = pd.concat([prjs, nprjs, nnprjs], join='inner', axis=1).index
fmeta = meta.loc[meta.BioProject.isin(finalindex)]

# run models
columns = fmeta.loc[fmeta.Disease != 'Healthy'].groupby(['BioProject','Disease']).count().index
shaps = pd.DataFrame(columns=columns)
aucrocs = pd.DataFrame(columns=columns)
prj = finalindex[0]
for prj in finalindex:
    ffmeta = fmeta.loc[fmeta.BioProject == prj]
    disease = ffmeta.loc[ffmeta.Disease != 'Healthy'].Disease.unique()[0]
    for disease in ffmeta.loc[ffmeta.Disease != 'Healthy'].Disease.unique():
        fffmeta = ffmeta.loc[(ffmeta.Disease == 'Healthy') | (ffmeta.Disease == disease)]
        strat = f.stratify(msptaxo, fffmeta, 'type')
        strat = f.upsample(strat)
        a, b, c = f.classifier(strat)
        shaps[(prj, disease)] = f.SHAP_bin(strat, a)

# filter shaps
fshaps = shaps.copy()
fshaps.index = fshaps.index.str.replace('\/.*','', regex=True)
fshaps.index = fshaps.index.str.replace('\(.*','', regex=True)
fshaps.index = fshaps.index.str.replace('\&.*','', regex=True)
fshaps = fshaps.loc[fshaps.abs().gt(0.013).any(axis=1)]

# countrymapping
countries = meta.groupby('BioProject').max()['Geography']
fshaps = fshaps.T.reset_index().set_index('BioProject').join(countries.to_frame()).set_index(['Geography','Disease']).T

f.setupplot()
f.clustermap(fshaps, metric='euclidean', figsize=(4,7))
f.savefig('shapscluster')
f.save(fshaps, 'shaps')
