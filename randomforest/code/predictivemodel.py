#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''

from sklearn.metrics import brier_score_loss
from scipy.spatial import distance
from scipy import stats
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skbio
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import plot_precision_recall_curve

taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

var = 'host_phenotype'
newmeta = meta.set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species'],how='outer').fillna(0).drop('species', axis=1).T

healthdf = msptaxo.join(newmeta[['health_status', var]]).fillna('H')
healthdf.loc[healthdf.health_status == 'H', 'host_phenotype'] = 'Healthy'
df = healthdf.copy()
#df = healthdf.loc[healthdf.health_status != 'H']
df.loc[df.host_phenotype == 'ME/CFS', 'host_phenotype'] = 'ME_CFS'
df.drop('health_status', axis=1, inplace=True)

#testdf = df.copy()
testdf = df.loc[(df.host_phenotype == 'Healthy') | (df.host_phenotype == 'LC')]
out = testdf.groupby('host_phenotype').filter(lambda x : len(x)>100)
testdf =  out.groupby('host_phenotype').sample(100)
classifier = RandomForestClassifier(n_estimators=1000)

X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
classifier.fit(X_train, y_train)
plot_roc_curve(classifier, X_test, y_test)
plot_precision_recall_curve(classifier, X_test, y_test)

fmtdf = pd.read_csv('../../../../data/FMT/downstream_data/mgs.csv', index_col=0).T
fmtmeta = pd.read_csv('../../../../data/FMT/FMT_analysis/data/newnewmetadata.csv').dropna(subset=['Sample ID'], axis=0).set_index('Sample ID')

Tfmtdf = fmtdf.T.join(taxo['species'], how='outer').fillna(0).drop('species', axis=1).T
mergd = Tfmtdf.join(fmtmeta['Type'], how='inner')
mergd.loc[mergd.Type != 'DONOR', 'Type'] = 'LC'
mergd.loc[mergd.Type == 'DONOR', 'Type'] = 'Healthy'
X = mergd.drop('Type', axis=1)
y = mergd.xs('Type', axis=1)
classifier.predict(X)

proba = pd.DataFrame(classifier.predict_proba(X), index=X.index, columns=y_test.unique())
proba['Type'] = y
probaj = proba.groupby('Type').mean()
probaj = proba.set_index('Type')
sns.clustermap(probaj, yticklabels=True)
plot_roc_curve(classifier, X, y)
plot_precision_recall_curve(classifier, X, y)
plot_confusion_matrix(classifier, X, y)
classifier.score(X,y)

