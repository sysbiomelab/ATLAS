#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy import stats
from scipy.spatial import distance
from scipy.stats import spearmanr
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import skbio

taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

var = 'host_phenotype'
newmeta = meta.set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
scaler = StandardScaler()
msptaxo.loc[:] = scaler.fit_transform(msptaxo)

crcmeta = meta.loc[meta.host_phenotype == 'CRC'].study_accession.unique() 
crcstudies = meta[meta.study_accession.isin(crcmeta)]
PRJDB4176 = crcstudies.loc[crcstudies.study_accession == 'PRJDB4176']
PRJEB6070  = crcstudies.loc[crcstudies.study_accession == 'PRJEB6070']

testdf = msptaxo.join(PRJDB4176.set_index('secondary_sample_accession')[var], how='inner')
ntestdf = msptaxo.join(PRJEB6070.set_index('secondary_sample_accession')[var], how='inner')
variables = ['CRC', 'Healthy']
ntestdf = ntestdf[ntestdf[var].isin(variables)]
testdf = testdf[testdf[var].isin(variables)]

classifier = RandomForestClassifier(n_estimators=500, n_jobs=-1,random_state=1)
X = testdf.drop(var, axis=1)
y = pd.get_dummies(testdf.xs(var, axis=1))['CRC']
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state = 1)
classifier.fit(X_train, y_train)
plot_confusion_matrix(classifier, X_test, y_test)
plot_roc_curve(classifier, X_test, y_test)
#plot_precision_recall_curve(classifier, X_test, y_test)
#classifier.score(X_test, y_test)
#degrees = 90
#plt.xticks(rotation=degrees)
plt.savefig('../results/inter.svg')
plt.show()

nX = ntestdf.drop(var, axis=1)
ny = pd.get_dummies(ntestdf.xs(var, axis=1))['CRC']
#ny = ntestdf.xs(var, axis=1)
classifier.score(nX, ny)

plot_confusion_matrix(classifier, nX, ny)
plot_roc_curve(classifier, nX, ny)
plot_precision_recall_curve(classifier, nX, ny)
plt.savefig('../results/intra.svg')
plt.show()

explainer = shap.Explainer(classifier)
shap_values = explainer(X)
#shap.summary_plot(shap_values[:,:,0], X.values, plot_type="bar", class_names= class_names, feature_names = X.columns)
