#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy.spatial import distance
from sklearn import metrics
from scipy import stats
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

crcmeta = meta.loc[meta.host_phenotype == 'CRC'].study_accession.unique() 
meta.loc[meta.host_phenotype == 'CRC'].study_accession.value_counts()
crcstudies = meta[meta.study_accession.isin(crcmeta)]

#new
PRJEB10878 = crcstudies.loc[crcstudies.study_accession == 'PRJEB10878']
#PRJEB10878 = PRJEB10878.loc[:, PRJEB10878.nunique() > 30]

PRJDB4176 = crcstudies.loc[crcstudies.study_accession == 'PRJDB4176']
PRJEB6070  = crcstudies.loc[crcstudies.study_accession == 'PRJEB6070']

var = 'host_phenotype'

#testdf = msp.T.join(PRJDB4176.set_index('secondary_sample_accession')[var], how='inner')
testdf = msp.T.join(PRJEB6070.set_index('secondary_sample_accession')[var], how='inner')
#ntestdf = msp.T.join(PRJEB6070.set_index('secondary_sample_accession')[var], how='inner')
ntestdf = msp.T.join(PRJEB10878.set_index('secondary_sample_accession')[var], how='inner')
variables = ['CRC', 'Healthy']
ntestdf = ntestdf[ntestdf[var].isin(variables)]
testdf = testdf[testdf[var].isin(variables)]
'''
# new go at randomly sampling from healthy pool
nX = ntestdf.drop(var, axis=1)
ny = ntestdf.xs(var, axis=1)
end = pd.DataFrame()
scores = []
for i in range(30):
    hsamples = meta.loc[meta['host_phenotype'] == 'Healthy'].sample(91, random_state=i)
    healthydf = msp.T.join(hsamples.set_index('secondary_sample_accession')[var], how='inner')
    healthy_testdf = testdf.loc[testdf.host_phenotype == 'CRC'].append(healthydf)
    classifier = RandomForestClassifier(random_state=i)
    X = healthy_testdf.drop(var, axis=1)
    y = healthy_testdf.xs(var, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=i)
    classifier.fit(X_train, y_train)
    explainer = shap.Explainer(classifier)
    shap_values = explainer(X)
    shapdf = pd.DataFrame(shap_values.values[0], index=X.columns, columns=y.unique())
    end[i] = shapdf['CRC']
    y_pred = classifier.predict(X_test)
    y_pred_proba = classifier.predict_proba(X_test)[:,1]
    scores.append(('self', roc_auc_score(y_test, y_pred_proba)))
    y_pred = classifier.predict(nX)
    y_pred_proba = classifier.predict_proba(nX)[:,1]
    scores.append(('other', roc_auc_score(ny, y_pred_proba)))
    #tshapdf = shapdf.join(taxo['species']).groupby('species').mean()
    tshapdf = end.join(taxo['species']).groupby('species').mean()
    ftshapdf = tshapdf.loc[np.abs(tshapdf > 0.008).any(axis=1)]
    #ftshapdf = tshapdf.loc[np.abs(tshapdf.CRC) > 0.008]
    #sns.heatmap(ftshapdf, cmap= 'vlag',square=True)

selfm = np.array(scores)[::2,1].astype(float).mean()
selfstd = np.array(scores)[::2,1].astype(float).std()
otherm = np.array(scores)[1::2,1].astype(float).mean()
otherstd = np.array(scores)[1::2,1].astype(float).std()

scores[1::2]
#ftshapdf = tshapdf.loc[np.abs(tshapdf.mean(axis=1)) > 0.005]
ftshapdf = tshapdf.loc[np.abs(stats.zscore(tshapdf.mean(axis=1))) > 2]
#ftshapdf = tshapdf.loc[:, (np.abs(tshapdf.apply(stats.zscore, axis=1)) > 3).any(axis=0)]
sns.heatmap(ftshapdf, cmap= 'vlag',square=True, yticklabels=True, center=0)
sns.swarmplot(data = ftshapdf.T, size=2)
sns.barplot(data = ftshapdf.T)
degrees = 90
plt.xticks(rotation=degrees)
plt.show()

'''
'''
nX = ntestdf.drop(var, axis=1)
ny = ntestdf.xs(var, axis=1)
hsamples = meta.loc[meta['host_phenotype'] == 'Healthy'].sample(91, random_state=1)
healthydf = msp.T.join(hsamples.set_index('secondary_sample_accession')[var], how='inner')
healthy_testdf = testdf.loc[testdf.host_phenotype == 'CRC'].append(healthydf)
classifier = RandomForestClassifier(random_state=1)
X = healthy_testdf.drop(var, axis=1)
y = healthy_testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)
classifier.fit(X_train, y_train)
plot_roc_curve(classifier, X_test, y_test, pos_label='CRC')
explainer = shap.Explainer(classifier)
shap_values = explainer(X)
shapdf = pd.DataFrame(shap_values.values[0], index=X.columns, columns=y.unique())
end[i] = shapdf['CRC']
#tshapdf = shapdf.join(taxo['species']).groupby('species').mean()
tshapdf = end.join(taxo['species']).groupby('species').mean()
ftshapdf = tshapdf.loc[np.abs(tshapdf > 0.008).any(axis=1)]
#ftshapdf = tshapdf.loc[np.abs(tshapdf.CRC) > 0.008]
#sns.heatmap(ftshapdf, cmap= 'vlag',square=True)
classifier.score(nX, ny)
plot_confusion_matrix(classifier, nX, ny)
plot_roc_curve(classifier, nX, ny, pos_label='CRC')
plot_precision_recall_curve(classifier, nX, ny)
'''


classifier = RandomForestClassifier()
X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
classifier.fit(X_train, y_train)
plot_confusion_matrix(classifier, X_test, y_test)
plot_roc_curve(classifier, X_test, y_test)
plot_precision_recall_curve(classifier, X_test, y_test)
classifier.score(X_test, y_test)
degrees = 90
plt.xticks(rotation=degrees)
plt.show()

nX = ntestdf.drop(var, axis=1)
ny = ntestdf.xs(var, axis=1)
classifier.score(nX, ny)

plot_confusion_matrix(classifier, nX, ny)
plot_roc_curve(classifier, nX, ny)
plot_precision_recall_curve(classifier, nX, ny)

explainer = shap.Explainer(classifier)
shap_values = explainer(X)
shapdf = pd.DataFrame(shap_values.values[0], index=X.columns, columns=y.unique())
tshapdf = shapdf.join(taxo['species']).groupby('species').mean()
ftshapdf = tshapdf.loc[np.abs(tshapdf.CRC) > 0.008]
sns.heatmap(ftshapdf, cmap= 'vlag',square=True)
