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
healthdf = msptaxo.join(newmeta['health_status'], how='inner').dropna()
healthydf = healthdf.loc[healthdf.health_status == 'H'].drop('health_status', axis=1).dropna()
healthydf[var] = 'Healthy'
df = healthdf.loc[healthdf.health_status != 'H'].drop('health_status', axis=1)
df = df.join(newmeta[var], how='inner').dropna()
df.host_phenotype.loc[df.host_phenotype == 'ME/CFS'] = 'ME_CFS'

excludedDiseases = [
'acute diarrhea',
'atherosclerosis',
'large adenoma',
'NGT',
'premature_734g',
'premature_2847g',
'Overweight',
'Obese',
'excluded controls',
'history of colorectal surgery',
'small adenoma',
'advanced adenoma',
'PD',
'GDM',
'UC',
'control',
'ob',
'BD',
'Underweight']
df = df.loc[~df.host_phenotype.isin(excludedDiseases)]

feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
shaps = pd.DataFrame()
plt.rcParams["figure.figsize"] = (7,7)

for i in df[var].unique():
    testdf = df.copy()
    diseaseddf = testdf.loc[testdf[var] == i]
    nondiseasedf=healthydf.sample(len(diseaseddf), random_state=1)
    testdf = pd.concat([nondiseasedf, diseaseddf])
    classifier = RandomForestClassifier(n_estimators=500, n_jobs=-1,random_state=1)
    X = testdf.drop(var, axis=1)
    y = pd.get_dummies(testdf.xs(var, axis=1))[i]
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state = 1)
    classifier.fit(X_train, y_train)
    #pickle.dump(classifier, open(i, 'wb'))
    #plot_roc_curve(classifier, X_test, y_test, pos_label=1)
    #plt.show()
    #plot_confusion_matrix(classifier, X_test, y_test, display_labels=['H', 'D'], colorbar=False, cmap='Reds')
    #plt.show()
    #plot_precision_recall_curve(classifier, X_test, y_test, pos_label=1)
    #plt.show()
    y_pred = classifier.predict(X_test)
    y_pred_proba = classifier.predict_proba(X_test)[:,1]
    fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba, pos_label=1)
    scores[i] = roc_auc_score(y_test, y_pred_proba)
    feature_imp[i] = pd.Series(classifier.feature_importances_,index=X.columns)
    explainer = shap.TreeExplainer(classifier)
    shaps_values = explainer(X)
    meanabsshap = pd.Series( np.abs(shaps_values[:, :, 0].values).mean(axis=0), index=X.columns)
    corrs = [spearmanr(shaps_values[:, x, 1].values, X.iloc[:,x])[0] for x in range(len(X.columns))]
    final = meanabsshap * np.sign(corrs)
    final.fillna(0, inplace=True)
    shaps[i] = final

plotdf = shaps[(np.abs(shaps) > 0.0125).any(axis=1)]
plotdf = plotdf[~plotdf.index.str.contains('unclassified')]
plotdf = plotdf.reindex(sorted(plotdf.columns), axis=1)

ax = sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm', center=0,xticklabels=True)
#plt.savefig('../results/shap2.svg')
plt.show()
