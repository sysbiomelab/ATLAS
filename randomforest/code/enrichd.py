#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from scipy.spatial import distance
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

enrich = pd.read_csv('../data/S3enrich.tsv', sep= '\t',index_col=0)
enrich.loc[enrich.Enriched == 'Depleted', enrich.select_dtypes(include=['number']).columns] *= -1
enrich.columns = enrich.columns.str.replace(r':.*', '', regex=True)
#enrich = enrich.reset_index().set_index(['msp', 'Unnamed', 'Enriched'])
#enrich.xs(enrich.Enriched == 'Enriched', level=0)
enrich = enrich.set_index('Unnamed').drop('Enriched', axis=1)
enrich = enrich[~enrich.index.str.contains('unclassified')]
enrich = enrich.T.groupby(enrich.columns).mean().T
plotdf = enrich[(np.abs(stats.zscore(enrich)) > 4.7).any(axis=1)]
sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm', center=0, z_score=True)

sunjaemeta= pd.read_csv('../data/metadata_tables.tsv', index_col=0, sep='\t')
taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

var = 'host_phenotype'
newmeta = meta.set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
healthdf = msptaxo.join(newmeta['health_status'], how='inner').dropna()
healthydf = healthdf.loc[healthdf.health_status == 'H'].drop('health_status', axis=1).dropna()
healthydf[var] = 'Healthy'
df = healthdf.loc[healthdf.health_status != 'H'].drop('health_status', axis=1)
df = df.join(newmeta[var], how='inner').dropna()
df.host_phenotype.loc[df.host_phenotype == 'ME/CFS'] = 'ME_CFS'

feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
shaps = pd.DataFrame()
plt.rcParams["figure.figsize"] = (3,3)
for i in df[var].unique():
    testdf = df.copy()
    diseaseddf = testdf.loc[testdf[var] == i]
    nondiseasedf=healthydf.sample(len(diseaseddf), random_state=1)
    testdf = pd.concat([nondiseasedf, diseaseddf])
    classifier = RandomForestClassifier(n_estimators=500, random_state=1)
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y)
    classifier.fit(X_train, y_train)
    #plot_roc_curve(classifier, X_test, y_test, pos_label=i)
    #plt.show()
    #plot_confusion_matrix(classifier, X_test, y_test, display_labels=['H', 'D'], colorbar=False, cmap='Reds')
    #plt.show()
    #plot_precision_recall_curve(classifier, X_test, y_test, pos_label=i)
    #plt.show()
    feature_imp[i] = pd.Series(classifier.feature_importances_,index=X.columns)
    explainer = shap.Explainer(classifier)
    shaps[i] = pd.Series(explainer(X).values.sum(axis=0)[:,0], index=X.columns)


#plotdf = feature_imp[(stats.zscore(feature_imp) > 4).any(axis=1)]
plotdf = feature_imp[(feature_imp > 0.02).any(axis=1)]
sns.clustermap(plotdf, yticklabels=True, cmap='Reds')
plotdf = shaps[(np.abs(stats.zscore(shaps)) > 2).any(axis=1)]
sns.clustermap(plotdf, yticklabels=True, cmap='coolwarm')

plt.show()
#plotdf.join(taxo.set_index('species')['gp']).set_index('gp').to_csv('shap2.txt',sep=' ')
