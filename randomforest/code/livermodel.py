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

livermeta = meta.loc[meta.host_phenotype == 'LC'].study_accession.unique() 
liverstudies = meta[meta.study_accession.isin(livermeta)]
PRJEB6337 = liverstudies.loc[liverstudies.study_accession == 'PRJEB6337']
PRJEB38481  = liverstudies.loc[liverstudies.study_accession == 'PRJEB38481']

var = 'host_phenotype'
testdf = msp.T.join(PRJEB6337.set_index('secondary_sample_accession')[var], how='inner')
ntestdf = msp.T.join(PRJEB38481.set_index('secondary_sample_accession')[var], how='inner')
variables = ['LC', 'Healthy']
ntestdf = ntestdf[ntestdf[var].isin(variables)]

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
#shap.plots.waterfall(shap_values.base_values, shap_values[0], X[0])
shap.summary_plot(shap_values, X.values, plot_type="bar", class_names= class_names, feature_names = X.columns)
shaps = pd.DataFrame(shap_values.values.sum(axis=0), index=X.columns, columns=y.unique())
from scipy import stats
plotdf = shaps[(np.abs(stats.zscore(shaps)) > 4).any(axis=1)]
