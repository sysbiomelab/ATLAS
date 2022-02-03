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
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skbio

sunjaemeta= pd.read_csv('../data/metadata_tables.tsv', index_col=0, sep='\t')
taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

diseasegroup = pd.read_csv('../data/cohort_colors.ob.txt')
diseasegroup.columns = range(3)
#split = diseasegroup[0].str.replace('.*:.*:', '')

split = diseasegroup[0].str.split(':',expand=True)
diseasegroup['disease'], diseasegroup['country'], diseasegroup['project'] = split[0], split[1], split[2]

mapping = diseasegroup[[2,'disease']].drop_duplicates()
mapping = mapping.append({'disease':'Healthy', 2:'Healthy'},ignore_index=True).set_index('disease')
newmeta = meta.set_index('host_phenotype').join(mapping[2]).reset_index().set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
var = 2
df = msptaxo.join(newmeta[[var, 'index']], how='inner').dropna()
df = df.loc[df[var] == 'Cancer'].drop(var,axis=1)
var = 'index'

feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
for i in df[var].unique():
    testdf = df.copy()
    testdf.loc[testdf[var] != i, var] = 'other'
    #nondiseasedf=testdf.loc[testdf[var] == 'other'].sample(1000)
    #diseaseddf = testdf.loc[testdf[var] == i]
    #testdf = pd.concat([nondiseasedf, diseaseddf])
    classifier = RandomForestClassifier()
    #classifier = RandomForestClassifier(class_weight={'other': 1, i: 50})
    #classifier = RandomForestClassifier(class_weight='balanced')
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, stratify=y)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    classifier.fit(X_train, y_train)
    feature_imp[i] = pd.Series(classifier.feature_importances_,index=X.columns)
    y_pred = classifier.predict(X_test)
    #classifier.score(X,y)
    #scores[i] = accuracy_score(y_test.values, y_pred)
    y_pred_proba = classifier.predict_proba(X_test)[:,1]
    fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba, pos_label=i)
    #scorecurve[i] = {'scores':roc_auc_score(y_test, y_pred_proba),'curves':[fpr, tpr]}
    scores[i] = roc_auc_score(y_test, y_pred_proba)
    curves[i] = [fpr, tpr]

#create ROC curve
for i,j in enumerate(curves):
    plt.plot(j[1], j[0], label=scores.index[i]+" AUC="+str(scores[i])[:6])
plt.plot([0, 1], [0, 1], "k--", lw=1)
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.show()

'''
sscores=scores.copy()
sfeature_imp = feature_imp.copy()

sns.heatmap(feature_imp.loc[(feature_imp > 0.012).T.any()], yticklabels=True, scale)

from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

testdf = df.copy()
classifier = RandomForestClassifier()
X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
classifier.fit(X_train, y_train)
plot_confusion_matrix(classifier, X_test, y_test)
#plot_roc_curve(classifier, X_test, y_test)
degrees = 90
plt.xticks(rotation=degrees)
plt.show()

fig, ax = plt.subplots(1,1)
ax = plot_roc_curve(classifier, X_test, y_test)
plot_roc_curve(classifier, X_test, y_test, ax=ax)

degrees = 90
plt.xticks(rotation=degrees)
plt.show()


final.sort_values(ascending=False, inplace=True)
sns.barplot(x=final, y=final.index)
plt.xlabel('Model Accuracy (%)')
plt.xlim(0,100)
plt.ylabel('Dataset used for MELD prediction')
plt.show()
plt.tight_layout()
plt.savefig('results/OvG_randomfores_MELD.pdf')

sns.clustermap(feature_imp.loc[(feature_imp > 0.01).T.any()], yticklabels=True)
