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

taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

diseasegroup = pd.read_csv('../data/cohort_colors.ob.txt')
diseasegroup.columns = range(3)
#split = diseasegroup[0].str.replace('.*:.*:', '')

split = diseasegroup[0].str.split(':',expand=True)
diseasegroup['disease'], diseasegroup['country'], diseasegroup['project'] = split[0], split[1], split[2]

var = 'host_phenotype'
#var = 2
mapping = diseasegroup[[2,'disease']].drop_duplicates()
mapping = mapping.append({'disease':'Healthy', 2:'Healthy'},ignore_index=True).set_index('disease')
#newmeta = meta.set_index(var).join(mapping[2]).set_index('secondary_sample_accession')
newmeta = meta.set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
healthdf = msptaxo.join(newmeta['health_status'], how='inner').dropna()
healthydf = healthdf.loc[healthdf.health_status == 'H'].drop('health_status', axis=1).dropna()
healthydf[var] = 'Healthy'
df = healthdf.loc[healthdf.health_status != 'H'].drop('health_status', axis=1)
df = df.join(newmeta[var], how='inner').dropna()

df.host_phenotype.loc[df.host_phenotype == 'ME/CFS'] = 'ME_CFS'
df = df.loc[df[var] != 'acute diarrhea'] 
df = df.loc[df[var] != 'atherosclerosis'] 
df = df.loc[df[var] != 'large adenoma'] 
df = df.loc[df[var] != 'history of colorectal surgery'] 
df = df.loc[df[var] != 'small adenoma'] 
df = df.loc[df[var] != 'PD'] 
df = df.loc[df[var] != 'GDM'] 
df = df.loc[df[var] != 'UC'] 
df = df.loc[df[var] != 'BD'] 
df = df.loc[df[var] != 'ob'] 
df = df.loc[df[var] != 'control'] 

scaler = StandardScaler()
df.loc[:, df.columns != var] = scaler.fit_transform(df.loc[:, df.columns != var])

feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
plt.rcParams["figure.figsize"] = (3,3)
for i in df[var].unique():
    testdf = df.copy()
    nondiseasedf=healthydf.sample(40, random_state=1)
    diseaseddf = testdf.loc[testdf[var] == i].sample(40, random_state=1)
    testdf = pd.concat([nondiseasedf, diseaseddf])
    classifier = RandomForestClassifier(n_estimators=500, random_state=1)
    #classifier = RandomForestClassifier(class_weight={'other': 1, i: 50})
    #classifier = RandomForestClassifier(class_weight='balanced')
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y)
    classifier.fit(X_train, y_train)
    #plot_roc_curve(classifier, X_test, y_test, pos_label=i)
    #plt.show()
    print(i)
    #plot_confusion_matrix(classifier, X_test, y_test, display_labels=['H', 'D'], colorbar=False, cmap='Reds')
    #plt.show()
    #plot_precision_recall_curve(classifier, X_test, y_test, pos_label=i)
    #plt.show()
    feature_imp[i] = pd.Series(classifier.feature_importances_,index=X.columns)
    y_pred = classifier.predict(X_test)
    y_pred_proba = classifier.predict_proba(X_test)[:,1]
    fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba, pos_label=i)
    #scores[i] = roc_auc_score(y_test, y_pred_proba)
    score = roc_auc_score(y_test, y_pred_proba)
    curves[i] = [fpr, tpr]
    #plt.title(f"{i} AUC={str(score)[:6]}")
    #plt.savefig(f"{i}.svg")

#left=0.1
#right=0.5
#bottom=0.1
#top=0.5
#create ROC curve
for i,j in enumerate(curves):
    plt.rcParams["figure.figsize"] = (5,5)
    plt.plot(j[0], j[1], label=scores.index[i]+" AUC="+str(scores[i])[:6])
    plot_roc_curve(classifier, X_test, y_test)
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    plt.show()
plt.plot([0, 1], [0, 1], "k--", lw=1)
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.show()

sns.clustermap(feature_imp.loc[(feature_imp > 0.02).T.any()], yticklabels=True)
plt.show()

# confusion matrix
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
plot_roc_curve(classifier, X_test, y_test)
plot_precision_recall_curve(classifier, X_test, y_test, pos_label='CRC')

degrees = 90
plt.xticks(rotation=degrees)
plt.show()

import shap
testdf = df.copy()
classifier = RandomForestClassifier(random_state=1)
X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y,stratify=y, random_state=1)
classifier.fit(X_train, y_train)
plot_confusion_matrix(classifier, X_test, y_test)
plt.xticks(rotation=90)
plt.show()
y_pred = classifier.predict(X_test)
y_pred_proba = classifier.predict_proba(X_test)
#fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
#scores[i] = roc_auc_score(y_test, y_pred_proba)
score = roc_auc_score(y_test, y_pred_proba, multi_class='ovo')
explainer = shap.Explainer(classifier)
shap_values = explainer(X)
#shap.plots.waterfall(shap_values.base_values, shap_values[0], X[0])
shap.summary_plot(shap_values, X.values, plot_type="bar", class_names= class_names, feature_names = X.columns)
shapdf = pd.DataFrame(shap_values[:,:,1].values, index=df.index, columns=df.loc[:, df.columns != var].columns)
#plotdf = shapdf.join(df[var]).groupby('host_phenotype').mean()
plotdf = shapdf.join(df[var])
sns.heatmap(plotdf, xticklabels=True)
plotdf.set_index(var).sum(axis=1).plot.bar()

from sklearn.multiclass import OneVsRestClassifier
model = RandomForestClassifier(random_state=1)
ovr = OneVsRestClassifier(model)
ovr.fit(X_train, y_train)
y_pred = ovr.predict(X_test)
y_pred_proba = ovr.predict_proba(X_test)
score = roc_auc_score(y_test, y_pred_proba, multi_class='ovr')
plot_confusion_matrix(ovr, X_test, y_test)

# roc curve for classes
fpr = {}
tpr = {}
thresh ={}

n_class = len(df[var].unique())

for i in range(n_class):    
    fpr[i], tpr[i], thresh[i] = roc_curve(y_test, y_pred_proba[:,i], pos_label=i)
    
# plotting    
plt.plot(fpr[0], tpr[0], linestyle='--',color='orange', label='Class 0 vs Rest')
plt.plot(fpr[1], tpr[1], linestyle='--',color='green', label='Class 1 vs Rest')
plt.plot(fpr[2], tpr[2], linestyle='--',color='blue', label='Class 2 vs Rest')
plt.plot(fpr[3], tpr[3], linestyle='--',color='yellow', label='Class 3 vs Rest')
plt.title('Multiclass ROC curve')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive rate')
plt.legend(loc='best')
plt.savefig('Multiclass ROC',dpi=300); 
