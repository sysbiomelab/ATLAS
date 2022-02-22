#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
from sklearn.multiclass import OneVsRestClassifier
from scipy.spatial import distance
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
from sklearn.metrics import average_precision_score
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import f1_score
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap
import skbio
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import plot_precision_recall_curve

taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv').set_index('secondary_sample_accession')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

var = 'host_phenotype'
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
df = msptaxo.join(meta[var], how='inner').dropna()

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

scaler = StandardScaler()
df.loc[:, df.columns != var] = scaler.fit_transform(df.loc[:, df.columns != var])

Y = label_binarize(df[var], classes=df[var].unique())

testdf = df.copy()
classifier = RandomForestClassifier(random_state=1)
X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, Y, stratify=y, random_state=1)
classifier.fit(X_train, y_train)
y_score = classifier.predict_proba(X_test)
y_pred = classifier.predict(X_test)
#ConfusionMatrixDisplay.from_estimator(classifier, X_test, y_test)
#PrecisionRecallDisplay.from_estimator(classifier, X_test, y_test)

# curves
fig, (roc, pr) = plt.subplots(1,2)
results = pd.DataFrame(index=['precision', 'recall', 'fpr', 'tpr', 'acc', 'aucroc', 'pr','f1'])
for i in range(len(df.host_phenotype.unique())):
    precision, recall, _ = precision_recall_curve(
        y_test[:, i],
        y_score[:, i])
    fpr, tpr, _ = roc_curve(
        y_test[:, i],
        y_score[:, i])
    accuracy = accuracy_score(y_test[i], y_pred[i]
    rocscore = roc_auc_score(y_test[i], y_score[i])
    prscore = average_precision_score(y_test[i], y_score[i])
    f1score = f1_score(y_test[i], y_pred[i])
    results[i] = precision, recall, fpr, tpr, rocscore, prscore, f1score
    pr.plot(results.loc['recall',i], results.loc['precision',i], lw=1, label='class {}'.format(df.host_phenotype.unique()[i]))
    roc.plot(results.loc['fpr',i], results.loc['tpr',i], lw=1, label='{}'.format(df.host_phenotype.unique()[i]))
results.columns = df.host_phenotype.unique()
results.T[['aucroc', 'pr', 'f1']]

plt.xlabel("recall")
plt.ylabel("precision")
plt.legend(loc="best")
plt.title("precision vs. recall curve")
plt.tight_layout()
#plt.xticks(rotation=90)
plt.show()
y_pred_proba = classifier.predict_proba(X_test)
print('AUCROC =', roc_auc_score(y_test[3], y_score[3]))
print('PR =', average_precision_score(y_test[3], y_score[3]))
#print('PR =', precision_recall_curve(y_test, y_pred_proba, multi_class='ovo'))
explainer = shap.Explainer(classifier)
shap_values = explainer(X)
#shap.plots.waterfall(shap_values.base_values, shap_values[0], X[0])
shap.summary_plot(shap_values, X.values, plot_type="bar", class_names= class_names, feature_names = X.columns)
shapdf = pd.DataFrame(shap_values[:,:,1].values, index=df.index, columns=df.loc[:, df.columns != var].columns)
#plotdf = shapdf.join(df[var]).groupby('host_phenotype').mean()
plotdf = shapdf.join(df[var])
sns.heatmap(plotdf, xticklabels=True)
plotdf.set_index(var).sum(axis=1).plot.bar()





model = RandomForestClassifier(random_state=1)
classifier = OneVsRestClassifier(model)
classifier.fit(X_train, y_train)
y_pred = ovr.predict(X_test)
y_pred_proba = ovr.predict_proba(X_test)
newscore = roc_auc_score(y_test, y_pred_proba, multi_class='ovr')
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
