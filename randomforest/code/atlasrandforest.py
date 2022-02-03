#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as PCA
from sklearn.preprocessing import StandardScaler
import skbio
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial import distance
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report)
import skbio

taxo = pd.read_csv('../../../../data/FMT/downstream_data/taxo.csv', index_col=0)
msp = pd.read_csv('../../../oldatlas/data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../../../oldatlas/data/unique_metadata.csv')

diseasegroup = pd.read_csv('../data/cohort_colors.ob.txt')
diseasegroup.columns = range(3)
diseasegroup[0] = diseasegroup[0].str.replace('.*:.*:', '')
diseasegroup.set_index(0, inplace=True)

meta = meta.set_index('study_accession').join(diseasegroup).set_index('secondary_sample_accession')

msptaxo = msp.join(taxo['species']).groupby('species').sum().T

var = 2
'''
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
pd.DataFrame(columns=df.columns, data=LabelEncoder().fit_transform(df.values.flatten()).reshape(df.shape))
samples_gutmeta['Type'] = LabelEncoder().fit_transform(samples_gutmeta.Type)
ohe = OneHotEncoder()
meta_ohe = ohe.fit_transform(df)
clf = RandomForestClassifier()
'''

df = msptaxo.join(meta[var], how='inner').dropna()

feature_imp = pd.DataFrame()
scores = pd.Series()
for i in df[var].unique():
    testdf = df.copy()
    testdf.loc[testdf[var] != i, var] = 'other'
    classifier = RandomForestClassifier()
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    classifier.fit(X_train, y_train)
    feature_imp[i] = pd.Series(classifier.feature_importances_,index=X.columns)
    y_pred = classifier.predict(X_test)
    #classifier.score(X,y)
    #scores[i] = accuracy_score(y_test.values, y_pred)
    y_pred_proba = classifier.predict_proba(X_test)[:,1]
    scores[i] = roc_auc_score(y_test, y_pred_proba)
    fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba, pos_label=i)
    plt.plot(fpr,tpr,label=str(i)+" AUC="+str(scores[i]))
    #roc_auc_score(y_test, y_pred)

'''
#define metrics
y_pred_proba = log_regression.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)

#create ROC curve
plt.plot(fpr,tpr,label="AUC="+str(scores[i]))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.show()
'''

sscores=scores.copy()
sfeature_imp = feature_imp.copy()

sns.heatmap(sfeature_imp.loc[(sfeature_imp > 0.012).T.any()], yticklabels=True)

from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

plot_confusion_matrix(classifier, X_test, y_test)
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
