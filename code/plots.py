#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
import functions
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
import pickle
import seaborn as sns
import shap
import skbio

taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0).T
allmeta = pd.read_csv('../metaAtlasFinish/final/basicMetadata.Theo.SL.2022.06.12.txt', sep='\t',  index_col=0)
hgmi = pd.read_csv('../metaAtlasFinish/finaltheo/hgmiedges.csv', index_col=1)
meta = hgmi.join(allmeta).reset_index().set_index('hgmi_id')

var = 'Disease'

shaps = pd.DataFrame()
#scores = pd.Series()
hgmaid = meta.index.unique()[0]
for hgmaid in meta.index.unique():
    df = msp.join(meta.loc[hgmaid].set_index('sample.ID')[var],how='inner').set_index(var)
    df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns).reset_index()
    model = RandomForestClassifier(
        bootstrap=True,
        ccp_alpha=0.0,
        class_weight='balanced_subsample',
        criterion='gini',
        max_depth=7,
        max_features='sqrt',
        max_leaf_nodes=None,
        max_samples=None,
        min_impurity_decrease=0.0005,
        min_samples_leaf=3,
        min_samples_split=10,
        min_weight_fraction_leaf=0.0,
        n_estimators=300,
        n_jobs=-1,
        oob_score=False,
        random_state=2,
        verbose=0,
        warm_start=False)
    X = df.drop(var, axis=1)
    y = pd.get_dummies(df.xs(var, axis=1))['Healthy']
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,random_state = 1)
    model.fit(X_train, y_train)
    #functions.crossval(model, X, y)
    y_pred = model.predict(X_test)
    explainer = shap.TreeExplainer(model)
    shaps_values = explainer(X)
    meanabsshap = pd.Series(np.abs(shaps_values.values).mean(axis=0)[:, 0], index=X.columns)
    corrs = [spearmanr(shaps_values.values[:, :, 0][:, x], X.iloc[:,x])[0] for x in range(len(X.columns))]
    final = meanabsshap * np.sign(corrs)
    final.fillna(0, inplace=True)
    shaps[hgmaid] = final
    #scores[hgmaid] = classification_report(y_true=y_test, y_pred=y_pred) + '\n\nAUCROC = ' + str(roc_auc_score(y_test, model.predict_proba(X_test)[:,1]))
    #with open(f"../results/{hgmaid}.txt", "w") as f: f.write(classification_report(y_true=y_test, y_pred=y_pred) + '\n\nAUCROC = ' + str(roc_auc_score(y_test, model.predict_proba(X_test)[:,1])))

#shaps = pd.read_csv('../results/shaps.csv', index_col=0)
shaps.columns = shaps.columns.astype(int)
#final.to_frame().join(taxo).sort_values(0).set_index("species").tail(20)[0].plot.barh()
#shaps.mean(axis=1).to_frame().join(taxo).sort_values(0).set_index("species").tail(20)[0].plot.barh()
a = shaps.join(taxo['species']).groupby('species').sum().T.join(meta.loc[meta.Disease != "Healthy"].groupby(level=0).first()['Disease']).set_index('Disease')
#means = means.loc[means[0].gt(0)]shaps.mean(axis=1).to_frame().join(taxo).sort_values(0).set_index("species")
#gcor = msp.T.join(taxo["species"]).groupby("species").sum().loc[ means.loc[means[0].gt(0)].index ].T.corr()
#lcor = msp.T.join(taxo["species"]).groupby("species").sum().loc[ means.loc[means[0].lt(0)].index ].T.corr()
#plt.tight_layout(); plt.show()

nc = means.join(pd.Series(clust, name="clust"))[["clust", 0]]

plotdf = shaps[(np.abs(shaps) > 0.01).any(axis=1)]

plotdf = plotdf.join(taxo['species']).set_index('species').T.join(meta.loc[meta.Disease != "Healthy"].groupby(level=0).first()['Disease']).set_index('Disease')
functions.clustermap(plotdf)


corr = msp.T.join(taxo["species"]).set_index("species").T.corr(method='spearman')
G = functions.network(corr, thresh = 0.5)
clust = functions.cluster(G)
functions.annotateplot(G, clust)
functions.clusterplot(G)
plt.show()

a = (
    shaps.join(taxo["species"])
    .groupby("species")
    .sum()
    .T.join(meta.loc[meta.Disease != "Healthy"].groupby(level=0).first()["Disease"])
    .set_index("Disease")
)
b = pd.concat([clust, a.loc['LC']], axis=1)
sns.stripplot(data=b, x="Group", y="LC", size=3, color="black")
sns.boxplot(data=b, x="Group", y="LC", showfliers=False, boxprops=dict(alpha=.5))
plt.show()


# Cluster correlation
b.dropna(inplace=True)
c = b.groupby('Group').mean()
renam = msp.T.join(taxo["species"]).set_index("species")
cluscor = renam.join(b['Group']).dropna().groupby('Group').mean().T.corr(method='spearman')

functions.clustermap(cluscor)
