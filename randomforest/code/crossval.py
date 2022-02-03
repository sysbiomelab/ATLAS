#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Random forest analysis of metadata
Theo Portlock
'''
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial import distance
import skbio
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

var = 'host_phenotype'
newmeta = meta.set_index('secondary_sample_accession')
msptaxo = msp.join(taxo['species']).groupby('species').sum().T
df = msptaxo.join(newmeta[var], how='inner').dropna()

df = df.loc[df[var] != 'excluded controls'] 
df = df.loc[df[var] != 'Healthy'] 
df = df.loc[df[var] != 'Obese'] 
df = df.loc[df[var] != 'acute diarrhea'] 
df = df.loc[df[var] != 'Overweight'] 
df = df.loc[df[var] != 'atherosclerosis'] 
df = df.loc[df[var] != 'large adenoma'] 
df = df.loc[df[var] != 'premature_734g'] 
df = df.loc[df[var] != 'premature_2847g'] 
df = df.loc[df[var] != 'Underweight'] 
df = df.loc[df[var] != 'GDM'] 
df = df.loc[df[var] != 'UC'] 
df = df.loc[df[var] != 'BD'] 
df = df.loc[df[var] != 'ob'] 
df = df.loc[df[var] != 'control'] 

'''
feature_imp = pd.DataFrame()
scorecurve = pd.DataFrame(columns=['scores', 'curves'])
scores = pd.Series()
curves = pd.Series()
for i in df[var].unique():
    testdf = df.copy()
    testdf.loc[testdf[var] != i, var] = 'other'
    nondiseasedf=testdf.loc[testdf[var] == 'other'].sample(1000)
    diseaseddf = testdf.loc[testdf[var] == i]
    testdf = pd.concat([nondiseasedf, diseaseddf])
    #classifier = RandomForestClassifier()
    #classifier = RandomForestClassifier(class_weight={'other': 1, i: 50})
    classifier = RandomForestClassifier(class_weight='balanced')
    X = testdf.drop(var, axis=1)
    y = testdf.xs(var, axis=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, stratify=y)
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

def evaluate(model, test_features, test_labels):
    predictions = model.score(test_features)
    return accuracy
'''

i = 'CRC'
reg = RandomForestClassifier()
testdf = df.copy()
testdf.loc[testdf[var] != i, var] = 'other'

groups = testdf.join(meta.set_index('secondary_sample_accession')['country']).country
X = testdf.drop(var, axis=1)
y = testdf.xs(var, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg.fit(X_train, y_train)
from sklearn.model_selection import cross_val_score
scores = cross_val_score(reg, X, y, cv=5, groups=groups)
np.sqrt(np.mean(np.absolute(scores)))

final = pd.Series()
final['Gut Species'] = runmodel(mspdf)
final['Oral Species'] = runmodel(omspdf)
final['Secondary Metabolites'] = runmodel(smdf)
final['Virulence Factors'] = runmodel(padf)
final['Protein Family'] = runmodel(pfdf)
final['Antimircobial Resistance'] = runmodel(cadf)
final['Cazymes'] = runmodel(czdf)
final.sort_values(ascending=False, inplace=True)
sns.barplot(x=final, y=final.index)
plt.xlabel('Model Accuracy (%)')
plt.xlim(0,100)
plt.ylabel('Dataset used for MELD prediction')
plt.show()
plt.tight_layout()
plt.savefig('results/OvG_randomfores_MELD.pdf')

'''
y_pred = clf.predict(X_test)

#This is jaccard
print("Accuracy RFC:",metrics.accuracy_score(y_test, y_pred))
print("Accuracy MSE:",metrics.mean_squared_error(y_test, y_pred))
print(metrics.confusion_matrix(y_test, y_pred))

evaluate(clf, X_test, y_test)

from sklearn.model_selection import GridSearchCV
# Create the parameter grid based on the results of random search
param_grid = {
    'bootstrap': [True],
    'max_depth': [80, 90, 100, 110],
    'max_features': [2, 3],
    'min_samples_leaf': [3, 4, 5],
    'min_samples_split': [8, 10, 12],
    'n_estimators': [100, 200, 300, 1000]
}
# Create a based model
rf = RandomForestRegressor()
# Instantiate the grid search model
grid_search = GridSearchCV(estimator = rf, param_grid = param_grid,
                          cv = 3, n_jobs = -1, verbose = 2)

grid_search.fit(X_train, y_train)
best_grid = grid_search.best_estimator_
grid_accuracy = evaluate(best_grid, X_test, y_test)
feature_imp = pd.Series(best_grid.feature_importances_,index=X.columns).sort_values(ascending=False)
sns.barplot(x=feature_imp, y=feature_imp.index)
plt.xlabel('Feature Importance Score')
#plt.ylabel('Metadata')
plt.savefig("meta_MELD_prediction.pdf")
'''

reg = RandomForestRegressor()
X = czdf.drop(variable, axis=1)
y = czdf.xs(variable, axis=1)
X_train, X_test, y_train, y_test = train_test_split(X, y)
reg.fit(X_train, y_train)
feature_imp = pd.Series(reg.feature_importances_,index=X.columns).sort_values(ascending=False)
sns.barplot(x=feature_imp, y=feature_imp.index)


'''
Ar_dist = distance.squareform(distance.pdist(czdf.drop('MELD', axis=1), metric="braycurtis"))
DM_dist = skbio.stats.distance.DistanceMatrix(Ar_dist)
PCoA = skbio.stats.ordination.pcoa(DM_dist)
plot = czdf.join(PCoA.samples[['PC1','PC2']].set_index(cadf.index))
sns.scatterplot(data=plot, x='PC1', y='PC2', hue='MELD', palette='coolwarm')

scaledData = StandardScaler().fit_transform(czdf.drop('MELD', axis=1))
pca = PCA(n_components=2)
result = pca.fit_transform(scaledData)
print(result)
czdf[['PC1', 'PC2']] = result
sns.scatterplot(data=czdf, x='PC1', y='PC2', hue='MELD', palette='coolwarm')

correlation = mspdf.corr(method='spearman')
columns = correlation.nlargest(10, 'MELD').index
correlation_map = np.corrcoef(mspdf[columns].values.T)
sns.set(font_scale=1.0)
heatmap = sns.heatmap(correlation_map, cbar=True, annot=True, square=True, fmt='.2f', yticklabels=columns.values, xticklabels=columns.values)
plt.show()

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
pipelines = []
pipelines.append(('ScaledLR', Pipeline([('Scaler', StandardScaler()),('LR',LinearRegression())])))
pipelines.append(('ScaledLASSO', Pipeline([('Scaler', StandardScaler()),('LASSO', Lasso())])))
pipelines.append(('ScaledEN', Pipeline([('Scaler', StandardScaler()),('EN', ElasticNet())])))
pipelines.append(('ScaledKNN', Pipeline([('Scaler', StandardScaler()),('KNN', KNeighborsRegressor())])))
pipelines.append(('ScaledCART', Pipeline([('Scaler', StandardScaler()),('CART', DecisionTreeRegressor())])))
pipelines.append(('ScaledGBM', Pipeline([('Scaler', StandardScaler()),('GBM', GradientBoostingRegressor())])))
pipelines.append(('ScaledRFR', Pipeline([('Scaler', StandardScaler()),('RFR', RandomForestRegressor())])))
results = []
names = []
for name, model in pipelines:
    kfold = KFold(n_splits=10)
    cv_results = cross_val_score(model, X_train, y_train, cv=kfold, scoring='neg_mean_squared_error')
    results.append(cv_results)
    names.append(name)
    msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
    print(msg)
scaler = StandardScaler().fit(X_train)
rescaledX = scaler.transform(X_train)
param_grid = dict(n_neighbors=np.array([2,3,4,5,6]))
model = KNeighborsRegressor()
kfold = KFold(n_splits=10)
grid = GridSearchCV(estimator=model, param_grid=param_grid, scoring='neg_mean_squared_error', cv=kfold)
grid_result = grid.fit(rescaledX, y_train)
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
from sklearn.metrics import mean_squared_error
scaler = StandardScaler().fit(X_train)
rescaled_X_train = scaler.transform(X_train)
model = KNeighborsRegressor(n_neighbors=5)
model.fit(rescaled_X_train, y_train)
# transform the validation dataset
rescaled_X_test = scaler.transform(X_test)
predictions = model.predict(rescaled_X_test)
print (mean_squared_error(y_test, predictions))
compare = pd.DataFrame({'Prediction': predictions, 'Test Data' : y_test})
mcompare = compare.reset_index().melt(id_vars=['index'])
