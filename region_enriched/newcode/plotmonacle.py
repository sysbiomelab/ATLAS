#!/usr/bin/env python3
import functions
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#meta = pd.read_csv('../../data/unique_metadata.csv').set_index('secondary_sample_accession')
data = pd.read_csv('../../data/plotdata.csv').set_index('Unnamed: 0')
data.columns = ['Component 1', 'Component 2'] 
meta = pd.read_csv('../../data/basicMetadata.Theo.SL.2022.08.18.csv', index_col=0)
#data = pd.read_csv('../results/monocleout.csv', index_col=0)
data.columns = ['Component 1', 'Component 2'] 
#enterotype = pd.read_csv('../data/enterotype.csv')

joined = data.join(meta)

#joined['Region'] = 'European'
#joined.loc[joined.westernised == 'NW', 'Region'] = 'Non-westernized'
#joined.loc[joined.country == 'JPN', 'Region'] = 'China/Japan/US'
#joined.loc[joined.country == 'CHN', 'Region'] = 'China/Japan/US'
#joined.loc[joined.country == 'USA', 'Region'] = 'China/Japan/US'
sns.scatterplot(data=joined,
        x=data.columns[0],
        y=data.columns[1],
        edgecolor=None,
        hue='enteroType',
        #palette='coolwarm',
        #hue_norm=(0,50),
        lw=0,
        s=5)
#sns.scatterplot(data=joined, x=data.columns[0], y=data.columns[1], edgecolor=None, hue='country', alpha=0.5)
#sns.scatterplot(data=joined, x=data.columns[0], y=data.columns[1], edgecolor=None, hue='Region', lw=0)

#functions.PCplot(joined,'westernised' , x='Component 1', y='Component 2')

import functions
plt.tight_layout()
plt.rcParams["figure.figsize"] = (7,6)
plt.savefig('../results/enteromonocle.pdf')

ax = sns.scatterplot(data=joined.groupby('country').mean(),
        x=data.columns[0],
        y=data.columns[1],
        edgecolor=None,
        hue=joined.groupby('country').first()['westernised'],
        alpha=0.5,
        size=joined.groupby('country').nunique()['sample_accession'],
        sizes=(0,2000))

for i in joined.groupby('country').mean().index.to_list():
    ax.annotate(i, xy=joined.groupby('country').mean().loc[i,['Component 1','Component 2']])

'''
df= joined.query('country != "SWE" and country != "DNK"')
sns.scatterplot(data=df,
        x=data.loc[df.index, 'Component 1'],
        y=data.loc[df.index, 'Component 2'],
        edgecolor=None,
        hue='gender',
        #palette='coolwarm',
        #hue_norm=(0,50),
        lw=0,
        s=5)
plt.tight_layout()
plt.rcParams["figure.figsize"] = (7,6)
plt.savefig('../results/filtgendermonocleplot.pdf')
