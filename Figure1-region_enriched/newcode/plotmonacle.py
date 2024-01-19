#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

meta = pd.read_csv('../data/unique_metadata.csv').set_index('secondary_sample_accession')
data = pd.read_csv('../data/plotdata.csv').set_index('Unnamed: 0')
data.columns = ['Component 1', 'Component 2'] 

joined = data.join(meta)

joined['Region'] = 'European'

joined.loc[joined.westernised == 'NW', 'Region'] = 'Non-westernized'

joined.loc[joined.country == 'JPN', 'Region'] = 'China/Japan/US'
joined.loc[joined.country == 'CHN', 'Region'] = 'China/Japan/US'
joined.loc[joined.country == 'USA', 'Region'] = 'China/Japan/US'

sns.scatterplot(data=joined, x=data.columns[0], y=data.columns[1], edgecolor=None, hue='Region', alpha=0.5)
plt.savefig('../results/Figure_1e.pdf')
