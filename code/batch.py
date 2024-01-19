#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Analysis of Batch effect on HGMA datasets
Theo Portlock
'''
import functions
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns

# Load species abundances (msp), metadata (meta), and taxonomical information (taxo)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0).T
meta = pd.read_csv('../metaAtlasFinish/final/basicMetadata.Theo.SL.2022.06.12.txt', sep='\t', index_col=0)
taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)

# Calculate genus level resolution
taxoMsp = msp.T.join(taxo['genus'], how='inner').groupby('genus').sum().T

# Compute beta diversity (Bray-Curtis)
beta = functions.beta(taxoMsp)


# Compute stats
from skbio.stats.distance import permanova
from skbio.stats.distance import anosim as permanova
from skbio import DistanceMatrix
## statistically different
permanova(DistanceMatrix(beta, beta.index), beta.join(meta)['BioProject'].astype('str').values)
## after shuffling group, there is no significance
newgroup = beta.join(meta)['BioProject'].to_frame()
newgroup['nbp'] = newgroup.sample(frac=1).values
permanova(DistanceMatrix(beta, beta.index), newgroup['nbp'].astype('str').values)


# Compute NMDS
nmds = functions.NMDS(beta)

# add batch information
batchdf = nmds.join(meta['BioProject'])

# Setup plotting parameters
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 12
matplotlib.rcParams['figure.figsize'] = 11.7,8.27

# Plot and save
sns.scatterplot(
    data=batchdf,
    x='NMDS1',
    y='NMDS2',
    hue='BioProject',
    size=0.1,
    linewidth=0,
)
plt.legend(
    title='BioProject',
    bbox_to_anchor=(1.001, 1),
    loc='upper left',
    fontsize='small',
    ncol=3
)
plt.tight_layout()
plt.savefig('../results/batchNMDS.svg')
plt.savefig('../results/batchNMDS.pdf')
plt.show()
