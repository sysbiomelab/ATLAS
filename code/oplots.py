#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Theo Portlock
'''
import functions
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

taxo = pd.read_csv('../data/gutTaxo.csv', index_col=0)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0)
meta = pd.read_csv('../metaAtlasFinish/finaltheo/samples.tsv')
meta.host_phenotype = meta.host_phenotype.fillna('Healthy')

taxoMsp = msp.join(taxo['genus']).groupby('genus').sum()
unclass = taxoMsp[taxoMsp.index.str.contains("unclassified")].sum()
taxoMsp = taxoMsp[~taxoMsp.index.str.contains("unclassified")].T
taxoMsp[-1] = unclass
df = taxoMsp.div(taxoMsp.sum(axis=1), axis=0)
df.to_csv('../results/genusAbund.csv')
