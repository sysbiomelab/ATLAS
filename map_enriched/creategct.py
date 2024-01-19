import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

genebag = pd.read_csv('../../../FMT/gutdataverse_files/IGC2.1990MSPs.tsv', sep='\t', index_col=0)
country_codes = pd.read_csv('../data/countrycodes.tsv', sep='\t',  index_col=1)
msp = pd.read_csv('../data/vect_atlas.csv', index_col=0).T

msp.iloc[0].to_frame().join(genebag["gene_name"]).groupby("gene_name").sum()

final = pd.DataFrame(index=genebag['gene_name'])
#final = pd.DataFrame()
for i in msp.columns: 
    final = final.join(msp[i].to_frame().join(genebag["gene_name"]).drop_duplicates().set_index("gene_name"))

card = pd.read_csv('../data/hs_10_4_igc2.CARD.tsv', sep='\t')
card_name = pd.read_csv('../data/aro.tsv', sep='\t', index_col=0)
card_name.index = card_name.index.str.replace('ARO:','')
card["gene_name"] = card.ORF_ID.str.split(expand=True)[0]
card["gene_name"] = card["gene_name"].str[0:-2]
card = card.set_index('gene_name')

