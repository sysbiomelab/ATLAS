import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection

#REALLY INTERESTING FILE

mspdf = pd.read_csv('../data/wellnessMgsMat.csv', index_col=0).T
mspdf.index = mspdf.index.str.replace('X',"")
mspdf.index = mspdf.index.str.replace('v',"")

taxaType='genus'

gmsp_samples = mspdf.T
gmsp_taxonomy = pd.read_csv("../data/gutTaxo.csv", index_col=0)
mspdf = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
mspdf.drop(list(mspdf.filter(regex = 'unclassified')), axis = 1, inplace = True)
metadf = pd.read_csv('../data/scapis_wellness_metabo_v1234_Sheet1.csv')
metadf = metadf.iloc[:,4:].set_index('Name').drop_duplicates().T

proteome = pd.read_csv('../data/proteomics.csv', index_col=0)

'''
correlationArray, uncorrectedPValueArray = spearmanr(df)
correlations = pd.DataFrame(
    correlationArray,
    index=df.columns,
    columns=df.columns)
uncorrectedPValues = pd.DataFrame(
    uncorrectedPValueArray,
    index=df.columns,
    columns=df.columns)

correlations = correlations.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)
uncorrectedPValues = uncorrectedPValues.dropna(thresh=3,axis=1).dropna(thresh=3, axis=0)

slicedCorrelations = correlations
slicedUncorrectedPValues = uncorrectedPValues

significantMatrix = pd.DataFrame(
    fdrcorrection(slicedUncorrectedPValues.values.flatten())[0].reshape(slicedUncorrectedPValues.shape),
    index = slicedUncorrectedPValues.index,
    columns = slicedUncorrectedPValues.columns)

edges = slicedCorrelations.stack().reset_index()
edges.columns = ['from','to','value']
edges.sort_values('value').tail(74).to_csv('newgo.csv',index=False)
edges.sort_values('value').tail(74)
edges.sort_values('value')

'''

pdf = df.xs(proteome.columns,axis=1)
mspd = df.xs(mspdf.columns,axis=1)
meta = df.xs(metadf.columns,axis=1)
#correlationArray, uncorrectedPValueArray = spearmanr(pdf, mspdf, axis=0)
#correlationArray, uncorrectedPValueArray = spearmanr(meta, proteome, axis=0)
c1, p1 = spearmanr(pdf, mspd, axis=0)
c2, p2 = spearmanr(pdf, meta, axis=0)
c3, p3 = spearmanr(mspd, meta, axis=0)

correlations1 = pd.DataFrame(
    c1,
    index=pdf.columns.append(mspd.columns),
    columns=pdf.columns.append(mspd.columns))
slicedCorrelations1 = correlations1.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

correlations2 = pd.DataFrame(
    c2,
    index=pdf.columns.append(meta.columns),
    columns=pdf.columns.append(meta.columns))
slicedCorrelations2 = correlations2.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

correlations3 = pd.DataFrame(
    c3,
    index=mspd.columns.append(meta.columns),
    columns=mspd.columns.append(meta.columns))
slicedCorrelations3 = correlations3.iloc[
        len(mspd.columns):,
        :len(mspd.columns)]

edges1 = slicedCorrelations1.stack().reset_index()
edges2 = slicedCorrelations2.stack().reset_index()
edges3 = slicedCorrelations3.stack().reset_index()
#edges = pd.concat([edges1,edges2,edges3])
#
#edges12 = pd.concat([edges1, edges2])
#edges13 = pd.concat([edges1, edges3])
#edges23 = pd.concat([edges2, edges3])
edges1.columns = ['from','to','value']
edges2.columns = ['from','to','value']
edges3.columns = ['from','to','value']

'''
edges1.drop_duplicates().sort_values('value').tail(20).to_csv('pdfMspd.csv',index=False)
edges2.drop_duplicates().sort_values('value').tail(20).to_csv('pdfMeta.csv',index=False)
edges3.drop_duplicates().sort_values('value').tail(20).to_csv('mspdMeta.csv',index=False)
'''

# for the pvalue number stuff
correlations1 = pd.DataFrame(
    c1,
    index=pdf.columns.append(mspd.columns),
    columns=pdf.columns.append(mspd.columns))
slicedCorrelations1 = correlations1.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

correlations2 = pd.DataFrame(
    c2,
    index=pdf.columns.append(meta.columns),
    columns=pdf.columns.append(meta.columns))
slicedCorrelations2 = correlations2.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

correlations3 = pd.DataFrame(
    c3,
    index=mspd.columns.append(meta.columns),
    columns=mspd.columns.append(meta.columns))
slicedCorrelations3 = correlations3.iloc[
        len(mspd.columns):,
        :len(mspd.columns)]



significantMatrix1 = pd.DataFrame(
    fdrcorrection(slicedCorrelations1.values.flatten())[0].reshape(slicedCorrelations1.shape),
    index = slicedCorrelations1.index,
    columns = slicedCorrelations1.columns)

significantMatrix2 = pd.DataFrame(
    fdrcorrection(slicedCorrelations2.values.flatten())[0].reshape(slicedCorrelations2.shape),
    index = slicedCorrelations2.index,
    columns = slicedCorrelations2.columns)

significantMatrix3 = pd.DataFrame(
    fdrcorrection(slicedCorrelations3.values.flatten())[0].reshape(slicedCorrelations3.shape),
    index = slicedCorrelations3.index,
    columns = slicedCorrelations3.columns)


edges1 = slicedCorrelations1.stack().reset_index()
edges2 = slicedCorrelations2.stack().reset_index()
edges3 = slicedCorrelations3.stack().reset_index()

sigdf = {}
sigdf[0] = {'from':'proteome', 'to':'microbiome', 'spearman':(slicedCorrelations1 > 0.3).sum().sum()}
sigdf[1] = {'from':'proteome', 'to':'metabolome', 'spearman':(slicedCorrelations2 > 0.3).sum().sum()}
sigdf[2] = {'from':'microbiome', 'to':'metabolome', 'spearman':(slicedCorrelations3 > 0.3).sum().sum()}
sigdf1 = pd.DataFrame(data = sigdf).T

sigdf1.to_csv('total_spearman.csv')

'''
