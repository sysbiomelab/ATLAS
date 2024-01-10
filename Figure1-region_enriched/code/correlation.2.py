import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection


mspdf = pd.read_csv('ATLAS/wellnessMgsMat.csv', index_col=0).T
mspdf.index = mspdf.index.str.replace('X',"")
mspdf.index = mspdf.index.str.replace('v',"")

taxaType='genus'

gmsp_samples = mspdf.T
gmsp_taxonomy = pd.read_csv("taxo.csv", index_col=0)
mspdf = gmsp_samples.join(gmsp_taxonomy[taxaType], how='inner').groupby(taxaType).sum().T
mspdf.drop(list(mspdf.filter(regex = 'unclassified')), axis = 1, inplace = True)

metadf = pd.read_csv('ATLAS/data/scapis_wellness_metabo_v1234_Sheet1.csv')
metadf = metadf.iloc[:,4:].set_index('Name').T
proteindf1 = pd.read_csv('ATLAS/data/olink.visit1.new.normalization.11.panels.txt', sep='\t', index_col=0)
proteindf2 = pd.read_csv('ATLAS/data/olink.visit2.new.normalization.11.panels.txt', sep='\t', index_col=0)
proteindf3 = pd.read_csv('ATLAS/data/olink.visit3.new.normalization.11.panels.txt', sep='\t', index_col=0)
proteindf4 = pd.read_csv('ATLAS/data/olink.visit4.new.normalization.11.panels.txt', sep='\t', index_col=0)

p1 = proteindf1.T.add_suffix('_1').T
p2 = proteindf2.T.add_suffix('_2').T
p3 = proteindf3.T.add_suffix('_3').T
p4 = proteindf4.T.add_suffix('_4').T

mergedp = pd.concat([p1,p2,p3,p4])

df = mergedp.join(metadf.join(mspdf))

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










pdf = df.xs(mergedp.columns,axis=1)
mspd = df.xs(mspdf.columns,axis=1)
meta = df.xs(metadf.columns,axis=1)

pdf.to_csv("prot.csv")
mspd.to_csv("msp_genus.csv")
meta.to_csv("mtb.csv")
#pair1=pdf.join(mspd, how="inner")
#pair2=pdf.join(meta, how="inner")
#pair3=mspd.join(meta, how="inner")

c1, p1 = spearmanr(pdf, mspd, axis=0)
c2, p2 = spearmanr(pdf, meta, axis=0)
c3, p3 = spearmanr(mspd, meta, axis=0)
#c1,p1=spearmanr(pair1)
#c2,p2=spearmanr(pair2)
#c3,p3=spearmanr(pair3)

##correlations and pvalues msp and proteins
correlations1 = pd.DataFrame(
    c1,
    index=pdf.columns.append(mspd.columns),
    columns=pdf.columns.append(mspd.columns))
slicedCorrelations1 = correlations1.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]
pval1 = pd.DataFrame(
    p1,
    index=pdf.columns.append(mspd.columns),
    columns=pdf.columns.append(mspd.columns))
sliced_pval1=pval1.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

#correlations metabolites and proteins
correlations2 = pd.DataFrame(
    c2,
    index=pdf.columns.append(meta.columns),
    columns=pdf.columns.append(meta.columns))
slicedCorrelations2= correlations2.iloc[
        len(pdf.columns):,
        :len(pdf.columns)]

#correlations metabolites and msp
correlations3 = pd.DataFrame(
    c3,
    index=mspd.columns.append(meta.columns),
    columns=mspd.columns.append(meta.columns))
slicedCorrelations3 = correlations3.iloc[
        len(mspd.columns):,
        :len(mspd.columns)]

#apply stack operation 

#fdr correction on pairs 1, commented
edges1 = slicedCorrelations1.stack().reset_index()
pv1=sliced_pval1.stack().reset_index()
#edges1[1] = pv1[0]
#edges1.dropna(inplace=True)
#edges1[2]=fdrcorrection(edges1[1])[0]

edges2 = slicedCorrelations2.stack().reset_index()
edges3 = slicedCorrelations3.stack().reset_index()

edges1.columns = ['from','to','value']
edges2.columns = ['from','to','value']
edges3.columns = ['from','to','value']

#save csv top 20
edges1.sort_values('value').tail(20).to_csv('pdfMspd.csv',index=False)
edges2.sort_values('value').tail(20).to_csv('pdfMeta.csv',index=False)
edges3.sort_values('value').tail(20).to_csv('mspdMeta.csv',index=False)

#save list of all correlations pairs
edges1[edges1['value']>0.4].sort_values('value').to_csv('pdfMspd_c4.csv',index=False)
edges1.sort_values('value').to_csv('pdfMspd_l.csv',index=False)
edges2.sort_values('value').to_csv('pdfMeta_l.csv',index=False)
edges3.sort_values('value').to_csv('mspdMeta_l.csv',index=False)

