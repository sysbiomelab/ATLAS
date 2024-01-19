import pandas as pd
import numpy as np
from scipy import stats
import math
import argparse
from collections import defaultdict, Counter
'''
disease_studies= ["Atherosclerosis_SW_id1","CVD_CN_id10", "GDM_CN_id12", "T2D_CN_id2", "T1D_FI_id21", "Cirrhosis_UK_id27",
                 "CFS_US_id34", "NAFLD_IT_id40", "NAFLD_ES_id40","CRC_JP_id45","CRC_IT_id46", "Behcet_CN_id52",
                 "CD_CN_id6", "Ankylosing_CN_id9","CRC_US_id11", "T2D_SW_id14", "IGT_SW_id14", "Obesity_DK_id16",
                 "IBD_ES_id17", "Cirrhosis_CN_id20", "Obesity_DK_id25", "NSCLC_FR_id26", "Melanoma_US_id3", "RA_CN_id31",
                 "T1D_LU_id35","RCC_FR_id41", "CRC_DE_id44", "T2D_ES_id53", "NAFLD_US_id7", "PD_DE_id8"]

exMatched = ["Melanoma_US_id3", "NAFLD_US_id7", "Cirrhosis_UK_id27", "NAFLD_ES_id40",
              "NAFLD_IT_id40", "T2D_ES_id53", "T1D_FI_id21"]
'''
def effectSize(caseVec, contVec, imputeVal=float(10) ):
    try:
        wtest = stats.mannwhitneyu(caseVec, contVec, alternative = "greater")
        N = len(caseVec) + len(contVec)
        Z = stats.norm.ppf(1-wtest.pvalue)
        es = (Z)/math.sqrt(N)
        if (es == float('inf')): es = imputeVal
        if (es == -float('inf')): es = - imputeVal
    except:
        #es=np.nan
        es=float(0)
    if es < 0 : es = float(0)
    #if es > 1 : es = float(1)
    return(es)
def effectSizeDown (caseVec, contVec, imputeVal=float(10)):
    try:
        wtest = stats.mannwhitneyu(caseVec, contVec, alternative = "less")
        N = len(caseVec) + len(contVec)
        Z = stats.norm.ppf(1-wtest.pvalue)
        es = (Z)/math.sqrt(N)
        if (es == float('inf')): es = imputeVal
        if (es == -float('inf')): es = - imputeVal
    except:
        #es=np.nan
        es=float(0)
    if es < 0 : es = float(0)
    #if es > 1 : es = float(1)
    return(es)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, metavar='FILE')
    parser.add_argument("--msptable", required=True, metavar='FILE')
    parser.add_argument("--cohorts", required=True, metavar='FILE')
    args=parser.parse_args()

    md=pd.read_csv(args.metadata, sep="\t")
    msp_tb=pd.read_csv(args.msptable, index_col=0, header=0)
    cohorts=pd.read_csv(args.msptable, header=None)
    cohorts=cohorts.rename(columns={0:'cohort', 1:'cohort_int',2:'category'})
    msp_tb=msp_tb*1e8

    #get list of study_id-disease
    #use only first samples (t1), not babies and with read coverage depth of 10M
    md=md.loc[ (md["timepoint"]=="t1") & (md["age_baby"].isna()) & (md["depth_10M"]=="ok") ] #only get samples at first timepoint and not baby age
    ##modify denmark ob in Phenotype with P in Health status
    md.loc[ md.host_phenotype=="ob",'health_status']="P"
    #get Disease samples (P) and build cohort list, same study, same country and same disease
    md_p=md.loc[ md["health_status"] =="P"  ] #Patho phenotye
    cohort_list=md_p[["host_phenotype","country","study_accession"]].agg(":".join, axis=1).unique()
    cohort_list.sort()

    # for each cohort in cohort list get sample_id list where
    es_ctrl1=[]
    es_ctrl2=[]
    #es_ctrl3=[]
    #cohort by cohort estimate effect size in the MSPs
    for i,cohort in enumerate(cohort_list):
        chr=i+1
        print(cohort,i)
        pheno=cohort.split(":")[0]
        country=cohort.split(":")[1]
        study=cohort.split(":")[2]

        #get samples_id from sample cohort (study and phenotype)
        md_cohort=md.loc[ (md["health_status"]=="P") & (md["study_accession"]==study) & (md["host_phenotype"]==pheno)]
        samples_cohort=list(md_cohort["secondary_sample_accession"].unique())
        #get control1, "matched" healthy samples, same study, same country
        md_ctrl=md.loc[ (md["health_status"]=="H") & (md["study_accession"]==study) & (md["country"]==country)]
        samples_ctrl1=list(md_ctrl["secondary_sample_accession"].unique())
        #get control2, "country" healthy samples, samecountry
        md_ctrl=md.loc[ (md["health_status"]=="H") & (md["country"]==country)]
        samples_ctrl2=list(md_ctrl["secondary_sample_accession"].unique())
        #get control3, "et" healthy samples in countries withitn the same set
        #print(len(samples_cohort),len(samples_ctrl1))
        if i==19:
            print(cohort)
            print(samples_ctrl1)
            print(samples_ctrl2)
        '''
        #healthy samples same group of countries. Group countries by enterotype
        if country in et1:
            md_ctrl=md.loc[ (md["health_status"] == "H") & (md["country"] in et1) ]
            samples_ctrl3=md_ctrl["sample_accession"]
        elif country in et2:
            md_ctrl=md.loc[ (md["health_status"] == "H") & (md["country"] in et2) ]
            samples_ctrl3=md_ctrl["samlpe_accession"]
        '''

        #estimate effect size for each MSP in cohort

        for j,msp in enumerate(msp_tb.index):
            bp=j+1
            #effect size control 1 mathced control samples
            msp_samples_ref=list(msp_tb.loc[msp,samples_cohort])
            msp_samples_ctrl1=list(msp_tb.loc[msp,samples_ctrl1])
            #estimate effect size only if more than 3 samples  and sum > 0
            if (len(msp_samples_ref) > 3 and len(msp_samples_ctrl1) > 3):
                if (sum(msp_samples_ref)>0) and (sum(msp_samples_ctrl1)>0):
                    es_upp=effectSize(msp_samples_ref, msp_samples_ctrl1)
                    es_down=effectSizeDown(msp_samples_ref, msp_samples_ctrl1)
                else:
                    es_upp=np.nan
                    es_down=np.nan
                #CHR BP P SNP for qqman R package
                #CHR cohort int id, BP msp int id,  P float effect size, SNP msp id
                es_row=[msp,bp,es_upp,es_down,cohort,chr,"matched"]
                es_ctrl1.append(es_row)
            #effect size control 2 healthy samples same country control
            msp_samples_ctrl2=np.array(msp_tb.loc[msp,samples_ctrl2])
            if (len(msp_samples_ref) > 3 and len(msp_samples_ctrl2) > 3):
                if sum(msp_samples_ref)>0 and sum(msp_samples_ctrl2)>0:
                    es_upp=effectSize(msp_samples_ref, msp_samples_ctrl2)
                    es_down=effectSizeDown(msp_samples_ref,msp_samples_ctrl2)
                else:
                    es_upp=np.nan
                    es_down=np.nan
                es_row=[msp,bp,es_upp,es_down,cohort,chr,"country"]
                es_ctrl2.append(es_row)
            #effect size control 2 healthy samples same country control
            '''
            #effect size control 1 mathced control samples
            list_ctrl=list(msp_tb.iloc[msp,samples_ctrl3])
            es_upp=effectSize(list_ref, list_ctrl1)
            es_row=[msp,es_upp,cohort,"ctrl3"]
            es_ctrl3.append(esrow)
            '''
    #transform results (msp effect size in cohort compared against control samples)
    #into dataframe for plotting into manhattan plots, ES=Effect Size
    #ctrl1 healthy, same study, same country
    df_es_ctrl1=pd.DataFrame(es_ctrl1, columns=["msp","msp_int","ES_up","ES_down","cohort","cohort_int","control"])
    df_es_ctrl1.to_csv("Effect_size.Matched_ctrl.tsv", sep="\t")
    #ctrl3 healthy same country
    df_es_ctrl2=pd.DataFrame(es_ctrl2, columns=["msp","msp_int","ES_up","ES_down","cohort","cohort_int","control"])
    df_qqman=df_es_ctrl2
    df_qqman.loc[df_qqman["ES_up"]>1, "ES_up"]=1
    df_qqman.loc[df_qqman["ES_down"]>1, "ES_down"]=1
    df_qqman.to_csv("Effect_size.Country_ctrl.tsv", sep="\t")

    ## filter up above 0.3, and down above 0.3
    ## count by msp mumber of samples found above, number of samples found below
    ES_sig_up=df_es_ctrl2[df_es_ctrl2["ES_up"]>0.3]
    #msp_ESup=ES_sig_up.groupby("msp")["ES_up"].sum() #.reset_index()
    ES_sig_dw=df_es_ctrl2[df_es_ctrl2["ES_down"]>0.3]
    #msp_ESdw=ES_sig_dw.groupby("msp")["ES_down"].sum()
    msp_up=Counter(ES_sig_up["msp"])
    msp_dw=Counter(ES_sig_dw["msp"])
    msp_list=set(msp_up.keys())
    msp_list=msp_list.union(set(msp_dw.keys()))


    ## return table msp cohort up cohort down, enriched and depleated
    volcano=[]
    for msp in msp_list:
        sm=msp_up[msp]+msp_dw[msp]
        dif=msp_up[msp]-msp_dw[msp]
        line=[msp,sm,dif]
        volcano.append(line)
    df_volc=pd.DataFrame(volcano, columns=["msp",'sum','dif'])
    df_volc.to_csv("Volcano_plot.Country.tsv", sep="\t")

    ######### What are the most common funcions amongst the Enriched/Depleated species in the Disease Cohorts?
    ######### Functional clusters in encriched/depleated MSP by disease
    #Get list of species enriched, and list of species depleated, with cohort information
    ES_sig_up=df_es_ctrl2[df_es_ctrl2["ES_up"]>0.3]
    ES_sig_up=pd.merge(ES_sig_up, cohorts, how='left', on=['cohort'])
    ES_sig_dw=df_es_ctrl2[df_es_ctrl2["ES_down"]>0.3]
    ES_sig_dw=pd.merge(ES_sig_dw, cohorts, how='left', on=['cohort'])
    cnt=Counter(cohorts["category"])
    for dis in cnt:
        if cnt[dis]>4:
            volcano=[]
            msp_up=Counter(ES_sig_up[ES_sig_up["category"]==dis]["msp"])
            print(msp_up)
            msp_dw=Counter(ES_sig_dw[ES_sig_dw["category"]==dis]["msp"])
            msp_list=set(msp_up.keys()).union(set(msp_dw.keys()))
            for msp in msp_list:
                if (msp_up[msp] > 1) or (msp_dw[msp] > 1):
                    sm=msp_up[msp]+msp_dw[msp]
                    dif=msp_up[msp]-msp_dw[msp]
                    line=[msp,sm,dif]
                    volcano.append(line)
            #df_volc=pd.DataFrame(volcano, columns=["msp",'sum','dif'])
            #df_volc.to_csv(dis+".disease.tsv", sep="\t")




    ###transform species into function_cluster, conserve cohort information
    #read csv with MSP to functional cluster_info

    #count number of times a function_cluster is enriched in a cohort

    #depleated by country or depleated by disease test


if __name__ == "__main__":
    main()
