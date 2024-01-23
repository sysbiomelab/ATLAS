

source("~/rscripts/2019_funcType/funcType.analysis.functions.r")
source("~/rscripts/2019_funcType/funcType.data.files.r")



metaRData = "~/sunjae.mAtlas.func.annotation/all.basic.clean.metadata.behcet.20190805.RData"
load(metaRData)

mgsRData = "~/sunjae.mAtlas.func.annotation/mgs.clean.updated.20190726.10M.RData"
load(mgsRData)

#funcModulesRData= "~/sunjae.mAtlas.func.annotation/funcModules.20190828.RData"
funcModulesRData= "~/sunjae.mAtlas.func.annotation/funcModules.20191227.RData"
load(funcModulesRData)

funcMatRData = "~/sunjae.mAtlas.func.annotation/funcMat.20190808.RData"
load(funcMatRData)

gssOssTabNewFile = "~/sunjae.mAtlas.func.annotation/MSP_set_class_tropism_20190523_short.txt"
gssOssTabNew = read.table(gssOssTabNewFile, header=T, sep="\t", stringsAsFactors = F)
rownames(gssOssTabNew) = gssOssTabNew$MSP
taxo = gssOssTabNew
gutTaxo = taxo

zMatUpdated = scale(t(mergeMatUpdated))
zMatUpdated = t(zMatUpdated)

getMspByModules <- function(funcModules, funcMat, cutRatio=0.75) {
  mspByModules = lapply(funcModules, function(currTerms){
    
    if (length(currTerms)<2) {
      currMspsByCurrTerms = funcMat[,colnames(funcMat) %in% currTerms]
    } else {
      currMat = funcMat[,colnames(funcMat) %in% currTerms]
      currMspsByCurrTerms = rowSums(currMat)
    }
    currMspsByCurrTerms = sort(currMspsByCurrTerms, decreasing = T)
    outMsps = currMspsByCurrTerms[currMspsByCurrTerms >= (length(currTerms)*cutRatio)]
    return(names(outMsps))
  })
  names(mspByModules) = names(funcModules)
  return(mspByModules)
}


mspByModules = getMspByModules(funcModules, funcMat)
mspCountsByModules = do.call(c, lapply(mspByModules,length) )
hist(mspCountsByModules)

getBSCorrelationsByModules <- function(mspByModules, zMat, pCut=0.01) {
  bsCorrMat <- do.call(rbind, lapply(mspByModules, function(currMsps) {
    inclusion = rep(0,dim(zMat)[1])
    names(inclusion) = rownames(zMat)
    inclusion[names(inclusion) %in% currMsps] = 1 
    
    
    corrVec = apply(zMat, 2, function(currCol){
      if (sum(currCol)==0) return(0)
      if (sum(inclusion)==0) return(0)
      corrOut = cor.test(currCol, inclusion)
      rho = corrOut$estimate
      pvalue = corrOut$p.value
      if (pvalue > pCut) rho=0
      return(rho)
    })
    
    
    return(corrVec)
    
  }))
  rownames(bsCorrMat) = names(mspByModules)
  return(bsCorrMat)
}
bsCorrMat <- getBSCorrelationsByModules(mspByModules, zMatUpdated)
#bsCorrMatRData = "~/sunjae.mAtlas.func.annotation/bsCorrMat.RData"
bsCorrMatRData = "~/sunjae.mAtlas.func.annotation/bsCorrMat.20191227.RData"
save(bsCorrMat, file=bsCorrMatRData)


