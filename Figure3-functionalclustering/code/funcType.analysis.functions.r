require(compiler)
require(reshape2)
require(igraph)

getEffectSize <- function(caseVec, contVec, imputeVal=10) {
  wtest = wilcox.test(caseVec, contVec, alternative = "greater")
  N = length(caseVec) + length(contVec)
  Z = qnorm(1-wtest$p.value)
  es = (Z)/sqrt(N)
  if (es == Inf) es = imputeVal
  if (es == -Inf) es = - imputeVal
  
  return(es)
}

getEffectSizeMat <- function(caseMat, contMat, imputeVal=10) {
  
  caseIds = colnames(caseMat)
  contIds = colnames(contMat)
  caseContMat = cbind(caseMat, contMat)
  
  caseInd = which(colnames(caseContMat) %in% caseIds)
  contInd = which(colnames(caseContMat) %in% contIds)
  
  
  effectSizes = apply(caseContMat, 1, function(currVec) {
    caseVec = currVec[caseInd]
    contVec = currVec[contInd]
    ES = getEffectSize(caseVec, contVec)
    return(ES)
  })
  
  return(effectSizes)
}


getEffectSizeOfWilcox <- function(caseVec, contVec, imputeVal=10) {
  if (mean(contVec)==0) return(0)
  
  #wtest = wilcox.test(caseVec, contVec, alternative = "greater")
  wtest = wilcox.test(contVec, caseVec, alternative = "greater")
  if (wtest$p.value == 1) return(0)
  
  N = length(caseVec) + length(contVec)
  #Z = qnorm(wtest$p.value)
  Z = qnorm(1-wtest$p.value)
  if (Z < 0) return(0)
  
  es = abs(Z)/sqrt(N)
  if (es == Inf) es=imputeVal
  
  return(es)
}
getEffectSizeUpOfWilcox <- function(caseVec, contVec, imputeVal=10) {
  if (mean(contVec)==0) return(0)
  
  wtest = wilcox.test(contVec, caseVec, alternative = "greater")
  if (wtest$p.value == 1) return(0)
  
  N = length(caseVec) + length(contVec)
  Z = qnorm(wtest$p.value)
  if (Z < 0) return(0)
  
  es = abs(Z)/sqrt(N)
  
  if (es == Inf) es=imputeVal
  return(es)
}
getEffectSizeDownOfWilcox <- function(caseVec, contVec, imputeVal=10) {
  if (mean(caseVec)==0) return(0)
  if (mean(caseVec)==0 & mean(contVec)==0) return(imputeVal)
  
  wtest = wilcox.test(contVec, caseVec, alternative = "less")
  if (wtest$p.value == 1) return(0)
  
  N = length(caseVec) + length(contVec)
  Z = qnorm(wtest$p.value)
  if (Z < 0) return(0)
  
  es = abs(Z)/sqrt(N)
  
  if (es == Inf) es=imputeVal
  return(es)
}

statisticsByCountrySpecific <- function(samplesByAllCountry, mergeMat, esCut=0.2, numCut=11) {
  
  targetCountries = names(samplesByAllCountry)
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaAllList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=3
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,samplesByAllCountry[[targetCountry]]]
    
    refCountries = names(samplesByAllCountry)
    indiRefCountries = names(samplesByAllCountry)
    
    refCountries = refCountries[refCountries != targetCountry]
    numRefCountries = length(refCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numRefCountries) {
        refCountry = refCountries[refInd]
        refMat = mergeMat[,samplesByAllCountry[[refCountry]]]
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      effectSizeMatAll = do.call(cbind, effectSizeList)
      rownames(effectSizeMatAll) = rownames(mergeMat)
    }
    
    allEffectiveAll = apply(effectSizeMatAll, 1, function(currRow) sum(currRow>=esCut))
    allTargetSpecies = rownames(effectSizeMatAll)[allEffectiveAll>=numCut]
    
    print(allTargetSpecies)
    print(getSpeciesName(allTargetSpecies,taxo))
    specificBacteriaAllList[[targetCountry]] = allTargetSpecies
    
  }
  return(specificBacteriaAllList)
}
statisticsByUrbanSpecific <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByNotIndiCountry)
  indiRefCountries = names(sampleListByIndiCountry)
  
  numTargetCountries = length(targetCountries)
  numIndiRefCountries = length(indiRefCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  
  specificBacteriaUrbanList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=4
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByNotIndiCountry[[targetCountry]]]
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      rownames(effectSizeMatIndi) = rownames(mergeMat)
      
    }
    
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=4]
    
    print(ruralTargetSpecies)
    print(getSpeciesName(ruralTargetSpecies, taxo))
    
    specificBacteriaUrbanList[[targetCountry]] = ruralTargetSpecies
    
  }
  return(specificBacteriaUrbanList)
}
statisticsByDiseaseSpecific <- function(sampleListByNotIndiCountry, samplesByDisease, mergeMat, esCut=0.2) {
  
  targetCountries = names(samplesByDisease)
  notIndiRefCountries = names(sampleListByNotIndiCountry)
  
  numTargetCountries = length(targetCountries)
  numNotIndiRefCountries = length(notIndiRefCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaDiseaseList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=4
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,samplesByDisease[[targetCountry]]]
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      rownames(effectSizeMatUrban) = rownames(mergeMat)
      
    }
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=8]
    
    print(urbanTargetSpecies)
    print(getSpeciesName(urbanTargetSpecies, taxo))
    
    specificBacteriaDiseaseList[[targetCountry]] = urbanTargetSpecies
    
  }
  return(specificBacteriaDiseaseList)
}


statisticsByIndiSpecific <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByIndiCountry)
  
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaIndiList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=1
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByIndiCountry[[targetCountry]]]
    
    notIndiRefCountries = names(sampleListByNotIndiCountry)
    indiRefCountries = names(sampleListByIndiCountry)
    
    notIndiRefCountries = notIndiRefCountries[notIndiRefCountries != targetCountry]
    numNotIndiRefCountries = length(notIndiRefCountries)
    
    indiRefCountries = indiRefCountries[indiRefCountries != targetCountry]
    numIndiRefCountries = length(indiRefCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatUrban = do.call(cbind, foldList)
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      
    }
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatIndi = do.call(cbind, foldList)
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      
    }
    
    rownames(foldChangeMatUrban) = rownames(mergeMat)
    rownames(effectSizeMatUrban) = rownames(mergeMat)
    
    rownames(foldChangeMatIndi) = rownames(mergeMat)
    rownames(effectSizeMatIndi) = rownames(mergeMat)
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=8]
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=3]
    
    overlapSpecies = urbanTargetSpecies[urbanTargetSpecies%in%ruralTargetSpecies]
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
    
    specificBacteriaIndiList[[targetCountry]] = overlapSpecies
    
    
  }
  return(specificBacteriaIndiList)
  
}
statisticsByUrbanCountry <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByNotIndiCountry)
  
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaUrbanList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=3
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByNotIndiCountry[[targetCountry]]]
    
    notIndiRefCountries = names(sampleListByNotIndiCountry)
    indiRefCountries = names(sampleListByIndiCountry)
    
    notIndiRefCountries = notIndiRefCountries[notIndiRefCountries != targetCountry]
    numNotIndiRefCountries = length(notIndiRefCountries)
    
    indiRefCountries = indiRefCountries[indiRefCountries != targetCountry]
    numIndiRefCountries = length(indiRefCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatUrban = do.call(cbind, foldList)
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      rownames(foldChangeMatUrban) = rownames(mergeMat)
      rownames(effectSizeMatUrban) = rownames(mergeMat)
      
    }
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatIndi = do.call(cbind, foldList)
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      rownames(foldChangeMatIndi) = rownames(mergeMat)
      rownames(effectSizeMatIndi) = rownames(mergeMat)
      
    }
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=7]
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=4]
    
    overlapSpecies = urbanTargetSpecies[urbanTargetSpecies%in%ruralTargetSpecies]
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
    specificBacteriaUrbanList[[targetCountry]] = overlapSpecies
    
  }
  return(specificBacteriaUrbanList)
}

statisticsByIndiCountry <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByIndiCountry)
  
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaIndiList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=1
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByIndiCountry[[targetCountry]]]
    
    notIndiRefCountries = names(sampleListByNotIndiCountry)
    indiRefCountries = names(sampleListByIndiCountry)
    
    notIndiRefCountries = notIndiRefCountries[notIndiRefCountries != targetCountry]
    numNotIndiRefCountries = length(notIndiRefCountries)
    
    indiRefCountries = indiRefCountries[indiRefCountries != targetCountry]
    numIndiRefCountries = length(indiRefCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatUrban = do.call(cbind, foldList)
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      
    }
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatIndi = do.call(cbind, foldList)
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      
    }
    
    rownames(foldChangeMatUrban) = rownames(mergeMat)
    rownames(effectSizeMatUrban) = rownames(mergeMat)
    
    rownames(foldChangeMatIndi) = rownames(mergeMat)
    rownames(effectSizeMatIndi) = rownames(mergeMat)
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=8]
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=3]
    
    overlapSpecies = urbanTargetSpecies[urbanTargetSpecies%in%ruralTargetSpecies]
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
    
    specificBacteriaIndiList[[targetCountry]] = overlapSpecies
    
    
  }
  return(specificBacteriaIndiList)
  
}
statisticsByDisease <- function(samplesByDisease, sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.2, numCut=17, allCut=7) {
  
  targetDiseases = names(samplesByDisease)
  numTargetDisease = length(targetDiseases)
  
  allCountries = names(sampleListByNotIndiCountry)
  numAllCountries = length(sampleListByNotIndiCountry)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaDiseaseList = list()
  for (diseaseInd in 1:numTargetDisease) {
    
    #diseaseInd=8
    targetDisease = targetDiseases[diseaseInd]
    print(targetDisease)
    targetMat = mergeMat[,samplesByDisease[[targetDisease]]]
    
    refDiseases = names(samplesByDisease)
    refDiseases = refDiseases[refDiseases != targetDisease]
    numRefDiseases = length(refDiseases)
    
    
    allMode = T
    if (allMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numAllCountries) {
        refCountry = allCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatAll = do.call(cbind, foldList)
      effectSizeMatAll = do.call(cbind, effectSizeList)
      rownames(foldChangeMatAll) = rownames(targetMat)
      rownames(effectSizeMatAll) = rownames(targetMat)
      
    }
    normalMode = F
    if (normalMode) {
      
      foldChanges = rowMeans(targetMat)/rowMeans(normalMat)
      names(foldChanges) = rownames(targetMat)
      names(foldChanges) = rownames(targetMat)
      foldChanges[is.nan(foldChanges)]=1
      foldChanges[foldChanges==0]=1e-1000
      foldChanges[foldChanges==Inf]=10^1000
      logFoldChanges = log10(foldChanges)
      
      foldChangeMatNormal = logFoldChanges
      
      effectSizeMatNormal = sapply(1:lenSpecies, function(ind){
        effectSize=getEffectSizeOfWilcox(normalMat[ind,],targetMat[ind,])
        return(effectSize)
      })
      names(foldChangeMatNormal) = rownames(targetMat)
      names(effectSizeMatNormal) = rownames(targetMat)
      
    }
    
    diseaseMode = T
    if (diseaseMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numRefDiseases) {
        refDisease = refDiseases[refInd]
        refMat = mergeMat[,samplesByDisease[[refDisease]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refDisease]] = logFoldChanges
        
        effectSizeList[[refDisease]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatDisease = do.call(cbind, foldList)
      effectSizeMatDisease = do.call(cbind, effectSizeList)
      rownames(foldChangeMatDisease) = rownames(mergeMat)
      rownames(effectSizeMatDisease) = rownames(mergeMat)
      
      
    }
    
    
    #normalTargetSpecies = names(effectSizeMatNormal[effectSizeMatNormal>=esCut])
    
    allEffectiveDisease = apply(effectSizeMatDisease, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveAll = apply(effectSizeMatAll, 1, function(currRow) sum(currRow>=esCut))
    
    diseaseTargetSpecies = rownames(effectSizeMatDisease)[allEffectiveDisease>=numCut]
    allTargetSpecies = rownames(effectSizeMatAll)[allEffectiveAll>=allCut]
    
    overlapSpecies = diseaseTargetSpecies[diseaseTargetSpecies %in% allTargetSpecies]
    
    specificBacteriaDiseaseList[[targetDisease]] = overlapSpecies
    table(diseaseTargetSpecies %in% allTargetSpecies)
    print(diseaseTargetSpecies)
    print(getSpeciesName(diseaseTargetSpecies,taxo))
    
    print(allTargetSpecies)
    print(getSpeciesName(allTargetSpecies,taxo))
    
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
  }
  return(specificBacteriaDiseaseList)
}
statisticsByDiseaseGeneral <- function(samplesByDisease, sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.2, numCut=17) {
  
  targetDiseases = names(samplesByDisease)
  numTargetDisease = length(targetDiseases)
  
  allCountries = names(sampleListByNotIndiCountry)
  numAllCountries = length(sampleListByNotIndiCountry)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaDiseaseList = list()
  for (diseaseInd in 1:numTargetDisease) {
    
    #diseaseInd=8
    targetDisease = targetDiseases[diseaseInd]
    print(targetDisease)
    targetMat = mergeMat[,samplesByDisease[[targetDisease]]]
    
    refDiseases = names(samplesByDisease)
    refDiseases = refDiseases[refDiseases != targetDisease]
    numRefDiseases = length(refDiseases)
    
    
    allMode = T
    if (allMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numAllCountries) {
        refCountry = allCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatAll = do.call(cbind, foldList)
      effectSizeMatAll = do.call(cbind, effectSizeList)
      rownames(foldChangeMatAll) = rownames(targetMat)
      rownames(effectSizeMatAll) = rownames(targetMat)
      
    }
    
    allEffectiveAll = apply(effectSizeMatAll, 1, function(currRow) sum(currRow>=esCut))
    allTargetSpecies = rownames(effectSizeMatAll)[allEffectiveAll>=allCut]
    specificBacteriaDiseaseList[[targetDisease]] = allTargetSpecies
    
    print(allTargetSpecies)
    print(getSpeciesName(allTargetSpecies,taxo))
    
  }
  return(specificBacteriaDiseaseList)
}

statisticsDownByUrbanCountry <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByNotIndiCountry)
  
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaUrbanList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=3
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByNotIndiCountry[[targetCountry]]]
    
    notIndiRefCountries = names(sampleListByNotIndiCountry)
    indiRefCountries = names(sampleListByIndiCountry)
    
    notIndiRefCountries = notIndiRefCountries[notIndiRefCountries != targetCountry]
    numNotIndiRefCountries = length(notIndiRefCountries)
    
    indiRefCountries = indiRefCountries[indiRefCountries != targetCountry]
    numIndiRefCountries = length(indiRefCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeDownOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatUrban = do.call(cbind, foldList)
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      rownames(foldChangeMatUrban) = rownames(mergeMat)
      rownames(effectSizeMatUrban) = rownames(mergeMat)
      
    }
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeDownOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatIndi = do.call(cbind, foldList)
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      rownames(foldChangeMatIndi) = rownames(mergeMat)
      rownames(effectSizeMatIndi) = rownames(mergeMat)
      
    }
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=7]
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=4]
    
    overlapSpecies = urbanTargetSpecies[urbanTargetSpecies%in%ruralTargetSpecies]
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
    specificBacteriaUrbanList[[targetCountry]] = overlapSpecies
    
  }
  return(specificBacteriaUrbanList)
}
statisticsDownByIndiCountry <- function(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2) {
  
  targetCountries = names(sampleListByIndiCountry)
  
  numTargetCountries = length(targetCountries)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaIndiList = list()
  for (countryInd in 1:numTargetCountries) {
    
    #countryInd=4
    targetCountry = targetCountries[countryInd]
    print(targetCountry)
    targetMat = mergeMat[,sampleListByIndiCountry[[targetCountry]]]
    
    notIndiRefCountries = names(sampleListByNotIndiCountry)
    indiRefCountries = names(sampleListByIndiCountry)
    
    notIndiRefCountries = notIndiRefCountries[notIndiRefCountries != targetCountry]
    numNotIndiRefCountries = length(notIndiRefCountries)
    
    indiRefCountries = indiRefCountries[indiRefCountries != targetCountry]
    numIndiRefCountries = length(indiRefCountries)
    
    urbanMode = T
    if (urbanMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numNotIndiRefCountries) {
        refCountry = notIndiRefCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeDownOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatUrban = do.call(cbind, foldList)
      effectSizeMatUrban = do.call(cbind, effectSizeList)
      
    }
    
    indiMode = T
    if (indiMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numIndiRefCountries) {
        refCountry = indiRefCountries[refInd]
        refMat = mergeMat[,sampleListByIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatIndi = do.call(cbind, foldList)
      effectSizeMatIndi = do.call(cbind, effectSizeList)
      
    }
    
    rownames(foldChangeMatUrban) = rownames(mergeMat)
    rownames(effectSizeMatUrban) = rownames(mergeMat)
    
    rownames(foldChangeMatIndi) = rownames(mergeMat)
    rownames(effectSizeMatIndi) = rownames(mergeMat)
    
    allEffectiveUrban = apply(effectSizeMatUrban, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveIndi = apply(effectSizeMatIndi, 1, function(currRow) sum(currRow>=esCut))
    
    urbanTargetSpecies = rownames(effectSizeMatUrban)[allEffectiveUrban>=8]
    ruralTargetSpecies = rownames(effectSizeMatIndi)[allEffectiveIndi>=3]
    
    overlapSpecies = urbanTargetSpecies[urbanTargetSpecies%in%ruralTargetSpecies]
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
    
    specificBacteriaIndiList[[targetCountry]] = overlapSpecies
    
    
  }
  return(specificBacteriaIndiList)
  
}
statisticsDownByDisease <- function(samplesByDisease, sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.2, numCut=17, allCut=8) {
  
  targetDiseases = names(samplesByDisease)
  numTargetDisease = length(targetDiseases)
  
  allCountries = names(sampleListByNotIndiCountry)
  numAllCountries = length(sampleListByNotIndiCountry)
  
  numMgs = dim(mergeMat)[1]
  lenSpecies = dim(mergeMat)[1]
  
  specificBacteriaDiseaseList = list()
  for (diseaseInd in 1:numTargetDisease) {
    
    #diseaseInd=1
    targetDisease = targetDiseases[diseaseInd]
    print(targetDisease)
    targetMat = mergeMat[,samplesByDisease[[targetDisease]]]
    
    refDiseases = names(samplesByDisease)
    refDiseases = refDiseases[refDiseases != targetDisease]
    numRefDiseases = length(refDiseases)
    
    
    allMode = T
    if (allMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numAllCountries) {
        refCountry = allCountries[refInd]
        refMat = mergeMat[,sampleListByNotIndiCountry[[refCountry]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refCountry]] = logFoldChanges
        
        effectSizeList[[refCountry]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeDownOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatAll = do.call(cbind, foldList)
      effectSizeMatAll = do.call(cbind, effectSizeList)
      rownames(foldChangeMatAll) = rownames(targetMat)
      rownames(effectSizeMatAll) = rownames(targetMat)
      
    }
    normalMode = F
    if (normalMode) {
      
      foldChanges = rowMeans(targetMat)/rowMeans(normalMat)
      names(foldChanges) = rownames(targetMat)
      names(foldChanges) = rownames(targetMat)
      foldChanges[is.nan(foldChanges)]=1
      foldChanges[foldChanges==0]=1e-1000
      foldChanges[foldChanges==Inf]=10^1000
      logFoldChanges = log10(foldChanges)
      
      foldChangeMatNormal = logFoldChanges
      
      effectSizeMatNormal = sapply(1:lenSpecies, function(ind){
        effectSize=getEffectSizeDownOfWilcox(normalMat[ind,],targetMat[ind,])
        return(effectSize)
      })
      names(foldChangeMatNormal) = rownames(targetMat)
      names(effectSizeMatNormal) = rownames(targetMat)
      
    }
    
    diseaseMode = T
    if (diseaseMode) {
      foldList = list()
      effectSizeList = list()
      
      for (refInd in 1:numRefDiseases) {
        refDisease = refDiseases[refInd]
        refMat = mergeMat[,samplesByDisease[[refDisease]]]
        foldChanges = rowMeans(targetMat)/rowMeans(refMat)
        names(foldChanges) = rownames(targetMat)
        foldChanges[is.nan(foldChanges)]=1
        foldChanges[foldChanges==0]=1e-1000
        foldChanges[foldChanges==Inf]=10^1000
        logFoldChanges = log10(foldChanges)
        
        foldList[[refDisease]] = logFoldChanges
        
        effectSizeList[[refDisease]] = sapply(1:lenSpecies, function(ind){
          effectSize=getEffectSizeDownOfWilcox(refMat[ind,],targetMat[ind,])
          return(effectSize)
        })
      }
      
      foldChangeMatDisease = do.call(cbind, foldList)
      effectSizeMatDisease = do.call(cbind, effectSizeList)
      
    }
    
    rownames(foldChangeMatDisease) = rownames(mergeMat)
    rownames(effectSizeMatDisease) = rownames(mergeMat)
    
    #normalTargetSpecies = names(effectSizeMatNormal[effectSizeMatNormal>=esCut])
    
    allEffectiveDisease = apply(effectSizeMatDisease, 1, function(currRow) sum(currRow>=esCut))
    allEffectiveAll = apply(effectSizeMatAll, 1, function(currRow) sum(currRow>=esCut))
    
    diseaseTargetSpecies = rownames(effectSizeMatDisease)[allEffectiveDisease>=numCut]
    allTargetSpecies = rownames(effectSizeMatAll)[allEffectiveAll>=allCut]
    
    overlapSpecies = diseaseTargetSpecies[diseaseTargetSpecies %in% allTargetSpecies]
    
    specificBacteriaDiseaseList[[targetDisease]] = overlapSpecies
    table(diseaseTargetSpecies %in% allTargetSpecies)
    print(diseaseTargetSpecies)
    print(getSpeciesName(diseaseTargetSpecies,taxo))
    
    print(allTargetSpecies)
    print(getSpeciesName(allTargetSpecies,taxo))
    
    print(overlapSpecies)
    print(getSpeciesName(overlapSpecies,taxo))
  }
  return(specificBacteriaDiseaseList)
}

getRelAbdOfTaxa <- function(targetMat, targetTaxo, taxo) {
  targetTaxoMat = getTaxaSumMat(targetMat, taxo, targetTaxo,T)
  targetTaxoMean = rowMeans(targetTaxoMat)
  relAbds = targetTaxoMean/sum(targetTaxoMean)
  return(relAbds)
}


getGenesFromStr <- function( str ) {
  if (is.na(str)) return("")
  if (is.null(str)) return("")
  if (str=="") return("")
  out = strsplit(str, split = ";")[[1]]
  return(out)
}
getGenesFromStr_cmp = cmpfun(getGenesFromStr)
pathWeighting <- function(koCatalog, koElementList) {
  koTab = table(koCatalog$ko)
  koFreq = as.numeric(koTab)
  names(koFreq) = names(koTab)
  koFreq = koFreq/sum(koFreq)
  
  weights = do.call(c,lapply(koElementList, function(kos) {
    return(sum(koFreq[kos],na.rm=T))
  }))
  return(weights)
}

collectVipToMat <- function(vipList) {
  allNames = unique(do.call(c,lapply(vipList, names)))
  numRow = length(allNames)
  numCol = length(vipList)
  vipMat = matrix(0,numRow,numCol)
  colnames(vipMat) = names(vipList)
  rownames(vipMat) = allNames
  
}

getVipScores <- function(datasetIDs, refIDs, targetMat, basicMetaMap) {
  labeledMatRef = targetMat[,colnames(targetMat) %in% refIDs]
  
  vipList = lapply(as.list(datasetIDs),function(datasetId){
    targetIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID == datasetId]
    labeledMatTarget = targetMat[,colnames(targetMat) %in% targetIDs]
    colnames(labeledMatTarget)=rep("target",dim(labeledMatTarget)[2])
    colnames(labeledMatRef)=rep("cont",dim(labeledMatRef)[2])
    labeledMatTargetComp = cbind(labeledMatTarget, labeledMatRef)
    matCompPlsDa = opls(t(labeledMatTargetComp), colnames(labeledMatTargetComp),plotL=F)
    
    vipKo = getVipVn(matCompPlsDa)
    vipKoAbv1 = vipKo[vipKo>1]
    return(vipKoAbv1)
  })
  return(vipList)
  
  
}
plotVipList <- function(vipList,marVec = c(5,10,4,1)) {
  par(mar=marVec)
  print(names(vipList))
  print("red gray")
  newVips = c(sort(vipList[[2]]),sort(vipList[[1]]))
  newCols = c(rep("#A9A9A966",length(vipList[[2]])),
              rep("#FF000088",length(vipList[[1]])))
  
  out=barplot(newVips, col=newCols,
              border=NA,
              las=1, horiz = T, cex.names = 0.7, cex.axis = 0.7)
  pre=out[length(vipList[[2]]),]
  pro=out[length(vipList[[2]])+1,]
  abline(h=(pre+pro)/2)
  par(mar=c(5.1,4.1,4.1,2.1))
}

### stefy paper functions ###
getVipList <- function(vips, currMat, currClasses) {
  ### presumed to have two classes ###
  
  uniqClasses = unique(currClasses)
  print(uniqClasses)
  samplesByClass = split(colnames(currMat), currClasses)
  samplesByClass = samplesByClass[uniqClasses]
  
  entities = rownames(currMat)
  
  class1Inds = sapply(entities, function(ent){
    meanClass1 = median(currMat[ent,samplesByClass[[1]]])
    meanClass2 = median(currMat[ent,samplesByClass[[2]]])
    if (meanClass1 > meanClass2) return(T)
    return(F)
  })
  class1Entities = entities[class1Inds]
  class2Entities = entities[!class1Inds]
  
  vipClass1 = vips[vips > 1 & names(vips) %in% class1Entities]
  vipClass2 = vips[vips > 1 & names(vips) %in% class2Entities]
  out = list(vipClass1, vipClass2)
  names(out)=uniqClasses
  return(out)
  
}
takeTopOfList <- function(vipList, topNum) {
  resList = lapply(vipList, function(x){
    x = sort(x,decreasing = T)
    if(ifelse(length(x)>topNum, T, F)) {
      x=x[1:topNum]
    }
    return(x)
  })
  return(resList)
}
labelVipList <- function(vipList,taxo) {
  out=lapply(vipList, function(x){
    names(x) = getSpeciesName(names(x),taxo)
    return(x)
  })
  return(out)
}

plotVolcano <- function(volcanoRes, labelMode = T, vBorder=2, hBorder=-log10(1e-3)) {
  require(ggplot2)
  fontSize=2.5
  if (labelMode) {
    ggOut = ggplot(volcanoRes, aes(x=lfc, y=sig, size=relAbd, label=label)) + 
      geom_jitter(aes(size = relAbd, colour = lfc)) + 
      scale_size(range=c(3,15)) + 
      scale_colour_gradientn(colours = c("#0000FF99", "#ffffff99","#FF000099")) + 
      theme(axis.ticks.length = unit(.2, "cm"), 
            axis.line = element_line(colour = "black"), 
            text = element_text(size=15), 
            legend.position="none", 
            panel.background=element_rect(fill="white")) + 
      xlab("")+ylab("") + ylim(c(-20,180))+
      geom_text_repel(size=fontSize,
                      segment.color="#80808088",
                      force=10) +
      geom_hline(yintercept = hBorder, colour="#808080AA", linetype='dashed') + 
      geom_vline(xintercept = vBorder, colour="#808080AA", linetype='dashed') +
      geom_vline(xintercept = -vBorder, colour="#808080AA", linetype='dashed') 
    print(ggOut)
    
  }
  if (!labelMode) {
    ggOut = ggplot(volcanoRes, aes(x=lfc, y=sig, size=relAbd, label=label)) + 
      geom_jitter(aes(size = relAbd, colour = lfc)) + 
      scale_size(range=c(3,15)) + 
      scale_colour_gradientn(colours = c("#0000FF99", "#ffffff99","#FF000099")) + 
      theme(axis.ticks.length = unit(.2, "cm"), 
            axis.line = element_line(colour = "black"), 
            text = element_text(size=15), 
            legend.position="none", 
            panel.background=element_rect(fill="white")) + 
      xlab("")+ylab("") + ylim(c(-20,180))+
      geom_hline(yintercept = hBorder, colour="#808080AA", linetype='dashed') + 
      geom_vline(xintercept = vBorder, colour="#808080AA", linetype='dashed') +
      geom_vline(xintercept = -vBorder, colour="#808080AA", linetype='dashed') 
    print(ggOut)
    
  }
  return(ggOut)
}


### stefy paper functions ###


attachLabels <- function(currMat, labelMap) {
  colLabels = labelMap$L1[match(colnames(currMat),labelMap$value)]
  colnames(currMat) = colLabels
  return(currMat)
} 

rowCVs <- function(mat) {
  rM = rowMeans(mat)
  rSD = rowSds(mat)
  rCV = rSD/rM*100
  names(rCV) = rownames(mat)
  return(rCV)
}
mergeMatrix <- function(x,y) {
  out = merge(x, y, by="row.names", all=T)
  rownames(out) = out[,1]
  out = as.matrix(out[,-1])
  out[is.na(out)]=0
  return(out)
}
getNonOverlapMat <-function(refMat, targetMat, debug=F) {
  if (debug) print(table(colnames(targetMat)%in%colnames(refMat)))
  targetMat = targetMat[,!colnames(targetMat)%in% colnames(refMat)]
  return(targetMat)
}
getNonOverlapVec <-function(refVec, targetVec) {
  print(table(names(targetVec)%in%names(refVec)))
  targetVec = targetVec[!names(targetVec)%in% names(refVec)]
  if (length(targetVec)==0)return(NULL)
  return(targetVec)
}

removeZeroLiset <- function(currList) {
  currList = lapply(currList, function(currSamples){
    return(currSamples[currSamples!=0])
  })
  return(currList)
}
removeNaList <- function(currList) {
  currList = lapply(currList, function(currSamples){
    return(currSamples[!is.na(currSamples)])
  })
  return(currList)
}
catchValues <- function(sampleList, mgsRichness) {
  listNames = names(sampleList)
  resOut = lapply(sampleList, function(samples) {
    samples = samples[samples %in% names(mgsRichness)]
    
    if (length(samples)>0) {
      out = mgsRichness[names(mgsRichness)%in%samples]
      out = out[!is.na(out)]
      return(out)
    }
    return(NULL)
  })
  return(resOut)
}
catchMatrix <- function(sampleList, mgsMat) {
  listNames = names(sampleList)
  resOut = lapply(sampleList, function(samples) {
    samples = samples[samples %in% colnames(mgsMat)]
    
    if (length(samples)>0) {
      out = mgsMat[,colnames(mgsMat)%in%samples]
      out = rowMeans(out)
      return(out)
    }
    return(NULL)
  })
  res = do.call(cbind,resOut)
  return(res)
}

getTaxaSumMat <- function(mgsVec, taxo, target, debug=F, naOmit=T) {
  numSamples = dim(mgsVec)[2]
  
  #checking msp with taxo
  checkField <- any(colnames(taxo)==target)
  if (!checkField) {
    print("please check target field")
    return(NULL)
  }
  checkMsps <- as.numeric(table(rownames(mgsVec) %in% rownames(taxo))["TRUE"])
  if (debug) {
    print("checking consistency of MSP names overlapped ...")
    print(paste(checkMsps/length(rownames(mgsVec))*100, "%"))
  }
  
  taxoTab = data.frame(id=rownames(taxo),target=taxo[,target],stringsAsFactors = F)
  taxoSplit = split(taxoTab$id, taxoTab$target)
  
  taxoSumMat = do.call(rbind, lapply(taxoSplit, function(msps){#
    msps = msps[msps%in%rownames(mgsVec)]
    if (length(msps)==0) return(rep(0,numSamples))
    if (length(msps)==1) return(mgsVec[rownames(mgsVec)%in%msps,])
    return(colSums(mgsVec[rownames(mgsVec)%in%msps,],na.rm = naOmit))
    
    
  }))
  rownames(taxoSumMat) = names(taxoSplit)
  if (debug) {
    print("dimension of output matrix...")
    print(dim(taxoSumMat))
  }
  
  return(taxoSumMat)
}

getTaxaAvgMat <- function(mgsVec, taxo, target, debug=F, naOmit=T) {
  numSamples = dim(mgsVec)[2]
  
  #checking msp with taxo
  checkField <- any(colnames(taxo)==target)
  if (!checkField) {
    print("please check target field")
    return(NULL)
  }
  checkMsps <- as.numeric(table(rownames(mgsVec) %in% rownames(taxo))["TRUE"])
  if (debug) {
    print("checking consistency of MSP names overlapped ...")
    print(paste(checkMsps/length(rownames(mgsVec))*100, "%"))
  }
  
  taxoTab = data.frame(id=rownames(taxo),target=taxo[,target],stringsAsFactors = F)
  taxoSplit = split(taxoTab$id, taxoTab$target)
  
  taxoAvgMat = do.call(rbind, lapply(taxoSplit, function(msps){#
    msps = msps[msps%in%rownames(mgsVec)]
    if (length(msps)==0) return(rep(0,numSamples))
    if (length(msps)==1) return(mgsVec[rownames(mgsVec)%in%msps,])
    return(colMeans(mgsVec[rownames(mgsVec)%in%msps,],na.rm = naOmit))
    
    
  }))
  rownames(taxoAvgMat) = names(taxoSplit)
  if (debug) {
    print("dimension of output matrix...")
    print(dim(taxoAvgMat))
  }
  
  return(taxoAvgMat)
}

parseTaxoBy <-function(taxo, parent, child) {
  pcUniq = unique(taxo[,c(parent,child)])
  pcList = split(pcUniq[,2],pcUniq[,1])
  return(pcList)
}

parseTaxo <- function(taxo) {
  #phylum -> class -> order -> family -> genus -> species
  
  pcUniq = unique(taxo[,c("phylum","class")])
  pcList = split(pcUniq[,2],pcUniq[,1])
  
  coUniq = unique(taxo[,c("class","order")])
  coList = split(coUniq[,2],coUniq[,1])
  
  ofUniq = unique(taxo[,c("order","family")])
  ofList = split(ofUniq[,2],ofUniq[,1])
  
  fgUniq = unique(taxo[,c("family","genus")])
  fgList = split(fgUniq[,2],fgUniq[,1])
  
  gsUniq = unique(taxo[,c("genus","species")])
  gsList = split(gsUniq[,2],gsUniq[,1])
  
  smUniq = unique(data.frame(taxo$species, rownames(taxo), stringsAsFactors = F))
  smList = split(smUniq[,2],smUniq[,1])
  
  resList = list(phylum=pcList, 
                 class=coList,
                 order=ofList,
                 family=fgList,
                 genus=gsList,
                 species=smList)
  
  return(resList)
  
}

getParts <- function(parsedTaxo, parentMean, childMean, target=NULL) {
  if (is.null(target)) hierarchy = parsedTaxo
  if (!is.null(target)) hierarchy = parsedTaxo[[target]]
  
  parts = lapply(seq_along(hierarchy), function(y,n,i) {
    parent= n[[i]]
    child=y[[i]]
    partVec = sapply(child,function(x) childMean[x] / as.numeric(parentMean[parent]) )
    names(partVec) = child
    return(as.list(partVec))
    
  } , y=hierarchy, n=names(hierarchy) )
  names(parts) = names(hierarchy)
  return(parts)
}

getCors <- function(parsedTaxo, parentMat, childMat, target=NULL) {
  
  
  if (is.null(target)) hierarchy = parsedTaxo
  if (!is.null(target)) hierarchy = parsedTaxo[[target]]
  
  parts = lapply(seq_along(hierarchy), function(y,n,i) {
    parent= n[[i]]
    child=y[[i]]
    corVec = sapply(child,function(x) cor(childMat[x,], parentMat[parent,],method="spearman") )
    names(corVec) = child
    return(as.list(corVec))
    
  } , y=hierarchy, n=names(hierarchy) )
  names(parts) = names(hierarchy)
  return(parts)
}

getPartLink <- function(taxo, mergeMat, exUnclassified = F) {
  parsedTaxo = parseTaxo(taxo)
  
  phylaMat = getTaxaSumMat(mergeMat,taxo,"phylum",T)
  classMat = getTaxaSumMat(mergeMat,taxo,"class",T)
  orderMat = getTaxaSumMat(mergeMat,taxo,"order",T)
  familyMat = getTaxaSumMat(mergeMat,taxo,"family",T)
  genusMat = getTaxaSumMat(mergeMat,taxo,"genus",T)
  speciesMat  = getTaxaSumMat(mergeMat,taxo,"species",T)
  
  phylaMean = rowMeans(phylaMat)
  classMean = rowMeans(classMat)
  orderMean = rowMeans(orderMat)
  familyMean = rowMeans(familyMat)
  genusMean = rowMeans(genusMat)
  speciesMean = rowMeans(speciesMat)
  
  pPart = phylaMean/sum(phylaMean)
  fPart = familyMean/sum(familyMean)
  
  pcPart = getParts(parsedTaxo, phylaMean, classMean, "phylum")
  coPart = getParts(parsedTaxo, classMean, orderMean, "class")
  ofPart = getParts(parsedTaxo, orderMean, familyMean, "order")
  fgPart = getParts(parsedTaxo, familyMean, genusMean, "family")
  gsPart = getParts(parsedTaxo, genusMean, speciesMean, "genus")
  
  pcMelt = melt(pcPart)
  coMelt = melt(coPart)
  ofMelt = melt(ofPart)
  fgMelt = melt(fgPart)
  gsMelt = melt(gsPart) 
  
  
  require(plotrix)
  ggplot(pcMelt) + geom_rect(aes(fill=L1,xmin=0,ymin=0,xmax=1,ymax=1)) +
    geom_rect(aes(fill=L2,xmin=0,ymin=0,xmax=1,ymax=1)) +
    theme(aspect.ratio=1) +
    coord_polar(theta="y")
  require(RColorBrewer)
  pie(1, radius=0.01, init.angle=90, col=c('white'), border = NA, labels='')
  floating.pie(0,0,sort(pPart),
               radius=5*iniR, 
               startpos=pi/2, col = brewer.pal(15,"RdYlBu"),
               border=NA)
  floating.pie(0,0,sort(pPart),  
               radius=3*iniR, 
               startpos=pi/2, 
               col="white",
               border=NA)
  
  pie(1, radius=0.01, init.angle=90, col=c('white'), border = NA, labels='')
  floating.pie(0,0,sort(pPart)[-c(11:15)],
               radius=5*iniR, 
               startpos=pi/2, col = brewer.pal(15,"RdYlBu"),
               border=NA)
  floating.pie(0,0,sort(pPart)[-c(11:15)],  
               radius=3*iniR, 
               startpos=pi/2, 
               col="white",
               border=NA)
  
  iniR=0.2
  pie(1, radius=0.01, init.angle=90, col=c('white'), border = NA, labels='')
  floating.pie(0,0,sort(fPart)[-c(82:90)],
               radius=5*iniR, 
               startpos=pi/2, col = brewer.pal(15,"RdYlBu"),
               border=NA)
  floating.pie(0,0,sort(fPart)[-c(82:90)],  
               radius=3*iniR, 
               startpos=pi/2, 
               col="white",
               border=NA)
  
  legend(0, 5*iniR, gsub("_"," ",names(pPart)[-1]), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
  
  aggregate(pcMelt$value, by=list(pcMelt$L1), FUN=sum)
  
  #col=as.character(colors[c('exons','NO','downstream','NO')]),
  
  
  
}

getTaxaStatPhylumFamily <- function(taxo, mergeMat) {
  taxoPhylumFamily = parseTaxoBy(taxo, "phylum", "family")
  
  phylaMat = getTaxaSumMat(mergeMat,taxo,"phylum",T)
  familyMat = getTaxaSumMat(mergeMat,taxo,"family",T)
  
  # classMat = getTaxaSumMat(mergeMat,taxo,"class",T)
  # orderMat = getTaxaSumMat(mergeMat,taxo,"order",T)
  # genusMat = getTaxaSumMat(mergeMat,taxo,"genus",T)
  # speciesMat  = getTaxaSumMat(mergeMat,taxo,"species",T)
  
  phylaMean = rowMeans(phylaMat)
  familyMean = rowMeans(familyMat)
  
  # classMean = rowMeans(classMat)
  # orderMean = rowMeans(orderMat)
  # genusMean = rowMeans(genusMat)
  # speciesMean = rowMeans(speciesMat)
  
  
  
  pfPart = getParts(taxoPhylumFamily, phylaMean, familyMean)
  pfCor = getCors(taxoPhylumFamily, phylaMat, familyMat)
  
  pfPartMelt = melt(pfPart)
  pfCorMelt = melt(pfCor)
  
  pfStat = cbind(pfPartMelt$value, pfCorMelt$value)
  rownames(pfStat) = pfPartMelt$L2
  colnames(pfStat) = c("","")
  pfStat[is.na(pfStat)]=0
  
  return(list(pfStat, pfPartMelt))
}
getTaxaStat <- function(taxo, mergeMat, exUnclassified = F ) {
  parsedTaxo = parseTaxo(taxo)
  
  phylaMat = getTaxaSumMat(mergeMat,taxo,"phylum",T)
  classMat = getTaxaSumMat(mergeMat,taxo,"class",T)
  orderMat = getTaxaSumMat(mergeMat,taxo,"order",T)
  familyMat = getTaxaSumMat(mergeMat,taxo,"family",T)
  genusMat = getTaxaSumMat(mergeMat,taxo,"genus",T)
  speciesMat  = getTaxaSumMat(mergeMat,taxo,"species",T)
  
  phylaMean = rowMeans(phylaMat)
  classMean = rowMeans(classMat)
  orderMean = rowMeans(orderMat)
  familyMean = rowMeans(familyMat)
  genusMean = rowMeans(genusMat)
  speciesMean = rowMeans(speciesMat)
  
  pPart = phylaMean/sum(phylaMean)
  
  pcPart = getParts(parsedTaxo, phylaMean, classMean, "phylum")
  coPart = getParts(parsedTaxo, classMean, orderMean, "class")
  ofPart = getParts(parsedTaxo, orderMean, familyMean, "order")
  fgPart = getParts(parsedTaxo, familyMean, genusMean, "family")
  gsPart = getParts(parsedTaxo, genusMean, speciesMean, "genus")
  
  pfPart = getParts(parsedTaxo, phylaMean, familyMean, "phylum")
  
  pcCor = getCors(parsedTaxo, phylaMat, classMat, "phylum")
  coCor = getCors(parsedTaxo, classMat, orderMat, "class")
  ofCor = getCors(parsedTaxo, orderMat, familyMat, "order")
  fgCor = getCors(parsedTaxo, familyMat, genusMat, "family")
  gsCor = getCors(parsedTaxo, genusMat, speciesMat, "genus")
  
  pcPartMelt = melt(pcPart)
  coPartMelt = melt(coPart)
  ofPartMelt = melt(ofPart)
  fgPartMelt = melt(fgPart)
  gsPartMelt = melt(gsPart)
  
  pcCorMelt = melt(pcCor)
  coCorMelt = melt(coCor)
  ofCorMelt = melt(ofCor)
  fgCorMelt = melt(fgCor)
  gsCorMelt = melt(gsCor)
  
  pcStat = cbind(pcPartMelt$value, pcCorMelt$value)
  coStat = cbind(coPartMelt$value, coCorMelt$value)
  ofStat = cbind(ofPartMelt$value, ofCorMelt$value)
  fgStat = cbind(fgPartMelt$value, fgCorMelt$value)
  gsStat = cbind(gsPartMelt$value, gsCorMelt$value)
  
  rownames(pcStat) = pcPartMelt$L2
  rownames(coStat) = coPartMelt$L2
  rownames(ofStat) = ofPartMelt$L2
  rownames(fgStat) = fgPartMelt$L2
  rownames(gsStat) = gsPartMelt$L2
  
  colnames(pcStat) = c("","")
  colnames(coStat) = c("","")
  colnames(ofStat) = c("","")
  colnames(fgStat) = c("","")
  colnames(gsStat) = c("","")
  
  
  pPart = cbind(pPart)
  colnames(pPart)=""
  
  taxoPhylumGenus = parseTaxoBy(taxo, "phylum", "genus")
  pgPart = getParts(taxoPhylumGenus, phylaMean, genusMean)
  pgCor = getCors(taxoPhylumGenus, phylaMat, genusMat)
  
  pgPartMelt = melt(pgPart)
  pgCorMelt = melt(pgCor)
  
  pgStat = cbind(pgPartMelt$value, pgCorMelt$value)
  rownames(pgStat) = pgPartMelt$L2
  colnames(pgStat) = c("","")
  pgStat[is.na(pgStat)]=0
  
  corrplot(pgStat[1:23,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.5,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  taxoPhylumFamily = parseTaxoBy(taxo, "phylum", "family")
  pfPart = getParts(taxoPhylumFamily, phylaMean, familyMean)
  pfCor = getCors(taxoPhylumFamily, phylaMat, familyMat)
  
  pfPartMelt = melt(pfPart)
  pfCorMelt = melt(pfCor)
  
  pfStat = cbind(pfPartMelt$value, pfCorMelt$value)
  rownames(pfStat) = pfPartMelt$L2
  colnames(pfStat) = c("","")
  pfStat[is.na(pfStat)]=0
  
  corrplot(pfStat[1:23,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.5,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  corrplot(pfStat[24:53,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.5,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  corrplot(pfStat[54:78,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.5,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  require(corrplot)
  corrplot(pfStat[79:90,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.5,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  
  
  corrplot(pcStat,
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.8,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  corrplot((pPart),
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.7,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  
  corrplot(coStat,
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.7,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
  ff = ofPartMelt[ofPartMelt$L1=="Bacteroidales","L2"]
  ff = ofPartMelt[ofPartMelt$L1=="Lactobacillales","L2"]
  ff = ofPartMelt[ofPartMelt$L1=="Enterobacterales","L2"]
  
  corrplot(ofStat[ff,],
           tl.col="black",
           tl.srt = 45,
           tl.cex = 0.7,
           cl.pos = "n",
           col=c("white", "black"),bg="lightblue")
  
}

getCorLink <- function(taxo, mergeMat, exUnclassified = F) { 
  parsedTaxo = parseTaxo(taxo)
  
  phylaMat = getTaxaSumMat(mergeMat,taxo,"phylum",T)
  classMat = getTaxaSumMat(mergeMat,taxo,"class",T)
  orderMat = getTaxaSumMat(mergeMat,taxo,"order",T)
  familyMat = getTaxaSumMat(mergeMat,taxo,"family",T)
  genusMat = getTaxaSumMat(mergeMat,taxo,"genus",T)
  speciesMat  = getTaxaSumMat(mergeMat,taxo,"species",T)
  
  pcCor = getCors(parsedTaxo, phylaMat, classMat, "phylum")
  coCor = getCors(parsedTaxo, classMat, orderMat, "class")
  ofCor = getCors(parsedTaxo, orderMat, familyMat, "order")
  fgCor = getCors(parsedTaxo, familyMat, genusMat, "family")
  gsCor = getCors(parsedTaxo, genusMat, speciesMat, "genus")
  
  pcCorMelt = melt(pcCor)
  coCorMelt = melt(coCor)
  ofCorMelt = melt(ofCor)
  fgCorMelt = melt(fgCor)
  gsCorMelt = melt(gsCor)
  
  
  corrplot(cbind(pcMelt$value, pcCorMelt$value))
  
  pcCor = lapply(seq_along(parsedTaxo$phylum), function(y,n,i) {
    parent= n[[i]]
    child=y[[i]]
    return(sapply(child,function(x) cor(phylaMat[parent,],classMat[x,],method="spearman")))
    
  } , y=parsedTaxo$phylum, n=names(parsedTaxo$phylum) )
  names(pcCor) = names(parsedTaxo$phylum)
  
  
  
}

getPhylaName <- function(mgsName, taxo) {
  phyla = taxo[match(mgsName, rownames(taxo)),"phylum"]
  return(phyla)
}
getFamilyName <- function(mgsName, taxo) {
  family = taxo[match(mgsName, rownames(taxo)),"family"]
  return(family)
}

getGenusName <- function(mgsName, taxo) {
  genus = taxo[match(mgsName, rownames(taxo)),"genus"]
  return(genus)
}

getSpeciesName <- function(mgsName, taxo) {
  species = taxo[match(mgsName, rownames(taxo)),"species"]
  return(species)
}

makeCorMat <- function(exprMat, mode="pearson") {
  if (mode == "pearson") {
    corMat = cor(t(exprMat))
  }
  if (mode == "spearman") {
    corMat = cor(t(exprMat), method = "spearman")
  }
  return(corMat)
}

makeCorTable <- function(exprMat, cutOff=0.99, mode="pearson", upperMode = F, self=F, debug=F) {
  minValue = -2
  
  if (mode == "pearson") {
    corMat = cor(t(exprMat))
  }
  if (mode == "incomplete") {
    corMat = cor(t(exprMat), use="pairwise.complete.obs")
  }
  if (mode == "spearman") {
    corMat = cor(t(exprMat), method = "spearman")
  }
  if (mode == "r2") {
    corMat = cor(t(exprMat))
    corMat = corMat^2
  }
  if (mode == "jcc") {
    corMat = getJaccMat(exprMat)
  }
  
  ### checking NA
  if (debug) print("na amount:")
  if (debug) print(table(is.na(corMat)))
  
  corMat[is.na(corMat)]=minValue
  corVec = corMat[lower.tri(corMat)]
  
  cutOffUpper = as.numeric( quantile(corVec, cutOff, na.rm=T) )
  cutOffLower = as.numeric( quantile(corVec[corVec!=minValue], 1-cutOff, na.rm=T) )
  cutOff = as.numeric( quantile(corMat[lower.tri(corMat)], cutOff, na.rm=T) )
  
  if(debug) print("quantile: ")
  if(debug) print(cutOffUpper)
  if(debug) print(cutOffLower)
  
  resTable = melt(corMat)[,1:3] 
  resTableNew = as.data.table(resTable)
  setkey(resTableNew, value)
  resTableNew = resTableNew[order(value, decreasing = T)]
  resTable = as.data.frame(resTableNew)
  
  if (!self) {
    resTable = resTable[resTable[,1]!=resTable[,2],]
  }
  
  if(debug) print("cut-off before")
  if(debug) print(dim(resTable))
  
  if (upperMode) resTable = resTable[resTable[,3]>=cutOff,]
  if (!upperMode) resTable = resTable[resTable[,3]>=cutOffUpper| resTable[,3]<=cutOffLower,]
  if(debug) print("cut-off on-the-way")
  if(debug) print(dim(resTable))
  
  resTable = resTable[!is.na(resTable[,3]),]
  resTable = resTable[resTable[,3]!=minValue,]
  
  if(debug) print("cut-off after")
  if(debug) print(dim(resTable))
  
  return(resTable)
}

makeCorNet <- function(corTable) {
  corNet = igraph::graph.data.frame(corTable[,1:2], directed=F)
  corNet = igraph::simplify(corNet, remove.multiple=T, remove.loops=T)
  return(corNet)
}

makeModuleList <- function(corNet, debug=F) {
  fc = cluster_walktrap(corNet)
  cluster = fc$membership
  geneCluster = data.frame(gene=V(corNet)$name, 
                           cluster=cluster,
                           stringsAsFactors=F)
  geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  geneClusterList = lapply(geneClusterList, "[[", "gene")
  geneClusterSizes = do.call(c, lapply(geneClusterList, length))
  
  if (debug) {
    print("... done")
    print(paste("... modularity:", as.character(modularity(fc))))
    print(paste("... no. clusters:", as.character(length(geneClusterList))))
    print(paste("... no. genes of max cluster:", as.character(sort(geneClusterSizes,T)[1])))
    
  }
  return(geneClusterList)
}

makeJaccNet <- function(jaccMat, cutJacc = 0.5) {
  jaccElementsAll = rownames(jaccMat)
  jaccTable = melt(jaccMat)
  jaccTable = jaccTable[jaccTable[,3] >= cutJacc & jaccTable[,1] != jaccTable[,2], ]
  jaccTable[,1]=as.character(jaccTable[,1])
  jaccTable[,2]=as.character(jaccTable[,2])
  
  jaccElementsLeft = unique(c(jaccTable[,1], jaccTable[,2]))
  jaccElementsExcluded = jaccElementsAll[!jaccElementsAll %in% jaccElementsLeft]
  
  jaccNet = igraph::graph.data.frame(jaccTable[,1:2], directed=F)
  jaccNet = igraph::simplify(jaccNet, remove.multiple=T, remove.loops=T)
  
  jaccNet = add_vertices(jaccNet, length(jaccElementsExcluded), attr=list(name=jaccElementsExcluded))
  
  return(jaccNet)
}


getClusterScore  <- function(corNet, currClusterGenes) {
  netGenes=V(corNet)$name
  currInd = match(currClusterGenes, netGenes)
  currGraph = subgraph(corNet, currInd)
  cc = transitivity(currGraph, "global")
  return(cc)
}

getClusterScores <- function(corNet, geneClusterList) {
  ccScores = lapply(geneClusterList, function(currClusterGenes){
    #currGraph = igraph::subgraph(corNet, currClusterGenes)
    currGraph = subgraph(corNet, currClusterGenes)
    cc = transitivity(currGraph, "global")
    return(cc)
  })
  ccScores = unlist(ccScores)
  names(ccScores) = names(geneClusterList)
  return(ccScores)
}

getModuleMatrix <- function(corNet, modulePruned, debug=F) {
  numEdges = ecount(corNet)
  adjMat = as_adj(corNet)
  degMat = getDegMat(corNet)
  
  numModules = length(modulePruned)
  if (debug) print(numModules)
  moduleMat = matrix(0, numModules, numModules)
  for (i in 1:(numModules-1)) {
    if (debug) print(i)
    for (j in (i+1):numModules) {
      moduleMat[i,j] = areModuleInteract(adjMat, degMat, numEdges, modulePruned[[i]], modulePruned[[j]])
    }
  }  
  moduleMat[lower.tri(moduleMat)] = t(moduleMat)[lower.tri(moduleMat)]
  rownames(moduleMat) = names(modulePruned)
  colnames(moduleMat) = names(modulePruned)
  return(moduleMat)
}

areModuleInteract <- function(adjMat, degMat, numEdges, genesA, genesB) {
  observedEdges = observedEdgesBtw(adjMat, genesA, genesB)
  expectedEdges = expectedEdgesBtw(degMat, numEdges, genesA, genesB)
  diff=(observedEdges-expectedEdges)/expectedEdges
  return(diff)
}
getModulePruned <- function(geneClusterList, cutCluster, debug=F) {
  geneClusterSizes = unlist(lapply(geneClusterList, length))
  names(geneClusterSizes) = names(geneClusterSizes)
  if (debug) print(table(geneClusterSizes>=cutCluster))
  geneClusterCut = names(geneClusterSizes[geneClusterSizes>=cutCluster])
  geneClusterList = geneClusterList[geneClusterCut]
  return(geneClusterList) 
}

getDegMat <- function(igraphNet) {
  #we regard network as undirected
  #degs = degree(igraphNet, mode="all")
  degs = igraph::degree(igraphNet, mode="total")
  names(degs) = V(igraphNet)$name
  degMat = as.matrix(degs) %*% t(as.matrix(degs))
  return(degMat)
}
expectedEdgesBtw <- function(degMat, numEdges, genesA, genesB) {
  degMatSel = degMat[genesA, genesB]/(2*numEdges)
  count = sum(degMatSel)
  return(count)
}
observedEdgesBtw <- function(adjMat, genesA, genesB) {
  adjMatSel = adjMat[genesA, genesB]
  count = sum(adjMatSel)
  return(count)
}
annotateModulesByCC <-function(corNet, moduleList,  cutCluster=3, cutCC = 0.5, debug=F) {
  
  modulePrunedList = getModulePruned(moduleList, cutCluster)
  moduleNames = names(modulePrunedList)
  
  # module matrix
  moduleMat = getModuleMatrix(corNet, modulePrunedList)
  #print(moduleMat)
  # first axis
  moduleCC = getClusterScores(corNet, modulePrunedList)
  highCCModuleNames = names(moduleCC[moduleCC >= cutCC])
  
  
  # summary of two axes
  summaryTypes = rep("None", length(modulePrunedList))
  names(summaryTypes) = moduleNames 
  summaryTypes[moduleNames %in% highCCModuleNames ] = "HighCC"
  summaryTypes[!moduleNames %in% highCCModuleNames ] = "LowCC"
  
  moduleCut=1
  moduleMelt = melt(moduleMat)
  moduleMelt = moduleMelt[moduleMelt[,3]>moduleCut,]
  moduleMelt[] = lapply(moduleMelt, as.character)
  modulesShown = unique.default(c(moduleMelt[,1], moduleMelt[,2]))
  moduleSizes = unlist(lapply(modulePrunedList, length))
  moduleAttr = data.frame(node=moduleNames, size=moduleSizes, summary = summaryTypes, cc = moduleCC, stringsAsFactors = F)
  
  modulesNotShown = moduleNames[!moduleNames %in% modulesShown]
  return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown = modulesNotShown ))
}


write.edge.cytoscape <- function(edgeTable, nodesNotShown, outFile) {
  lenTable=dim(edgeTable)[1]
  lenNodesNotShown = length
  print(head(edgeTable))
  sink(outFile, append=F)
  cat("source")
  cat("\t")
  cat("target")
  cat("\n")
  for (i in 1:lenTable) {
    # ints = unlist(edgeTable[i,])
    cat(edgeTable[i,1])
    cat("\t")
    cat(edgeTable[i,2])
    cat("\n")
  }
  for (node in nodesNotShown) {
    cat(node)
    cat("\n")
  }
  sink() 
}
write.node.cytoscape <- function(nodeTable, outFile) {
  write.table(nodeTable, outFile, row.names = F, quote = F, sep="\t")
}
writeModuleNetworkAnnot <- function(moduleOut, nodeFile, edgeFile) {
  nodeTable = moduleOut[["nodeTable"]]
  edgeTable = moduleOut[["edgeTable"]]
  nodesNotShown = moduleOut[["nodesNotShown"]]
  write.node.cytoscape(nodeTable, nodeFile)
  write.edge.cytoscape(edgeTable, nodesNotShown, edgeFile)
}
