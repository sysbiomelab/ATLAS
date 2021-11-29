#save.image(file="workspace.2019.08.13.RData")
#load(file="workspace.2019.08.13.RData")

#source("C:\\Data/comparative.analysis.healthy.sweden/microbiome_atlas_comparative/check.mgs.functions.r")
source("check.mgs.functions.r")

loadLibrary=T
if (loadLibrary) {
  ###lmer test
  library(MuMIn)
  
  
  require(dplyr)
  #require(h20)
  require(betapart)
  require(viridis)
  require(colourpicker)
  require(circlize)
  require(markovchain)
  require(diagram)
  require(ggtern)
  require(survminer)
  require(survival)
  require(qqman)  
  require(igraph)
  require(momr)
  require(monocle)
  require(ppcor)
  library(ggplot2)
  library(zoo)
  require(reshape2)
  require(data.table)
  require(matrixStats)
  require(gplots)
  require(lme4)
  require(ade4)
  require(pls)
  require(car)
  require(lmerTest)
  require(corrplot)
  require(vegan)
  require(ape)
  require(Rtsne)
  require(RColorBrewer)
  require(randomForest)
  require(ropls)
  require(BiocManager)
  require(plotrix)
  #require(pROC)
  require(ROCR)
  require(networkD3)
  require(tidyverse)
  require(RColorBrewer)
  require(ggrepel)
  require(prabclus)
  require(dbscan)  
  require(mclust)
  require(DirichletMultinomial)
  require(densityClust)
  
  #install("C:/Users/SunjaeLee/Downloads/rPython/")
  #install.packages("http://genome.crg.es/~didac/ggsunburst/ggsunburst_0.0.10.tar.gz", repos=NULL, type="source")
  
  setwd("C:\\Data/comparative.analysis.healthy.sweden/")
  
  gutMsp1992Data = "C:\\Data/comparative.analysis.healthy.sweden/hs_10.4_1992_MSP_freeze2_20180905.RData"
  mgsRData = "C:\\Data/comparative.analysis.healthy.sweden/mgs.all.10M.RData"
  mgsCleanRData = "C:\\Data/comparative.analysis.healthy.sweden/mgs.clean.10M.RData"
  mgsClean2RData = "C:\\Data/comparative.analysis.healthy.sweden/mgs.clean.madagascar.10M.RData"
  mgsCleanUpdatedRData = "C:\\Data/comparative.analysis.healthy.sweden/mgs.clean.updated.20190501.10M.RData"
  
  basicMetaCleanRData = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.20190109.RData"
  basicMetaCleanRData2 = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.20190114.RData"
  basicMetaCleanRData3 = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.updated.20190501.RData"
  #basicMetaCleanRDataUpdated = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.madagascar.20190125.RData"
  #commented at 2020-06-19
  basicMetaCleanRDataUpdated2 = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.behcet.20190805.RData"
  
  hqDataFile="C:/Data/comparative.analysis.healthy.sweden/high.qualtiy.metadata.tables.txt"
  hqDataRData="C:/Data/comparative.analysis.healthy.sweden/high.qualtiy.metadata.tables.RData"
  load(hqDataRData)
  #hqDataTab = read.table(hqDataFile, header=T, sep="\t", stringsAsFactors = F)
  datasetOverviewFile = "C:/Data/comparative.analysis.healthy.sweden/atlas.dataset.overview.txt"
  datasetOverviewTab = read.delim(datasetOverviewFile, header=T, sep="\t", stringsAsFactors = F)
  keggPathDescTabRData = "C:/Data/comparative.analysis.healthy.sweden/keggPathDescTab.RData"
  #moduleDescTabRData = "C:/Data/comparative.analysis.healthy.sweden/moduleDescTab.RData"
  moduleDescTabRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/KEGG module/moduleDescTab.RData" 
  
  phageomeData = "C://Data/comparative.analysis.healthy.sweden/mergePhageMatFinal.RData"
  
  load(hqDataRData)
  load(phageomeData)
  
  load(keggPathDescTabRData)
  load(moduleDescTabRData)
  moduleDescTab$FullName = paste(moduleDescTab$moduleId, moduleDescTab$moduleName, sep = ":")
  
  exModules = moduleDescTab$FullName[moduleDescTab$system =="Module set"]
  moduleDescTab = moduleDescTab[moduleDescTab$system !="Module set",]
  
  #keggPathDescTab
  #load("colon.mgs.med.vec.fecies.RData")
  
  # pp = getTaxaAvgMat(mgs_med_vec,taxo,"phylum")
  # pp2 = rowMedians(pp)
  # names(pp2)=rownames(pp)
}

loadFunc = T
if (loadFunc) {
  
  #### general functions ####
  getAbSubj <- function(str) {
    subj = strsplit(str, split = "_")[[1]][1]
    return(subj)
  }
  getAbSubjs <- function(strVec) {
    return(sapply(strVec, getAbSubj))
  }
  
  getVisitIds <- function(str) {
    visitId = strsplit(str, split = "_")[[1]][1]
    return(visitId)
  }
  
  getFields <- function(str) {
    str = gsub("visit1_","",str)
    str = gsub("visit2_","",str)
    str = gsub("visit3_","",str)
    str = gsub("visit4_","",str)
    return(str)
  }
  
  getRelAbd <- function(targetMat) {
    relAbds = sweep(targetMat, MARGIN=2, colSums(targetMat), FUN="/")
    return(relAbds)
  }
  jacc <- function(x,y) {
    prod = sum(x*y)
    sumx = sum(x)
    sumy = sum(y)
    uni  = sumx + sumy - prod
    return(prod/uni)
  }
  
  checkOverlap <- function(xList, yList) {
    xLen = length(xList)
    yLen = length(yList)
    xy = sum((xList %in% yList))
    print(c(xLen, yLen, xy, xy/(xLen+yLen-xy)))
  }
  # contourMat = pcoaMgsMat
  # target1 = newCaseByEnteroList$ETF
  # target2 = newContByEnteroList$ETF
  # 
  getAnosimBoth <- function(mgsMat, target1, target2) {
    targetAll = unique(c(target1, target2))
    vegOut=vegdist(t(mgsMat[,targetAll]),"bray")
    targetClasses = rep("class1",length(targetAll))
    targetClasses[targetAll %in% target2] = "class2"
    set.seed(1)
    anosimOut = anosim(vegOut, targetClasses)
    print(anosimOut)
    return(anosimOut)
  }
  getPermanovaBoth <- function(mgsMat, target1, target2) {
    targetAll = unique(c(target1, target2))
    vegOut=vegdist(t(mgsMat[,targetAll]),"bray")
    targetClasses = rep("class1",length(targetAll))
    targetClasses[targetAll %in% target2] = "class2"
    set.seed(1)
    permanovaOut = adonis2(vegOut ~ targetClasses, permutations = 999, method="bray")
    print(permanovaOut)
    return(permanovaOut)
  }
  
  
  plotContourBoth <- function(contourMat, target1, target2) {
    targetAll = unique(c(target1, target2))
    k <- 9
    my.cols.1 <- (brewer.pal(k, "PuRd"))
    my.cols.2 <- (brewer.pal(k, "Blues"))
    
    
    targetInd1 = rownames(contourMat)%in% target1
    targetInd2 = rownames(contourMat)%in% target2
    
    targetMat1 = contourMat[targetInd1,]
    targetMat2 = contourMat[targetInd2,]
    targetMatAll = rbind(targetMat1, targetMat2)
    
    print(summary(out))
    z.1 <- kde2d(targetMat1[,1], targetMat1[,2], n=50)
    z.2 <- kde2d(targetMat2[,1], targetMat2[,2], n=50)
    
    plot(contourMat[targetAll,], xlab="PCoA1", ylab="PCoA2", col="#80808033", pch=19, cex=.4, xlim=c(-1,1),ylim=c(-1,1))
    contour(z.1, drawlabels=FALSE, nlevels=k, col=my.cols.1, add=TRUE)
    contour(z.2, drawlabels=FALSE, nlevels=k, col=my.cols.2, add=TRUE)
    
    abline(h=0, v=0, col="black")
    return(list(z.1, z.2))
    
  }
  
  plotContour <- function(contourMat) {
    k <- 11
    my.cols <- rev(brewer.pal(k, "RdYlBu"))
    
    z <- kde2d(contourMat[,1], contourMat[,2], n=50)
    
    plot(contourMat, xlab="X label", ylab="Y label", col="#80808033", pch=19, cex=.4, xlim=c(-1,1),ylim=c(-1,1))
    contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
    
    abline(h=0, v=0, col="black")
    return(z)
    
  }
  
  getNumbersFromWellnessVisits <- function(visits) {
    visits = gsub("v","",visits)
    visits = as.numeric(visits)
    return(visits)
  }
  
  
  getChiPval <- function(chisq, df) {
    return(1 - pchisq(chisq,df))
  }
  
  writeJaccNetToCytoscape <- function(samples, jaccMat, metaTab, nodeFile, edgeFile, jaccCut=0.5) {
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector = c(col_vector,col_vector)
    
    currJaccMat = jaccMat[samples, samples]
    currJaccTab = melt(currJaccMat, variable.factor=F)
    currJaccTab$Var1 = as.character(currJaccTab$Var1)
    currJaccTab$Var2 = as.character(currJaccTab$Var2)
    currJaccTab = currJaccTab[currJaccTab$value >= jaccCut,]
    currJaccTab = currJaccTab[currJaccTab$Var1 != currJaccTab$Var2,]
    
    nodes = unique(c(currJaccTab$Var1, currJaccTab$Var2))
    subjects= metaTab$metadata.ID[match(samples,metaTab$sample.ID)]
    subtypes= metaTab$subtype[match(samples,metaTab$sample.ID)]
    types= metaTab$type[match(samples,metaTab$sample.ID)]
    
    nodesNotShown = samples[!samples %in% nodes]
    nodeTable = data.frame(node=samples,
                           subject=subjects,
                           type = types,
                           subtype = subtypes)
    edgeTable = currJaccTab 
    
    outList = list(nodeTable = nodeTable, 
                   edgeTable = edgeTable, 
                   nodesNotShown = nodesNotShown )
    
    writeModuleNetworkAnnot(outList, nodeFile, edgeFile)
    return(outList)
  }
  
  copyLowerToUpper <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  getJaccMat <- function(mat) {
    
    jaccMat = matrix(0,dim(mat)[2],dim(mat)[2])
    rownames(jaccMat) = colnames(mat)
    colnames(jaccMat) = colnames(mat)
    
    combs = t(combn(ncol(mat),2))
    jaccs = apply(combs, 1, function(x) jacc(mat[,x[1]], mat[,x[2]]))
    jaccMat[combs]=jaccs
    diag(jaccMat)=1
    jaccMat = copyLowerToUpper(jaccMat)
    return(jaccMat)
  }
  
  getJaccFrom <- function(strs1, strs2) {
    len1 = length(strs1)
    len2 = length(strs2)
    if (len1 == 0 | len2 ==0 ) return(NA)
    lenO = sum(strs1 %in% strs2)
    lenU = len1 + len2 - lenO
    return(lenO/lenU)
  }
  
  #### disease network ####
  makeSpeciesDiseaseNet <- function(speciesTab) {
    diseaseNet = igraph::graph.data.frame(speciesTab[,1:2], directed=F)
    diseaseNet = set_edge_attr(diseaseNet, "weight", value= speciesTab[,3])
    diseaseNet = igraph::simplify(diseaseNet, remove.multiple=T, remove.loops=T)
    return(diseaseNet)
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
  
  
  makeModuleList <- function(jaccNet, weightMode = F, debug=F) {
    fc = cluster_walktrap(jaccNet)
    if (weightMode) fc = cluster_walktrap(jaccNet, weights = E(jaccNet)$weight)
    
    #fc = cluster_edge_betweenness(jaccNet)
    #fc = cluster_louvain(jaccNet)
    
    cluster = fc$membership
    geneCluster = data.frame(gene=V(jaccNet)$name, 
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
  
  
  getClusterScore  <- function(jaccNet, currClusterGenes) {
    netGenes=V(jaccNet)$name
    currInd = match(currClusterGenes, netGenes)
    currGraph = subgraph(jaccNet, currInd)
    cc = transitivity(currGraph, "global")
    return(cc)
  }
  getClusterScores <- function(jaccNet, geneClusterList) {
    ccScores = lapply(geneClusterList, function(currClusterGenes){
      #currGraph = igraph::subgraph(corNet, currClusterGenes)
      currGraph = subgraph(jaccNet, currClusterGenes)
      cc = transitivity(currGraph, "global")
      if (is.nan(cc)) cc = 0
      return(cc)
    })
    ccScores = unlist(ccScores)
    names(ccScores) = names(geneClusterList)
    return(ccScores)
  }
  
  getModuleMatrix <- function(jaccNet, modulePruned, debug=F) {
    numEdges = ecount(jaccNet)
    adjMat = as_adj(jaccNet)
    degMat = getDegMat(jaccNet)
    
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
  
  prunListByNetworkGenesAndCut <- function(modList, networkGenes, clusterCut) {
    outList = lapply(modList, function(modGenes){
      modGenes = modGenes[modGenes %in% networkGenes]
    })
    outListLens = lapply(outList, length)
    outList = outList[outListLens>=clusterCut]
    return(outList)
    
  }
  areModuleInteract <- function(adjMat, degMat, numEdges, genesA, genesB) {
    observedEdges = observedEdgesBtw(adjMat, genesA, genesB)
    expectedEdges = expectedEdgesBtw(degMat, numEdges, genesA, genesB)
    diff=(observedEdges-expectedEdges)/expectedEdges
    return(diff)
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
  getModulePruned <- function(geneClusterList, cutCluster, debug=F) {
    geneClusterSizes = unlist(lapply(geneClusterList, length))
    names(geneClusterSizes) = names(geneClusterSizes)
    if (debug) print(table(geneClusterSizes>=cutCluster))
    geneClusterCut = names(geneClusterSizes[geneClusterSizes>=cutCluster])
    geneClusterList = geneClusterList[geneClusterCut]
    return(geneClusterList) 
  }
  getInteractionsFromModules <- function( coexpTable, moduleList, targets, debug=T ) {
    moduleNames = names(moduleList)
    if (!all(targets %in% moduleNames)) {
      print("... some targets were not shown in modulePruned")
      print("... please check it again")
      return(NULL)
    }
    
    targetModules = moduleList[targets]
    targetModuleGenes = unique.default(do.call(c, targetModules))
    targetModuleMaps = rep(NA, length(targetModuleGenes))
    for (i in 1:length(targets)) {
      #currTarget = targets[1]
      currTarget = targets[i]
      currTargetGenes = unlist(moduleList[[currTarget]])
      targetModuleMaps[targetModuleGenes %in% currTargetGenes] = currTarget
    }
    
    targetInd = coexpTable[,1]%in%targetModuleGenes & coexpTable[,2]%in%targetModuleGenes
    targetEdgeTable = coexpTable[targetInd, 1:3]
    colnames(targetEdgeTable) = c("Source", "Target","Corr")
    
    targetNodeTable = data.frame(node=targetModuleGenes,
                                 symbol=getSymbols(targetModuleGenes),
                                 moduleType=targetModuleMaps,
                                 stringsAsFactors=F)
    
    
    return(list(nodeTable = targetNodeTable, edgeTable = targetEdgeTable))
  }
  
  annotateModulesByCC <-function(jaccNet, moduleList,  cutCluster=5, cutCC = 0.5, debug=F) {
    
    modulePrunedList = getModulePruned(moduleList, cutCluster)
    moduleNames = names(modulePrunedList)
    
    # module matrix
    moduleMat = getModuleMatrix(jaccNet, modulePrunedList)
    print(moduleMat)
    # first axis
    moduleCC = getClusterScores(jaccNet, modulePrunedList)
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
    moduleMelt = moduleMelt[!is.na(moduleMelt$value),]
    modulesShown = unique.default(c(moduleMelt[,1], moduleMelt[,2]))
    moduleSizes = unlist(lapply(modulePrunedList, length))
    moduleAttr = data.frame(node=moduleNames, size=moduleSizes, summary = summaryTypes, cc = moduleCC, stringsAsFactors = F)
    
    modulesNotShown = moduleNames[!moduleNames %in% modulesShown]
    return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown = modulesNotShown ))
  }
  
  getEnrichMentGsByList <- function(refGeneSet, targetGenes, debug=F) {
    allRefGenes = unique.default(do.call(c, refGeneSet))
    if (debug) print(table(targetGenes %in% allRefGenes))
    targetGenes = targetGenes[targetGenes %in% allRefGenes]
    
    hyperPs = lapply(refGeneSet, function(currRefGenes) {
      hyperP = getHyperP( targetGenes, currRefGenes, allRefGenes, F ) 
      return(hyperP)
    })
    hyperPs = unlist(hyperPs, use.names = F)
    names(hyperPs) = names(refGeneSet)
    return(hyperPs)
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
  
  
  
  #### Enterotyping ####
  
  getSubectIds <- function(strs) {
    newStrs = gsub("_v[1-4]","",strs)
    return(newStrs)
  }
  getEntCount <- function(enterotypes) {
    etf <- sum(enterotypes == "ET-Firmicutes")
    etb <- sum(enterotypes == "ET-Bacteroides")
    etp <- sum(enterotypes == "ET-Prevotella")
    return(c(etf, etb, etp))
  }
  getEntCountFrom <- function(basicMetaMap, cut0, cut1, by="GeneRichness") {
    entCount = getEntCount(basicMetaMap$enteroType[basicMetaMap[,by] <= cut1 & basicMetaMap[,by] > cut0])
    return(entCount)
  }
  getSurvivalRatesFrom <- function(survivalStats, cut0, cut1, by="mgs.abd") {
    survRates = survivalStats$surv.rates[survivalStats[,by] <= cut1 & survivalStats[,by] > cut0]
    return(survRates)
  }
  
  getSelectedStatFrom <- function(dataTab, cut0, cut1, by, target) {
    selected = dataTab[dataTab[,by] <= cut1 & dataTab[,by] > cut0, target]
    return(selected)
  }
  
  
  getMarkovStatFrom <- function(markovStatTab, cut0, cut1, by="meanAbd", target="inflow") {
    targetRates = markovStatTab[markovStatTab[,by] <= cut1 & markovStatTab[,by] > cut0, target]
    return(targetRates)
  }
  getBinTab <- function(basicMetaMap, by="GeneRichness") {
    richnessMin = min(basicMetaMap[,by])
    richness01 = quantile(basicMetaMap[,by], 0.1)
    richness02 = quantile(basicMetaMap[,by], 0.2)
    richness03 = quantile(basicMetaMap[,by], 0.3)
    richness04 = quantile(basicMetaMap[,by], 0.4)
    richness05 = quantile(basicMetaMap[,by], 0.5)
    richness06 = quantile(basicMetaMap[,by], 0.6)
    richness07 = quantile(basicMetaMap[,by], 0.7)
    richness08 = quantile(basicMetaMap[,by], 0.8)
    richness09 = quantile(basicMetaMap[,by], 0.9)
    richnessMax = max(basicMetaMap[,by])
    
    entTab00 = getEntCountFrom(basicMetaMap, richnessMin-1, richness01, by)
    entTab01 = getEntCountFrom(basicMetaMap, richness01, richness02, by)
    entTab02 = getEntCountFrom(basicMetaMap, richness02, richness03, by)
    entTab03 = getEntCountFrom(basicMetaMap, richness03, richness04, by)
    entTab04 = getEntCountFrom(basicMetaMap, richness04, richness05, by)
    entTab05 = getEntCountFrom(basicMetaMap, richness05, richness06, by)
    entTab06 = getEntCountFrom(basicMetaMap, richness06, richness07, by)
    entTab07 = getEntCountFrom(basicMetaMap, richness07, richness08, by)
    entTab08 = getEntCountFrom(basicMetaMap, richness08, richness09, by)
    entTab09 = getEntCountFrom(basicMetaMap, richness09, richnessMax+1, by)
    
    entTabByRichnes = rbind(entTab00,
                            entTab01,
                            entTab02,
                            entTab03,
                            entTab04,
                            entTab05,
                            entTab06,
                            entTab07,
                            entTab08,
                            entTab09)
    entTabByRichnes = sweep(entTabByRichnes, MARGIN=1, rowSums(entTabByRichnes)/100, FUN="/")
    return(entTabByRichnes)
  }
  getBinTabForSurvival <- function(survivalStats, by="mgs.abd") {
    abdMin = min(survivalStats[,by])
    abd02 = quantile(survivalStats[,by], 0.2)
    abd04 = quantile(survivalStats[,by], 0.4)
    abd06 = quantile(survivalStats[,by], 0.6)
    abd08 = quantile(survivalStats[,by], 0.8)
    abdMax = max(survivalStats[,by])
    
    survTab00 = getSurvivalRatesFrom(survivalStats, abdMin-1, abd02, by)
    survTab02 = getSurvivalRatesFrom(survivalStats, abd02, abd04, by)
    survTab04 = getSurvivalRatesFrom(survivalStats, abd04, abd06, by)
    survTab06 = getSurvivalRatesFrom(survivalStats, abd06, abd08, by)
    survTab08 = getSurvivalRatesFrom(survivalStats, abd08, abdMax+1, by)
    
    survListByAbd = list(per20=survTab00,
                         per40=survTab02,
                         per60=survTab04,
                         per80=survTab06,
                         per100=survTab08)
    
    return(survListByAbd)
  }
  getBinTabForSurvival2 <- function(survivalStats, by="mgs.abd") {
    abdMin = min(survivalStats[,by])
    abd02 = quantile(survivalStats[,by], 0.25)
    abd04 = quantile(survivalStats[,by], 0.5)
    abd06 = quantile(survivalStats[,by], 0.75)
    abdMax = max(survivalStats[,by])
    
    print(abdMin)
    print(abd02)
    print(abd04)
    print(abd06)
    print(abdMax)
    
    survTab00 = getSurvivalRatesFrom(survivalStats, abdMin-1, abd02, by)
    survTab02 = getSurvivalRatesFrom(survivalStats, abd02, abd04, by)
    survTab04 = getSurvivalRatesFrom(survivalStats, abd04, abd06, by)
    survTab06 = getSurvivalRatesFrom(survivalStats, abd06, abdMax+1, by)
    
    survListByAbd = list(per25=survTab00,
                         per50=survTab02,
                         per75=survTab04,
                         per100=survTab06)
    
    return(survListByAbd)
  }
  getRichnessBinByOutflowMean <- function(dataTab, by="outflow", target="richness") {
    qMin = min(outflowTab[,by])
    q01 = quantile(outflowTab[,by], 0.1)
    q02 = quantile(outflowTab[,by], 0.2)
    q03 = quantile(outflowTab[,by], 0.3)
    q04 = quantile(outflowTab[,by], 0.4)
    q05 = quantile(outflowTab[,by], 0.5)
    q06 = quantile(outflowTab[,by], 0.6)
    q07 = quantile(outflowTab[,by], 0.7)
    q08 = quantile(outflowTab[,by], 0.8)
    q09 = quantile(outflowTab[,by], 0.9)
    qMax = max(outflowTab[,by])
    
    selected00 = getSelectedStatFrom(outflowTab, qMin-1, q01, by, target)
    selected01 = getSelectedStatFrom(outflowTab, q01, q02, by, target)
    selected02 = getSelectedStatFrom(outflowTab, q02, q03, by, target)
    selected03 = getSelectedStatFrom(outflowTab, q03, q04, by, target)
    selected04 = getSelectedStatFrom(outflowTab, q04, q05, by, target)
    selected05 = getSelectedStatFrom(outflowTab, q05, q06, by, target)
    selected06 = getSelectedStatFrom(outflowTab, q06, q07, by, target)
    selected07 = getSelectedStatFrom(outflowTab, q07, q08, by, target)
    selected08 = getSelectedStatFrom(outflowTab, q08, q09, by, target)
    selected09 = getSelectedStatFrom(outflowTab, q09, qMax+1, by, target)
    
    selctedBins = list(per10=selected00,
                       per20=selected01,
                       per30=selected02,
                       per40=selected03,
                       per50=selected04,
                       per60=selected05,
                       per70=selected06,
                       per80=selected07,
                       per90=selected08,
                       per100=selected09)
    
    return(selctedBins)
  }
  
  makeKeggQuery <- function(pathList, mspKeggList, pathId, mspId, viewMode = F) {
    koList = mspKeggList[[mspId]]
    prefix1 = "https://www.kegg.jp/kegg-bin/show_pathway?map="
    prefix2 = "&multi_query=7165+red,blue"
    suffices = unlist(sapply(koList, function(ko) {
      suffix1 = "%0d%0a"
      suffix2 = "+red"
      concats = paste(suffix1, ko, suffix2, sep="")
      return(concats)
    }))
    suffices = paste(suffices, collapse = "")
    finalUrl = paste(prefix1, pathId, prefix2, suffices, sep = "")
    if (viewMode) browseURL(finalUrl)
    return(finalUrl)  
  }
  
  makeCurrTab <- function(currRow) {
    currTab = data.frame(case=currRow[1:3], cont=currRow[4:6])
    rownames(currTab) = c("ET-Firmicutes", "ET-Bacteroides", "ET-Prevotella")
    return(currTab)
  }
  grepEnteroTypesFrom <- function(sampleList, basicMetaMap) {
    newEnteroList = lapply(sampleList, function(currSamples){
      enterotypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% currSamples]
      return(enterotypes)
    })
    return(newEnteroList)
  }
  
  grepGeneRichnesssFrom <- function(sampleList, basicMetaMap) {
    newRichnessList = lapply(sampleList, function(currSamples){
      richness = basicMetaMap$GeneRichness[basicMetaMap$sample.ID %in% currSamples]
      return(richness)
    })
    return(newRichnessList)
  }
  
  grepGeneRichnesssFromNew <- function(geoList, basicMetaMap) {
    newRichnessList = lapply(geoList, function(cohortList){
      return(lapply(cohortList, function(currSamples)  return(basicMetaMap$GeneRichness[basicMetaMap$sample.ID %in% currSamples])))
    })
    return(newRichnessList)
  }
  
  getEntTab <- function(enteroList) {
    entGeoTab = do.call(rbind,lapply(enteroList,  function(enterotypes){
      entCounts = getEntCount(enterotypes)
      return(entCounts)
    }))
    colnames(entGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
    return(entGeoTab)
  }
  
  getEntGeoPlot <- function(enteroList, ordered=T, margins=c(8,4,4,1)) {
    entGeoTab = do.call(rbind,lapply(enteroList,  function(enterotypes){
      entCounts = getEntCount(enterotypes)
      return(entCounts)
    }))
    colnames(entGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
    
    entGeoTab = sweep(entGeoTab, MARGIN=1, rowSums(entGeoTab)/100, FUN="/")
    if (ordered) entGeoTab = entGeoTab[order(entGeoTab[,2]),]
    if (ordered) entGeoTab = entGeoTab[order(entGeoTab[,1]),]
    par(mar=margins)
    
    etfCol = "#77e2c799"
    etbCol = "#bee24599"
    etpCol = "#f35b38"
    
    tradVec = c(255,127,0,120)/255
    #tradVec = c(255,127,0,200)/255
    tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
    
    cjpVec = c(211,211,211,120)/255
    #cjpVec = c(211,211,211,200)/255
    cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
    
    eurVec = c(117,112,179,120)/255
    #eurVec = c(117,112,179,200)/255
    eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
    
    etpCol = tradCol
    etbCol = cjpCol
    etfCol = eurCol
    
    barplot(t(entGeoTab),#[c(2,4,3,5,6,7,8,9),]
            col=c(etfCol,etbCol,etpCol), las=2)
    return(entGeoTab)
  }
  
  #### meta-analysis ####
  checkMinEffectSizeOfFeatures <- function(datasetList, abdMat) {
    
    datasetSeqs = seq_along(datasetList)
    features = rownames(abdMat)
    
    minEffectSizeList = list()
    for (i in datasetSeqs) {
      samplesI = datasetList[[i]]
      minEffectSize = rep(10, length(features))
      names(minEffectSize) = features
      
      for (j in datasetSeqs) {
        samplesJ = datasetList[[j]]
        if (i!=j) {
          print(c(i,j))
          currCaseMat = abdMat[, colnames(abdMat) %in% samplesI]
          currContMat = abdMat[, colnames(abdMat) %in% samplesJ]
          
          currEsVec <- getEffectSizeMat(currCaseMat, currContMat)
          minEffectSize <- rowMins(cbind(minEffectSize, currEsVec))
          
        }
      }
      minEffectSizeList[[i]] = minEffectSize
      
    }
    minEffectSizeMat = do.call(cbind, minEffectSizeList)
    colnames(minEffectSizeMat) = names(datasetList)
    rownames(minEffectSizeMat) = rownames(abdMat)
    return(minEffectSizeMat)
    
  }
  checkMaxEffectSizeOfFeatures <- function(datasetList, abdMat) {
    
    datasetSeqs = seq_along(datasetList)
    features = rownames(abdMat)
    
    maxEffectSizeList = list()
    for (i in datasetSeqs) {
      samplesI = datasetList[[i]]
      maxEffectSize = rep(0, length(features))
      names(maxEffectSize) = features
      
      for (j in datasetSeqs) {
        samplesJ = datasetList[[j]]
        if (i!=j) {
          print(c(i,j))
          currCaseMat = abdMat[, colnames(abdMat) %in% samplesI]
          currContMat = abdMat[, colnames(abdMat) %in% samplesJ]
          
          currEsVec <- getEffectSizeMat(currCaseMat, currContMat)
          maxEffectSize <- rowMaxs(cbind(maxEffectSize, currEsVec))
          
        }
      }
      maxEffectSizeList[[i]] = maxEffectSize
      
    }
    maxEffectSizeMat = do.call(cbind, maxEffectSizeList)
    colnames(maxEffectSizeMat) = names(datasetList)
    rownames(maxEffectSizeMat) = rownames(abdMat)
    return(maxEffectSizeMat)
    
  }
  checkBestFeaturesInES <- function(datasetList, abdMat, esCut=0.8) {
    
    datasetSeqs = seq_along(datasetList)
    features = rownames(abdMat)
    
    featureList = list()
    for (i in datasetSeqs) {
      samplesI = datasetList[[i]]
      featureCounts = rep(0, length(features))
      names(featureCounts) = features
      
      for (j in datasetSeqs) {
        samplesJ = datasetList[[j]]
        if (i!=j) {
          print(c(i,j))
          currCaseMat = abdMat[, colnames(abdMat) %in% samplesI]
          currContMat = abdMat[, colnames(abdMat) %in% samplesJ]
          
          currEsVec <- getEffectSizeMat(currCaseMat, currContMat)
          currEsVec <- currEsVec >= esCut
          featureCounts <- featureCounts + currEsVec
          
        }
      }
      featureList[[i]] = featureCounts
      
      
    }
    featureMat = do.call(cbind, featureList)
    colnames(featureMat) = names(datasetList)
    return(featureMat)
    
  }
  getAvgAbdMatBy <- function(abdMat, datasetList) {
    
    avgMat = do.call(cbind, lapply(datasetList, function(dsSamples) {
      currMat = abdMat[, colnames(abdMat) %in% dsSamples]
      return(rowMeans(currMat))
    }))
    
    rownames(avgMat) = rownames(abdMat)
    colnames(avgMat) = names(datasetList)
    
    return(avgMat)
  }
  
  #### volcano plots ####
  getVolcanoStat <- function(targetMat, hgcSamples, lgcSamples, taxo, speciesMode=F) {
    targetMean = rowMeans(targetMat)
    targetPart  = targetMean/sum(targetMean)
    
    
    hgcMat = targetMat[,hgcSamples]
    lgcMat = targetMat[,lgcSamples]
    
    relAbdMat = getRelAbd(targetMat)
    relAbdHgcMat = relAbdMat[,hgcSamples]
    relAbdLgcMat = relAbdMat[,lgcSamples]
    
    allMean = rowMeans(relAbdMat)
    hgcMean = rowMeans(relAbdHgcMat)
    lgcMean = rowMeans(relAbdLgcMat)
    
    lenSpecies = dim(hgcMat)[1]
    pValues = sapply(1:lenSpecies, function(ind) {
      mergedVec = c(hgcMat[ind,],lgcMat[ind,])
      mergedFac = c(rep("HGC",length(hgcMat[ind,])),
                    rep("LGC",length(lgcMat[ind,])))
      
      wilcoxOut = wilcox.test(mergedVec ~ mergedFac,)
      return(wilcoxOut$p.value)
    })
    pValues[is.nan(pValues)]=1
    adjPs = p.adjust(pValues,method="BH")
    logAdjPs = -log10(adjPs)
    logAdjPs[logAdjPs>300] = 300
    
    lfcs = sapply(1:lenSpecies,function(ind){
      hmeans = mean(hgcMat[ind,])
      lmeans = mean(lgcMat[ind,]) 
      return(hmeans/lmeans)
    })
    names(lfcs)=rownames(hgcMat)
    lfcs[is.nan(lfcs)]=1
    lfcs[lfcs==0]=2^-15
    lfcs[is.infinite(lfcs)]=2^15
    lfcs = log2(lfcs)
    
    statsTab = data.frame(lfc=lfcs, 
                          pvalue=pValues,
                          qvalue=adjPs,
                          sig=logAdjPs, 
                          relAbd=allMean, 
                          relAbd_HGC = hgcMean,
                          relAbd_LGC = lgcMean,
                          label=names(lfcs), 
                          stringsAsFactors = F)
    
    if (speciesMode) {
      statsTab = data.frame(lfc=lfcs, 
                            pvalue=pValues,
                            qvalue=adjPs,
                            sig=logAdjPs,
                            relAbd=allMean,
                            relAbd_HGC = hgcMean,
                            relAbd_LGC = lgcMean,
                            label=getSpeciesName(names(lfcs),taxo), 
                            stringsAsFactors = F)
      
    }
    
    return(statsTab)
    
    
  }
  getVolcanoStatPaired <- function(targetMat, hgcSamples, lgcSamples, taxo, speciesMode=F) {
    targetMean = rowMeans(targetMat)
    targetPart  = targetMean/sum(targetMean)
    
    
    hgcMat = targetMat[,match(hgcSamples, colnames(targetMat))]
    lgcMat = targetMat[,match(lgcSamples, colnames(targetMat))]
    
    relAbdMat = getRelAbd(targetMat)
    relAbdHgcMat = relAbdMat[,match(hgcSamples, colnames(targetMat))]
    relAbdLgcMat = relAbdMat[,match(lgcSamples, colnames(targetMat))]
    
    allMean = rowMeans(relAbdMat)
    hgcMean = rowMeans(relAbdHgcMat)
    lgcMean = rowMeans(relAbdLgcMat)
    
    lenSpecies = dim(hgcMat)[1]
    pValues = sapply(1:lenSpecies, function(ind) {
      mergedVec = c(hgcMat[ind,],lgcMat[ind,])
      mergedFac = c(rep("HGC",length(hgcMat[ind,])),
                    rep("LGC",length(lgcMat[ind,])))
      
      wilcoxOut = wilcox.test(mergedVec ~ mergedFac, paird=T)
      return(wilcoxOut$p.value)
    })
    pValues[is.nan(pValues)]=1
    adjPs = p.adjust(pValues,method="BH")
    logAdjPs = -log10(adjPs)
    logAdjPs[logAdjPs>300] = 300
    
    lfcs = sapply(1:lenSpecies,function(ind){
      hmeans = mean(hgcMat[ind,])
      lmeans = mean(lgcMat[ind,]) 
      return(hmeans/lmeans)
    })
    names(lfcs)=rownames(hgcMat)
    lfcs[is.nan(lfcs)]=1
    lfcs[lfcs==0]=2^-15
    lfcs[is.infinite(lfcs)]=2^15
    lfcs = log2(lfcs)
    
    statsTab = data.frame(lfc=lfcs, 
                          pvalue=pValues,
                          qvalue=adjPs,
                          sig=logAdjPs, 
                          relAbd=allMean, 
                          relAbd_HGC = hgcMean,
                          relAbd_LGC = lgcMean,
                          label=names(lfcs), 
                          stringsAsFactors = F)
    
    if (speciesMode) {
      statsTab = data.frame(lfc=lfcs, 
                            pvalue=pValues,
                            qvalue=adjPs,
                            sig=logAdjPs,
                            relAbd=allMean,
                            relAbd_HGC = hgcMean,
                            relAbd_LGC = lgcMean,
                            label=getSpeciesName(names(lfcs),taxo), 
                            stringsAsFactors = F)
      
    }
    
    return(statsTab)
    
    
  }
  
  
  
  
  #### disease signatures ####
  
  plsrDiseaseSignature <- function(sampleTab, featureTab) {
    healthInd <- sampleTab$type == "health"
    diseaseInd <- sampleTab$type == "disease"
    
    currTab = sampleTab[,c("type","geography","sequencer")]
    
    lapply(as.list(colnames(featureTab)), function(currFeature) {
      dataTab = cbind(currTab, feature=featureTab[,currFeature])
      healthTab = dataTab[healthInd,]
      diseaseTab = dataTab[diseaseInd,]
      dataTab = rbind(healthTab, diseaseTab)
      
      out <- tryCatch({
        plsrOut = plsr(feature ~ type + geography, ncomp=15, data=dataTab, validation = "LOO")
        summary(plsrOut)
        plot(RMSEP(plsrOut), legendpos = "topright")
        plot(plsrOut, ncomp = 15, asp = 1, line = TRUE)
        plot(plsrOut, plottype = "scores", comps = 1:4)
        
        coef(plsrOut)
        
      }, error = function(cond){
        message(cond)
        print(currFeature)
        return(c(NA,NA))
      })
    })
    statOut = do.call(rbind,lapply(as.list(colnames(featureTab)), function(currFeature){
      
      dataTab = cbind(currTab, feature=featureTab[,currFeature])
      healthTab = dataTab[healthInd,]
      diseaseTab = dataTab[diseaseInd,]
      dataTab = rbind(healthTab, diseaseTab)
      
      out <-tryCatch({
        lmerOut = plsr(feature ~   type + (1|geography), data=dataTab, na.action = na.exclude)
        
        beta = summary(lmerOut)$coefficients["typehealth","Estimate"]
        betaP = summary(lmerOut)$coefficients["typehealth","Pr(>|t|)"]
        geoCoeff = coef(lmerOut)$geography[,1]
        geoCoeffNames = rownames(coef(lmerOut)$geography)
        names(geoCoeff) = geoCoeffNames
        
        return(c(beta,betaP))
      },
      error = function(cond) {
        message(cond)
        print(currFeature)
        return(c(NA,NA))
      })
    }))
    rownames(statOut) = colnames(featureTab)
    colnames(statOut) = c("estimate", "pvalue")
    
    
    
    
  }
  lmerDiseaseSignature <- function(sampleTab, featureTab) {
    healthInd <- sampleTab$type == "health"
    diseaseInd <- sampleTab$type == "disease"
    
    currTab = sampleTab[,c("type","geography","sequencer")]
    geoOut = lapply(as.list(colnames(featureTab)), function(currFeature){
      
      dataTab = cbind(currTab, feature=featureTab[,currFeature])
      healthTab = dataTab[healthInd,]
      diseaseTab = dataTab[diseaseInd,]
      dataTab = rbind(healthTab, diseaseTab)
      out <-tryCatch({
        lmerOut = lmer(feature ~   type + (1|geography), data=dataTab, na.action = na.exclude)
        
        beta = summary(lmerOut)$coefficients["typehealth","Estimate"]
        betaP = summary(lmerOut)$coefficients["typehealth","Pr(>|t|)"]
        geoCoeff = coef(lmerOut)$geography[,1]
        geoCoeffNames = rownames(coef(lmerOut)$geography)
        names(geoCoeff) = geoCoeffNames
        
        return(geoCoeff)
      },
      error = function(cond) {
        message(cond)
        print(currFeature)
        return(rep(NA,18))
      })
    })
    geoMat = do.call(cbind, geoOut)
    colnames(geoMat) = colnames(featureTab)
    geoMat[is.na(geoMat)] = 0
    
    statOut = do.call(rbind,lapply(as.list(colnames(featureTab)), function(currFeature){
      
      dataTab = cbind(currTab, feature=featureTab[,currFeature])
      healthTab = dataTab[healthInd,]
      diseaseTab = dataTab[diseaseInd,]
      dataTab = rbind(healthTab, diseaseTab)
      
      out <-tryCatch({
        lmerOut = lmer(feature ~   type + (1|geography), data=dataTab, na.action = na.exclude)
        
        beta = summary(lmerOut)$coefficients["typehealth","Estimate"]
        betaP = summary(lmerOut)$coefficients["typehealth","Pr(>|t|)"]
        geoCoeff = coef(lmerOut)$geography[,1]
        geoCoeffNames = rownames(coef(lmerOut)$geography)
        names(geoCoeff) = geoCoeffNames
        
        return(c(beta,betaP))
      },
      error = function(cond) {
        message(cond)
        print(currFeature)
        return(c(NA,NA))
      })
    }))
    rownames(statOut) = colnames(featureTab)
    colnames(statOut) = c("estimate", "pvalue")
    
    
    return(list(geography=geoMat, disease=statOut))
  }
  
  
  
  
  #### pie plots ####
  getRelAbdOfTaxa <- function(targetMat, targetTaxo, taxo) {
    targetTaxoMat = getTaxaSumMat(targetMat, taxo, targetTaxo,T)
    targetTaxoMean = rowMeans(targetTaxoMat)
    relAbds = targetTaxoMean/sum(targetTaxoMean)
    return(relAbds)
  }
  
  getRelAbdWithOthersImputed <- function(currRelAbd, topList){
    newRelAbd = currRelAbd[match(topList, names(currRelAbd))]
    relAbdOthers = sum(currRelAbd[!names(currRelAbd)%in%topList]) 
    newRelAbd = c(newRelAbd, relAbdOthers)
    names(newRelAbd) = c(topList, "others")
    return(newRelAbd)
  }
  
  #### geography signals ####
  showPerBarPlotByCountry <- function(samplesByCountry, target, mergeMat, cutOff=1e-7) {
    perList <- lapply(samplesByCountry, function(currSamples){
      currVec = mergeMat[target,currSamples]
      present = sum(currVec >= cutOff)/length(currVec)
      absent = sum(currVec < cutOff)/length(currVec)
      return(c(present, absent))
    })
    perMat <- do.call(cbind, perList)
    
    par(mar=c(5,10,4,2))
    barplot(perMat, horiz = T, las=1)
    par(mar=c(5.1,4.1,4.1,2.1))
    return(perMat)
  }
  showBoxPlotsByCountry <- function(samplesByCountry, target, mergeMat) {
    abdList <- lapply(samplesByCountry, function(currSamples){
      return(mergeMat[target,currSamples])
    })
    par(mar=c(5,10,4,2))
    boxplot(abdList,horizontal = T,las=1, ylim=c(0,1e-5))
    par(mar=c(5.1,4.1,4.1,2.1))
    return(abdList)
  }
  
  #### hidden cluster analysis ####
  getKeyClusterFromFirst <- function(targetMat, eps=0.3, pts=10) {
    set.seed(1)
    #cl = hdbscan(targetMat, minPts = pts)
    cl = dbscan(targetMat, eps=eps, minPts = pts)
    #plot(cl, show_flat = TRUE)
    #plot(cl)
    print(sort(table(cl$cluster)))
    if (all(cl$cluster=="0")) return(NULL)
    
    keyCluster=names(sort(table(cl$cluster),decreasing = T))
    keyCluster = keyCluster[keyCluster!="0"]
    keyCluster = keyCluster[1]
    print(keyCluster)
    keyClusterInd = which(cl$cluster==keyCluster)
    keyClusterSamples = rownames(targetMat)[keyClusterInd]
    return(keyClusterSamples)
  }
  getKeyClusterFromSecond <- function(targetMat, pts=10) {
    set.seed(1)
    cl = hdbscan(targetMat, minPts = pts)
    plot(cl, show_flat = TRUE)
    plot(cl)
    print(sort(table(cl$cluster)))
    
    keyCluster=names(sort(table(cl$cluster),decreasing = T))[2]
    print(keyCluster)
    keyClusterInd = which(cl$cluster==keyCluster)
    keyClusterSamples = rownames(targetMat)[keyClusterInd]
    return(keyClusterSamples)
  }
  
  
  
  #### functional clusters ####
  getGenesFromVfTerms <- function(vfTerms, patricVfDescTab) {
    genes = patricVfDescTab$Gene[match(vfTerms, patricVfDescTab$RefSeq.Locus.Tag)]
    return(genes)
  }
  getProductsFromVfTerms <- function(vfTerms, patricVfDescTab) {
    products = patricVfDescTab$Product[match(vfTerms, patricVfDescTab$RefSeq.Locus.Tag)]
    return(products)
  }
  getClassesFromVfTerms <- function(vfTerms, patricVfDescTab) {
    classes = patricVfDescTab$Classification[match(vfTerms, patricVfDescTab$RefSeq.Locus.Tag)]
    return(classes)
  }
  grepMspsByKOs <- function(targetKOs, gutKoBestMat, taxo, cutRatio=0.75, debug=F) {
    if (debug) print(length(targetKOs))
    if (length(targetKOs)==1) {
      currMspsByKos = (gutKoBestMat[,targetKOs])
    } else {
      currMspsByKos = rowSums(gutKoBestMat[,targetKOs])
    }
    currMspsByKos = sort(currMspsByKos, decreasing = T)
    currMspsByKosTab = data.frame(msp=names(currMspsByKos), 
                                  species=getSpeciesName(names(currMspsByKos), taxo), 
                                  count = currMspsByKos, stringsAsFactors = F)
    
    cutMspsByKosTab = currMspsByKosTab[currMspsByKosTab$count >= (length(targetKOs)*cutRatio),]
    return(list(msps=cutMspsByKosTab$msp, wholeData=currMspsByKosTab))
    
  }
  
  grepBestMspsAssociatedWithModules <- function(funcModules, funcMat,debug=F) {
    
    mspList = lapply(funcModules, function(currFuncs) {
      lenCurrFuncs = length(currFuncs)
      if (length(currFuncs)==1) {
        currMspsByFuncs = funcMat[,currFuncs]
      } else {
        currMspsByFuncs = rowSums(funcMat[,currFuncs])
      }
      maxRatio = max(currMspsByFuncs)
      currMspsByFuncs = currMspsByFuncs[currMspsByFuncs==maxRatio]
      return(names(currMspsByFuncs))
    })
    return(mspList)
  }
  
  grepBestScoresMspsAssociatedWithModules <- function(funcModules, funcMat,debug=F) {
    
    scoreList = lapply(funcModules, function(currFuncs) {
      lenCurrFuncs = length(currFuncs)
      if (length(currFuncs)==1) {
        currMspsByFuncs = funcMat[,currFuncs]
      } else {
        currMspsByFuncs = rowSums(funcMat[,currFuncs])
      }
      maxRatio = max(currMspsByFuncs)/lenCurrFuncs
      return(maxRatio)
    })
    return(scoreList)
  }
  grepMspsAssociatedWithModules <- function(funcModules, funcMat, cutRatio=0.75, debug=F) {
    
    mspList = lapply(funcModules, function(currFuncs) {
      lenCurrFuncs = length(currFuncs)
      if (length(currFuncs)==1) {
        currMspsByFuncs = funcMat[,currFuncs]
      } else {
        currMspsByFuncs = rowSums(funcMat[,currFuncs])
      }
      currMspsByFuncs = sort(currMspsByFuncs, decreasing = T)
      currMspsByFuncs = currMspsByFuncs[currMspsByFuncs>=(lenCurrFuncs*cutRatio)]
      return(names(currMspsByFuncs))
    })
    return(mspList)
  }
  
  
  
  grepMspsByFuncModules <- function(funcKoModules, gutKoBestMat, taxo, targetModule, cutRatio=0.75, debug=F) {
    targetKOs = funcKoModules[[targetModule]]
    if (debug) print(length(targetKOs))
    if (length(targetKOs)==1) {
      currMspsByKos = (gutKoBestMat[,targetKOs])
    } else {
      currMspsByKos = rowSums(gutKoBestMat[,targetKOs])
    }
    currMspsByKos = sort(currMspsByKos, decreasing = T)
    currMspsByKosTab = data.frame(msp=names(currMspsByKos), 
                                  species=getSpeciesName(names(currMspsByKos), taxo), 
                                  count = currMspsByKos, stringsAsFactors = F)
    
    cutMspsByKosTab = currMspsByKosTab[currMspsByKosTab$count >= (length(targetKOs)*cutRatio),]
    return(list(msps=cutMspsByKosTab$msp, wholeData=currMspsByKosTab))
    
  }
  checkFractionOfGenus <- function(gutTaxo, targetMsps, plotMode = F) {
    gutGenusCount = table(gutTaxo$genus)
    gutGenusCountVec = as.numeric(gutGenusCount)
    names(gutGenusCountVec) = names(gutGenusCount)
    
    targetGenera = gutTaxo$genus[gutTaxo$MSP %in% targetMsps]
    targetGenusCount = table(targetGenera)
    targetGenusCountVec = as.numeric(targetGenusCount)
    names(targetGenusCountVec) = names(targetGenusCount)
    
    matchedGutGenusCountVec = gutGenusCountVec[match(names(targetGenusCountVec),names(gutGenusCountVec))]
    outVec = targetGenusCountVec*100/matchedGutGenusCountVec
    if (plotMode) {
      par(mar=c(4,7,3,1))
      showVec = sort(outVec)
      showVec = showVec[!grepl("unclassified", names(showVec))]
      barplot(showVec, xlab="found species (%)",
              col=colorRampPalette(c("#0000FF", "white", "#808080"))(length(outVec)), 
              cex.axis = 0.6, cex.names = 0.6,
              horiz = T, las=1)
      par(mar=c(5,4,4,1))
    }
    return(outVec)
    
  }
  findModulesBy <- function(funcModules, targetTerms) {
    areContained = do.call(c, lapply(funcModules, function(x) any(x %in% targetTerms) ))
    return(names(funcModules)[areContained])
  }
  
  
  #### function comparisons ####
  getChisqStat <- function(refFuncMat, targetFuncMat, funcs) {
    
    pvalues = sapply(funcs, function(currFunc) {
      refCurrFunc = refFuncMat[currFunc,]
      targetCurrFunc = targetFuncMat[currFunc,]
      
      refPresent = sum(refCurrFunc==1)
      refAbsent = sum(refCurrFunc==0)
      
      targetPresent = sum(targetCurrFunc==1)
      targetAbsent = sum(targetCurrFunc==0)
      
      contigTab = cbind(refFlow=c(present=refPresent, absent=refAbsent),
                        targetFlow=c(present=targetPresent, absent=targetAbsent))
      pvalue=chisq.test(contigTab)$p.value
      return(pvalue)
    }, simplify = F)
    pvalues = unlist(pvalues)
    
    oddsRatios = sapply(funcs, function(currFunc) {
      refCurrFunc = refFuncMat[currFunc,]
      targetCurrFunc = targetFuncMat[currFunc,]
      
      refPresent = sum(refCurrFunc==1)
      refAbsent = sum(refCurrFunc==0)
      
      targetPresent = sum(targetCurrFunc==1)
      targetAbsent = sum(targetCurrFunc==0)
      
      
      oddsRatio = (targetPresent*refAbsent)/(targetAbsent*refPresent)
      if (refPresent*targetAbsent == 0) return(Inf) 
      return(oddsRatio)
    }, simplify = F)
    oddsRatios = unlist(oddsRatios)
    
    chisqStatTab = data.frame(term=names(pvalues),
                              pvalue=pvalues,
                              oddsRatio=oddsRatios,
                              logP = -log10(pvalues),
                              logOdds = log(oddsRatios))
    
    return(chisqStatTab)
  }
  getFractionMat <- function(fMat, bMat, pMat, funcs) {
    ratios = sapply(funcs, function(currFunc) {
      fFunc = fMat[currFunc,]
      bFunc = bMat[currFunc,]
      pFunc = pMat[currFunc,]
      
      fPresent = sum(fFunc==1)
      fAbsent = sum(fFunc==0)
      
      bPresent = sum(bFunc==1)
      bAbsent = sum(bFunc==0)
      
      pPresent = sum(pFunc==1)
      pAbsent = sum(pFunc==0)
      
      fRatio = fPresent/(fPresent + fAbsent)
      bRatio = bPresent/(bPresent + bAbsent)
      pRatio = pPresent/(pPresent + pAbsent)
      
      return(c(fRatio, bRatio, pRatio))
    }, simplify = F)
    ratioMat = do.call(rbind, ratios)
    rownames(ratioMat) = funcs
    colnames(ratioMat) = c("F","B","P")
    return(ratioMat)      
    
  }
  getVarScoreMat <- function(funcMat, funcs, etfSpecies, etbSpecies, etpSpecies) {
    
    selFuncMat = funcMat[,colnames(funcMat) %in% funcs]
    varMat = apply(selFuncMat, 2, function(currFunc) {
      fClass = rep("no",dim(selFuncMat)[1])
      bClass = rep("no",dim(selFuncMat)[1])
      pClass = rep("no",dim(selFuncMat)[1])
      fClass[rownames(selFuncMat)%in% etfSpecies]="yes"
      bClass[rownames(selFuncMat)%in% etbSpecies]="yes"
      pClass[rownames(selFuncMat)%in% etpSpecies]="yes"
      
      currData = data.frame(fClass= fClass, bClass=bClass, pClass=pClass, func=currFunc)
      coefs = summary(lm(func~fClass+bClass+pClass, data=currData))$coefficients[,1]
      fstats = summary(lm(func~fClass+bClass+pClass, data=currData))$fstatistic 
      pvalue = pf(fstats[1], fstats[2], fstats[3], lower.tail = F) 
      
      fRatio = mean(currFunc[fClass=="yes"]) / mean(currFunc[fClass=="no"])
      bRatio = mean(currFunc[bClass=="yes"]) / mean(currFunc[bClass=="no"])
      pRatio = mean(currFunc[pClass=="yes"]) / mean(currFunc[pClass=="no"])
      
      fBeta = coefs[2]
      bBeta = coefs[3]
      pBeta = coefs[4]
      
      fBeta = fBeta^2
      bBeta = bBeta^2
      pBeta = pBeta^2
      
      sumBeta = fBeta+bBeta+pBeta
      vars = c(fBeta, bBeta, pBeta)/sumBeta
      out = c(vars, fRatio, bRatio, pRatio, pvalue)
      names(out) = c("fBeta", "bBeta", "pBeta", "fRatio", "bRatio", "pRatio", "pvalue")
      return(out)
      
    })
    
    varMat = t(varMat)
    varMat = as.data.frame(varMat)
    varMat$term = rownames(varMat)
    return(varMat)
  }
  getVarScoreMatTwo <- function(funcMat, funcs, inSpecies, outSpecies) {
    
    selFuncMat = funcMat[,colnames(funcMat) %in% funcs]
    varMat = apply(selFuncMat, 2, function(currFunc) {
      fClass = rep("no",dim(selFuncMat)[1])
      bClass = rep("no",dim(selFuncMat)[1])
      fClass[rownames(selFuncMat)%in% inSpecies]="yes"
      bClass[rownames(selFuncMat)%in% outSpecies]="yes"
      
      currData = data.frame(fClass= fClass, bClass=bClass,  func=currFunc)
      coefs = summary(lm(func~fClass+bClass , data=currData))$coefficients[,1]
      fstats = summary(lm(func~fClass+bClass , data=currData))$fstatistic 
      pvalue = pf(fstats[1], fstats[2], fstats[3], lower.tail = F) 
      
      fRatio = mean(currFunc[fClass=="yes"]) / mean(currFunc[fClass=="no"])
      bRatio = mean(currFunc[bClass=="yes"]) / mean(currFunc[bClass=="no"])
      
      fBeta = coefs[2]
      bBeta = coefs[3]
      
      fBeta = fBeta^2
      bBeta = bBeta^2
      
      sumBeta = fBeta+bBeta 
      vars = c(fBeta, bBeta )/sumBeta
      out = c(vars, fRatio, bRatio,   pvalue)
      names(out) = c("inBeta", "outBeta",  "inRatio", "outRatio",  "pvalue")
      return(out)
      
    })
    
    varMat = t(varMat)
    varMat = as.data.frame(varMat)
    varMat$term = rownames(varMat)
    return(varMat)
  }
  
  selectVarScoreMat <- function(varMat, pCut=1e-2) {
    pvalueInd = varMat$pvalue < pCut
    
    betaMaxInds = apply(varMat, 1, function(currRow){
      return(which.max(currRow[1:3]))
    })
    
    ratioMaxInds = apply(varMat, 1, function(currRow){
      return(which.max(currRow[4:6]))
    })
    
    consistentInds <- betaMaxInds == ratioMaxInds
    
    
    sigMat = varMat[pvalueInd & consistentInds, ] 
    #sigMat = sigMat[fInd | bInd | pInd, ] 
    return(sigMat)
  }
  selectVarScoreMatTwo <- function(varMat, pCut=1e-2) {
    pvalueInd = varMat$pvalue < pCut
    
    betaMaxInds = apply(varMat, 1, function(currRow){
      return(which.max(currRow[1:2]))
    })
    
    ratioMaxInds = apply(varMat, 1, function(currRow){
      return(which.max(currRow[3:4]))
    })
    
    consistentInds <- betaMaxInds == ratioMaxInds
    
    
    sigMat = varMat[pvalueInd & consistentInds, ] 
    #sigMat = sigMat[fInd | bInd | pInd, ] 
    return(sigMat)
  }
  
  
  
  
}


load("current.example.MGS.other.omics.RData")

repeatResult = T
if (repeatResult) {
  
  mixedModelMetaboMgsTabRData = "J://Deposit//Project//2018_microbiome_atlas//microbiome.atlas.Rproj//mixedModelMetaboMgsTab.2021.1024.RData"
  mixedModelMetaboFamilyTabRData = "J://Deposit//Project//2018_microbiome_atlas//microbiome.atlas.Rproj//mixedModelMetaboFamilyTab.2021.1026.RData"
  load(mixedModelMetaboMgsTabRData)
  load(mixedModelMetaboFamilyTabRData)
  
  
  circosMode = T
  if (circosMode) {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    #selectedTab = mixedModelMetaboFamilyTab[mixedModelMetaboFamilyTab$exp.var > 0.1, ]
    selectedTab = mixedModelMetaboMgsTab[mixedModelMetaboMgsTab$exp.var > 0.25, ]
    
    # explained variance > 10%
    selectedTabPos = selectedTab[selectedTab$tMgs > 0,]
    # positive relations
    selectedTabNeg = selectedTab[selectedTab$tMgs < 0,]
    # negative relations
    
    circosInput = data.frame(from=getGenusName(selectedTabPos$mgs, taxo),
                             to=selectedTabPos$metabo,
                             value=(abs(selectedTabPos$exp.var)),
                             stringsAsFactors = F)
    
    
    circosInput = circosInput[!grepl("unclassified",circosInput$from),]
    
    elements = (unique(c(circosInput$from, circosInput$to)))
    lenElements = length(elements)
    
    circosGridCols = (col_vector)[1:length(elements)]
    names(circosGridCols) = elements
    
    #circos.par(track.margin=c(0,0)) 
    set.seed(1)
    circos.initialize(elements, xlim=c(0,1))
    chordDiagram(circosInput)
    chordDiagram(circosInput, annotationTrack = "grid", preAllocateTracks = 1, grid.col = circosGridCols)
    circos.trackPlotRegion(track.index = 1,  panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(mean(xlim), ylim[1] + .1, cex = 0.6,  sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      circos.axis(h = "top", labels = F, major.tick = F, 
                  sector.index = sector.name, track.index = 2)
    }, bg.border = NA)
    circos.clear() 
    
    
    
  }
  
  
}


########################################3
# below needed when we generate results


loadWellnessClinData = T
if (loadWellnessClinData) {
  wellnessClinData = "J://Deposit//Project//2017.Wellness//clinical.meta.data//20170822//wellness.meta.data.extended.99.subjects.170626.txt"
  wellnessClinTab = read.delim(wellnessClinData, sep="\t", stringsAsFactors = F)
  colnames(wellnessClinTab)[1] = "subject_id"
  
  
  ###ongoing
  meltClinMatrix <- function(wellnessClinTab) {
    defaultFields = c("subject_id","scapis_id","SUBJID","birthdate","height","imt.cca","blood.type","atherosclerosis","age.visit1")
    exFields = c("type","fasting","datetime","directanalysis","biobanksamples","faeces","other_smoking","other_smoking_hours","other_antibiotics","other_antibiotic_details","other_health_change","other_health_change_details")
    defaultTab = wellnessClinTab[,defaultFields]
    #View(defaultTab[defaultTab$atherosclerosis==2,])
    
    clinVisitTab = wellnessClinTab[,!colnames(wellnessClinTab) %in% defaultFields]
    rownames(clinVisitTab) = paste("X", wellnessClinTab$subject_id, sep="")
    
    clinVisitExMatch = sapply(colnames(clinVisitTab), function(currColumn){
      matched = sapply(exFields, function(exField) {
        return(grepl(exField, currColumn))
      })
      matched = any(matched)
      return(matched)
    })
    
    clinVisitTab = clinVisitTab[,!clinVisitExMatch]
    
    clinVisitIds = sapply(colnames(clinVisitTab), getVisitIds)
    clinVisitIdsUnique = unique(clinVisitIds)
    clinVisitFields = sapply(colnames(clinVisitTab), getFields)
    clinVisitFieldsUnique = unique(clinVisitFields)
    
    
    #subject_visit by clinical parameters
    
    mergedClinVisitTab = do.call(rbind, sapply(clinVisitFieldsUnique, function(currField) {
      
      data = clinVisitTab[,clinVisitFields==currField]
      visitInfo = clinVisitIds[clinVisitFields==currField]
      visitInfo = gsub("isit","",visitInfo)
      colnames(data) = visitInfo
      data = as.matrix(data)
      melted = melt(data)
      melted$field = currField
      melted$Var3 = paste(melted$Var1, melted$Var2, sep = "_")
      return(melted)
    }, simplify = F))
    
    mergedClinVisitMat = acast(mergedClinVisitTab, Var3 ~ field, value.var = "value")
    dim(mergedClinVisitMat)
    
    return(mergedClinVisitMat)
  }
  
  clinicalMatrix = meltClinMatrix(wellnessClinTab)
  colnames(clinicalMatrix) = gsub("analysis_code_","",colnames(clinicalMatrix))
  colnames(clinicalMatrix) = gsub("_result","",colnames(clinicalMatrix))
  #cbind(colnames(clinicalMatrix))
  
}

loadMode = T
if (loadMode) {
  funcModulesRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation//funcModules.20191227.RData"
  load(funcModulesRData)
  
  koDescMapFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation\\ko_desc_mapping.20190603.RData"
  load(koDescMapFile)
  koDescMap = koDescMap[!duplicated(koDescMap),]
  
  
  pcoaAllUpdatedData = "C:/Data/comparative.analysis.healthy.sweden/pcoaAllData.updated.20190805.RData"
  load(pcoaAllUpdatedData)
  
  pcoaMat = pcoaOut$vectors 
  pcoaGenusMat = pcoaGenusOut$vectors
  pcoaFamilyMat = pcoaFamilyOut$vectors
  pcoaClassMat = pcoaClassOut$vectors
  pcoaOrderMat = pcoaOrderOut$vectors
  dim(pcoaOut$vectors)
  
  wellnessMgsDist=vegdist(t(wellnessMgsMat2), method="bray")
  wellnessPcoaOut=pcoa(wellnessMgsDist)
  wellnessPcoaMat = wellnessPcoaOut$vectors 
  
  gutMsp1992Data = "C:\\Data/comparative.analysis.healthy.sweden/hs_10.4_1992_MSP_freeze2_20180905.RData"
  load(gutMsp1992Data)
  rm(MSP_gc)
  rm(MSP_id)
  gutMsps = rownames(taxo)
  #rm(taxo)
  
  
  suppressMsps = c("msp_1974","msp_2153")
  suppressMsps2 = c("msp_1974","msp_2153", "msp_0864")
  gutMsps = gutMsps[!gutMsps %in% suppressMsps2]
  
  taxo = taxo[!rownames(taxo) %in% suppressMsps,]
  gutTaxo = taxo[gutMsps,]
  
  mgsInfoFile="C:\\Data/comparative.analysis.healthy.sweden/mgs.files.info.txt"
  mgsInfo = read.table(mgsInfoFile, header = T, sep = "\t",stringsAsFactors = F)
  
  gssSpecies = gutTaxo$MSP[gutTaxo$MSP %in% taxo$MSP[taxo$type_sh=="gss"]]
  ossSpecies = gutTaxo$MSP[gutTaxo$MSP %in% taxo$MSP[taxo$type_sh=="oss"]]
  gutOnlySpecies = gutTaxo$MSP[gutTaxo$MSP %in% taxo$MSP[taxo$type_sh=="gut"]]
  oralOnlySpecies = gutTaxo$MSP[gutTaxo$MSP %in% taxo$MSP[taxo$type_sh=="oral"]]
  
  mgsCleanUpdatedRData2 = "C:\\Data/comparative.analysis.healthy.sweden/mgs.clean.updated.20190726.10M.RData"
  load(mgsCleanUpdatedRData2)
  
  mergeMatUpdated = mergeMatUpdated[!rownames(mergeMatUpdated)%in%suppressMsps,]
  
  basicMetaCleanRDataUpdated2 = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.behcet.20190805.RData"
  load(basicMetaCleanRDataUpdated2)
  
  
}

loadWellnessMgsData = T
if (loadWellnessMgsData) {
  
  wellnessMetadata = basicMetaMapUpdated[basicMetaMapUpdated$dataset.ID=="id28",]
  
  mgsMat = mergeMatUpdated
  metaTab = basicMetaMapUpdated
  
  wellnessMetadata = metaTab[metaTab$dataset.ID=="id28",]
  wellnessSubjs = (wellnessMetadata$metadata.ID)
  missingSubjs1 = names(table(wellnessSubjs)[table(wellnessSubjs)==3])
  exceptionSubjs = c("X3452","X3834","X3918","X4215")
  missingSubjs2 = c("X3913","X3535","X3016")
  missingSubjs1 = missingSubjs1[!missingSubjs1 %in% exceptionSubjs]
  
  wellnessMetadata = wellnessMetadata[!wellnessMetadata$metadata.ID %in% missingSubjs1, ]
  wellnessMetadata = wellnessMetadata[!wellnessMetadata$metadata.ID %in% missingSubjs2, ]
  wellnessMetadata = wellnessMetadata[!wellnessMetadata$metadata.ID %in% exceptionSubjs, ]
  wellnessMgsMat = mergeMatUpdated[,match(wellnessMetadata$sample.ID, colnames(mergeMatUpdated))]
  wellnessMgsMatDigit = (wellnessMgsMat>0) * 1
  wellnessMgsMeans = rowMeans(wellnessMgsMat)
  wellnessSamples = wellnessMetadata$sample.ID
  wellnessSubjects = gsub("_v[1-4]","",wellnessSamples)
  wellnessVisits = gsub("X[0-9]+_","",wellnessSamples)
  wellnessSamplesBySubjects = split(wellnessSamples, wellnessSubjects)
  wellnessUniqSubjs = unique(wellnessMetadata$metadata.ID)
  wellnessClinicalMatrix = t(clinicalMatrix)[,match(wellnessSamples, rownames(clinicalMatrix))]
  
}


# wellnessMetaboDist=vegdist(t(wellnessMetabolome), method="bray")
# wellnessMetaboPcoaOut=pcoa(wellnessMetaboDist)
# wellnessMetaboPcoaMat = wellnessMetaboPcoaOut$vectors 
# 
# plot(wellnessMetaboPcoaMat[,1:2])
# 


loadWellnessMetabolomeData = T
if (loadWellnessMetabolomeData) {
  metabolomeFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.multiomics\\scapis_wellness_metabo_v1234.txt"
  metabolomeTab = read.delim(metabolomeFile, sep="\t")
  
  View(metabolomeTab)
  
  wellnessMetabolomeMap = metabolomeTab[,1:5]
  wellnessMetabolomeMap$NewId = gsub("Meta", "", metabolomeTab$variable)
  wellnessMetabolome = metabolomeTab[,-c(1:5)]
  
  metabolomeNames = paste(wellnessMetabolomeMap$Name, wellnessMetabolomeMap$NewId, sep = "_")
  #metabolomeNames[1:5]
  
  rownames(wellnessMetabolome) = metabolomeNames
  colnames(wellnessMetabolome) = gsub("_","_v", colnames(wellnessMetabolome))
  wellnessMetabolome = wellnessMetabolome[,!is.na(colSums(wellnessMetabolome))]
  
  #dim(wellnessMetabolome)
  
  
}



generateResult = F
if (generateResult) {
  
  wellnessMgsMatWithMetabo = wellnessMgsMat[,colnames(wellnessMgsMat) %in% colnames(wellnessMetabolome)]
  wellnessMetabolomeWithMgs = wellnessMetabolome[,match(colnames(wellnessMgsMatWithMetabo), colnames(wellnessMetabolome))]
  
  wellnessMgsMatWithMetabo = wellnessMgsMatWithMetabo[rowSums(wellnessMgsMatWithMetabo)!=0,]
  wellnessGenusMat = getTaxaSumMat(wellnessMgsMatWithMetabo, taxo, "genus" )
  wellnessGenusMat = wellnessGenusMat[rowSums(wellnessGenusMat)!=0,]
  wellnessFamilyMat = getTaxaSumMat(wellnessMgsMatWithMetabo, taxo, "family" )
  wellnessFamilyMat = wellnessFamilyMat[rowSums(wellnessFamilyMat)!=0,]
  
  
  
  mixedEffectMode = T
  if (mixedEffectMode) {
    wellnessMetabolomeWithMgsT = as.data.frame(t(wellnessMetabolomeWithMgs))
    wellnessMetabolomeWithMgsT$sample = rownames(wellnessMetabolomeWithMgsT)
    wellnessMetabolomeWithMgsT$subject =  gsub("_v[1-4]", "", wellnessMetabolomeWithMgsT$sample)
    wellnessMetabolomeWithMgsT$subject =  factor(wellnessMetabolomeWithMgsT$subject)
    wellnessMetabolomeWithMgsT$visit = gsub("^X[0-9]+_", "", wellnessMetabolomeWithMgsT$sample)
    metaboFeatures = rownames(wellnessMetabolomeWithMgs)
    
    wellnessMgsMatWithMetaboT = as.data.frame(t(wellnessMgsMatWithMetabo))
    wellnessMgsMatWithMetaboT$sample = rownames(wellnessMgsMatWithMetaboT)
    wellnessMgsMatWithMetaboT$subject =  gsub("_v[1-4]", "", wellnessMgsMatWithMetaboT$sample)
    wellnessMgsMatWithMetaboT$subject =  factor(wellnessMgsMatWithMetaboT$subject)
    wellnessMgsMatWithMetaboT$visit = gsub("^X[0-9]+_", "", wellnessMgsMatWithMetaboT$sample)
    mgsFeatures = rownames(wellnessMgsMatWithMetabo)
    
    wellnessGenusMatWithMetaboT = as.data.frame(t(wellnessGenusMat))
    wellnessGenusMatWithMetaboT$sample = rownames(wellnessGenusMatWithMetaboT)
    wellnessGenusMatWithMetaboT$subject =  gsub("_v[1-4]", "", wellnessGenusMatWithMetaboT$sample)
    wellnessGenusMatWithMetaboT$subject =  factor(wellnessGenusMatWithMetaboT$subject)
    wellnessGenusMatWithMetaboT$visit = gsub("^X[0-9]+_", "", wellnessGenusMatWithMetaboT$sample)
    genusFeatures = rownames(wellnessGenusMat)
    
    wellnessFamilyMatWithMetaboT = as.data.frame(t(wellnessFamilyMat))
    wellnessFamilyMatWithMetaboT$sample = rownames(wellnessFamilyMatWithMetaboT)
    wellnessFamilyMatWithMetaboT$subject =  gsub("_v[1-4]", "", wellnessFamilyMatWithMetaboT$sample)
    wellnessFamilyMatWithMetaboT$subject =  factor(wellnessFamilyMatWithMetaboT$subject)
    wellnessFamilyMatWithMetaboT$visit = gsub("^X[0-9]+_", "", wellnessFamilyMatWithMetaboT$sample)
    familyFeatures = rownames(wellnessFamilyMat)
    
    
    
    ###checking models ####
    
    ###### warning~! it takes a day ######
    mixedModelMetaboMgsTab = do.call(rbind,sapply(metaboFeatures, function(currMetabo) {
      print(paste(which(metaboFeatures == currMetabo), date(), sep="|||||"))
      mixedModelResult = do.call(rbind.data.frame, sapply(mgsFeatures, function(currMgs) {
        #print(which(mgsFeatures == currMgs))
        currData = data.frame(metabolite = wellnessMetabolomeWithMgsT[,currMetabo],
                              mgs = wellnessMgsMatWithMetaboT[,currMgs],
                              subj = wellnessMetabolomeWithMgsT$subject,
                              stringsAsFactors = F)
        lmerMetaboMgs = lmer(metabolite ~ mgs + (1|subj), data = currData)
        if (dim(summary(lmerMetaboMgs)$coefficients)[1]<2) {
          return(list(pMgs=NA,  tMgs = NA, exp.var = exp.vNAar))
          return(result)
        }
        #print(summary(lmerMetaboMgs))
        #print(which(mgsFeatures == currMgs))
        varNames = rownames(summary(lmerMetaboMgs)$coefficients)
        exp.var = r.squaredGLMM(lmerMetaboMgs)[1]
        
        tMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Estimate"], NA)
        pMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Pr(>|t|)"], NA)
        
        return(list(pMgs=pMgs,   tMgs = tMgs, exp.var = exp.var, mgs = currMgs, metabo = currMetabo))
        
      }, simplify=F))
      return(mixedModelResult)
    }, simplify = F))
    mixedModelMetaboMgsTab$species = getSpeciesName(mixedModelMetaboMgsTab$mgs, taxo)
    mixedModelMetaboMgsTabRData = "J://Deposit//Project//2018_microbiome_atlas//microbiome.atlas.Rproj//mixedModelMetaboMgsTab.2021.1024.RData"
    save(mixedModelMetaboMgsTab, file=mixedModelMetaboMgsTabRData)
    
    
    ###### warning~! it takes a half day  ######
    mixedModelMetaboFamilyTab = do.call(rbind,sapply(metaboFeatures, function(currMetabo) {
      print(paste(which(metaboFeatures == currMetabo), date(), sep="|||||"))
      mixedModelResult = do.call(rbind.data.frame, sapply(familyFeatures, function(currMgs) {
        #print(which(mgsFeatures == currMgs))
        currData = data.frame(metabolite = wellnessMetabolomeWithMgsT[,currMetabo],
                              mgs = wellnessFamilyMatWithMetaboT[,currMgs],
                              subj = wellnessMetabolomeWithMgsT$subject,
                              stringsAsFactors = F)
        lmerMetaboMgs = lmer(metabolite ~ mgs + (1|subj), data = currData)
        if (dim(summary(lmerMetaboMgs)$coefficients)[1]<2) {
          return(list(pMgs=NA,  tMgs = NA, exp.var = exp.vNAar))
          return(result)
        }
        #print(summary(lmerMetaboMgs))
        #print(which(mgsFeatures == currMgs))
        varNames = rownames(summary(lmerMetaboMgs)$coefficients)
        exp.var = r.squaredGLMM(lmerMetaboMgs)[1]
        
        tMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Estimate"], NA)
        pMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Pr(>|t|)"], NA)
        
        return(list(pMgs=pMgs,   tMgs = tMgs, exp.var = exp.var, mgs = currMgs, metabo = currMetabo))
        
      }, simplify=F))
      return(mixedModelResult)
    }, simplify = F))
    mixedModelMetaboFamilyTabRData = "J://Deposit//Project//2018_microbiome_atlas//microbiome.atlas.Rproj//mixedModelMetaboFamilyTab.2021.1026.RData"
    save(mixedModelMetaboFamilyTab, file=mixedModelMetaboFamilyTabRData)
    
    
  } 
  
  
}
 