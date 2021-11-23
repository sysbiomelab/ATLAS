require(ggplot2)
require(ggrepel)
require(reshape2)
getSpeciesName <- function(mgsName, taxo) {
  species = taxo[match(mgsName, rownames(taxo)),"species"]
  return(species)
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

getEnrichMentGsByGs <- function(refGeneSet, targetGeneSet) {
  allRefGenes = unique.default(do.call(c, refGeneSet))
  allTargetGenes = unique.default(do.call(c, targetGeneSet))

  hyperRes = lapply(refGeneSet, function(currRefGenes) {
    hyperPs = lapply(targetGeneSet, function(currTargetGenes){
      currTargetGenes = currTargetGenes[currTargetGenes %in% allRefGenes]
      hyperP = getHyperP( currTargetGenes, currRefGenes, allRefGenes )
    })
    return(hyperPs)
  })
}
getHyperP <-function( targetGenes, refGenes, allGenes ) {
  targetGenes = targetGenes[targetGenes %in% allGenes]

  numAll = length(allGenes)
  numOverlap = length( refGenes[refGenes %in% targetGenes] )
  numTrue = length(refGenes)
  numPositive = length(targetGenes)

  out = 1-phyper(numOverlap-1, numTrue, numAll - numTrue, numPositive)
  return(out)
}
getVectorFromTable <- function(tabObj) {
  tabNames = names(tabObj)
  tabValues = as.numeric(tabObj)
  names(tabValues) = tabNames
  return(tabValues)
}

#dinput data with msp ES above 0.3, numer of cohort where it appears and number of cohorts is enriched minus number of cohorts is depleated

msp_vol<-read.csv("Volcano_plot.Country.tsv", sep="\t", row.names=1 )
taxo<-read.csv("../taxo.csv",row.names=1)
msp_vol$species<- getSpeciesName(msp_vol$msp, taxo)
msp_vol$species<- sapply(strsplit(as.character(msp_vol$species),'/'), "[", 1)
msp_vol$species<- sapply(strsplit(as.character(msp_vol$species),'&'), "[", 1)
#pdf(file= "VolcanoPlot.pdf", width=6, height=4 )
jpeg(f="MSPinDisease.volcanoplot.jpg", w=4*300, h=3*300,pointsize=100)

set.seed(2)
pos <- position_jitter(width = 0.3, seed = 1)
ggplot(msp_vol[abs(msp_vol$dif)<=5,], aes(x=dif, y=sum)) + #, label=label
  geom_jitter(shape = 21, colour = "black", width = 0.6, height = 0.6, size=5, aes(fill = dif)) +
  geom_jitter( data = msp_vol[abs(msp_vol$dif)>5,], mapping = aes(x = dif, y=sum, fill=dif), shape = 21, size = 5, position = pos)+
  geom_text_repel(data = msp_vol[abs(msp_vol$dif)>5,],
                  aes(x=dif,y=sum,label=paste(species,msp, sep="\n")),
                  size=8, segment.color = 'grey', position=pos, max.overlaps=50)+

  #scale_fill_gradient2(low = "blue", mid="white",  high = "red") +
  scale_fill_gradient2(low = "red", mid="white",  high = "blue") +
  scale_x_continuous(breaks=seq(-10,10,4),minor_breaks = seq(-10, 10, 1))+
  scale_y_continuous(breaks=seq(0,10,4),minor_breaks = seq(0, 10, 1))+

  geom_hline(yintercept = 3, colour="gray") +
  geom_vline(xintercept = 2, colour="gray") +
  geom_vline(xintercept = -2, colour="gray") +
  geom_abline(slope=1)+
  geom_abline(slope=-1) +
  xlab("Number of times Enriched -  Number of times Depleted ") + ylab("Number of Cohorts with Effect Size above 0.3")+
  theme(panel.grid.major = element_line(size = 0.5, linetype = "dotted", colour = "gray"),
     panel.grid.minor = element_line(size = 0.5, linetype = "dotted", colour = "gray"),
     legend.position = "none",
     axis.title=element_text(size=30),
     panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"))
dev.off()

#contains "funcMat" var matrix with gene names as colums, msp as rows and 0/1 presence absence
load("sunjae/funcMat.20190808.RData")
#contains "funcModules" var list of modules and their gene name definitions
load("sunjae/funcModules.20191227.RData")
#build mspByMod list
mspByModList =grepMspsAssociatedWithModules(funcModules = funcModules,
                                               funcMat = funcMat,
                                               cutRatio = 0.75)
# Optional subset selected modules and  filter modules with minumum number of MSP
#names(assocFuncMsps) = names(funcModuleSel)
#assocFuncMspNums = do.call(c, lapply(assocFuncMsps, length))
funcModuleNums = do.call(c, lapply(mspByModList, length))
#mspCut = 10
#assocFuncMsps = assocFuncMsps[assocFuncMspNums>=mspCut]

#read Effect size dataframe of disease cohorts
es_cohort<-read.csv("Effect_size.Country_ctrl.tsv" ,row.names=1, sep="\t")
#Build lists depletedMspByCohortWCountryControl and enrichedMspByCohortsWCountryControl
#list of enriched species with effectsize >0.3 by cohort
efup<-es_cohort[which(es_cohort$ES_up>0.3),]
efdw<-es_cohort[which(es_cohort$ES_down>0.3),]
#convert into list
efdw_list<-split(efdw$msp,efdw$cohort)
efup_list<-split(efup$msp,efup$cohort)

#estimate hypergeometric prob,
#compare set of MSPs for func cluster against set of es>0.3 MSPs by cohort
pCut=1e-3
  enrichSigDownDiseaseFuncMat = getEnrichMentGsByGs(mspByModList, efdw_list)
  enrichSigDownDiseaseFuncMelt = melt(enrichSigDownDiseaseFuncMat)
  enrichSigDownDiseaseFuncMeltSig=enrichSigDownDiseaseFuncMelt[enrichSigDownDiseaseFuncMelt$value<=pCut, ]
  #funcDownDiseaseCounts = getVectorFromTable(table(enrichSigDownDiseaseFuncMeltSig$Var2))
  funcDownDiseaseCounts = getVectorFromTable(table(enrichSigDownDiseaseFuncMeltSig$L1))

  enrichSigUpDiseaseFuncMat = getEnrichMentGsByGs(mspByModList, efup_list)
  enrichSigUpDiseaseFuncMelt = melt(enrichSigUpDiseaseFuncMat)
  enrichSigUpDiseaseFuncMeltSig=enrichSigUpDiseaseFuncMelt[enrichSigUpDiseaseFuncMelt$value<=pCut, ]
  #funcUpDiseaseCounts = getVectorFromTable(table(enrichSigUpDiseaseFuncMeltSig$Var2))
  funcUpDiseaseCounts = getVectorFromTable(table(enrichSigUpDiseaseFuncMeltSig$L1))
  #get names of functional clusters
  allFuncNames = unique(c(names(funcDownDiseaseCounts),
                         names(funcUpDiseaseCounts)))
  #build dataframe matching eriched and depleated sets
  funcClusterStatTab = data.frame(name=allFuncNames,
                                 down=funcDownDiseaseCounts[match(allFuncNames, names(funcDownDiseaseCounts))],
                                 up=funcUpDiseaseCounts[match(allFuncNames, names(funcUpDiseaseCounts))],
                                 size=funcModuleNums[match(allFuncNames, names(funcModuleNums))],
                                 stringsAsFactors = F)
  #change NA values to 0
  funcClusterStatTab[is.na(funcClusterStatTab)]=0
#estimate coord for volcano plot
funcClusterStatTab$diff =  funcClusterStatTab$up - funcClusterStatTab$down
funcClusterStatTab$sum = funcClusterStatTab$up + funcClusterStatTab$down

set.seed(1)
ggplot(funcClusterStatTab2, aes(x=diff, y=sum)) +
  geom_jitter(shape = 21, colour = "#00000055", width = 0.5, height = 0.5, size=2, aes(fill = outflowSig)) +
  #scale_fill_manual(values=c("#80808033", "#6a0dad88")) +
  #geom_jitter(shape = 21, colour = "black", width = 0.5, height = 0.5, size=2, aes(fill = diff)) +
  #scale_fill_gradient2(low = "blue", mid="white",  high = "red") +
  geom_hline(yintercept = 3, colour="gray") +
  geom_vline(xintercept = 2, colour="gray") +
  geom_vline(xintercept = -2, colour="gray") +
  geom_abline(slope=1)+
  geom_abline(slope=-1) +
  xlab("|Enrichment| - |Depletion|") + ylab("|Enrichment| + |Depletion|")+
  scale_y_continuous(breaks=c(2,4,6,8))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5,
                                        linetype = "solid"))
