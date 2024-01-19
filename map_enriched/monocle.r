require(igraph)
require(monocle)

### very important ###
trajectoryMode = T
if (trajectoryMode) {
  
  #wholeDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/whole.trajectory.cds.RData"
  #normDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/norm.trajectory.cds.RData"
  #normDataTrajectory2RData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/norm.trajectory2.cds.RData"
  #exceptEurDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/except.eur.trajectory.cds.RData"
  #exceptCjpDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/except.cjp.trajectory.cds.RData"
  #exceptTradDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/except.trad.trajectory.cds.RData"
  #prunedDataTrajectoryRData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData/pruned.trajectory.cds.RData"
  wholeDataTrajectoryRData = '../data/whole.trajectory.cds.RData'
  normDataTrajectoryRData = '../data/norm.trajectory.cds.RData'
  #normDataTrajectory2RData = '../data/norm.trajectory2.cds.RData'
  normDataTrajectory2RData = '../data/new.norm.trajectory2.cds.RData'
  exceptEurDataTrajectoryRData = '../data/except.eur.trajectory.cds.RData'
  exceptCjpDataTrajectoryRData = '../data/except.cjp.trajectory.cds.RData'
  exceptTradDataTrajectoryRData = '../data/except.trad.trajectory.cds.RData'
#theo
#for annotateddataframe
  library(Biobase)
  source("check.mgs.functions.r")
  load('../data/vect_atlas.RData')
  mergeMatUpdated <- vect_atlas
  load('../../../FMT/downstream_data/hs_10.4_1992_MSP_freeze2_20180905.RData')
  basicMetaMapUpdated = read.csv('../data/unique_metadata.csv')

#old data
  load('../data/all.basic.clean.metadata.behcet.20190805.RData')
#srr55... is the sample.id
#dataset.id id53
  
  loadData = T
  if (loadData) {
    newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
    newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
    newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
    
    cjpSamples = unlist(newNormalSamplesByGeo[c("China","Japan", "US")])
    eurSamples = unlist(newNormalSamplesByGeo[europeanCountries])
    eurFiSamples = unlist(newNormalSamplesByGeo[c(europeanCountries, "Finland")])
    trdSamples = unlist(newNormalSamplesByGeo[c(traditionalCountries,"Thai")])
  }
  
  normDataMode2 = F
  if (normDataMode2) {
    
    
    
    generateMode = T
    if (generateMode) {
      
      allNormMode = T
      if (allNormMode) {
        #normAllMat = mergeMatUpdated*10^9
        normAllMat = vect_atlas*10^9
        #normAllMat = normAllMat[, colnames(normAllMat) %in% newNormalAllIds_wo_FI]
	sharedIDs = intersect(basicMetaMapUpdated$secondary_sample_accession, colnames(normAllMat))
	#sharedIDs = intersect(basicMetaMapUpdated$sample.ID, colnames(normAllMat))
        normAllMat = normAllMat[, sharedIDs]
	#basicMetaMapUpdated = basicMetaMapUpdated[sharedIDs, ]
	basicMetaMapUpdated =  basicMetaMapUpdated[ which(basicMetaMapUpdated$secondary_sample_accession==sharedIDs), ]
        #basicMetaMapUpdated[, basicMetaMapUpdated.secondary_sample_accession %in% colnames(normAllMat)]
        #min(normAllMat[normAllMat!=0]) --> 3.8194
        
       # trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
       #                                  Well=basicMetaMapUpdated$dataset.ID,
       #                                  Hours=0,
       #                                  Media="NONE")
# theo
        trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$secondary_sample_accession,
                                         Well=basicMetaMapUpdated$study_accession,
                                         Hours=0,
                                         Media="NONE")
		
        
        #rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
#theo
        rownames(trajectorySampleTab) = basicMetaMapUpdated$secondary_sample_accession
        #trajectorySampleTab = trajectorySampleTab[rownames(trajectorySampleTab) %in% newNormalAllIds_wo_FI,]
        trajectorySampleTab = trajectorySampleTab[sharedIDs,]
        
        trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                       biotype="protein_coding",
                                       num_cells_expressed = 1,
                                       use_for_ordering=FALSE)
        
        rownames(trajectoryTaxaTab) = rownames(taxo)
        trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]
        
        pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
        fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
        
        cds <- newCellDataSet(as.matrix(c(normAllMat), 
                              phenoData = pd, 
                              featureData = fd,
                              lowerDetectionLimit = 0.1)
        
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        
        cds <- reduceDimension(cds, 
                               max_components = 2,
                               method = 'DDRTree')
        
        cds <- orderCells(cds)
        
        save(cds, file=normDataTrajectory2RData)
      }
      
      
      
      
    }
    
    allNormMode  = T
    if (allNormMode) {
      getPercentGeoByState <- function(targetSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) {
        
        
        currCounts = do.call(c, lapply(newNormalSamplesByGeo, function(currSamples) {
          counts = sum(currSamples %in% targetSamples )
          return(counts)
        }))
        
        return(currCounts*100/newNormalSampleCountsByGeo)
      }
      
      firstBranchSamples = unlist(samplesByState[c("1","2","9","8")])
      midBranchSamples = unlist(samplesByState[c("3")])
      lastBranchSamples = unlist(samplesByState[c("4","5","6","7")])
      
      
    }
  }
  
  normDataMode = T
  if (normDataMode) {
    
    pcaTsneMode = F
    if (pcaTsneMode) {
      
      loadData = T
      if (loadData) {
        newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
        newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
        newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
        
        cjpSamples = unlist(newNormalSamplesByGeo[c("China","Japan", "US")])
        eurSamples = unlist(newNormalSamplesByGeo[europeanCountries])
        trdSamples = unlist(newNormalSamplesByGeo[traditionalCountries])
        
        length(cjpSamples)
        length(eurSamples)
        length(trdSamples)
        
        normAllMat = mergeMatUpdated*10^9
        normAllMat = normAllMat[, colnames(normAllMat) %in% newNormalAllIds]
      }
      
      dimReductionMode = T
      if (dimReductionMode) {
        normPcaOut = prcomp(t(normAllMat))
        eigs <- normPcaOut$sdev^2
        cumsum(eigs[1:20]) / sum(eigs)
        normPcaMat = normPcaOut$x[,1:11]
      }
      
      tsneMode = T
      if (tsneMode) {
        set.seed(1)
        tsneOut = Rtsne(normPcaMat, dims = 2, perplexity = 100)
        normTsneMat = tsneOut$Y
        rownames(normTsneMat) = rownames(normPcaMat)
        
        # tsneOut = Rtsne(t(normAllMat), dims = 2, perplexity = 50)
        # normTsneMat = tsneOut$Y
        # rownames(normTsneMat) = colnames(normAllMat)
        
        tsneCols = rep("gray",dim(normTsneMat)[1])
        names(tsneCols) = rownames(normTsneMat)
        tsneCols[cjpSamples] = cjpCol
        tsneCols[eurSamples] = eurCol
        tsneCols[trdSamples] = tradCol
        
        plot(normTsneMat, col=tsneCols, pch=16,xlab="tSNE1",ylab="tSNE2", xaxt="n",yaxt="n" )
      }
      
      
      
      
    }
    
    generateMode = F
    if (generateMode) {
      normAllMat = mergeMatUpdated*10^9
      normAllMat = normAllMat[, colnames(normAllMat) %in% newNormalAllIds]
      
      trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
                                       Well=basicMetaMapUpdated$dataset.ID,
                                       Hours=0,
                                       Media="NONE")
      
      rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
      trajectorySampleTab = trajectorySampleTab[rownames(trajectorySampleTab) %in% newNormalAllIds,]
      
      trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                     biotype="protein_coding",
                                     num_cells_expressed = 1,
                                     use_for_ordering=FALSE)
      
      rownames(trajectoryTaxaTab) = rownames(taxo)
      trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]
      
      pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
      fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
      
      cds <- newCellDataSet(normAllMat, 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit = 0.1)
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      cds <- reduceDimension(cds, 
                             max_components = 2,
                             method = 'DDRTree')
      
      cds <- orderCells(cds)
      
      save(cds, file=normDataTrajectoryRData)
      
    }
    
    loadMode = T
    if (loadMode) load(normDataTrajectoryRData)
    #save(cds, file=normDataTrajectoryRData)
    
    checkTrajectory = T
    if (checkTrajectory) {
      
      mergeMatUpdatedScaled = t(scale(t(mergeMatUpdatedScaled)))
      mergeMatInflowScaled  = colSums(mergeMatUpdatedScaled[rownames(mergeMatUpdatedScaled)%in%inflowSpecies,])
      mergeMatOutflowScaled  = colSums(mergeMatUpdatedScaled[rownames(mergeMatUpdatedScaled)%in%outflowSpecies,])
      
      basicMetaMapUpdated$inflowScaled = mergeMatInflowScaled[match(basicMetaMapUpdated$sample.ID , names(mergeMatInflowScaled))]
      basicMetaMapUpdated$outflowScaled = mergeMatOutflowScaled[match(basicMetaMapUpdated$sample.ID , names(mergeMatOutflowScaled))]
      
      inflowByEnterotype = split(basicMetaMapUpdated$inflowScaled, basicMetaMapUpdated$enteroType)
      outflowByEnterotype = split(basicMetaMapUpdated$outflowScaled, basicMetaMapUpdated$enteroType)
      
      par(mar=c(10,4,4,1))
      boxplot(inflowByEnterotype, las=2)
      boxplot(outflowByEnterotype, las=2)
      
      
      getTrajectoryData <- function(cds, newNormalMetadata) {
        requireNamespace("igraph")
        gene_short_name <- NA
        sample_name <- NA
        sample_state <- pData(cds)$State
        data_dim_1 <- NA
        data_dim_2 <- NA
        lib_info_with_pseudo <- pData(cds)
        x = 1; y = 2
        data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
        #data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
        #  select_(data_dim_1 = x, data_dim_2 = y) %>% rownames_to_column("sample_name") %>% 
        #  mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
        #                                       rownames_to_column("sample_name"), by = "sample_name")
        #output = data_df[,c("data_dim_1","data_dim_2")]
        #rownames(output) = data_df$sample_name
        #output$geo = newNormalMetadata$Geography[match(rownames(output),newNormalMetadata$sample.ID)]
        return(output)
      }
      getAvgPointsByGeo <- function(trajectoryData, newNormalSamplesByGeo, plotMode = T, labelMode = T) {
        newNormalSamplesByGeo = newNormalSamplesByGeo[names(newNormalSamplesByGeo) !="Finland"]
        geoPoints = do.call(rbind,lapply(newNormalSamplesByGeo, function(currSamples){
          currData = trajectoryData[rownames(trajectoryData) %in% currSamples,]
          return(colMeans(currData[,1:2]))
        }))
        geoSds = do.call(c,lapply(newNormalSamplesByGeo, function(currSamples){
          currData = trajectoryData[rownames(trajectoryData) %in% currSamples,]
          print(length(currSamples))
          #print(dim(currData))
          #print(colSds(currData))
          currSds = c(sd(currData[,1]),sd(currData[,2]))
          currSqSds = sqrt(sum(currSds^2))*2
          return(currSqSds)
        }))
        geoPointTab = data.frame(geoPoints, geo=rownames(geoPoints), sd=geoSds)
        if (plotMode) {
          cols = c("Tanzania"=tradCol,"Peru"=tradCol,"Madagascar"=tradCol,"Mongo"=tradCol,"Thai"=tradCol,"Fiji"=tradCol, "India"=tradCol,"Finland"=eurCol,"China"=cjpCol,"Japan"=cjpCol,"US"=cjpCol,"Italy"=eurCol,"Spain"=eurCol,"Luxembourg"=eurCol,"Germany"=eurCol,"UK"=eurCol,"Denmark"=eurCol,"Sweden"=eurCol)
          
          if (labelMode) {
            p = ggplot(geoPointTab, aes(x=data_dim_1, y=data_dim_2, color=factor(geo), label=geo, size=sd)) + 
              geom_point() + 
              ggrepel::geom_text_repel(size=3.5, color="black") +
              scale_color_manual(values=cols) + 
              xlab("Component 1")+ylab("Component 2")+ 
              theme(panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_blank(),
                    strip.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    strip.text.y = element_blank(),
                    legend.position="none", 
                    panel.border = element_rect(colour = "black", fill=NA, size=.5))
            print(p)
          }
          if (!labelMode) {
            p = ggplot(geoPointTab, aes(x=data_dim_1, y=data_dim_2, color=factor(geo), label=geo, size=sd)) + 
              geom_point() + 
              scale_color_manual(values=cols) + 
              xlab("Component 1")+ylab("Component 2")+ 
              theme(panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_blank(),
                    strip.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    strip.text.y = element_blank(),
                    legend.position="none", 
                    panel.border = element_rect(colour = "black", fill=NA, size=.5))
            print(p)
          }
          
          
        }
        return(geoPointTab)
      }
      
      
      
      cols = c("Tanzania"=tradCol,"Peru"=tradCol,"Madagascar"=tradCol,"Mongo"=tradCol,"Thai"=tradCol,"Fiji"=tradCol, "India"=tradCol,"Finland"=eurCol,"China"=cjpCol,"Japan"=cjpCol,"US"=cjpCol,"Italy"=eurCol,"Spain"=eurCol,"Luxembourg"=eurCol,"Germany"=eurCol,"UK"=eurCol,"Denmark"=eurCol,"Sweden"=eurCol)
      
      trajectoryData = getTrajectoryData(cds, newNormalMetadata)
      
      getAvgPointsByGeo(trajectoryData, newNormalSamplesByGeo, T, F)
      
      trajectoryData2 = getTrajectoryData(cds, basicMetaMapUpdated)
#PROBABLY STOP HERE
      
      trajectoryData$Age = basicMetaMapUpdated$Age[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      trajectoryData$BMI = basicMetaMapUpdated$BMI[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      trajectoryData$Sequencer = basicMetaMapUpdated$Sequencer[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      trajectoryData$Gender = basicMetaMapUpdated$Gender[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      trajectoryData$Enterotype = basicMetaMapUpdated$enteroType[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      trajectoryData$inflow = mergeMatInflowScaled[match(rownames(trajectoryData), names(mergeMatInflowScaled))]
      trajectoryData$outflow = mergeMatOutflowScaled[match(rownames(trajectoryData), names(mergeMatOutflowScaled))]
      trajectoryData$richness = basicMetaMapUpdated$GeneRichness[match(rownames(trajectoryData), basicMetaMapUpdated$sample.ID)]
      
      
      ggplot(trajectoryData, aes(x=data_dim_1, y=data_dim_2, color=factor(geo), label=geo)) + 
        geom_point() + 
        scale_color_manual(values=cols) + 
        xlab("Component 1")+ylab("Component 2")+ 
        theme(panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text.x = element_blank(),
              axis.text.y = element_blank(),
              strip.text.y = element_blank(),
              legend.position="none", 
              panel.border = element_rect(colour = "black", fill=NA, size=.5))
      
      
      #trajectory.inflow.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$richness),], aes(x=data_dim_1, y=data_dim_2, color=richness)) + 
        geom_point(alpha=0.5) + scale_colour_gradient2(low="blue",mid = "gray",high = "red", midpoint=6e5)+
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      
      #trajectory.inflow.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$inflow),], aes(x=data_dim_1, y=data_dim_2, color=inflow)) + 
        geom_point(alpha=0.5) + scale_colour_gradient2(low="blue",mid = "gray",high = "red", midpoint=20)+
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      
      #trajectory.outflow.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$outflow),], aes(x=data_dim_1, y=data_dim_2, color=outflow)) + 
        geom_point(alpha=0.5) + scale_colour_gradient2(low="blue",mid = "gray",high = "red", midpoint=100)+
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      #trajectory.age.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$Age),], aes(x=data_dim_1, y=data_dim_2, color=Age)) + 
        geom_point(alpha=0.5) + scale_colour_gradient2(low="blue",mid = "gray",high = "red", midpoint=50)+
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      #trajectory.bmi.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$BMI),], aes(x=data_dim_1, y=data_dim_2, color=BMI)) + 
        geom_point(alpha=0.5) + scale_colour_gradient2(low="blue",mid = "gray",high = "red", midpoint=20)+
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      #trajectory.sequencer.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$Sequencer),], aes(x=data_dim_1, y=data_dim_2, color=factor(Sequencer))) + 
        geom_point(alpha=0.3)+ scale_color_brewer(palette="Dark2") +
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      
      #trajectory.sequencer.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$Sequencer),], aes(x=data_dim_1, y=data_dim_2, color=factor(Sequencer))) + 
        geom_point(alpha=0.3)+ scale_color_brewer(palette="Dark2") +
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      #trajectory.gender.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$Gender),], aes(x=data_dim_1, y=data_dim_2, color=factor(Gender))) + 
        geom_point(alpha=0.3)+ scale_color_brewer(palette="Dark2") +
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      #trajectory.enterotype.mapping
      ggplot(trajectoryData[!is.na(trajectoryData$Enterotype),], aes(x=data_dim_1, y=data_dim_2, color=factor(Enterotype))) + 
        geom_point(alpha=0.3)+ scale_color_brewer(palette="Dark2") +
        xlab("Component 1")+ylab("Component 2") + theme_bw()
      
      
      
      plot_cell_trajectory(cds, color_by = "Pseudotime")
      plot_cell_trajectory(cds, color_by = "State") +
        scale_color_manual(breaks = c("1", "2", "3","4","5"), 
                           values=c(tradCol, cjpCol, cjpCol, cjpCol, eurCol)) + 
        xlab("Component 1")+ylab("Component 2")+
        theme(axis.text.x=element_blank(), #legend.position = "none",
              
              axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      
      
      brachInfo = pData(cds)
      samplesByState = split(as.character(brachInfo$Library), brachInfo$State)
      
      targetData = newNormalAllList$Finland_id21
      lapply(samplesByState, function(currSamples)any(currSamples %in% targetData))
      
      targetState = 5
      currStateInfo = (basicMetaMapUpdated[ basicMetaMapUpdated$sample.ID %in% rownames(brachInfo[brachInfo$State==targetState,]),])
      mean(currStateInfo$GeneRichness)
      sort(table(currStateInfo$Geography))
      
      newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
      newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
      newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
      
      getPercentGeoByState <- function(targetSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) {
        
        
        currCounts = do.call(c, lapply(newNormalSamplesByGeo, function(currSamples) {
          counts = sum(currSamples %in% targetSamples )
          return(counts)
        }))
        
        return(currCounts*100/newNormalSampleCountsByGeo)
      }
      
      firstBranchSamples = samplesByState[["1"]]
      midBranchSamples = unlist(samplesByState[c("2","3","4")])
      lastBranchSamples = samplesByState[["5"]]
      
      
      View(basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% firstBranchSamples, ])
      View(basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% midBranchSamples, ])
      
      
      firstBranchStat = getPercentGeoByState(firstBranchSamples,
                                             newNormalSamplesByGeo,
                                             newNormalSampleCountsByGeo)
      
      midBranchStat = getPercentGeoByState(midBranchSamples,#unlist(samplesByState[c("2","3","4")])
                                           newNormalSamplesByGeo,
                                           newNormalSampleCountsByGeo)
      
      lastBranchStat = getPercentGeoByState(lastBranchSamples,#samplesByState[["5"]]
                                            newNormalSamplesByGeo,
                                            newNormalSampleCountsByGeo)
      sort(getPercentGeoByState(samplesByState[["7"]],
                                newNormalSamplesByGeo,
                                newNormalSampleCountsByGeo))
      
      
      ordered0 = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                   "Japan","China","US", "Finland",
                   "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
      
      ordered = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                  "Japan","China","US", 
                  "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
      
      
      tradVec = c(255,127,0,120)/255
      tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
      
      cjpVec = c(211,211,211,120)/255
      cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
      
      eurVec = c(117,112,179,120)/255
      eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
      
      
      orderedCols = c(rep(tradCol, 7),
                      rep(cjpCol, 3), #4
                      rep(eurCol, 7))
      par(mfrow=c(4,1))
      par(mar=c(1,4,1,1))
      barplot(firstBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
      par(mar=c(1,4,1,1))
      barplot(midBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
      par(mar=c(1,4,1,1))
      barplot(lastBranchStat[ordered], col=orderedCols, las=2)
      #barplot 
      
      ### check age distrubtion by 
      
      agesByNewNormalGeo <- lapply(newNormalSamplesByGeo, function(currSamples){
        ages = basicMetaMapUpdated$Age[basicMetaMapUpdated$sample.ID %in% currSamples]
        return(ages)
      })
      
      par(mfrow=c(2,1))
      par(mar=c(3,4,1,1))
      boxplot(agesByNewNormalGeo[ordered0], las=2, ylab="Age")  
      
      generateGeoDiffStatForRezaPouyan = T
      if (generateGeoDiffStatForRezaPouyan) {
        firstBranchSamples = samplesByState[["1"]]
        midBranchSamples = unlist(samplesByState[c("2","3","4")])
        lastBranchSamples = samplesByState[["5"]]
        
        firstSampleTab = data.frame(sample=firstBranchSamples, type="TraditionalCluster")
        midSampleTab = data.frame(sample=midBranchSamples, type="UsChinaJapanCluster")
        lastSampleTab = data.frame(sample=lastBranchSamples, type="EuropeanCluster")
        allBranchSampleTab = rbind(firstSampleTab, 
                                   midSampleTab,
                                   lastSampleTab)
        write.table(allBranchSampleTab, "C://Data//comparative.analysis.healthy.sweden//countryClusterSampleTab.txt", sep="\t")
        
        volcano.stats.mgs.first.vs.mid = getVolcanoStat(mergeMatUpdated, midBranchSamples, firstBranchSamples, taxo, T)
        volcano.stats.mgs.mid.vs.last = getVolcanoStat(mergeMatUpdated, lastBranchSamples, midBranchSamples, taxo, T)
        volcano.stats.mgs.first.vs.last = getVolcanoStat(mergeMatUpdated, lastBranchSamples, firstBranchSamples, taxo, T)
        
        statsMgsForFirstVsMidNewFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.first.vs.mid.new.txt"
        statsMgsForMidVsLastNewFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.mid.vs.last.new.txt"
        statsMgsForFirstVsLastNewFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.first.vs.last.new.txt"
        
        write.table(volcano.stats.mgs.first.vs.mid, statsMgsForFirstVsMidNewFile, sep="\t")
        write.table(volcano.stats.mgs.mid.vs.last, statsMgsForMidVsLastNewFile, sep="\t")
        write.table(volcano.stats.mgs.first.vs.last, statsMgsForFirstVsLastNewFile, sep="\t")
        
      }
      
      generateDetailGeoDiffStatForReza = T
      if (generateDetailGeoDiffStatForReza) {
        
        newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
        newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
        newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
        
        europeCountries = c("Denmark", "Germany", "Italy", "Japan", "Luxembourg", "Spain", "Sweden", "UK")
        tradCountries = c("Fiji", "India", "Madagascar", "Mongo", "Peru", "Tanzania", "Thai")
        geoCompList = list(CN = newNormalSamplesByGeo$China, 
                           EU = unlist(newNormalSamplesByGeo[europeCountries]), 
                           US = newNormalSamplesByGeo$US, 
                           JP = newNormalSamplesByGeo$Japan, 
                           TR = unlist(newNormalSamplesByGeo[tradCountries]))
        
        geoCompMelt = melt(geoCompList)
        colnames(geoCompMelt) = c("sample.ID","Country_code")
        
        geoCompMeltFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats//geo.diff.stats.samples.20190925.txt"
        write.table(geoCompMelt, geoCompMeltFile, sep="\t", quote = F, row.names = F)
        
        geoCompLen = length(geoCompList)
        geoSeqs = seq_along(geoCompList)
        
        resList = list()
        ind=1
        iSeqs = seq(1, geoCompLen-1)
        for (i in iSeqs) {
          iGeoSamples = geoCompList[[i]]
          iGeo = names(geoCompList)[i]
          
          jSeqs = seq(i+1, geoCompLen)
          for (j in jSeqs) {
            if (i!=j) 
              jGeoSamples = geoCompList[[j]]
            jGeo = names(geoCompList)[j]
            
            ijVolcanoStats = getVolcanoStat(mergeMatUpdated, jGeoSamples, iGeoSamples, taxo, T)
            colnames(ijVolcanoStats)= c("lfc_target_over_ref", "pvalue", "qvalue", "sig", "relAbd", "relAbd_target", "relAbd_reference", "label")
            ijVolcanoStats$msp = rownames(ijVolcanoStats)
            ijVolcanoStats$target = iGeo
            ijVolcanoStats$reference = jGeo
            resList[[ind]] = ijVolcanoStats
            ind = ind+1
          }
        }
        resTab = do.call(rbind, resList)
        geoCompWilcoxStatsFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats//geo.diff.stats.wilcoxon.output.20190925.txt"
        write.table(resTab, geoCompWilcoxStatsFile, sep="\t", quote = F, row.names = F)
        
      }
      
      generatePouyanStat = F
      if (generatePouyanStat) {
        
        firstSamples = newNormalSamplesByGeo[ordered[1:7]]
        midSamples = newNormalSamplesByGeo[ordered[8:10]]
        lastSamples = newNormalSamplesByGeo[ordered[11:17]]
        
        firstSamples = unlist(firstSamples)
        midSamples = unlist(midSamples)
        lastSamples = unlist(lastSamples)
        
        volcano.stats.mgs.first.vs.mid = getVolcanoStat(mergeMatUpdated, midSamples, firstSamples, taxo, T)
        volcano.stats.mgs.mid.vs.last = getVolcanoStat(mergeMatUpdated, lastSamples, midSamples, taxo, T)
        volcano.stats.mgs.first.vs.last = getVolcanoStat(mergeMatUpdated, lastSamples, firstSamples, taxo, T)
        
        
        statsMgsForFirstVsMidFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.first.vs.mid.txt"
        statsMgsForMidVsLastFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.mid.vs.last.txt"
        statsMgsForFirstVsLastFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.geo.difference.stats\\richness.volcano.stats.mgs.first.vs.last.txt"
        
        write.table(volcano.stats.mgs.first.vs.mid, statsMgsForFirstVsMidFile, sep="\t")
        write.table(volcano.stats.mgs.mid.vs.last, statsMgsForMidVsLastFile, sep="\t")
        write.table(volcano.stats.mgs.first.vs.last, statsMgsForFirstVsLastFile, sep="\t")
        
      }
    }
    
    checkWithRichnessSignature = T
    if (checkWithRichnessSignature) {
      
      loadRichnessNewCutMode = T
      if (loadRichnessNewCutMode) {
        adjpCut = 1e-3
        statsMgsForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.normal.txt"
        statsMgsForNorm = read.table(statsMgsForNormSamples, sep="\t", header=T, stringsAsFactors = F)
        statsMgsForNormSel = statsMgsForNorm[statsMgsForNorm$qvalue <= adjpCut,]
        # statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 0,]
        # statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -0,]
        statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 2,]
        statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -2,]
        
        hgcSpecies = rownames(statsMgsForNormSelHgc)
        lgcSpecies = rownames(statsMgsForNormSelLgc)
        
        length(hgcSpecies)
        length(lgcSpecies)
        
        if (F) {
          mergeScaleMat =  t(scale(t(mergeMatUpdated)))
          hgcSignatureSum = colSums(mergeScaleMat[hgcSpecies,])/sqrt(length(hgcSpecies))
          lgcSignatureSum = colSums(mergeScaleMat[lgcSpecies,])/sqrt(length(lgcSpecies))
          
          gssLen = sum(rownames(mergeScaleMat) %in% gssSpecies)
          ossLen = sum(rownames(mergeScaleMat) %in% ossSpecies)
          
          
          gssSignatureSum = colSums(mergeScaleMat[rownames(mergeScaleMat) %in% gssSpecies,])/sqrt(length(gssSpecies))
          ossSignatureSum = colSums(mergeScaleMat[rownames(mergeScaleMat) %in% ossSpecies,])/sqrt(length(ossSpecies))
          
          totalSignatures = hgcSignatureSum - lgcSignatureSum
          trajectoryData2$richnessSignature = totalSignatures[match(rownames(trajectoryData2), names(totalSignatures))]
          
          gutOralSignatures = ossSignatureSum - gssSignatureSum
          trajectoryData2$gutOralSignature = gutOralSignatures[match(rownames(trajectoryData2), names(gutOralSignatures))]
          
          oralSignature = ossSignatureSum
          trajectoryData2$oralSignature = oralSignature[match(rownames(trajectoryData2), names(oralSignature))]
          
          
        }
        normMgsMat = mergeMatUpdated[, colnames(mergeMatUpdated) %in% newNormalAllIds]
        normMgsMat = normMgsMat[rowSums(normMgsMat)!=0,]
        normMgsScaleMat = t(scale(t(normMgsMat)))
        
        
        hgcSpecies = hgcSpecies[hgcSpecies %in% rownames(normMgsScaleMat)]
        lgcSpecies = lgcSpecies[lgcSpecies %in% rownames(normMgsScaleMat)]
        length(hgcSpecies)
        length(lgcSpecies)
        
        
        hgcSignatureSum = colSums(normMgsScaleMat[hgcSpecies,])/sqrt(length(hgcSpecies))
        lgcSignatureSum = colSums(normMgsScaleMat[lgcSpecies,])/sqrt(length(lgcSpecies))
        
        gssLen = sum(rownames(normMgsScaleMat) %in% gssSpecies)
        ossLen = sum(rownames(normMgsScaleMat) %in% ossSpecies)
        
        
        gssSignatureSum = colSums(normMgsScaleMat[rownames(normMgsScaleMat) %in% gssSpecies,])/sqrt(length(gssSpecies))
        ossSignatureSum = colSums(normMgsScaleMat[rownames(normMgsScaleMat) %in% ossSpecies,])/sqrt(length(ossSpecies))
        
        totalSignatures = hgcSignatureSum - lgcSignatureSum
        trajectoryData$richnessSignature = totalSignatures[match(rownames(trajectoryData), names(totalSignatures))]
        
        gutOralSignatures = ossSignatureSum - gssSignatureSum
        trajectoryData$gutOralSignature = gutOralSignatures[match(rownames(trajectoryData), names(gutOralSignatures))]
        
        oralSignature = ossSignatureSum
        trajectoryData$oralSignature = oralSignature[match(rownames(trajectoryData), names(oralSignature))]
        
        library(scales)
        ggplot(trajectoryData, aes(x=data_dim_1, y=data_dim_2, color=richnessSignature, label=geo)) + 
          geom_point() + 
          scale_color_gradient2(low = "#2c528c88", mid = "#ffffff77",high = "#ff000088") + 
          xlab("Component 1")+ylab("Component 2")+ 
          theme(panel.background = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text.x = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_blank(), #legend.position="none", 
                panel.border = element_rect(colour = "black", fill=NA, size=.5))
        
        ggplot(trajectoryData, aes(x=data_dim_1, y=data_dim_2, color=gutOralSignatures, label=geo)) + 
          geom_point(aes(size=gutOralSignature)) + 
          scale_color_gradient2(low = "#2c528c33", mid = "#ffffff33", midpoint = 5,high = "#ff0000") + 
          #scale_color_gradient(low =  "#ffffff33", high = "#ff0000") + 
          xlab("Component 1")+ylab("Component 2")+ 
          theme(panel.background = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text.x = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_blank(), #legend.position="none", 
                panel.border = element_rect(colour = "black", fill=NA, size=.5))
        
        ggplot(trajectoryData, aes(x=data_dim_1, y=data_dim_2, color=gutOralSignatures, label=geo)) + 
          geom_point() + 
          #scale_color_gradient2(low = "#2c528c88", mid = "#ffffff77",high = "#ff000088") + 
          scale_color_gradient(low = "#ffffff77", high = "#ff0000") + 
          xlab("Component 1")+ylab("Component 2")+ 
          theme(panel.background = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text.x = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_blank(), #legend.position="none", 
                panel.border = element_rect(colour = "black", fill=NA, size=.5))
        
      }
      
    }
    
    checkEnteroType = F
    if (checkEnteroType) {
      stateEnteroTypeList = grepEnteroTypesFrom(samplesByState, basicMetaMapUpdated)
      getEntGeoPlot(stateEnteroTypeList)
    }
    
    diffMode = F
    if (diffMode) {
      diff_test_rest_RData = "j://Deposit/Project/2018_microbiome_atlas/trajectoryData//diff.norm.by.pseudotime.RData"
      diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
      save(diff_test_res, file=diff_test_rest_RData)
      
      diffOut = diff_test_res[,c("gene_short_name", "pval", "qval")]
      diffOut = diffOut[order(diffOut[,"pval"]),]
      diffOutClassified = diffOut[!grepl("unclassified", diffOut$gene_short_name),]
      
      head(diffOutClassified)
      
      sig_genes <- rownames(diffOutClassified[1:5,])
      
      sig_genes = c("msp_0905c",
                    "msp_1169c",
                    "msp_0656",
                    "msp_0832",
                    "msp_0187",
                    "msp_0209",
                    "msp_0224")
      sig_genes = c("msp_1341",
                    "msp_0483",
                    "msp_0956",
                    "msp_0811",
                    "msp_0006",
                    "msp_0792")
      
      sig_genes = c("msp_0216",
                    "msp_0368c",
                    "msp_1458",
                    "msp_0983",
                    "msp_0197")
      
      
      sig_genes = c("msp_0174",
                    "msp_0458",
                    "msp_0108",
                    "msp_0005",
                    "msp_0501")
      
      sig_genes = c("msp_0053",
                    "msp_1268",
                    "msp_0916",
                    "msp_0170",
                    "msp_0025")
      
      sig_genes = c("msp_0042",
                    "msp_1088",
                    "msp_0220",
                    "msp_0041",
                    "msp_0338",
                    "msp_0313")
      
      
      
      cds_subset <- cds[sig_genes,]
      plot_genes_in_pseudotime(cds_subset, color_by = "State")
      
      
    }
    
  }
  
  exceptEurMode = T
  if (exceptEurMode) {
    
    generateMode = T
    if (generateMode) {
      exceptEurSamples = newNormalAllIds[!newNormalAllIds %in% eurFiSamples]
      normAllMat = mergeMatUpdated*10^9
      normAllMat = normAllMat[, colnames(normAllMat) %in% exceptEurSamples]
      
      trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
                                       Well=basicMetaMapUpdated$dataset.ID,
                                       Hours=0,
                                       Media="NONE")
      
      rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
      trajectorySampleTab = trajectorySampleTab[rownames(trajectorySampleTab) %in% exceptEurSamples,]
      
      trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                     biotype="protein_coding",
                                     num_cells_expressed = 1,
                                     use_for_ordering=FALSE)
      
      rownames(trajectoryTaxaTab) = rownames(taxo)
      trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]
      
      pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
      fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
      
      cds <- newCellDataSet(normAllMat, 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit = 0.1)
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      cds <- reduceDimension(cds, 
                             max_components = 2,
                             method = 'DDRTree')
      
      cds <- orderCells(cds)
      
      save(cds, file=exceptEurDataTrajectoryRData)
    }
    
    loadMode = T
    if (loadMode) {
      load(exceptEurDataTrajectoryRData)
    }
    
    checkTrajectory = T
    if (checkTrajectory) {
      plot_cell_trajectory(cds, color_by = "Pseudotime")
      plot_cell_trajectory(cds, color_by = "State")
      
      
      brachInfo = pData(cds)
      samplesByState = split(as.character(brachInfo$Library), brachInfo$State)
      
      newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
      newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
      newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
      
      getPercentGeoByState <- function(targetSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) {
        currCounts = do.call(c, lapply(newNormalSamplesByGeo, function(currSamples) {
          counts = sum(currSamples %in% targetSamples )
          return(counts)
        }))
        
        return(currCounts*100/newNormalSampleCountsByGeo)
      }
      
      geoByStateMat = do.call(rbind,lapply(samplesByState, function(currSamples) getPercentGeoByState(currSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) ))
      View(geoByStateMat)
      
      barplotMode = T
      if (barplotMode) {
        firstSeqs = c("1","2","3","4","5","6","7")
        midSeqs = c("8")
        lastSeqs = c("9")
        firstBranchSamples = unlist(samplesByState[firstSeqs])
        midBranchSamples = unlist(samplesByState[midSeqs])
        lastBranchSamples = unlist(samplesByState[lastSeqs])
        
        firstBranchStat = getPercentGeoByState(firstBranchSamples,
                                               newNormalSamplesByGeo,
                                               newNormalSampleCountsByGeo)
        
        midBranchStat = getPercentGeoByState(midBranchSamples,
                                             newNormalSamplesByGeo,
                                             newNormalSampleCountsByGeo)
        
        lastBranchStat = getPercentGeoByState(lastBranchSamples,
                                              newNormalSamplesByGeo,
                                              newNormalSampleCountsByGeo)
        
        ordered = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                    "Japan","China","US", 
                    "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
        
        
        tradVec = c(255,127,0,120)/255
        tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
        
        cjpVec = c(211,211,211,120)/255
        cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
        
        eurVec = c(117,112,179,120)/255
        eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
        
        
        orderedCols = c(rep(tradCol, 7),
                        rep(cjpCol, 3), #4
                        rep(eurCol, 7))
        par(mfrow=c(4,1))
        par(mar=c(1,4,1,1))
        barplot(firstBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(midBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(lastBranchStat[ordered], col=orderedCols, las=2)
        #barplot 
        
      }
      
    }
  }
  
  exceptCjpMode = T
  if (exceptCjpMode) {
    
    generateMode = F
    if (generateMode) {
      exceptCjpSamples = newNormalAllIds[!newNormalAllIds %in% cjpSamples]
      normAllMat = mergeMatUpdated*10^9
      normAllMat = normAllMat[, colnames(normAllMat) %in% exceptCjpSamples]
      
      trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
                                       Well=basicMetaMapUpdated$dataset.ID,
                                       Hours=0,
                                       Media="NONE")
      
      rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
      trajectorySampleTab = trajectorySampleTab[rownames(trajectorySampleTab) %in% exceptCjpSamples,]
      
      trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                     biotype="protein_coding",
                                     num_cells_expressed = 1,
                                     use_for_ordering=FALSE)
      
      rownames(trajectoryTaxaTab) = rownames(taxo)
      trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]
      
      pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
      fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
      
      cds <- newCellDataSet(normAllMat, 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit = 0.1)
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      cds <- reduceDimension(cds, 
                             max_components = 2,
                             method = 'DDRTree')
      
      cds <- orderCells(cds)
      
      save(cds, file=exceptCjpDataTrajectoryRData)
      
    }
    
    loadMode = T
    if (loadMode) {
      load(exceptCjpDataTrajectoryRData)
    }
    
    checkTrajectory = T
    if (checkTrajectory) {
      plot_cell_trajectory(cds, color_by = "Pseudotime")
      plot_cell_trajectory(cds, color_by = "State")
      
      brachInfo = pData(cds)
      samplesByState = split(as.character(brachInfo$Library), brachInfo$State)
      
      newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
      newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
      newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
      
      getPercentGeoByState <- function(targetSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) {
        currCounts = do.call(c, lapply(newNormalSamplesByGeo, function(currSamples) {
          counts = sum(currSamples %in% targetSamples )
          return(counts)
        }))
        
        return(currCounts*100/newNormalSampleCountsByGeo)
      }
      
      geoByStateMat = do.call(rbind,lapply(samplesByState, function(currSamples) getPercentGeoByState(currSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) ))
      View(geoByStateMat)
      
      barplotMode = T
      if (barplotMode) {
        firstSeqs = c("1")
        midSeqs = c("2","3","4","5","6","7","8")
        lastSeqs = c("9","10","11")
        firstBranchSamples = unlist(samplesByState[firstSeqs])
        midBranchSamples = unlist(samplesByState[midSeqs])
        lastBranchSamples = unlist(samplesByState[lastSeqs])
        
        firstBranchStat = getPercentGeoByState(firstBranchSamples,
                                               newNormalSamplesByGeo,
                                               newNormalSampleCountsByGeo)
        
        midBranchStat = getPercentGeoByState(midBranchSamples,
                                             newNormalSamplesByGeo,
                                             newNormalSampleCountsByGeo)
        
        lastBranchStat = getPercentGeoByState(lastBranchSamples,
                                              newNormalSamplesByGeo,
                                              newNormalSampleCountsByGeo)
        
        ordered = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                    "Japan","China","US", 
                    "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
        
        
        tradVec = c(255,127,0,120)/255
        tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
        
        cjpVec = c(211,211,211,120)/255
        cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
        
        eurVec = c(117,112,179,120)/255
        eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
        
        
        orderedCols = c(rep(tradCol, 7),
                        rep(cjpCol, 3), #4
                        rep(eurCol, 7))
        par(mfrow=c(4,1))
        par(mar=c(1,4,1,1))
        barplot(firstBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(midBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(lastBranchStat[ordered], col=orderedCols, las=2)
        #barplot 
        
      }
      
    }
  }
  
  exceptTradMode = T
  if (exceptTradMode) {
    
    generateMode = F
    if (generateMode) {
      exceptTrdSamples = newNormalAllIds[!newNormalAllIds %in% trdSamples]
      normAllMat = mergeMatUpdated*10^9
      normAllMat = normAllMat[, colnames(normAllMat) %in% exceptTrdSamples]
      
      trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
                                       Well=basicMetaMapUpdated$dataset.ID,
                                       Hours=0,
                                       Media="NONE")
      
      rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
      trajectorySampleTab = trajectorySampleTab[rownames(trajectorySampleTab) %in% exceptTrdSamples,]
      
      trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                     biotype="protein_coding",
                                     num_cells_expressed = 1,
                                     use_for_ordering=FALSE)
      
      rownames(trajectoryTaxaTab) = rownames(taxo)
      trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]
      
      pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
      fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
      
      cds <- newCellDataSet(normAllMat, 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit = 0.1)
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      cds <- reduceDimension(cds, 
                             max_components = 2,
                             method = 'DDRTree')
      
      cds <- orderCells(cds)
      
      save(cds, file=exceptTradDataTrajectoryRData)
      
    }
    
    loadMode = T
    if (loadMode) {
      load(exceptTradDataTrajectoryRData)
    }
    
    checkTrajectory = T
    if (checkTrajectory) {
      plot_cell_trajectory(cds, color_by = "Pseudotime")
      plot_cell_trajectory(cds, color_by = "State")
      
      brachInfo = pData(cds)
      samplesByState = split(as.character(brachInfo$Library), brachInfo$State)
      
      newNormalMetadata = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds, ]
      newNormalSamplesByGeo = split(newNormalMetadata$sample.ID, newNormalMetadata$Geography)
      newNormalSampleCountsByGeo = do.call(c, lapply(newNormalSamplesByGeo, length))
      
      getPercentGeoByState <- function(targetSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) {
        currCounts = do.call(c, lapply(newNormalSamplesByGeo, function(currSamples) {
          counts = sum(currSamples %in% targetSamples )
          return(counts)
        }))
        
        return(currCounts*100/newNormalSampleCountsByGeo)
      }
      
      geoByStateMat = do.call(rbind,lapply(samplesByState, function(currSamples) getPercentGeoByState(currSamples, newNormalSamplesByGeo, newNormalSampleCountsByGeo) ))
      View(geoByStateMat)
      
      barplotMode = T
      if (barplotMode) {
        firstSeqs = c("7","8")
        midSeqs = c("2","3","4","5","6")
        lastSeqs = c("1","9")
        firstBranchSamples = unlist(samplesByState[firstSeqs])
        midBranchSamples = unlist(samplesByState[midSeqs])
        lastBranchSamples = unlist(samplesByState[lastSeqs])
        
        firstBranchStat = getPercentGeoByState(firstBranchSamples,
                                               newNormalSamplesByGeo,
                                               newNormalSampleCountsByGeo)
        
        midBranchStat = getPercentGeoByState(midBranchSamples,
                                             newNormalSamplesByGeo,
                                             newNormalSampleCountsByGeo)
        
        lastBranchStat = getPercentGeoByState(lastBranchSamples,
                                              newNormalSamplesByGeo,
                                              newNormalSampleCountsByGeo)
        
        ordered = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                    "Japan","China","US", 
                    "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
        
        
        tradVec = c(255,127,0,120)/255
        tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
        
        cjpVec = c(211,211,211,120)/255
        cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
        
        eurVec = c(117,112,179,120)/255
        eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
        
        
        orderedCols = c(rep(tradCol, 7),
                        rep(cjpCol, 3), #4
                        rep(eurCol, 7))
        par(mfrow=c(4,1))
        par(mar=c(1,4,1,1))
        barplot(firstBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(midBranchStat[ordered], col=orderedCols, las=2, axisnames = F)
        par(mar=c(1,4,1,1))
        barplot(lastBranchStat[ordered], col=orderedCols, las=2)
        #barplot 
        
      }
      
    }
    
  }
  
  wholeDataMode = F
  if (wholeDataMode) {
    
    generateMode = F
    if (generateMode) {
      trajectorySampleTab = data.frame(Library=basicMetaMapUpdated$sample.ID,
                                       Well=basicMetaMapUpdated$dataset.ID,
                                       Hours=0,
                                       Media="NONE")
      
      rownames(trajectorySampleTab) = basicMetaMapUpdated$sample.ID
      
      trajectoryTaxaTab = data.frame(gene_short_name = taxo$species, 
                                     biotype="protein_coding",
                                     num_cells_expressed = 1,
                                     use_for_ordering=FALSE)
      
      rownames(trajectoryTaxaTab) = rownames(taxo)
      trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(mergeMatUpdated), rownames(taxo)),]
      
      pd <- new("AnnotatedDataFrame", data = trajectorySampleTab)
      fd <- new("AnnotatedDataFrame", data = trajectoryTaxaTab)
      cds <- newCellDataSet(mergeMatUpdated*10^9, 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit = 0.1)
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      #cds <- detectGenes(cds, min_expr = 0.1)
      #expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 1))
      
      
      cds <- reduceDimension(cds, 
                             max_components = 2,
                             method = 'DDRTree')
      
      cds <- orderCells(cds)
    }
    
    loadMode = T
    if (loadMode) load(wholeDataTrajectoryRData)
    
    #save(cds, file=wholeDataTrajectoryRData)
    
    plot_cell_trajectory(cds, color_by = "Pseudotime")
    plot_cell_trajectory(cds, color_by = "State")
    
    
    plot_cell_trajectory(cds, color_by = "State") +
      facet_wrap(~State, nrow = 1)
    
    
    
    brachInfo = pData(cds)
    samplesByState = split(as.character(brachInfo$Library), brachInfo$State)
    babyData = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id38"]
    
    targetData = diseaseUniqSubjectList$T2D_ES_id53
    lapply(samplesByState, function(currSamples)any(currSamples %in% targetData))
    
    targetState = 2
    currStateInfo = (basicMetaMapUpdated[ basicMetaMapUpdated$sample.ID %in% rownames(brachInfo[brachInfo$State==targetState,]),])
    mean(currStateInfo$GeneRichness)
    
    sort(table(currStateInfo$Geography))
    
    
    # plot_genes_branched_pseudotime(lung[lung_genes,],
    #                                branch_point = 1,
    #                                color_by = "Time",
    #                                ncol = 1)
    
  }
  
  prunedDataMode = T
  if (prunedDataMode) {
    
  }
  
}

checkDysbiosisGenus = F
if (checkDysbiosisGenus) {
  
  corrGenusMat = cor(genusMatUpdated)
  
  boxplot(as.numeric(corrGenusMat[newDiseaseIds, newWesternIds]))
  boxplot(as.numeric(corrGenusMat[newWesternIds, newWesternIds]))
  boxplot(as.numeric(corrGenusMat[newDiseaseIds, newDiseaseIds]))
  
  corrList = list(withinIndustrial = as.numeric(corrGenusMat[newWesternIds, newWesternIds]),
                  withinTraditional = as.numeric(corrGenusMat[newNonWesternIds, newNonWesternIds]),
                  withinDiseases = as.numeric(corrGenusMat[newDiseaseIds, newDiseaseIds]))
  
  
  mean(as.numeric(corrGenusMat[newDiseaseIds, newWesternIds]))
  mean(as.numeric(corrGenusMat[newWesternIds, newWesternIds]))
  mean(as.numeric(corrGenusMat[newWesternIds, newNonWesternIds]))
  mean(as.numeric(corrGenusMat[newNonWesternIds, newNonWesternIds]))
  mean(as.numeric(corrGenusMat[newDiseaseIds, newDiseaseIds]))
  
}

### very important ###
analyzeFunctionCluster = T
if (analyzeFunctionCluster) {
  
  # funcModNetNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcModeNet.nodeTable.20190828.txt"
  # funcModNetAttachedNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcModeNet.attached.nodeTable.20190902.txt"
  # funcModNetEdgeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcModeNet.edgeTable.20190828.txt"
  # 
  # funcGemModNetNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcGemModeNet.nodeTable.20190828.txt"
  # funcGemModNetAttachedNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcGemModeNet.attached.nodeTable.20190828.txt"
  # funcGemModNetEdgeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcGemModeNet.edgeTable.20190828.txt"
  # 
  funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
  load(funcMatRData)
  
  funcModNetNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation//funcModeNet.nodeTable.20191227.txt"
  funcModNetAttachedNodeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation\\funcModeNet.attached.nodeTable.20191227.txt"
  funcModNetEdgeTableFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation\\funcModeNet.edgeTable.20191227.txt"
  
  #funcModulesRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcModules.20190828.RData"
  funcModulesRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functional.annotation//funcModules.20191227.RData"
  #funcGemModulesRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\funcGemModules.20190828.RData"
  
  checkMaxNumberOfMatchedSpecies <- function(funcModules, funcMat, suppressMsps2) {
    funcMat = funcMat[ !rownames(funcMat) %in% suppressMsps2, ]
    maxNums = lapply(funcModules, function(items){
      if(length(items)==1) {
        currMat = funcMat[,colnames(funcMat)%in% items]
        rowCount = sum(currMat>0)
        return(rowCount)
        
      }
      currMat = funcMat[,colnames(funcMat)%in% items]
      rowCounts = rowSums(currMat>0) 
      return(max(rowCounts))
    })
    return(maxNums)
  }
  maxNums = checkMaxNumberOfMatchedSpecies(funcModules, funcMat, suppressMsps2)
  smallFcs = names(maxNums[maxNums<3])
  
  length(funcModules)
  
  loadAnnots = T
  if (loadAnnots) { # seven types of functional/phenotypic annotations
    
    ### virulence factors ###
    vfMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/vfBestMat.RData"
    load(vfMatRData) #vfMatDigit
    vfTerms = colnames(vfMatDigit)
    patricVfDescTabFile = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/PATRIC_sp_gene.txt"
    patricVfDescTab = read.table(patricVfDescTabFile, header=T, sep="\t", stringsAsFactors = F)
    
    ### JGI phenotype ###
    jgiMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/jgiMat.RData"
    load(jgiMatRData)
    jgiTerms = colnames(jgiMat)
    
    ### reactions ###
    gemRData = "J://Deposit/Project/2018_microbiome_atlas//atlas.GEM.model.related/GEM.from.Reza.amazing/absentPresentInModels.RData"
    load(gemRData)
    reactionTerms = colnames(gemMat)
    
    ### antibiotic resistance ###
    mustardMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/mustard.mat.RData"
    load(mustardMatRData)
    arTerms = colnames(mustardMat)
    
    ### antismash - secondary metabolism ###
    antismashMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/msp1992.igc2.antismashMat.RData"
    load(antismashMatRData)
    antismashMat[1:3,1:3]
    antismashMelt = melt(antismashMat)
    head(antismashMelt)
    antismashMelt = antismashMelt[antismashMelt$value!=0,]
    dim(antismashMelt)
    
    writeMode = F
    if (writeMode) write.table(antismashMelt, "antismash.list.txt", sep="\t")
    antismashTerms = colnames(antismashMat)
    
    ### Cazy  ### 
    newGutCazyMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/CAZY.INRA/igc2.new.cazy.mat.RData"
    load(newGutCazyMatRData)
    colnames(cazyMat) = paste("cazy", colnames(cazyMat), sep=".")
    cazyTerms = colnames(cazyMat)
    
    ### PFAM ###
    igc2PfamMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PFAM/igc2PfamMat.RData"
    load(igc2PfamMatRData)
    pfamTerms = colnames(pfamMat)
    
    ### KEGG ###
    gutKoBestMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/gutKoBestMat.1990.RData"
    load(gutKoBestMatRData)
    koTerms = colnames(gutKoBestMat)
    
    ### KEGG pathway ###
    # already loaded from above scripts
    # data: bacteriaKoElementList
    
    ### KEGG modules ###
    keggModuleDescTabRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/KEGG module/moduleDescTab.RData"
    koModuleListRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/KEGG module/koModuleList.20190219.RData"
    load(keggModuleDescTabRData) #moduleDescTab
    load(koModuleListRData)
    microbeModuleTab = moduleDescTab[moduleDescTab$isBacteriaModule|moduleDescTab$isArchaeaModule,]
    nonMicrobeModuleTab = moduleDescTab[!(moduleDescTab$isBacteriaModule|moduleDescTab$isArchaeaModule),]
    nonMicrobeModuleTab2 = moduleDescTab[grepl("eukaryotes",moduleDescTab$moduleName),]
    
    microbeModules = paste(microbeModuleTab$moduleId, microbeModuleTab$moduleName, sep = ":")
    nonMicrobeModules = paste(nonMicrobeModuleTab$moduleId, nonMicrobeModuleTab$moduleName, sep = ":")
    nonMicrobeModules2 = paste(nonMicrobeModuleTab2$moduleId, nonMicrobeModuleTab2$moduleName, sep = ":")
    nonMicrobeModules = unique(c(nonMicrobeModules,nonMicrobeModules2))
    
    
    ### antismash ### spreaded
    antismashEnrichVec = getEnrichMentGsByList(funcModules, antismashTerms)
    antismashEnrichVec[antismashEnrichVec<=0.01]
    length(antismashEnrichVec[antismashEnrichVec<1])
    antismashClusters = names(antismashEnrichVec[antismashEnrichVec<1])
    
    ### virulence factors ###
    vfEnrichVec = getEnrichMentGsByList(funcModules, vfTerms)
    sort( vfEnrichVec[vfEnrichVec<=0.01])
    sort( vfEnrichVec[vfEnrichVec<1])
    length( sort( vfEnrichVec[vfEnrichVec<1]))
    vfClusters = names(vfEnrichVec[vfEnrichVec<1])
    
    ### antibiotic resistance ### spreaded
    arEnrichVec = getEnrichMentGsByList(funcModules, arTerms)
    sort( arEnrichVec[arEnrichVec<=0.01])
    length(sort( arEnrichVec[arEnrichVec<1]))
    arClusters = names(arEnrichVec[arEnrichVec<1])
    
    ### phenotype - JGI ### spreaded, but funcModules-149 with anaerobe, rod-shaped and non-sporulating
    jgiEnrichVec = getEnrichMentGsByList(funcModules, jgiTerms)
    sort( jgiEnrichVec[jgiEnrichVec<=0.01])
    length(sort( jgiEnrichVec[jgiEnrichVec<1]))
    jgiClusters = names(jgiEnrichVec[jgiEnrichVec<1])
    
    ### cazy ### mostly spreaded
    cazyEnrichVec = getEnrichMentGsByList(funcModules, cazyTerms)
    length(sort( cazyEnrichVec[cazyEnrichVec<1]))
    #sort( cazyEnrichVec[cazyEnrichVec<=0.01])
    cazyClusters = names(cazyEnrichVec[cazyEnrichVec<1])
    
    ### pfam ###
    pfamEnrichVec = getEnrichMentGsByList(funcModules, pfamTerms) # it contains "NA" values --> should be removed
    length(sort( pfamEnrichVec[pfamEnrichVec<1]))
    #sort( pfamEnrichVec[pfamEnrichVec<=0.01])
    pfamClusters = names(pfamEnrichVec[pfamEnrichVec<1])
    
    ### KEGG ###
    keggEnrichVec = getEnrichMentGsByList(funcModules, koTerms)
    length(sort( keggEnrichVec[keggEnrichVec<1]))
    #sort( keggEnrichVec[keggEnrichVec<=0.01])
    keggClusters = names(keggEnrichVec[keggEnrichVec<1])
    
    #pie(table(funcModuleNums==1))
    #median(funcModuleNums[funcModuleNums!=1])
    
    
    
  }
  
  loadFuncModule = T
  if (loadFuncModule) {
    
    load(funcModulesRData)
    #load(funcGemModulesRData)
    
    funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
    load(funcMatRData)
    
    funcModuleNums = do.call(c, lapply(funcModules, length))
    funcModuleSel = funcModules[funcModuleNums>=3]
    signletons = names(funcModuleNums[funcModuleNums==1])
    singletons = names(funcModuleNums[funcModuleNums==1])
    
    table(funcModuleNums>1)
    table(funcModuleNums>=3)
    
    bileAcidKeggs =  c("K01442","K00076", "K23231", "K22604", "K22605", "K22606", "K22607", "K15868", "K15871", "K15869", "K15870", "K15872", "K15873", "K15874", "K07007")
    
    checkCountMode = F
    if (checkCountMode) {
      enrichedClusters = c(4153, 4105, 339,315, 20 , 18,16)
      names(enrichedClusters) = c("KEGG",
                                  "PFAM",
                                  "virulence",
                                  "CAZyme",
                                  "phenotype",
                                  "secondary met.", 
                                  "antibiotic res.")
      
      par(mar=c(10,4,4,1))
      barplot(enrichedClusters, las=2, log="y",
              col=colorRampPalette(c("gray","white", "blue"))(8))
    }
    
  }
  
  checkFunctionTableS4 = T
  if (checkFunctionTableS4) {
    tabS4FuncFile = "C:\\Data\\comparative.analysis.healthy.sweden\\functions.in.table.s4.txt"
    tabS4FuncList = read.table(tabS4FuncFile, sep="\t", stringsAsFactors = F, header=F)[,]
    tabS4ClassList = rep(NA, length(tabS4FuncList))
    
    tabS4ClassList[tabS4FuncList %in% vfTerms] = "Virluence(PATRIC)"
    tabS4ClassList[tabS4FuncList %in% jgiTerms] = "phenotype(JGI)"
    tabS4ClassList[tabS4FuncList %in% koTerms] = "KEGG orthology"
    tabS4ClassList[tabS4FuncList %in% pfamTerms] = "protein domain(PFAM)"
    tabS4ClassList[tabS4FuncList %in% cazyTerms] = "CAZyme"
    tabS4ClassList[tabS4FuncList %in% arTerms] = "AMR(Mustard)"
    tabS4ClassList[tabS4FuncList %in% antismashTerms] = "Secondary metabolites (Anti-SMASH)"
    # write.table(data.frame(func=tabS4FuncList, class=tabS4ClassList),
    #             file="functions.class.in.table.s4.txt", 
    #             sep="\t")
    # 
  }
  
  checkTableS10 = T
  if (checkTableS10) {
    clusterName = paste("CL-",names(funcModules),sep="")
    
    mspByModList =grepMspsAssociatedWithModules(funcModules = funcModules,
                                                funcMat = funcMat[!rownames(funcMat) %in% suppressMsps2, ], 
                                                cutRatio = 0.75)
    mspLenByModList = unlist(lapply(mspByModList,length))
    clusterMsps = unlist(lapply(mspByModList, function(x) paste(x, collapse = ";")))
    clusterSpecies = unlist(lapply(mspByModList, function(x) paste(getSpeciesName(x,taxo), collapse = ";")))
    
    clusterTerms = unlist(lapply(funcModules, function(x) paste(x, collapse = ";")))
    
    koCheck = unlist(lapply(funcModules, function(x) any(x %in% koTerms)))
    pfamCheck = unlist(lapply(funcModules, function(x) any(x %in% pfamTerms)))
    vfCheck = unlist(lapply(funcModules, function(x) any(x %in% vfTerms)))
    cazyCheck = unlist(lapply(funcModules, function(x) any(x %in% cazyTerms)))
    antismashCheck = unlist(lapply(funcModules, function(x) any(x %in% antismashTerms)))
    arCheck = unlist(lapply(funcModules, function(x) any(x %in% arTerms)))
    jgiCheck = unlist(lapply(funcModules, function(x) any(x %in% jgiTerms)))
    
    
    
    koTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% koTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    
    igc2PfamDescRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PFAM/igc2PfamDesc.RData"
    load(igc2PfamDescRData)
    pfamDescTab = unique(igc2PfamDesc[,c("pfam_id","pfam_name", "pfam_desc")])
    
    getDescFromPfamTerms <- function(pfamTerms, pfamDescTab) {
      desc = pfamDescTab$pfam_desc[match(pfamTerms, pfamDescTab$pfam_name)]
      return(desc)
    }
    
    pfamTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% pfamTerms ] 
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    vfTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% vfTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    cazyTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% cazyTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    arTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% arTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    antismashTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% antismashTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    
    jgiTermsInMod = unlist(lapply(funcModules, function(x) {
      currTerms = x[x %in% jgiTerms]
      if (length(currTerms)==0) return("-")
      return(paste(currTerms, collapse = ";"))
    }))
    
    #koDesc
    #pfamDesc
    #vfProducts
    
    keggDescs = do.call(c, lapply(funcModules, function(x) {
      descs = getDescFromKoTerms(koTerms[koTerms %in% x], koDescMap)
      descs = descs[descs!=""]
      descs = descs[!is.na(descs)]
      str = getStrCatBy(descs)
      descs = unique(getItemsFromStr(str))
      out = getStrCatBy(descs)
      if (out=="") out = "-"
      return(out)
    }) )
    
    pfamDescs = do.call(c, lapply(funcModules, function(x) {
      descs = getDescFromPfamTerms(pfamTerms[pfamTerms %in% x], pfamDescTab)
      descs = descs[descs!=""]
      descs = descs[!is.na(descs)]
      str = getStrCatBy(descs)
      descs = unique(getItemsFromStr(str))
      out = getStrCatBy(descs)
      if (out=="") out = "-"
      return(out)
    }) )
    
    vfProducts = do.call(c, lapply(funcModules, function(x) {
      prods = getProductsFromVfTerms(vfTerms[vfTerms %in% x],patricVfDescTab)
      prods = prods[prods!=""]
      prods = prods[!is.na(prods)]
      str = getStrCatBy(prods)
      prods = unique(getItemsFromStr(str))
      out = getStrCatBy(prods)
      if (out=="") out = "-"
      return(out)
    }) )
    
    vfClasses = do.call(c, lapply(matchedFunc, function(x) {
      clss = getClassesFromVfTerms(vfTerms[vfTerms %in% x],patricVfDescTab)
      clss = clss[clss!=""]
      clss = clss[!is.na(clss)]
      str = getStrCatBy(clss)
      clss = unique(getItemsFromStr(str))
      out = getStrCatBy(clss)
      if (out=="") out = "-"
      return(out)
    }) )
    
    loadEnrichmentMode = T
    if (loadEnrichmentMode) {
      keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20200102.RData"
      #keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20190905.RData"
      keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20200102.RData"
      #keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20190905.RData"
      load(keggEnrichPathMatRData)
      load(keggEnrichModuleMatRData)
      
      keggEnrichPathMelt = melt(keggEnrichPathMat)
      keggEnrichPathMelt = keggEnrichPathMelt[keggEnrichPathMelt$value<=0.01,]
      keggEnrichPathMelt$Var1 = as.character(keggEnrichPathMelt$Var1)
      colnames(keggEnrichPathMelt) = c("pathway","funcModule","pvalue")
      exPaths = c("ko04626:Plant-pathogen interaction",
                  "ko05120:Epithelial cell signaling in Helicobacter pylori infection",
                  "ko04621:NOD-like receptor signaling pathway")
      keggEnrichPathMelt = keggEnrichPathMelt[!keggEnrichPathMelt$pathway %in% exPaths,]
      keggEnrichPathMelt$pathway[keggEnrichPathMelt$pathway=="ko05111:Vibrio cholerae pathogenic cycle"] <- "ko05111:Biofilm formation - Vibrio cholerae"
      
      keggEnrichModuleMelt = melt(keggEnrichModuleMat)
      keggEnrichModuleMelt = keggEnrichModuleMelt[keggEnrichModuleMelt$value<=0.01,]
      keggEnrichModuleMelt$Var1 = as.character(keggEnrichModuleMelt$Var1)
      keggEnrichModuleMelt = keggEnrichModuleMelt[keggEnrichModuleMelt$Var1 %in% microbeModules,]
      keggEnrichModuleMelt = keggEnrichModuleMelt[!keggEnrichModuleMelt$Var1 %in% nonMicrobeModules,]
      colnames(keggEnrichModuleMelt) = c("module","funcModule","pvalue")
    }
    
    keggEnModules = do.call(c, sapply(as.numeric(names(funcModules)), function(x) {
      currEnriched <- keggEnrichModuleMelt[keggEnrichModuleMelt$funcModule == x,]
      if (dim(currEnriched)[1] == 0 ) return("")
      return(getStrCatBy(currEnriched$module))
    }, simplify = F) )
    keggEnModSubSystems = do.call(c, sapply(as.numeric(names(funcModules)), function(x) {
      currEnriched <- keggEnrichModuleMelt[keggEnrichModuleMelt$funcModule == x,]
      if (dim(currEnriched)[1] == 0 ) return("")
      currModules = currEnriched$module
      currModuleIds = getIdsFromSplit(currModules)
      currSubsystems = moduleDescTab$subsystem[moduleDescTab$moduleId %in%currModuleIds]
      currSubsystems = unique(currSubsystems)
      return(getStrCatBy(currSubsystems))
    }, simplify = F) )
    keggEnPathways = do.call(c, sapply(as.numeric(names(funcModules)), function(x) {
      currEnriched <- keggEnrichPathMelt[keggEnrichPathMelt$funcModule == x,]
      if (dim(currEnriched)[1] == 0 ) return("")
      return(getStrCatBy(currEnriched$pathway))
    }, simplify = F) )
    
    normalClusters = names(funcModules)[!names(funcModules) %in% c(smallFcs, signletons)]
    normalClustersMod = paste("CL-",(normalClusters),sep="")
    View(normalClustersMod)
    
    tabS10 = data.frame(ID=clusterName,
                        size=funcModuleNums,
                        numEnrichedSpecies = mspLenByModList,
                        sigletons = ifelse(funcModuleNums==1, "Yes","No"), 
                        enrichedMsp = clusterMsps,
                        enrichedSpecies = clusterSpecies,
                        keggModules = keggEnModules,
                        keggSubsystems = keggEnModSubSystems,
                        keggPathwys = keggEnPathways,
                        koTerm = koTermsInMod,
                        pfamTerm = pfamTermsInMod,
                        virulenceTerm=vfTermsInMod,
                        cazyTerm=cazyTermsInMod,
                        amrTerm = arTermsInMod,
                        antismashTerm = antismashTermsInMod,
                        phenotypeTerm = jgiTermsInMod,
                        koTermDesc = keggDescs,
                        pfamTermDesc = pfamDescs,
                        virulenceDesc = vfProducts,
                        allTerms = clusterTerms)
    tabS10File = "C://Data//comparative.analysis.healthy.sweden//tabS10.txt"
    #write.table(tabS10, tabS10File, sep="\t", quote = F, row.names = F)
    
    tabS10File = "C://Data//comparative.analysis.healthy.sweden//tabS10.csv"
    #write.csv(tabS10, tabS10File, row.names = F)
    
  }
  
  generateMode = F
  if (generateMode) {
    ### KEGG pathways ###
    funcKoModules <- lapply(funcModules, function(x) x[x%in% koTerms])
    funcKoModules <- funcKoModules[do.call(c,lapply(funcKoModules, length))!=0]
    keggEnrichPathMat = getEnrichMentGsByGs(funcKoModules, bacteriaKoElementList)
    keggEnrichPathMelt = melt(keggEnrichPathMat)
    keggEnrichPathMelt = keggEnrichPathMelt[keggEnrichPathMelt$value<=0.05,]
    
    #keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20190905.RData"
    keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20200102.RData"
    save(keggEnrichPathMat, file=keggEnrichPathMatRData)
    
    ### KEGG modules ###
    funcKoModules <- lapply(funcModules, function(x) x[x%in% koTerms])
    funcKoModules <- funcKoModules[do.call(c,lapply(funcKoModules, length))!=0]
    keggEnrichModuleMat = getEnrichMentGsByGs(funcKoModules, koModuleList[names(koModuleList) %in% microbeModules])
    keggEnrichModuleMelt = melt(keggEnrichModuleMat)
    keggEnrichModuleMelt = keggEnrichModuleMelt[keggEnrichModuleMelt$value<=0.05,]
    
    #keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20190905.RData"
    keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20200102.RData"
    save(keggEnrichModuleMat, file=keggEnrichModuleMatRData)
    
  }
  
  
  
  attachMode = T
  if (attachMode) {
    # funcModules$`230`
    # 
    # outs <- do.call(c,lapply(funcModules, function(x) any(x %in% vfTerms)))
    # outTab = data.frame(names(funcModules),outs)
    # View(outTab[outTab[,2],])
    
    # termTest = do.call(c,lapply(funcModules, function(x)any(x%in%jgiTerms)))
    # termTest[termTest]
    # funcModules[["80"]]
    loadEnrichmentMode = T
    if (loadEnrichmentMode) {
      keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20200102.RData"
      #keggEnrichPathMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichPathMat.20190905.RData"
      keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20200102.RData"
      #keggEnrichModuleMatRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\functypeData\\keggEnrichModuleMat.20190905.RData"
      load(keggEnrichPathMatRData)
      load(keggEnrichModuleMatRData)
      
      keggEnrichPathMelt = melt(keggEnrichPathMat)
      keggEnrichPathMelt = keggEnrichPathMelt[keggEnrichPathMelt$value<=0.01,]
      keggEnrichPathMelt$Var1 = as.character(keggEnrichPathMelt$Var1)
      colnames(keggEnrichPathMelt) = c("pathway","funcModule","pvalue")
      exPaths = c("ko04626:Plant-pathogen interaction",
                  "ko05120:Epithelial cell signaling in Helicobacter pylori infection",
                  "ko04621:NOD-like receptor signaling pathway")
      keggEnrichPathMelt = keggEnrichPathMelt[!keggEnrichPathMelt$pathway %in% exPaths,]
      keggEnrichPathMelt$pathway[keggEnrichPathMelt$pathway=="ko05111:Vibrio cholerae pathogenic cycle"] <- "ko05111:Biofilm formation - Vibrio cholerae"
      
      keggEnrichModuleMelt = melt(keggEnrichModuleMat)
      keggEnrichModuleMelt = keggEnrichModuleMelt[keggEnrichModuleMelt$value<=0.01,]
      keggEnrichModuleMelt$Var1 = as.character(keggEnrichModuleMelt$Var1)
      keggEnrichModuleMelt = keggEnrichModuleMelt[keggEnrichModuleMelt$Var1 %in% microbeModules,]
      keggEnrichModuleMelt = keggEnrichModuleMelt[!keggEnrichModuleMelt$Var1 %in% nonMicrobeModules,]
      colnames(keggEnrichModuleMelt) = c("module","funcModule","pvalue")
    }
    
    ### attaching information ###
    funcModNetNodeTable = read.table(funcModNetNodeTableFile, stringsAsFactors = F, header=T, sep="\t")
    attachInfo <- function(funcModNetNodeTable,  funcModules,  antismashTerms, vfTerms, arTerms, jgiTerms, cazyTerms, pfamTerms, koTerms, koDescMap, moduleDescTab, patricVfDescTab) {
      
      
      nodes = as.character(funcModNetNodeTable$node)
      matchedFunc = funcModules[match(nodes, names(funcModules))]
      
      antismashCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% antismashTerms)))
      antismashCheck = antismashCheck*1
      antismashIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(antismashTerms[antismashTerms %in% x])) )
      
      vfCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% vfTerms)))
      vfCheck = vfCheck*1
      vfIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(vfTerms[vfTerms %in% x])) )
      vfGenes = do.call(c, lapply(matchedFunc, function(x) {
        
        genes = getGenesFromVfTerms(vfTerms[vfTerms %in% x],patricVfDescTab)
        genes = genes[genes!=""]
        genes = genes[!is.na(genes)]
        str = getStrCatBy(genes)
        genes = unique(getItemsFromStr(str))
        return(getStrCatBy(genes))
        
      }) )
      vfProducts = do.call(c, lapply(matchedFunc, function(x) {
        
        prods = getProductsFromVfTerms(vfTerms[vfTerms %in% x],patricVfDescTab)
        prods = prods[prods!=""]
        prods = prods[!is.na(prods)]
        str = getStrCatBy(prods)
        prods = unique(getItemsFromStr(str))
        return(getStrCatBy(prods))
        
      }) )
      vfClasses = do.call(c, lapply(matchedFunc, function(x) {
        
        clss = getClassesFromVfTerms(vfTerms[vfTerms %in% x],patricVfDescTab)
        clss = clss[clss!=""]
        clss = clss[!is.na(clss)]
        str = getStrCatBy(clss)
        clss = unique(getItemsFromStr(str))
        return(getStrCatBy( clss ))
        
      }) )
      
      arCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% arTerms)))
      arCheck = arCheck*1
      arIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(arTerms[arTerms %in% x])) )
      
      jgiCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% jgiTerms)))
      jgiCheck = jgiCheck*1
      jgiIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(jgiTerms[jgiTerms %in% x])) )
      
      cazyCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% cazyTerms)))
      cazyCheck = cazyCheck*1
      cazyIds = do.call(c, lapply(matchedFunc, function(x) {
        str = getStrCatBy(cazyTerms[cazyTerms %in% x])
        str = gsub("cazy.","",str)
        return(str)
      } ) )
      
      pfamCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% pfamTerms)))
      pfamCheck = pfamCheck*1
      pfamIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(pfamTerms[pfamTerms %in% x])) )
      
      keggCheck = do.call(c, lapply(matchedFunc, function(x) any(x %in% koTerms)))
      keggCheck = keggCheck*1
      keggIds = do.call(c, lapply(matchedFunc, function(x) getStrCatBy(koTerms[koTerms %in% x])) )
      keggGenes = do.call(c, lapply(matchedFunc, function(x) {
        genes = getGeneFromKoTerms(koTerms[koTerms %in% x], koDescMap)
        genes = genes[genes!=""]
        genes = genes[!is.na(genes)]
        str = getStrCatBy(genes)
        genes = unique(getItemsFromStr(str))
        return(getStrCatBy(genes))
      }) )
      
      keggDescs = do.call(c, lapply(matchedFunc, function(x) {
        descs = getDescFromKoTerms(koTerms[koTerms %in% x], koDescMap)
        descs = descs[descs!=""]
        descs = descs[!is.na(descs)]
        str = getStrCatBy(descs)
        descs = unique(getItemsFromStr(str))
        return(getStrCatBy(descs))
      }) )
      
      ### KEGG modules and pathway
      keggEnModules = do.call(c, lapply(as.list(nodes), function(x) {
        currEnriched <- keggEnrichModuleMelt[keggEnrichModuleMelt$funcModule == x,]
        if (dim(currEnriched)[1] == 0 ) return("")
        return(getStrCatBy(currEnriched$module))
      }) )
      keggEnModSubSystems = do.call(c, lapply(as.list(nodes), function(x) {
        currEnriched <- keggEnrichModuleMelt[keggEnrichModuleMelt$funcModule == x,]
        if (dim(currEnriched)[1] == 0 ) return("")
        currModules = currEnriched$module
        currModuleIds = getIdsFromSplit(currModules)
        currSubsystems = moduleDescTab$subsystem[moduleDescTab$moduleId %in%currModuleIds]
        currSubsystems = unique(currSubsystems)
        return(getStrCatBy(currSubsystems))
      }) )
      keggEnPathways = do.call(c, lapply(as.list(nodes), function(x) {
        currEnriched <- keggEnrichPathMelt[keggEnrichPathMelt$funcModule == x,]
        if (dim(currEnriched)[1] == 0 ) return("")
        return(getStrCatBy(currEnriched$pathway))
      }) )
      funcModNetAttachNodeTable = funcModNetNodeTable
      funcModNetAttachNodeTable$hasAntismash = antismashCheck
      funcModNetAttachNodeTable$hasVirulence = vfCheck
      funcModNetAttachNodeTable$hasAR = arCheck
      funcModNetAttachNodeTable$hasJgi = jgiCheck
      funcModNetAttachNodeTable$hasCazy = cazyCheck
      funcModNetAttachNodeTable$hasPfam = pfamCheck
      funcModNetAttachNodeTable$hasKegg = keggCheck
      
      funcModNetAttachNodeTable$antismahsTerms = antismashIds
      funcModNetAttachNodeTable$virulenceTerms = vfIds
      funcModNetAttachNodeTable$virulenceGenes = vfGenes
      funcModNetAttachNodeTable$virulenceProducts = vfProducts
      funcModNetAttachNodeTable$virulenceClasses = vfClasses
      funcModNetAttachNodeTable$ARTerms = arIds
      funcModNetAttachNodeTable$jgiTerms = jgiIds
      funcModNetAttachNodeTable$cazyTerms = cazyIds
      funcModNetAttachNodeTable$pfamTerms = pfamIds
      funcModNetAttachNodeTable$keggTerms = keggIds
      funcModNetAttachNodeTable$keggDescs = keggDescs
      funcModNetAttachNodeTable$keggGenes = keggGenes
      funcModNetAttachNodeTable$keggModules = keggEnModules
      funcModNetAttachNodeTable$keggSubsystems = keggEnModSubSystems
      funcModNetAttachNodeTable$keggPathway = keggEnPathways
      
      return(funcModNetAttachNodeTable)
      
    }
    
    
    funcModNetAttachNodeTable = attachInfo(funcModNetNodeTable, 
                                           funcModules, 
                                           antismashTerms, 
                                           vfTerms, 
                                           arTerms, 
                                           jgiTerms,
                                           cazyTerms,
                                           pfamTerms, 
                                           koTerms,
                                           koDescMap,
                                           moduleDescTab,
                                           patricVfDescTab)
    
    
    
    
    writeMode = F
    if (writeMode) write.table(funcModNetAttachNodeTable, funcModNetAttachedNodeTableFile, row.names = F, quote = F, sep="\t")
    
  } 
  
  countMode = T
  if (countMode) {
    pfamCheck = unlist(lapply(funcModules, function(x) any(x %in% pfamTerms)))
    koCheck = unlist(lapply(funcModules, function(x) any(x %in% koTerms)))
    antismashCheck = unlist(lapply(funcModules, function(x) any(x %in% antismashTerms)))
    cazyCheck = unlist(lapply(funcModules, function(x) any(x %in% cazyTerms)))
    vfCheck = unlist(lapply(funcModules, function(x) any(x %in% vfTerms)))
    arCheck = unlist(lapply(funcModules, function(x) any(x %in% arTerms)))
    jgiCheck = unlist(lapply(funcModules, function(x) any(x %in% jgiTerms)))
    
    table(pfamCheck)
    table(koCheck)
    table(antismashCheck)
    table(cazyCheck)
    table(vfCheck)
    table(jgiCheck)
    
    pfamClusters = names(funcModules)[pfamCheck]
    koClusters = names(funcModules)[koCheck]
    cazyClusters = names(funcModules)[cazyCheck]
    antismashClusters = names(funcModules)[antismashCheck]
    vfClusters = names(funcModules)[vfCheck]
    jgiClusters = names(funcModules)[jgiCheck]
    arClusters = names(funcModules)[arCheck]
    
    table(pfamClusters %in% singletons)
    table(koClusters %in% singletons)
    table(cazyClusters %in% singletons)
    table(antismashClusters %in% singletons)
    table(vfClusters %in% singletons)
    table(jgiClusters %in% singletons)
    table(arClusters %in% singletons)
    
    fcStatTab = data.frame(stringsAsFactors=FALSE,
                           Number = c(2901, 1204, 2957, 1197, 153, 162, 9, 9, 249, 92, 14, 6, 14, 2),
                           Category = c("PFAM", "PFAM", "KEGG", "KEGG", "CAZy", "CAZy", "antiSMASH",
                                        "antiSMASH", "virulence", "virulence", "JGI-GOLD",
                                        "JGI-GOLD", "AMR", "AMR"),
                           Type = c("Singleton", "non-singleton", "Singleton", "non-singleton",
                                    "Singleton", "non-singleton", "Singleton",
                                    "non-singleton", "Singleton", "non-singleton", "Singleton",
                                    "non-singleton", "Singleton", "non-singleton")
    )
    fcStatTab$Category = factor(fcStatTab$Category, levels = c("KEGG", "PFAM", "virulence","CAZy", "JGI-GOLD","antiSMASH","AMR"))
    
    ggplot(fcStatTab[fcStatTab$Category%in%c("KEGG","PFAM","virulence", "CAZy"),], aes(x=Category, fill=Type, y=Number)) +geom_bar(stat="identity")
    ggplot(fcStatTab[!fcStatTab$Category%in%c("KEGG","PFAM","virulence","CAZy"),], aes(x=Category, fill=Type, y=Number)) +geom_bar(stat="identity")
    
    
  }
  
  keggModuleVsFunctionalCluster = T
  if (keggModuleVsFunctionalCluster) {
    
    funcJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//funcJaccMat.20191227.RData"
    load(funcJaccMatRData)
    
    mucinCazymes = c("cazy.GT27", "cazy.GT14", "cazy.GT11", "cazy.GT10", "cazy.GH89", "cazy.GH84", "cazy.GH33", "cazy.GH20", "cazy.GH29", "cazy.GH18", "cazy.GH123", "cazy.GH109", "cazy.CBM51", "cazy.CBM50", "cazy.CBM32")
    storageCazymes = c("cazy.CBM34", "cazy.CBM48", "cazy.CBM20", "cazy.GH13", "cazy.GH133", "cazy.GH32", "cazy.GH36", "cazy.GH37", "cazy.GH77", "cazy.GH57", "cazy.GT101", "cazy.GT26", "cazy.GT3", "cazy.GT35", "cazy.GT5")
    pectinCazymes = c("cazy.CBM77","cazy.CE8","cazy.GH105","cazy.GH28", "cazy.PL1")
    
    koMicrobeModuleList = koModuleList[microbeModules]
    koMicrobeModuleNums = unlist(lapply(koMicrobeModuleList, length))
    koMicrobeModuleList = koMicrobeModuleList[koMicrobeModuleNums>=5]
    
    moduleJaccList = lapply(koMicrobeModuleList, function(kos) {
      currJaccMat = funcJaccMat[rownames(funcJaccMat)%in% kos,rownames(funcJaccMat)%in% kos]
      currJaccs = currJaccMat[lower.tri(currJaccMat)]
      return(currJaccs)
    })
    
    pathJaccList = lapply(bacteriaKoElementList, function(kos) {
      currJaccMat = funcJaccMat[rownames(funcJaccMat)%in% kos,rownames(funcJaccMat)%in% kos]
      currJaccs = currJaccMat[lower.tri(currJaccMat)]
      return(currJaccs)
    })
    
    fcJaccList = lapply(funcModules, function(kos) {
      currJaccMat = funcJaccMat[rownames(funcJaccMat)%in% kos,rownames(funcJaccMat)%in% kos]
      currJaccs = currJaccMat[lower.tri(currJaccMat)]
      return(currJaccs)
    })
    
    
    normalClusters = names(funcModules)[!names(funcModules) %in% c(smallFcs, signletons)]
    normFuncModules = funcModules[normalClusters]
    
    exModules = c("M00161:Photosystem II ", "M00163:Photosystem I ")
    
    moduleJaccMeans = unlist(lapply(moduleJaccList[!names(moduleJaccList) %in% exModules], mean))
    pathJaccMeans = unlist(lapply(pathJaccList, mean))
    fcJaccMeans = unlist(lapply(fcJaccList[normalClusters], mean))
    
    jaccCompList = list(module=moduleJaccMeans,
                        pathway=pathJaccMeans,
                        fc=fcJaccMeans)
    jaccCompMelt = melt(jaccCompList)
    jaccCompMelt[,2]=as.character(jaccCompMelt[,2])
    colnames(jaccCompMelt) = c("jacc","category")
    jaccCompMelt$category[jaccCompMelt$category=="fc"] = "Functional cluster"
    jaccCompMelt$category[jaccCompMelt$category=="module"] = "KEGG module"
    jaccCompMelt$category[jaccCompMelt$category=="pathway"] = "KEGG pathway"
    jaccCompMelt$category = factor(jaccCompMelt$category, levels = c("KEGG pathway", "KEGG module", "Functional cluster"))
    
    p=ggplot(jaccCompMelt) + geom_boxplot(aes(x=category,y=jacc)) + 
      scale_fill_manual(values = c("white", "white", "white")) +
      xlab("") + ylab("") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 60, hjust=1),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    #ggsave("jacc.kegg.functional.cluster.pdf", p, width=2.35, heigth=3.29, units = "in")
    
    
    boxplot(jaccCompList, las=1)
    
    
    jaccOutList = lapply(koMicrobeModuleList, function(kos) {
      currJaccMat = funcJaccMat[rownames(funcJaccMat)%in% kos,rownames(funcJaccMat)%in% kos]
      currJaccs = currJaccMat[lower.tri(currJaccMat)]
      return(currJaccs)
    })
    jaccOutNums = unlist(lapply(jaccOutList, length))
    jaccOutMedians = unlist(lapply(jaccOutList, median))
    jaccOutMeans = unlist(lapply(jaccOutList, mean))
    par(mar=c(5,20,4,1))
    boxplot(jaccOutList[names(sort(jaccOutMedians,T))[1:20]], cex=0.8, horizontal=T,las=1)
    
    hist(jaccOutMedians)
    hist(jaccOutMeans)
    
    targetModule = koMicrobeModuleList$`M00126:Tetrahydrofolate biosynthesis, GTP => THF `
    targetModule = koMicrobeModuleList$`M00161:Photosystem II `
    tt = unlist(lapply(funcModules, function(x) any(x %in% targetModule)))
    tt[tt]
    targetModule[targetModule %in% funcModules$`49`]
    
    boxplot(jaccOutList$`M00126:Tetrahydrofolate biosynthesis, GTP => THF `)
    
    View(jaccOutMedians)
    
    targetModule = koModuleList$`M00023:Tryptophan biosynthesis, chorismate => tryptophan `
    View(funcJaccMat[rownames(funcJaccMat) %in% targetModule, colnames(funcJaccMat) %in% targetModule])
    
    
    targetList = list(LCarnintine=funcModules$`143`,
                      Cholin=funcModules$`4071`,
                      TMAO=funcModules$`236`)
    
    targetList = list(LCarnintine=funcModules$`143`,
                      Butyrate=funcModules$`332`,
                      Acetogenesis=c(funcModules$`361`, funcModules$`354`),
                      #AcetateWL_2=funcModules$`354`,
                      BCAA=funcModules$`12`, # folate, BCAA
                      #Methan=funcModules$`20`,
                      mucin_sialidase=funcModules$`303`,
                      #Cholin=funcModules$`4071`,
                      TMAO=funcModules$`236`,
                      Propioante1=funcModules$`4`
                      #,Propioante2=funcModules$`10`
                      
    )
    targetNames = c("L-carnitine", "Butyrate", "Acetogenesis",  "BCAA", "Host-mucin(sialidase)", "Choline", "TMAO", "Propionate")
    
    # Bile1=funcModules$`12`,
    # Bile2=funcModules$`2908`,
    # Bile3=funcModules$`4534`
    
    mspByTargetModList =grepMspsAssociatedWithModules(funcModules = targetList,
                                                      funcMat = funcMat[!rownames(funcMat) %in% suppressMsps2, ], 
                                                      cutRatio = 0.75)
    
    
    
    
    
    
    
    getCircosInputForFunctionCluster <- function(mspByTargetModList, targetNames,  taxo, mode = NULL, cutRatio=0.75) {
      names(mspByTargetModList) = targetNames
      outTab = melt(mspByTargetModList)
      rownames(outTab) = paste(outTab[,1], outTab[,2])
      #return(outTab)
      outTab$species = getSpeciesName(outTab[,1], taxo)
      outTab$genus = getGenusName(outTab[,1], taxo)
      outTab$phylum = getPhylaName(outTab[,1], taxo)
      if (!is.null(mode)) {
        outTab = outTab[,c("L1",mode)] # genus, species, phylum
        colnames(outTab)=c("from","to")
        outTab = outTab[outTab$to != "unclassified",]
        outTab = outTab[!grepl("unclassified",outTab$to),]
      }
      
      return(outTab)
      
    }
    
    circosOut = getCircosInputForFunctionCluster(mspByTargetModList, targetNames, taxo , "phylum", 0.75)
    View(circosOut)
    if (F) {
      
      elements = (unique(c(circosOut$from, circosOut$to)))
      lenElements = length(elements)
      
      circosGridCols = (col_vector)[1:length(elements)]
      names(circosGridCols) = elements
      
      #circos.par(track.margin=c(0,0)) 
      set.seed(1)
      #circos.initialize(elements, xlim=c(0,1))
      #chordDiagram(circosOut)
      chordDiagram(circosOut, annotationTrack = "grid", preAllocateTracks = 1, grid.col = circosGridCols)
      circos.trackPlotRegion(track.index = 1,  panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + .1, cex = 0.5,  sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
        circos.axis(h = "top", labels = F, major.tick = F, 
                    sector.index = sector.name, track.index = 2)
      }, bg.border = NA)
      circos.clear() 
      
    }    
    
    
    funcModuleNums[1:5]
    nums = sort(mspByModNums[!names(mspByModNums) %in% signletons & !names(mspByModNums) %in% smallFcs],T)
    barplot(nums[nums>5])
    numStatTab = data.frame(clusterName = names(nums),
                            numSpecies=nums,
                            clusterSize=funcModuleNums[match(names(nums), names(funcModuleNums))],
                            inflow=(inflowsOfMspByModList)[match(names(nums), names(inflowsOfMspByModList))],
                            outflow=(outflowsOfMspByModList)[match(names(nums), names((outflowsOfMspByModList)))],
                            stringsAsFactors = F)
    numStatTab = numStatTab[numStatTab$numSpecies!=0, ]
    
    plot(numStatTab$inflow, numStatTab$numSpecies)
    plot(numStatTab$outflow, numStatTab$numSpecies)
    
    plot(numStatTab$inflow, numStatTab$clusterSize, log="y")
    plot(numStatTab$outflow, numStatTab$clusterSize, log="y")
    
    ggplot(numStatTab, aes(x=numSpecies, y=clusterSize)) + 
      geom_point(aes(fill=numSpecies, size=clusterSize), shape=21) +
      xlab("Number of enriched species") + ylab("Cluster size") +
      xlim(-10,650)+
      scale_size(range = c(2, 10))+ 
      coord_trans(y="log10")+
      #scale_fill_gradient2(mid="#ffffff88",  high = "#6a0dad88") +
      scale_fill_gradient2(mid="#ffffff88",  high = "#964B00AA") +
      theme(panel.grid.minor = element_line(colour="#80808044", size=0.5),
            panel.grid.major = element_line(colour="#80808044", size=0.5),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
    inflowCut=quantile(inflow, 0.9)
    outflowCut=quantile(outflow, 0.9)
    
    ggplot(numStatTab, aes(x=numSpecies, y=clusterSize)) + 
      #geom_point(aes(fill=numSpecies, size=clusterSize), shape=21) +
      geom_point(aes(fill=inflow >= inflowCut, size=clusterSize), shape=21) +
      xlab("Number of enriched species") + ylab("Cluster size") +
      xlim(-10,650)+
      scale_size(range = c(2, 10))+ 
      coord_trans(y="log10")+
      #scale_fill_gradient2(mid="#ffffff88",  high = "#6a0dad88") +
      scale_fill_manual(name = 'inflow > top-10%', values = setNames(c('blue','#ffffff11'),c(T, F))) +
      theme(panel.grid.minor = element_line(colour="#80808044", size=0.5),
            panel.grid.major = element_line(colour="#80808044", size=0.5),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
    ggplot(numStatTab, aes(x=numSpecies, y=clusterSize)) + 
      #geom_point(aes(fill=numSpecies, size=clusterSize), shape=21) +
      geom_point(aes(fill= outflow >= outflowCut, size=clusterSize), shape=21) +
      xlab("Number of enriched species") + ylab("Cluster size") +
      xlim(-10,650)+
      scale_size(range = c(2, 10))+ 
      coord_trans(y="log10")+
      #scale_fill_gradient2(mid="#ffffff88",  high = "#6a0dad88") +
      scale_fill_manual(name = 'outflow > top-10%', values = setNames(c('red','#ffffff11'),c(T, F))) +
      theme(panel.grid.minor = element_line(colour="#80808044", size=0.5),
            panel.grid.major = element_line(colour="#80808044", size=0.5),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
    plot(nums, funcModuleNums[match(names(nums), names(funcModuleNums))], 
         xlab="number of enriched species", 
         ylab="cluster size", log="y")
    
    
  }
  
  hgclgcList = list(hgc=hgcSpecies, 
                    lgc=lgcSpecies)
  
  hgclgcFuncAssocMat = getEnrichMentGsByGs(hgclgcList, mspByModList)
  hgclgcFuncAssocTab = data.frame(name=rownames(hgclgcFuncAssocMat), hgclgcFuncAssocMat)
  
  hgclgcFuncAssocMat = getEnrichMentGsByGs(hgclgcList, mspByTargetModList)
  
  
  getFuncMspMapping = T
  if (getFuncMspMapping) {
    mspByModList =grepMspsAssociatedWithModules(funcModules = funcModules,
                                                funcMat = funcMat[!rownames(funcMat) %in% suppressMsps2, ], 
                                                cutRatio = 0.75)
    inflowOutflowMode = F
    if (inflowOutflowMode) {
      inflowsOfMspByModList = do.call(c, lapply(mspByModList[normalClusters], function(x) mean(inflow[names(inflow)%in% x])))
      outflowsOfMspByModList = do.call(c, lapply(mspByModList[normalClusters], function(x) mean(outflow[names(outflow)%in% x])))
      
    }
    
    checkInflowOutflowEnrichment = T
    if (checkInflowOutflowEnrichment) {
      
      inflowEnrichment = getEnrichMentGsByList(mspByModList, inflowSpecies)
      outflowEnrichment = getEnrichMentGsByList(mspByModList, outflowSpecies)
      enrichmentTab = data.frame(id = names(mspByModList),
                                 inflowP = inflowEnrichment, 
                                 outflowP = outflowEnrichment, 
                                 stringsAsFactors = F)
      enrichmentNormTab = enrichmentTab[enrichmentTab$id %in% normalClusters,]
      
      
    }
    
    inflowOutflowNormalClusters = data.frame(name=names(inflowsOfMspByModList),
                                             inflow=(inflowsOfMspByModList), 
                                             outflow=(outflowsOfMspByModList),
                                             stringsAsFactors = F)
    mspByModNums = do.call(c,lapply(mspByModList, length))
    
    checkSpecifcity<- function(msps, taxo, target) {
      targetTaxa = taxo[taxo$MSP %in% msps,target]
      targetTaxa = targetTaxa[targetTaxa!="unclassified"]
      targetTaxa = targetTaxa[!grepl("unclassified",targetTaxa)]
      targetTaxaUniq = unique(targetTaxa)
      if (length(targetTaxaUniq)==1) return(targetTaxaUniq)
      return(NA)
    }
    
    checkSpecificies = unlist(lapply(mspByModList, function(x) checkSpecifcity(x, taxo, "family")))
    sort(table(checkSpecificies))
    checkSpecificies = checkSpecificies[!is.na(checkSpecificies)]
    checkSpecificies[checkSpecificies=="Prevotellaceae"]
    checkSpecificies[checkSpecificies=="Spirochaetaceae"]
    checkSpecificies[checkSpecificies=="Enterococcaceae"]
    
    checkSpecificies = unlist(lapply(mspByModList, function(x) checkSpecifcity(x, taxo, "class")))
    sort(table(checkSpecificies))
    checkSpecificies = checkSpecificies[!is.na(checkSpecificies)]
    checkSpecificies[checkSpecificies=="Verrucomicrobiae"]
    
    
    checkSpecificies = unlist(lapply(mspByModList, function(x) checkSpecifcity(x, taxo, "genus")))
    table(checkSpecificies)
    checkSpecificies = checkSpecificies[!is.na(checkSpecificies)]
    checkSpecificies[checkSpecificies=="Prevotella"]
    checkSpecificies[checkSpecificies=="Bacteroides"]
    checkSpecificies[checkSpecificies=="Veillonella"]
    checkSpecificies[checkSpecificies=="Streptococcus"]
    checkSpecificies[checkSpecificies=="Sphaerochaeta"]
    
    #mspByModList =grepMspsAssociatedWithModules(funcModules = funcModules,funcMat = funcMat, cutRatio = 0.95)
    mspNamesByModList = lapply(mspByModList, function(x) getSpeciesName(x, taxo))
    
    bestMspByModList = grepBestMspsAssociatedWithModules(funcModules = funcModules, funcMat = funcMat)
    bestScoresMspByModList = grepBestScoresMspsAssociatedWithModules(funcModules = funcModules, funcMat = funcMat)
    bestScoresMspByModList = unlist(bestScoresMspByModList)
    bestMspNamesByModList = lapply(bestMspByModList, function(x) getSpeciesName(x, taxo))
    
    mspByModSelList =grepMspsAssociatedWithModules(funcModules = funcModuleSel,funcMat = funcMat, cutRatio = 0.75)
    mspBySecMetList = grepMspsAssociatedWithModules(funcModules = funcModules[antismashClusters],funcMat = funcMat, cutRatio = 0.75)
    
    mspByModMelt = melt(mspByModList)
    mspByModMelt$value = as.character(mspByModMelt$value)
    mspByModMelt$L1 = as.character(mspByModMelt$L1)
    moduleByMspList = split(mspByModMelt[,2], mspByModMelt[,1])
    moduleNumsByMsp = do.call(c, lapply(moduleByMspList, length))
    
    sort(moduleNumsByMsp)[1:5]
    sort(moduleNumsByMsp, T)[1:5]
    
    crisprModules = c("1440","20","251","315","339","411","415","443","445","48","57","633","90","93")
    competenceModules = c("101","105","12","128","14","383","606","63","64","8")
    
    
    if (F) {
      funcModules$`332`
      funcModules$`143`
      funcModules$`236`
      funcModules$`4071`
      funcModules$`3206`
      funcModules$`354`
      funcModules$`332`
      
      K08298 
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K01687"))) # 
      
      tt = unlist(lapply(funcModules, function(x) any(x=="cazy.GH33"))) #303, sialidase
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K11264"))) #propanoyl-coA
      tt = unlist(lapply(funcModules, function(x) any(x=="K18426"))) #propanoyl-coA
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K01847"))) #propanoyl-coA
      tt = unlist(lapply(funcModules, function(x) any(x=="K01026"))) #propanoyl-coA
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K13922"))) #propionaldehyde
      
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K15023"))) #361, acetogen
      tt = unlist(lapply(funcModules, function(x) any(x=="K00198"))) #354, acetogen
      tt = unlist(lapply(funcModules, function(x) any(x=="K01067"))) #4, pyruvate acetate  
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K08299"))) #143, L-carnitine degradation
      tt = unlist(lapply(funcModules, function(x) any(x=="K08298"))) #3206
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K18277")))# None
      tt = unlist(lapply(funcModules, function(x) any(x=="K22444")))# None # l-carnitine --> TMA
      tt = unlist(lapply(funcModules, function(x) any(x=="K07811")))#236 # TMA --> TMAO
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K20038"))) #4071 # choline trimethylamine-lyase; choline --> TMA
      tt = unlist(lapply(funcModules, function(x) any(x=="Rieske")))#1086 # carnitine monooxygenase motif; carnitine --> TMA  
      tt = unlist(lapply(funcModules, function(x) any(x=="Rieske_2")))#10 # carnitine monooxygenase motif; carnitine --> TMA  
      tt = unlist(lapply(funcModules, function(x) any(x=="Ring_hydroxyl_A")))#10 # carnitine monooxygenase motif; carnitine --> TMA  
      tt = unlist(lapply(funcModules, function(x) any(x=="K01034")))
      tt = unlist(lapply(funcModules, function(x) any(x=="K20811"))) #inulin, inulosucrase
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K01667"))) #tryptophan --> indole
      tt = unlist(lapply(funcModules, function(x) any(x=="K01905"))) #  acetate
      names(tt[tt])
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K00627"))) #  acetate
      names(tt[tt])
      
      tt = unlist(lapply(funcModules, function(x) any(x=="K20626"))) #  propionate lactate-pathway
      tt = unlist(lapply(funcModules, function(x) any(x=="K00248"))) #  propionate lactate-pathway
      
      K00837
      
      tt = unlist(lapply(funcModules, function(x) any(x %in% mucinCazymes))) #12, 117, 294, 303,  ...
      tt = unlist(lapply(funcModules, function(x) any(x %in% mgePfamTerms))) #10, 12, 48, 51, ...
      tt = unlist(lapply(funcModules, function(x) any(x %in% bileAcidKeggs))) #12, 2908, 4534
      
      bileAcidKeggs[bileAcidKeggs %in% funcModules[["12"]] ]
      bileAcidKeggs[bileAcidKeggs %in% funcModules[["2908"]] ]
      bileAcidKeggs[bileAcidKeggs %in% funcModules[["4534"]] ]
      
      names(tt[tt])
      funcModules[names(tt[tt])]
      
      sort(bestMspNamesByModList$`12`)
      sort(bestMspNamesByModList$`117`)
      
      sort(bestMspNamesByModList$`53`)
      sort(bestMspNamesByModList$`132`)
      
      sort(bestMspNamesByModList$`143`)
      sort(bestMspNamesByModList$`332`)
      sort(bestMspNamesByModList$`236`)
      sort(bestMspNamesByModList$`1086`)
      sort(bestMspNamesByModList$`4071`) # 
      sort(bestMspNamesByModList$`303`) # 
      sort(bestMspNamesByModList$`4479`) # 
      sort(bestMspNamesByModList$`2288`) # acetate
      sort(bestMspNamesByModList$`2283`) # indole
      sort(bestMspNamesByModList$`1321`)
      
      sort(mspNamesByModList$`12`) # mucin-degrader??
      sort(mspNamesByModList$`143`)
      sort(mspNamesByModList$`332`)
      sort(mspNamesByModList$`236`)
      sort(mspNamesByModList$`3206`)
      sort(mspNamesByModList$`4071`)
      sort(mspNamesByModList$`1086`)
      
      sort(mspNamesByModList$`494`)
      sort(mspNamesByModList$`4`)
      
      funcModules$`152`
      funcModules$`236`
      funcModules$`2283`
      funcModules$`2288`
      funcModules$`1322`
    }
    
    
    
  }
  
  associateModuleMspMode = T
  if (associateModuleMspMode) {
    
    
    loadPair = T
    if (loadPair) {
      diseaseWithCountryPairTab = rbind(c("Atherosclerosis_SW_id1", "Sweden"),
                                        c("CVD_CN_id10", "China"),
                                        c("GDM_CN_id12", "China"),
                                        c("T2D_CN_id2", "China"),
                                        c("T1D_FI_id21", "Finland"),
                                        c("Cirrhosis_UK_id27", "UK"),
                                        c("CFS_US_id34", "US"),
                                        c("NAFLD_IT_id40", "Italy"),
                                        c("NAFLD_ES_id40", "Spain"),
                                        c("CRC_JP_id45", "Japan"),
                                        c("CRC_IT_id46", "Italy"),
                                        c("Behcet_CN_id52", "China"),
                                        c("CD_CN_id6", "China"),
                                        c("Ankylosing_CN_id9", "China"),
                                        c("CRC_US_id11", "US"),
                                        c("T2D_SW_id14", "Sweden"),
                                        c("IGT_SW_id14", "Sweden"),
                                        c("Obesity_DK_id16", "Denmark"),
                                        c("IBD_ES_id17", "Spain"),
                                        c("Cirrhosis_CN_id20", "China"),
                                        c("Obesity_DK_id25", "Denmark"),
                                        c("NSCLC_FR_id26", ""), 
                                        c("Melanoma_US_id3", "US"),
                                        c("RA_CN_id31", "China"),
                                        c("T1D_LU_id35", "Luxembourg"),
                                        c("RCC_FR_id41", ""),
                                        c("CRC_DE_id44", "Germany"),
                                        c("T2D_ES_id53", "Spain"),
                                        c("NAFLD_US_id7", "US"),
                                        c("PD_DE_id8", "Germany"))
      
      diseaseWithCountryPairTab2 = rbind(c("Atherosclerosis_SW_id1", "Sweden"),
                                         c("CVD_CN_id10", "China"),
                                         c("GDM_CN_id12", "China"),
                                         c("T2D_CN_id2", "China"),
                                         c("T1D_FI_id21", "Finland"),
                                         c("Cirrhosis_UK_id27", "UK"),
                                         c("CFS_US_id34", "US"),
                                         c("NAFLD_IT_id40", "Italy"),
                                         c("NAFLD_ES_id40", "Spain"),
                                         c("CRC_JP_id45", "Japan"),
                                         c("CRC_IT_id46", "Italy"),
                                         c("Behcet_CN_id52", "China"),
                                         c("CD_CN_id6", "China"),
                                         c("Ankylosing_CN_id9", "China"),
                                         c("CRC_US_id11", "US"),
                                         c("T2D_SW_id14", "Sweden"),
                                         c("IGT_SW_id14", "Sweden"),
                                         c("Obesity_DK_id16", "Denmark"),
                                         c("IBD_ES_id17", "Spain"),
                                         c("Cirrhosis_CN_id20", "China"),
                                         c("Obesity_DK_id25", "Denmark"),
                                         c("NSCLC_FR_id26", "Germany"), 
                                         c("Melanoma_US_id3", "US"),
                                         c("RA_CN_id31", "China"),
                                         c("T1D_LU_id35", "Luxembourg"),
                                         c("RCC_FR_id41", "Germany"),
                                         c("CRC_DE_id44", "Germany"),
                                         c("T2D_ES_id53", "Spain"),
                                         c("NAFLD_US_id7", "US"),
                                         c("PD_DE_id8", "Germany"))
      
      diseaseWithEtPairTab = rbind(c("Atherosclerosis_SW_id1", "etFirmicutes"),
                                   c("CVD_CN_id10", "etBacteroides"),
                                   c("GDM_CN_id12", "etBacteroides"),
                                   c("T2D_CN_id2", "etBacteroides"),
                                   c("T1D_FI_id21", "etFirmicutes"),
                                   c("Cirrhosis_UK_id27", "etFirmicutes"),
                                   c("CFS_US_id34", "etBacteroides"),
                                   c("NAFLD_IT_id40", "etFirmicutes"),
                                   c("NAFLD_ES_id40", "etFirmicutes"),
                                   c("CRC_JP_id45", "etBacteroides"),
                                   c("CRC_IT_id46", "etFirmicutes"),
                                   c("Behcet_CN_id52", "etBacteroides"),
                                   c("CD_CN_id6", "etBacteroides"),
                                   c("Ankylosing_CN_id9", "etBacteroides"),
                                   c("CRC_US_id11", "etBacteroides"),
                                   c("T2D_SW_id14", "etFirmicutes"),
                                   c("IGT_SW_id14", "etFirmicutes"),
                                   c("Obesity_DK_id16", "etFirmicutes"),
                                   c("IBD_ES_id17", "etFirmicutes"),
                                   c("Cirrhosis_CN_id20", "etBacteroides"),
                                   c("Obesity_DK_id25", "etFirmicutes"),
                                   c("NSCLC_FR_id26", "etFirmicutes"), 
                                   c("Melanoma_US_id3", "etBacteroides"),
                                   c("RA_CN_id31", "etBacteroides"),
                                   c("T1D_LU_id35", "etFirmicutes"),
                                   c("RCC_FR_id41", "etFirmicutes"),
                                   c("CRC_DE_id44", "etFirmicutes"),
                                   c("T2D_ES_id53", "etFirmicutes"),
                                   c("NAFLD_US_id7", "etBacteroides"),
                                   c("PD_DE_id8", "etFirmicutes"))
      
      diseaseWithMatchedPairTab = rbind(c("Atherosclerosis_SW_id1", "Sweden_id1"),
                                        c("CVD_CN_id10", "China_id10"),
                                        c("GDM_CN_id12", "China_id12"),
                                        c("T2D_CN_id2", "China_id2"),
                                        c("T1D_FI_id21", "Finland_id21"),
                                        c("Cirrhosis_UK_id27", ""),
                                        c("CFS_US_id34", "US_id34"),
                                        c("NAFLD_IT_id40", ""),
                                        c("NAFLD_ES_id40", ""),
                                        c("CRC_JP_id45", "Japan_id45"),
                                        c("CRC_IT_id46", "Italy_id46"),
                                        c("Behcet_CN_id52", "China_id52"),
                                        c("CD_CN_id6", ""),
                                        c("Ankylosing_CN_id9", "China_id9"),
                                        c("CRC_US_id11", "US_id11"),
                                        c("T2D_SW_id14", "Sweden_id14"),
                                        c("IGT_SW_id14", "Sweden_id14"),
                                        c("Obesity_DK_id16", "Denmark_id16"),
                                        c("IBD_ES_id17", "Spain_id17"),
                                        c("Cirrhosis_CN_id20", "China_id20"),
                                        c("Obesity_DK_id25", "Denmark_id25"),
                                        c("NSCLC_FR_id26", ""), 
                                        c("Melanoma_US_id3", ""),
                                        c("RA_CN_id31", "China_id31"),
                                        c("T1D_LU_id35", "Luxembourg_id35"),
                                        c("RCC_FR_id41", ""),
                                        c("CRC_DE_id44", "Germany_id44"),
                                        c("T2D_ES_id53", ""),
                                        c("NAFLD_US_id7", ""),
                                        c("PD_DE_id8", "Germany_id8"))
      diseaseWithCountryPairTab = data.frame(diseaseWithCountryPairTab, stringsAsFactors = F)
      diseaseWithMatchedPairTab = data.frame(diseaseWithMatchedPairTab, stringsAsFactors = F)
      diseaseWithEtPairTab = data.frame(diseaseWithEtPairTab, stringsAsFactors = F)
      
      colnames(diseaseWithCountryPairTab) = c("reference","target")
      colnames(diseaseWithMatchedPairTab) = c("reference","target")
      colnames(diseaseWithEtPairTab) = c("reference","target")
      
    }
    
    loadDiseaseSigMode = T
    if (loadDiseaseSigMode) {
      
      allPairMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairMelt.disease.signatures.RData"
      load(allPairMeltRData)
      allPairDownMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairDownMelt.disease.signatures.RData"
      load(allPairDownMeltRData)
      
      cohortOrderFile = "C://Data//comparative.analysis.healthy.sweden//disease.cohort.order.txt"
      cohortOrders = read.table(cohortOrderFile, stringsAsFactors = F, header=F)[,]
      names(cohortOrders)=1:30
      
      allPairMeltImputed = allPairMelt
      #allPairMeltImputed = allPairMeltImputed[!allPairMeltImputed$disease_cohort %in% c("NSCLC_FR_id26", "RCC_FR_id41", "T1D_FI_id21"),]
      allPairMeltImputed$effect_size[allPairMeltImputed$effect_size < 0] = 0
      allPairMeltImputed$effect_size[allPairMeltImputed$effect_size > 1] = 1
      allPairMeltImputed$effect_size[is.na(allPairMeltImputed$effect_size)] = 0
      
      manhattanMode = T
      if (manhattanMode) {
        uniqMsp = unique(allPairMeltImputed$msp) 
        uniqMspId = 1:length(uniqMsp)
        names(uniqMspId) =uniqMsp
        
        allPairMeltImputed$CHR = names(cohortOrders)[match(allPairMeltImputed$disease_cohort, cohortOrders)]
        allPairMeltImputed$CHR = as.numeric(allPairMeltImputed$CHR)
        #allPairMeltImputed$CHR = as.numeric(gsub(".*_id", "", allPairMeltImputed$disease_cohort))
        allPairMeltImputed$BP = uniqMspId[match(allPairMeltImputed$msp, names(uniqMspId))]
        allPairMeltImputed$P = allPairMeltImputed$effect_size
        allPairMeltImputed$SNP = allPairMeltImputed$msp
        
      }
      
      allPairMeltImputedMatched = allPairMeltImputed[allPairMeltImputed$control=="matched",]
      allPairMeltImputedCountry = allPairMeltImputed[allPairMeltImputed$control=="country",]
      allPairMeltImputedCluster = allPairMeltImputed[allPairMeltImputed$control=="cluster",]
      
      allPairMeltImputedMatchedSig = allPairMeltImputedMatched[allPairMeltImputedMatched$effect_size>=0.3,]
      allPairMeltImputedCountrySig = allPairMeltImputedCountry[allPairMeltImputedCountry$effect_size>=0.3,]
      allPairMeltImputedClusterSig = allPairMeltImputedCluster[allPairMeltImputedCluster$effect_size>=0.3,]
      
      
      allPairDownMeltImputed = allPairDownMelt
      #allPairDownMeltImputed = allPairDownMeltImputed[!allPairDownMeltImputed$disease_cohort %in% c("NSCLC_FR_id26", "RCC_FR_id41", "T1D_FI_id21"),]
      allPairDownMeltImputed$effect_size[allPairDownMeltImputed$effect_size < 0] = 0
      allPairDownMeltImputed$effect_size[allPairDownMeltImputed$effect_size > 1] = 1
      allPairDownMeltImputed$effect_size[is.na(allPairDownMeltImputed$effect_size)] = 0
      
      manhattanMode = T
      if (manhattanMode) {
        uniqMsp = unique(allPairDownMeltImputed$msp) 
        uniqMspId = 1:length(uniqMsp)
        names(uniqMspId) =uniqMsp
        allPairDownMeltImputed$CHR = names(cohortOrders)[match(allPairDownMeltImputed$disease_cohort, cohortOrders)]
        allPairDownMeltImputed$CHR = as.numeric(allPairDownMeltImputed$CHR)
        #allPairDownMeltImputed$CHR = as.numeric(gsub(".*_id", "", allPairDownMeltImputed$disease_cohort))
        allPairDownMeltImputed$BP = uniqMspId[match(allPairDownMeltImputed$msp, names(uniqMspId))]
        allPairDownMeltImputed$P = allPairDownMeltImputed$effect_size
        allPairDownMeltImputed$SNP = allPairDownMeltImputed$msp
      }
      
      allPairDownMeltImputedMatched = allPairDownMeltImputed[allPairDownMeltImputed$control=="matched",]
      allPairDownMeltImputedCountry = allPairDownMeltImputed[allPairDownMeltImputed$control=="country",]
      allPairDownMeltImputedCluster = allPairDownMeltImputed[allPairDownMeltImputed$control=="cluster",]
      
      allPairDownMeltImputedMatchedSig = allPairDownMeltImputedMatched[allPairDownMeltImputedMatched$effect_size>=0.3,]
      allPairDownMeltImputedCountrySig = allPairDownMeltImputedCountry[allPairDownMeltImputedCountry$effect_size>=0.3,]
      allPairDownMeltImputedClusterSig = allPairDownMeltImputedCluster[allPairDownMeltImputedCluster$effect_size>=0.3,]
      
      depletedMspByCohortsWCountryControl = split(allPairDownMeltImputedCountrySig$msp,
                                                  allPairDownMeltImputedCountrySig$disease_cohort)
      
      enrichedMspByCohortsWCountryControl = split(allPairMeltImputedCountrySig$msp,
                                                  allPairMeltImputedCountrySig$disease_cohort)
      
    }
    
    loadDiseaseSignatures = F
    if (loadDiseaseSignatures) {
      
      allPairMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairMelt.disease.signatures.RData"
      load(allPairMeltRData)
      allPairDownMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairDownMelt.disease.signatures.RData"
      load(allPairDownMeltRData)
      
      esCut = 0.3
      
      allPairMeltSel = allPairMelt[allPairMelt$effect_size>=esCut,]
      allPairMeltSel = allPairMeltSel[!is.na(allPairMeltSel$effect_size),]
      #allPairMeltSel = allPairMeltSel[allPairMeltSel$disease_cohort != "T1D_FI_id21",]
      diseaseUpSpecies = unique(allPairMeltSel$msp)
      
      allPairDownMeltSel = allPairDownMelt[allPairDownMelt$effect_size>=esCut,]
      allPairDownMeltSel = allPairDownMeltSel[!is.na(allPairDownMeltSel$effect_size),]
      diseaseDownSpecies = unique(allPairDownMeltSel$msp)
      
      diseaseSigList = list(up=diseaseUpSpecies, down=diseaseDownSpecies)
    }
    
    checkAssociations = T
    if (checkAssociations) {
      
      withoutSelection = T
      if (withoutSelection) {
        pCut=1e-3
        enrichSigDownDiseaseFuncMat = getEnrichMentGsByGs(mspByModList, depletedMspByCohortsWCountryControl)
        enrichSigDownDiseaseFuncMelt = melt(enrichSigDownDiseaseFuncMat)
        enrichSigDownDiseaseFuncMeltSig=enrichSigDownDiseaseFuncMelt[enrichSigDownDiseaseFuncMelt$value<=pCut, ]
        funcDownDiseaseCounts = getVectorFromTable(table(enrichSigDownDiseaseFuncMeltSig$Var2))
        
        enrichSigUpDiseaseFuncMat = getEnrichMentGsByGs(mspByModList, enrichedMspByCohortsWCountryControl)
        enrichSigUpDiseaseFuncMelt = melt(enrichSigUpDiseaseFuncMat)
        enrichSigUpDiseaseFuncMeltSig=enrichSigUpDiseaseFuncMelt[enrichSigUpDiseaseFuncMelt$value<=pCut, ]
        funcUpDiseaseCounts = getVectorFromTable(table(enrichSigUpDiseaseFuncMeltSig$Var2))
        
        allFuncNames = unique(c(names(funcDownDiseaseCounts),
                                names(funcUpDiseaseCounts)))
        
        funcClusterStatTab = data.frame(name=allFuncNames,
                                        down=funcDownDiseaseCounts[match(allFuncNames, names(funcDownDiseaseCounts))],
                                        up=funcUpDiseaseCounts[match(allFuncNames, names(funcUpDiseaseCounts))],
                                        size=funcModuleNums[match(allFuncNames, names(funcModuleNums))],
                                        stringsAsFactors = F)
        
        funcClusterStatTab[is.na(funcClusterStatTab)]=0
        #funcClusterStatTab$diff = funcClusterStatTab$down - funcClusterStatTab$up
        funcClusterStatTab$diff =  funcClusterStatTab$up - funcClusterStatTab$down
        
        funcClusterStatTab$sum = funcClusterStatTab$up + funcClusterStatTab$down
        funcClusterStatTab$fcInFewSpecies = F
        funcClusterStatTab$fcInFewSpecies[funcClusterStatTab$name %in% smallFcs] = T
        
        signletons = names(funcModuleNums[funcModuleNums==1])
        funcClusterStatTab$singleton = F
        funcClusterStatTab$singleton[funcClusterStatTab$name %in% signletons] = T
        funcClusterStatTab2 = funcClusterStatTab
        funcClusterStatTab2 = funcClusterStatTab2[!funcClusterStatTab2$fcInFewSpecies,]
        funcClusterStatTab2 = funcClusterStatTab2[!funcClusterStatTab2$singleton,]
        funcClusterStatTab2$inflowEnrichmentP = enrichmentNormTab$inflowP[match(funcClusterStatTab2$name, enrichmentNormTab$id)]
        funcClusterStatTab2$outflowEnrichmentP = enrichmentNormTab$outflowP[match(funcClusterStatTab2$name, enrichmentNormTab$id)]
        funcClusterStatTab2$inflowSig = funcClusterStatTab2$inflowEnrichmentP <= 0.01
        funcClusterStatTab2$outflowSig = funcClusterStatTab2$outflowEnrichmentP <= 0.01
        
        
        
        set.seed(1)
        ggplot(funcClusterStatTab2, aes(x=diff, y=sum)) + 
          geom_jitter(shape = 21, colour = "#00000055", width = 0.5, height = 0.5, size=2, aes(fill = outflowSig)) + 
          scale_fill_manual(values=c("#80808033", "#6a0dad88")) +
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
        
        View(funcClusterStatTab[funcClusterStatTab$name %in% antismashClusters,])
        View(funcClusterStatTab[funcClusterStatTab$name %in% vfClusters,])
        View(funcClusterStatTab[funcClusterStatTab$name %in% cazyClusters,])
        
        table(funcClusterStatTab2$outflowEnrichmentP[funcClusterStatTab2$diff >= 2 & funcClusterStatTab2$sum >= 3] <=  1e-3)
        table(funcClusterStatTab2$inflowEnrichmentP[funcClusterStatTab2$diff >= 2& funcClusterStatTab2$sum >= 3] <= 1e-3)
        
        
        table(funcClusterStatTab2$outflowEnrichmentP[funcClusterStatTab2$diff <= -2 & funcClusterStatTab2$sum >= 3] <=  1e-3)
        table(funcClusterStatTab2$inflowEnrichmentP[funcClusterStatTab2$diff <= -2 & funcClusterStatTab2$sum >= 3] <= 1e-3)
        
        
        
        
        set.seed(1)
        ggplot(funcClusterStatTab2, aes(x=inflowEnrichmentP, y=outflowEnrichmentP)) + 
          geom_point(shape = 21, aes(fill=diff, size=size)) + 
          scale_fill_gradient2(low = "blue", mid="white",  high = "red") +
          scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10') +
          scale_size(range = c(2, 10))+ 
          xlab("Inflow enrichment") + ylab("Outflow enrichment")+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position = "none",
                panel.background = element_rect(fill = "white", 
                                                colour = "black", 
                                                size = 0.5, 
                                                linetype = "solid"))
        
        
        
        ggplot(numStatTab, aes(x=numSpecies, y=clusterSize)) + 
          geom_point(aes(fill=numSpecies, size=clusterSize), shape=21) +
          xlab("Number of enriched species") + ylab("Cluster size") +
          xlim(-10,650)+
          scale_size(range = c(2, 10))+ 
          coord_trans(y="log10")+
          #scale_fill_gradient2(mid="#ffffff88",  high = "#6a0dad88") +#964B00
          scale_fill_gradient2(mid="#ffffff88",  high = "#964B0088") +
          theme(panel.grid.minor = element_line(colour="#80808044", size=0.5),
                panel.grid.major = element_line(colour="#80808044", size=0.5),
                legend.position = "none",
                panel.background = element_rect(fill = "white", 
                                                colour = "black", 
                                                size = 0.5, 
                                                linetype = "solid")) 
        
        
        
      }
      
      withSelection = T
      if (withSelection) {
        pCut=1e-3
        enrichSigDownDiseaseFuncMat = getEnrichMentGsByGs(mspByModSelList, depletedMspByCohortsWCountryControl)
        enrichSigDownDiseaseFuncMat = enrichSigDownDiseaseFuncMat[,colnames(enrichSigDownDiseaseFuncMat) %in% names(funcModuleSel)]
        enrichSigDownDiseaseFuncMelt = melt(enrichSigDownDiseaseFuncMat)
        enrichSigDownDiseaseFuncMeltSig=enrichSigDownDiseaseFuncMelt[enrichSigDownDiseaseFuncMelt$value<=pCut, ]
        #lapply(funcModules[names(sort(table(enrichSigDownDiseaseFuncMeltSig$Var2), decreasing = T)[1:10])], function(x) getDescFromKoTerms(x, koDescMap))
        funcDownDiseaseCounts = getVectorFromTable(table(enrichSigDownDiseaseFuncMeltSig$Var2))
        
        enrichSigUpDiseaseFuncMat = getEnrichMentGsByGs(mspByModSelList, enrichedMspByCohortsWCountryControl)
        enrichSigUpDiseaseFuncMat = enrichSigUpDiseaseFuncMat[,colnames(enrichSigUpDiseaseFuncMat) %in% names(funcModuleSel)]
        enrichSigUpDiseaseFuncMelt = melt(enrichSigUpDiseaseFuncMat)
        enrichSigUpDiseaseFuncMeltSig=enrichSigUpDiseaseFuncMelt[enrichSigUpDiseaseFuncMelt$value<=pCut, ]
        funcUpDiseaseCounts = getVectorFromTable(table(enrichSigUpDiseaseFuncMeltSig$Var2))
        
        allFuncNames = unique(c(names(funcDownDiseaseCounts),
                                names(funcUpDiseaseCounts)))
        
        funcClusterStatTab = data.frame(name=allFuncNames,
                                        down=funcDownDiseaseCounts[match(allFuncNames, names(funcDownDiseaseCounts))],
                                        up=funcUpDiseaseCounts[match(allFuncNames, names(funcUpDiseaseCounts))])
        funcClusterStatTab[is.na(funcClusterStatTab)]=0
        funcClusterStatTab$diff = funcClusterStatTab$down - funcClusterStatTab$up
        funcClusterStatTab$sum = funcClusterStatTab$up + funcClusterStatTab$down
        
        set.seed(1)
        ggplot(funcClusterStatTab, aes(x=diff, y=sum)) + 
          geom_jitter(shape = 21, colour = "black", width = 0.5, height = 0.5, size=2, aes(fill = diff)) + 
          scale_fill_gradient2(low = "blue", mid="white",  high = "red") +
          geom_hline(yintercept = 3, colour="gray") +
          geom_vline(xintercept = 2, colour="gray") +
          geom_vline(xintercept = -2, colour="gray") +
          geom_abline(slope=1)+
          geom_abline(slope=-1) +
          xlab("|Depletion| - |Enrichment|") + ylab("|Depletion| + |Enrichment|")+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                legend.position = "none",
                panel.background = element_rect(fill = "white", 
                                                colour = "black", 
                                                size = 0.5, 
                                                linetype = "solid"))
        
        #ggplot(funcClusterStatTab, aes(x=diff,y=sum)) + geom_point()
        
        sort(table(enrichSigUpDiseaseFuncMeltSig$Var2), decreasing = T)[1:10]
        hist(table(enrichSigUpDiseaseFuncMeltSig$Var2))
        lapply(funcModules[names(sort(table(enrichSigUpDiseaseFuncMeltSig$Var2), decreasing = T)[1:10])], function(x) getDescFromKoTerms(x, koDescMap))
        
      }
      
      ### secondary only ###
      
      pCut=1e-3
      enrichSigDownDiseaseFuncMat = getEnrichMentGsByGs(mspBySecMetList, depletedMspByCohortsWCountryControl)
      enrichSigDownDiseaseFuncMelt = melt(enrichSigDownDiseaseFuncMat)
      enrichSigDownDiseaseFuncMeltSig=enrichSigDownDiseaseFuncMelt[enrichSigDownDiseaseFuncMelt$value<=pCut, ]
      #lapply(funcModules[names(sort(table(enrichSigDownDiseaseFuncMeltSig$Var2), decreasing = T)[1:10])], function(x) getDescFromKoTerms(x, koDescMap))
      funcDownDiseaseCounts = getVectorFromTable(table(enrichSigDownDiseaseFuncMeltSig$Var2))
      
      enrichSigUpDiseaseFuncMat = getEnrichMentGsByGs(mspBySecMetList, enrichedMspByCohortsWCountryControl)
      enrichSigUpDiseaseFuncMelt = melt(enrichSigUpDiseaseFuncMat)
      enrichSigUpDiseaseFuncMeltSig=enrichSigUpDiseaseFuncMelt[enrichSigUpDiseaseFuncMelt$value<=pCut, ]
      funcUpDiseaseCounts = getVectorFromTable(table(enrichSigUpDiseaseFuncMeltSig$Var2))
      
      allFuncNames = unique(c(names(funcDownDiseaseCounts),
                              names(funcUpDiseaseCounts)))
      
      funcClusterStatTab = data.frame(name=allFuncNames,
                                      down=funcDownDiseaseCounts[match(allFuncNames, names(funcDownDiseaseCounts))],
                                      up=funcUpDiseaseCounts[match(allFuncNames, names(funcUpDiseaseCounts))])
      funcClusterStatTab[is.na(funcClusterStatTab)]=0
      funcClusterStatTab$diff = funcClusterStatTab$down - funcClusterStatTab$up
      funcClusterStatTab$sum = funcClusterStatTab$up + funcClusterStatTab$down
      
      ggplot(funcClusterStatTab, aes(x=diff, y=sum)) + 
        geom_jitter(shape = 21, colour = "black", width = 0.5, height = 0.5, size=2, aes(fill = diff)) + 
        scale_fill_gradient2(low = "blue", mid="white",  high = "red") +
        geom_hline(yintercept = 3, colour="gray") +
        geom_vline(xintercept = 2, colour="gray") +
        geom_vline(xintercept = -2, colour="gray") +
        geom_abline(slope=1)+
        geom_abline(slope=-1) +
        xlab("|Depletion| - |Enrichment|") + ylab("|Depletion| + |Enrichment|")+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      #ggplot(funcClusterStatTab, aes(x=diff,y=sum)) + geom_point()
      
      sort(table(enrichSigUpDiseaseFuncMeltSig$Var2), decreasing = T)[1:10]
      hist(table(enrichSigUpDiseaseFuncMeltSig$Var2))
      lapply(funcModules[names(sort(table(enrichSigUpDiseaseFuncMeltSig$Var2), decreasing = T)[1:10])], function(x) getDescFromKoTerms(x, koDescMap))
      
      
    }
    
    loadRichnessSignatures = T
    if (loadRichnessSignatures) {
      statsMgsForAllSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.all.txt"
      statsMgsForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.normal.txt"
      statsMgsForDiseaseSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.disease.txt"
      statsMgsForIndustrialSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.industrial.txt"
      statsMgsForTraditionalSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.traditional.txt"
      
      adjpCut = 1e-3
      fcCut = 1
      statsMgsForNorm = read.table(statsMgsForNormSamples, sep="\t", header=T, stringsAsFactors = F)
      statsMgsForNormSel = statsMgsForNorm[statsMgsForNorm$qvalue <= adjpCut,]
      #statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 0,]
      #statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < 0,]
      statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > fcCut,]
      statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < (-1)*fcCut,]
      
      
      quantCut = 0.5
      hgcCut = quantile(statsMgsForNormSelHgc$relAbd_HGC, quantCut)
      lgcCut = quantile(statsMgsForNormSelLgc$relAbd_LGC, quantCut)
      
      hgcSpecies = rownames(statsMgsForNormSelHgc[statsMgsForNormSelHgc$relAbd_HGC >= hgcCut,])
      lgcSpecies = rownames(statsMgsForNormSelLgc[statsMgsForNormSelLgc$relAbd_LGC >= lgcCut,])
      
      richnessSpecies = c(hgcSpecies, lgcSpecies)
      
      ###
      statsMgsForNormSelHgcNew = statsMgsForNormSelHgc[order(statsMgsForNormSelHgc$qvalue),]
      statsMgsForNormSelLgcNew = statsMgsForNormSelLgc[order(statsMgsForNormSelLgc$qvalue),]
      hgcTop25Species = rownames(statsMgsForNormSelHgcNew[1:25,])
      lgcTop25Species = rownames(statsMgsForNormSelLgcNew[1:25,])
      
      richnessTop25Species = c(hgcTop25Species, lgcTop25Species)
      
      hgcLgcList = list(hgc=hgcSpecies, lgc=lgcSpecies)
      hgcLgcTop25List = list(hgc=hgcTop25Species, lgc=lgcTop25Species)
      
    }
    
    loadGeoSignatures = T
    if (loadGeoSignatures) {
      
      
      
    }
    
    tt = unlist(lapply(funcModules, function(x) any(x == "K03416")))
    tt[tt]
    
    
    getFuncClusterInfo = T
    if (getFuncClusterInfo) {
      secretionClusters = c("8","120","425","250")
      mucinCazymes = c("cazy.GT27", "cazy.GT14", "cazy.GT11", "cazy.GT10", "cazy.GH89", "cazy.GH84", "cazy.GH33", "cazy.GH20", "cazy.GH29", "cazy.GH18", "cazy.GH123", "cazy.GH109", "cazy.CBM51", "cazy.CBM50", "cazy.CBM32")
      storageCazymes = c("cazy.CBM34", "cazy.CBM48", "cazy.CBM20", "cazy.GH13", "cazy.GH133", "cazy.GH32", "cazy.GH36", "cazy.GH37", "cazy.GH77", "cazy.GH57", "cazy.GT101", "cazy.GT26", "cazy.GT3", "cazy.GT35", "cazy.GT5")
      pectinCazymes = c("cazy.CBM77","cazy.CE8","cazy.GH105","cazy.GH28", "cazy.PL1")
      
      mucinInd=do.call(c,lapply(funcModules, function(x) any(x %in% mucinCazymes)))
      mucinClusters = names(mucinInd[mucinInd])
      
      storageInd=do.call(c,lapply(funcModules, function(x) any(x %in% storageCazymes)))
      storageClusters = names(storageInd[storageInd])
      
      pectinInd=do.call(c,lapply(funcModules, function(x) any(x %in% pectinCazymes)))
      pectinClusters = names(pectinInd[pectinInd])
      
      mucinClusterInfo = lapply(funcModules[mucinClusters], function(terms){
        termStr = paste(terms[terms %in% mucinCazymes], collapse =";")
        return(termStr)
      })
      mucinClusterInfoList = unlist(mucinClusterInfo)
      
      storageClusterInfo = lapply(funcModules[storageClusters], function(terms){
        termStr = paste(terms[terms %in% storageCazymes], collapse =";")
        return(termStr)
      })
      storageClusterInfoList = unlist(storageClusterInfo)
      
      pectinClusterInfo = lapply(funcModules[pectinClusters], function(terms){
        termStr = paste(terms[terms %in% pectinCazymes], collapse =";")
        return(termStr)
      })
      pectinClusterInfoList = unlist(pectinClusterInfo)
      
      antismashClusterInfo = lapply(funcModules[antismashClusters], function(terms){
        termStr = paste(terms[terms %in% antismashTerms], collapse =";")
        return(termStr)
      })
      antismashClusterInfoList = unlist(antismashClusterInfo)
      
      vfClusterInfo = lapply(funcModules[vfClusters], function(terms){
        terms = terms[terms %in% vfTerms]
        if (length(terms)==0) return(NA)
        if (length(terms)!=0) {
          terms=getClassesFromVfTerms(terms, patricVfDescTab)
          terms=unique(terms)
        }
        #if (terms =="NA") return(NA)
        termStr = paste(terms, collapse =";")
        return(termStr)
      })
      vfClusterInfoList = unlist(vfClusterInfo)
      
      cazyClusterInfo = lapply(funcModules[cazyClusters], function(terms){
        termStr = paste(terms[terms %in% cazyTerms], collapse =";")
        return(termStr)
      })
      cazyClusterInfoList = unlist(cazyClusterInfo)
      
    }
    
    analyzHgcLgcMode = T
    if (analyzHgcLgcMode) {
      moduleHgcLgcEnrichMat = getEnrichMentGsByGs(mspByModList, hgcLgcList)
      moduleHgcLgcEnrichMelt = melt(moduleHgcLgcEnrichMat)
      moduleHgcLgcEnrichMelt$Var2 = as.character(moduleHgcLgcEnrichMelt$Var2)
      moduleHgcLgcEnrichMelt = moduleHgcLgcEnrichMelt[moduleHgcLgcEnrichMelt$value<=0.05,]
      
      attachClusterInfoMode = F
      if (attachClusterInfoMode) {
        antismashMatched = do.call(c,sapply(moduleHgcLgcEnrichMelt$Var2, function(ind) {
          if (any(antismashClusters == ind)) {
            return(antismashClusterInfoList[ind])
          }
          return(NA)
        }, simplify = F))
        moduleHgcLgcEnrichMelt$antismash = antismashMatched
        
        mucinMatched = do.call(c,sapply(moduleHgcLgcEnrichMelt$Var2, function(ind) {
          if (any(mucinClusters == ind)) {
            return(mucinClusterInfoList[ind])
          }
          return(NA)
        }, simplify = F))
        moduleHgcLgcEnrichMelt$mucin = mucinMatched
        
        storageMatched = do.call(c,sapply(moduleHgcLgcEnrichMelt$Var2, function(ind) {
          if (any(storageClusters == ind)) {
            return(storageClusterInfoList[ind])
          }
          return(NA)
        }, simplify = F))
        moduleHgcLgcEnrichMelt$storage = storageMatched
        
        pectinMatched = do.call(c,sapply(moduleHgcLgcEnrichMelt$Var2, function(ind) {
          if (any(pectinClusters == ind)) {
            return(pectinClusterInfoList[ind])
          }
          return(NA)
        }, simplify = F))
        moduleHgcLgcEnrichMelt$pectin = pectinMatched
        
        
        vfMatched = do.call(c,sapply(moduleHgcLgcEnrichMelt$Var2, function(ind) {
          if (any(vfClusters == ind)) {
            return(vfClusterInfoList[ind])
          }
          return(NA)
        }, simplify = F))
        moduleHgcLgcEnrichMelt$vf = vfMatched
        
        moduleHgcLgcEnrichMeltFile = "C://Data//comparative.analysis.healthy.sweden//moduleHgcLgcEnrichMelt.txt"
        write.table(moduleHgcLgcEnrichMelt, file=moduleHgcLgcEnrichMeltFile, sep="\t")
        
        rm(moduleHgcLgcEnrichMelt)
      }
      
      if (F) {
        
        #antismashClusters
        View(moduleHgcLgcEnrichMelt[moduleHgcLgcEnrichMelt$Var2 %in% antismashClusters,])
        View(moduleHgcLgcEnrichMelt[moduleHgcLgcEnrichMelt$Var2 %in% mucinClusters,])
        View(moduleHgcLgcEnrichMelt[moduleHgcLgcEnrichMelt$Var2 %in% secretionClusters,])
        
        tempCheck = moduleHgcLgcEnrichMelt[moduleHgcLgcEnrichMelt$Var2 %in% names(antismashEnrichVec[antismashEnrichVec<=0.1]), ]
        View(tempCheck)
        tempCheck$value = -log10(tempCheck$value)
        ggplot(tempCheck, aes(x=Var1, y=value, group=Var2)) +geom_bar(stat="identity")
        
      }
      
    }
    
    moduleDiseaseSigEnrichMat = getEnrichMentGsByGs(mspByModList, diseaseSigList)
    moduleDiseaseSigEnrichMelt = melt(moduleDiseaseSigEnrichMat)
    moduleDiseaseSigEnrichMelt = moduleDiseaseSigEnrichMelt[moduleDiseaseSigEnrichMelt$value<=0.05,]
    
    ### mucin-degrading ###
    if (F) {
      View(taxo[taxo$MSP %in% mspByModList$`518`, ]) # gt20
      View(taxo[taxo$MSP %in% mspByModList$`511`, ]) # gt10
      
      View(taxo[taxo$MSP %in% mspByModList$`1356`, ]) # gt27
      View(taxo[taxo$MSP %in% mspByModList$`435`, ]) # gt11
      
      View(taxo[taxo$MSP %in% mspByModList$`294`, ]) # gh89
      View(taxo[taxo$MSP %in% mspByModList$`508`, ]) # gh84
      View(taxo[taxo$MSP %in% mspByModList$`303`, ]) # gh33
    }
    
    vfClusters = names(sort(vfEnrichVec)[1:5])
    
    vfClusterMspList = mspByModList[vfClusters]
    vfClusterMspMelt = melt(vfClusterMspList)
    vfClusterMspMelt$count = 1
    vfClusterMspMat = acast(vfClusterMspMelt, value~ L1, value.var = "count")
    vfClusterMspMat[is.na(vfClusterMspMat)] = 0
    rownames(vfClusterMspMat) = getSpeciesName(rownames(vfClusterMspMat), taxo)
    dim(vfClusterMspMat)
    
    heatmap.2(vfClusterMspMat, trace="none")
    heatmap.2(vfClusterMspMat,
              margins = c(4,12), 
              key=F, lhei=c(z,20),
              trace="none", 
              cexRow = 0.7, cexCol = 0.9,
              col=colorRampPalette(c("white","red"))(256),
              sepwidth = c(0.01,0.01),
              sepcolor = "gray",
              colsep = 0:10,
              rowsep = 0:56,
              srtCol = 45)
    
    tt=do.call(c,lapply(funcModules, function(x) any(x %in% mucinCazymes)))
    
    mucinClusters = names(tt[tt])
    mucinClusterMspList = mspByModList[mucinClusters]
    mucinClusterMspMelt = melt(mucinClusterMspList)
    mucinClusterMspMelt$count = 1
    mucinClusterMspMat = acast(mucinClusterMspMelt, value~ L1, value.var = "count")
    dim(mucinClusterMspMat)
    mucinClusterMspMat[is.na(mucinClusterMspMat)] = 0
    mucinClusterMspMatRowCounts = rowSums(mucinClusterMspMat)
    mucinClusterMspMatRowCountsWithName = mucinClusterMspMatRowCounts
    names(mucinClusterMspMatRowCountsWithName) = getSpeciesName(names(mucinClusterMspMatRowCounts),taxo)
    mucinClusterMspMatRowCountsWithName = sort(mucinClusterMspMatRowCountsWithName)
    
    par(mar=c(5,25,4,1))
    barplot(rev(rev(mucinClusterMspMatRowCountsWithName[mucinClusterMspMatRowCountsWithName>=10])), horiz = T, las=1, cex.axis = 0.7, cex.names = 0.7)
    
    tt=do.call(c,lapply(funcModules, function(x) any(x %in% storageCazymes)))
    storageClusters = names(tt[tt])
    storageClusterMspList = mspByModList[storageClusters]
    storageClusterMspMelt = melt(storageClusterMspList)
    storageClusterMspMelt$count = 1
    storageClusterMspMat = acast(storageClusterMspMelt, value~ L1, value.var = "count")
    dim(storageClusterMspMat)
    storageClusterMspMat[is.na(storageClusterMspMat)] = 0
    storageClusterMspMatRowCounts = rowSums(storageClusterMspMat)
    storageClusterMspMatRowCountsWithName = storageClusterMspMatRowCounts
    names(storageClusterMspMatRowCountsWithName) = getSpeciesName(names(storageClusterMspMatRowCounts),taxo)
    storageClusterMspMatRowCountsWithName = sort(storageClusterMspMatRowCountsWithName)
    
    par(mar=c(5,25,4,1))
    barplot(rev(rev(storageClusterMspMatRowCountsWithName[storageClusterMspMatRowCountsWithName>=10])), horiz = T, las=1, cex.axis = 0.7, cex.names = 0.7)
    
    modLengths = do.call(c, lapply(funcModules, length))
    modLengths = sort(modLengths)
    
    
  }
  
  checkKoMode = F
  if (checkKoMode) {
    targetModules = findModulesBy(funcModules, c("K00132","K04072","K04073","K18366"))
    targetModules = findModulesBy(funcModules, c("K00100"))
    targetModules = findModulesBy(funcModules, c("K00634"))
    targetModules = findModulesBy(funcModules, c("K01034"))
    targetModules = findModulesBy(funcModules, c("K01035"))
    targetModules = findModulesBy(funcModules, c("K19709"))
    targetModules = findModulesBy(funcModules, c("K01034","K01035","K19709"))
    targetModules
    
    getDescFromKoTerms(funcKoModules[["132"]], koDescMap)
    unlist(funcKoModules[targetModules])
    
    
  }
  
  checkModuleMode = T
  if (checkModuleMode) {
    butryCluster132 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "132")
    View(butryCluster132$wholeData)
    
    
    butryCluster332 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "332")
    View(butryCluster332$wholeData)
    
    butryCluster621 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "621")
    View(butryCluster621$wholeData)
    
    butryCluster1241 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "1241")
    View(butryCluster1241$wholeData)
    
    butryCluster1539 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "1539")
    View(butryCluster1539$wholeData)
    
    butryCluster3404 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "3404")
    View(butryCluster3404$wholeData)
    
    
    butryCluster_132_1539_3404 = grepMspsByKOs(unlist(funcKoModules[targetModules]), gutKoBestMat, taxo)
    View(butryCluster_132_1539_3404$wholeData)
    
    macroCluster4930 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "4930")
    View(macroCluster4930$wholeData)
    
    macroCluster4969 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "4969")
    View(macroCluster4969$wholeData)
    
    
    
    acetateCluster361 = grepMspsByFuncModules(funcModules, gutKoBestMat, taxo, "361")
    acetateCluster354 = grepMspsByFuncModules(funcModules, gutKoBestMat, taxo, "354")
    
    propionateCluster551 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "551")
    propionateCluster620 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "620")
    propionateCluster609 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "609")
    propionateCluster100 = grepMspsByFuncModules(funcKoModules, gutKoBestMat, taxo, "100")
    
    check = checkFractionOfGenus(gutTaxo, butryCluster332$msps, T)
    check2 = checkFractionOfGenus(gutTaxo, acetateCluster354$msps, T)
    check3 = checkFractionOfGenus(gutTaxo, butryCluster1241$msps, T)
    
    zz =do.call(c, lapply(funcGemModules, function(x) any(x %in% funcModules[["332"]])))
    
    View(acetateCluster354$wholeData[acetateCluster354$msps,])   
  }
  
  associateSamplesWithCluster = T
  if (associateSamplesWithCluster) {
    
    generateMode = F
    if (generateMode) {
      funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
      load(funcMatRData)
      
      if (F) {
        naTests = unlist(lapply(funcModules, function(x) any(x %in% "NA")))
        funcModules = funcModules[!names(funcModules) %in%c("1766","7749")]
        save(funcModules, file = funcModulesRData)
        
      }
      
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
      bsCorrMatRData = "j://Deposit/Project/2018_microbiome_atlas//functypeData//bsCorrMat.RData"
      save(bsCorrMat, file=bsCorrMatRData)
    }
    
    loadMode = T
    if (loadMode) {
      
      bsCorrMatRData = "j://Deposit/Project/2018_microbiome_atlas//functypeData//bsCorrMat.RData"
      load(bsCorrMatRData)
      
      length(funcModules)
      dim(bsCorrMat)
      rownames(bsCorrMat) = names(funcModules)
      
      bsCorrMaxRows = rowMaxs(bsCorrMat)
      names(bsCorrMaxRows)= rownames(bsCorrMat)
      
      bsCorrMelt = melt(bsCorrMat)
      bsCorrMeltCut = bsCorrMelt[bsCorrMelt$value>=0.3,]
      bsCorrMeltCut$Var1 = as.character(bsCorrMeltCut$Var1)
      bsCorrMeltCut$Var2 = as.character(bsCorrMeltCut$Var2)
      sigSamplesByCluster = split(bsCorrMeltCut$Var2, bsCorrMeltCut$Var1)
      
      colVec = rep("#80808011", 5224)
      names(colVec) = rownames(pcoaGenusMat)
      colVec[(sigSamplesByCluster[["367"]])] = "#FF0000"
      
      plot(pcoaGenusMat[,1], 
           pcoaGenusMat[,2], col=colVec, pch=16)
      
      
      
    }
    
  }
  
  ### attach annotations
  
}

addEnteroType = T
if (addEnteroType) {
  
  loadEnteroLists = T
  if (loadEnteroLists) {
    
    
    babyList = list(baby.finland = newNormalList$Finland_id21,
                    baby.T1D = newDiseaseList$T1D_FI,
                    preterm = newPerturbCaseList$Preterm_US)
    
    
    newCaseContByGeoList= list(Sweden_D=c(newDiseaseList$Atherosclerosis_SE, 
                                          newDiseaseList$T2D_SE),
                               Sweden_N=c(newNormalList$Sweden_id28,
                                          newNormalList$Sweden_id1,
                                          newNormalList$Sweden_id14),
                               Denmark_D = c(newDiseaseList$Obesity_DK_16,
                                             newDiseaseList$Obesity_DK_25),
                               Denmark_N=c(newNormalList$Denmark_id16, 
                                           newNormalList$Denmark_id25),
                               UK_D=c(newDiseaseList$Cirrhosis_GB),
                               UK_N=c(newNormalList$UK_id29,
                                      newNormalList$UK_id39),
                               Spain_D=c(newDiseaseList$IBD_ES,
                                         newDiseaseList$T2D_ES),
                               Spain_N=c(newNormalList$Spain_id17),
                               Italy_D=c(newDiseaseList$NAFLD_IT,
                                         newDiseaseList$CRC_IT),
                               Italy_N=c(newNormalList$Italy_id37,
                                         newNormalList$Italy_id46),
                               Germany_D=c(newDiseaseList$CRC_DE,
                                           newDiseaseList$Parkinson_DE),
                               Germany_N=c(newNormalList$Germany_id44,
                                           newNormalList$Germany_id8),
                               Luxembourg_D=c(newDiseaseList$T1D_LU),
                               Luxembourg_N=c(newNormalList$Luxembourg_id35),
                               Japan_D=c(newDiseaseList$CRC_JP),
                               Japan_N=c(newNormalList$Japan_id45),
                               US_D=c(newDiseaseList$CRC_US,
                                      newDiseaseList$Melanoma_US,
                                      newDiseaseList$CFS_US,
                                      newDiseaseList$NAFLD_US),
                               US_N=c(newNormalList$US_id11,
                                      newNormalList$US_id32,
                                      newNormalList$US_id34,
                                      newNormalList$US_id36,
                                      newNormalList$US_id43),
                               China_D=c(newDiseaseList$CVD_CN,
                                         newDiseaseList$GDM_CN,
                                         newDiseaseList$T2D_CN,
                                         newDiseaseList$Cirrhosis_CN,
                                         newDiseaseList$Rheumatoid_CN,
                                         newDiseaseList$Crohn_CN,
                                         newDiseaseList$Ankylosing_CN,
                                         newDiseaseList$Behcet_CN),
                               China_N=c(newNormalList$China_id10,
                                         newNormalList$China_id12,
                                         newNormalList$China_id2,
                                         newNormalList$China_id20,
                                         newNormalList$China_id31,
                                         newNormalList$China_id52,
                                         newNormalList$China_id6,
                                         newNormalList$China_id9))
    
    newCaseContByEnteroList = grepEnteroTypesFrom(rev(newCaseContByGeoList), basicMetaMapUpdated) 
    targeCountry="Sweden"
    getEntTab(newCaseContByGeoList[paste(targeCountry, c("_N","_D"), sep="")])
    entGeoTab = getEntGeoPlot(newCaseContByEnteroList, F)
    chisq.test(entGeoTab[paste(targeCountry, c("_N","_D"), sep=""),])
    chisq.test(entGeoTab[c("Sweden_N","Sweden_D"),])
    
    newNormalEnteroList = grepEnteroTypesFrom(c(newNormalList,newNonWesternList), basicMetaMapUpdated) 
    newDiseaseEnteroList = grepEnteroTypesFrom(newDiseaseList, basicMetaMapUpdated) 
    newPerturbationEnteroList = grepEnteroTypesFrom(newPerturbCaseList, basicMetaMapUpdated) 
    babyEnteroList = grepEnteroTypesFrom(babyList, basicMetaMapUpdated) 
    
    newGeoEnteroList = grepEnteroTypesFrom(newGeoList, basicMetaMapUpdated) 
    newContEnteroList = grepEnteroTypesFrom(newContByEnteroList, basicMetaMapUpdated) 
    newCaseEnteroList = grepEnteroTypesFrom(newCaseByEnteroList, basicMetaMapUpdated) 
    diseaseContComp = list(ETB_norm = newContEnteroList$ETB,
                           ETB_case = newCaseEnteroList$ETB,
                           ETF_norm = newContEnteroList$ETF,
                           ETF_case = newCaseEnteroList$ETF)
  }
  
  
  overallCaseContComparison = T
  if (overallCaseContComparison) {
    entGeoTab = getEntGeoPlot(diseaseContComp)
  }
  
  normalByGeoComparison = T
  if (normalByGeoComparison) {
  }
  
  caseContByGeoComparison = F
  if (caseContByGeoComparison) {
    names(newGeoList)
    names(newDiseaseList)
  }
  
  
  thaiList = list(thai=normUniqCountryList$Thailand,
                  shortImmigrant=thaiShort,
                  longImmigrant=thaiLong,
                  us=normUniqCountryList$US)
  
  entGeoTab = getEntGeoPlot(grepEnteroTypesFrom(thaiList, basicMetaMapUpdated))
  
  entGeoTab = getEntGeoPlot(newContEnteroList)
  entGeoTab = getEntGeoPlot(newGeoEnteroList)
  entGeoTab = getEntGeoPlot(newGeoEnteroList) #[c(1,2,3,4,5,6,7,8,9,11,10,12:17)]
  entGeoTab = getEntGeoPlot(newGeoEnteroList) #[c(1,2,4)]
  entGeoTab = getEntGeoPlot(newNormalEnteroList)
  entGeoTab = getEntGeoPlot(newDiseaseEnteroList, T, c(10,4,4,1))
  entGeoTab = getEntGeoPlot(newPerturbationEnteroList)
  entGeoTab = getEntGeoPlot(babyEnteroList, ordered = F)
  
  etTypesByRichnessMode = T
  if (etTypesByRichnessMode) {
    diseaseMetaData = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newDiseaseIds,]
    etfMetaData = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newEtfCountriesSamples,]
    etbMetaData = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newEtbCountriesSamples,]
    etpMetaData = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newEtpCountriesSamples,]
    
    entTabByRichness = getBinTab(basicMetaMapUpdated)
    #entTabByRichness = getBinTab(diseaseMetaData)
    entTabByRichness = getBinTab(etfMetaData)
    
    etfCol = "#77e2c799"
    etbCol = "#bee24599"
    etpCol = "#f35b38"
    
    tradVec = c(255,127,0,120)/255
    tradVec = c(255,127,0,200)/255
    tradCol = rgb(tradVec[1],tradVec[2],tradVec[3],tradVec[4])
    
    cjpVec = c(211,211,211,120)/255
    cjpVec = c(211,211,211,200)/255
    cjpCol = rgb(cjpVec[1],cjpVec[2],cjpVec[3],cjpVec[4])
    
    eurVec = c(117,112,179,120)/255
    eurVec = c(117,112,179,200)/255
    eurCol = rgb(eurVec[1],eurVec[2],eurVec[3],eurVec[4])
    
    etpCol = tradCol
    etbCol = cjpCol
    etfCol = eurCol
    
    
    
    colnames(entTabByRichness) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
    rownames(entTabByRichness) = NULL
    #par(mar=c(10,10,10,10))
    par(mar=c(5,4,4,1))
    barplot(t(entTabByRichness),
            col=c(etfCol,etbCol,etpCol), las=2)
    legend("topright", 
           bty="n",
           inset = c(1.2,0),xpd=T,cex=0.7,
           legend=c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella"),
           fill=c(etfCol,etbCol,etpCol))
  }
  
  etTypesByBmiMode = T
  if (etTypesByBmiMode) {
    
    View(basicMetaMapUpdated)
    
    imputedInd <- basicMetaMapUpdated$Imputation=="Imputed"
    nonImputedInd <- basicMetaMapUpdated$Imputation=="None"
    bmiInds <- basicMetaMapUpdated$WhatImputed  %in% c("AGX","AXX","XXX")
    bmiNaInds <- is.na(basicMetaMapUpdated$BMI)
    
    noImputedMetaMap = basicMetaMapUpdated[!bmiNaInds & nonImputedInd,]
    bmiOnlyMetaMap = basicMetaMapUpdated[!bmiNaInds & bmiInds & imputedInd, ]
    
    bmiAllMetaMap = rbind(noImputedMetaMap, bmiOnlyMetaMap)
    
    bmiEtfAllMetaMap = bmiAllMetaMap[bmiAllMetaMap$sample.ID %in% etfSamples, ]
    bmiEtbAllMetaMap = bmiAllMetaMap[bmiAllMetaMap$sample.ID %in% etbSamples, ]
    bmiEtpAllMetaMap = bmiAllMetaMap[bmiAllMetaMap$sample.ID %in% etpSamples, ]
    
    entTabByBMI = getBinTab(bmiEtpAllMetaMap, by="BMI")
    
    colnames(entTabByBMI) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
    rownames(entTabByBMI) = NULL
    
    par(mar=c(5,4,4,1))
    barplot(t(entTabByBMI),
            col=c("#77e2c799","#bee24599","#f35b38"), las=2)
    legend("topright", 
           bty="n",
           inset = c(1.2,0),xpd=T,cex=0.7,
           legend=c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella"),
           fill=c("#77e2c799","#bee24599","#f35b38"))
    
    
  }
  
  
  
}

checkEnteroType = T
if (checkEnteroType) {
  
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
  
  
  ### MAIN FIGURE FOR ET-types by geography ###
  geographyWithEnterotypeMode = T
  if (geographyWithEnterotypeMode) {
    
    genusMode = T
    if (genusMode) {
      
      notIndiGeos = basicMetaMapUpdated$Geography[match(newNormalIds, basicMetaMapUpdated$sample.ID)]
      notIndiEnterotype = basicMetaMapUpdated$enteroType[match(newNormalIds, basicMetaMapUpdated$sample.ID)]
      
      notIndiEnterptypeByGeos = split(notIndiEnterotype,notIndiGeos)
      entGeoTab = do.call(rbind,lapply(notIndiEnterptypeByGeos,  function(enterotypes){
        entCounts = getEntCount(enterotypes)
        return(entCounts)
      }))
      colnames(entGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
      entGeoTab = sweep(entGeoTab, MARGIN=1, rowSums(entGeoTab)/100, FUN="/")
      entGeoTab = entGeoTab[order(entGeoTab[,1]),]
      par(mar=c(8,3,2,1))
      seqs = c("Finland","Japan","China","US",
               "Luxembourg","Germany","Italy","Spain","Denmark","UK","Sweden")
      barplot(t(entGeoTab[seqs,]), #2:11, #t(entGeoTab[c(2,4,3,5,6,7,8,9),])
              col=c(etfCol,etbCol,etpCol), las=2)
      
      indiGeos = basicMetaMapUpdated$Geography[match(newNonWesternIds, basicMetaMapUpdated$sample.ID)]
      indiEnterotype = basicMetaMapUpdated$enteroType[match(newNonWesternIds, basicMetaMapUpdated$sample.ID)]
      
      indiEnterptypeByGeos = split(indiEnterotype,indiGeos)
      indiEntGeoTab = do.call(rbind,lapply(indiEnterptypeByGeos,  function(enterotypes){
        entCounts = getEntCount(enterotypes)
        return(entCounts)
      }))
      colnames(indiEntGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
      indiEntGeoTab = sweep(indiEntGeoTab, MARGIN=1, rowSums(indiEntGeoTab)/100, FUN="/")
      indiEntGeoTab = indiEntGeoTab[order(indiEntGeoTab[,2]),]
      par(mar=c(8,3,2,1))
      barplot(t(indiEntGeoTab),
              col=c(etfCol,etbCol,etpCol), las=2)
      allEntGeoTab = rbind(indiEntGeoTab, entGeoTab[seqs,])
      par(mar=c(8,3,2,1))
      barplot(t(allEntGeoTab),
              col=c(etfCol,etbCol,etpCol), las=2)
      
    }
    
    mgsMode = F
    if (mgsMode) {
      notIndiGeos = basicMetaMapUpdated$Geography[match(newNormalIds, basicMetaMapUpdated$sample.ID)]
      notIndiEnterotype = basicMetaMapUpdated$enteroTypeM[match(newNormalIds, basicMetaMapUpdated$sample.ID)]
      
      notIndiEnterptypeByGeos = split(notIndiEnterotype,notIndiGeos)
      entGeoTab = do.call(rbind,lapply(notIndiEnterptypeByGeos,  function(enterotypes){
        entCounts = getEntCount(enterotypes)
        return(entCounts)
      }))
      colnames(entGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
      entGeoTab = sweep(entGeoTab, MARGIN=1, rowSums(entGeoTab)/100, FUN="/")
      entGeoTab = entGeoTab[order(entGeoTab[,1]),]
      par(mar=c(8,4,4,1))
      barplot(t(entGeoTab[2:11,]), #t(entGeoTab[c(2,4,3,5,6,7,8,9),])
              col=c("#77e2c799","#bee24599","#f35b38"), las=2)
      
      indiGeos = basicMetaMapUpdated$Geography[match(newNonWesternIds, basicMetaMapUpdated$sample.ID)]
      indiEnterotype = basicMetaMapUpdated$enteroTypeM[match(newNonWesternIds, basicMetaMapUpdated$sample.ID)]
      
      indiEnterptypeByGeos = split(indiEnterotype,indiGeos)
      indiEntGeoTab = do.call(rbind,lapply(indiEnterptypeByGeos,  function(enterotypes){
        entCounts = getEntCount(enterotypes)
        return(entCounts)
      }))
      colnames(indiEntGeoTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
      indiEntGeoTab = sweep(indiEntGeoTab, MARGIN=1, rowSums(indiEntGeoTab)/100, FUN="/")
      indiEntGeoTab = indiEntGeoTab[order(indiEntGeoTab[,2]),]
      par(mar=c(8,4,4,1))
      barplot(t(indiEntGeoTab),
              col=c("#77e2c799","#bee24599","#f35b38"), las=2)
      
      
    }
    
  }
  
  checkNextTime = F
  if (checkNextTime) {
    
    ### SUPPLEMENTARY FIGURES ###
    geographyAndDiseaseByEnterotypeMode = F
    if (geographyAndDiseaseByEnterotypeMode) {
      
      etfNormalSamples = etFirmicutesCountriesSamples[etFirmicutesCountriesSamples  %in% normalIDs]
      etbNormalSamples = etBacteroidesCountriesSamples[etBacteroidesCountriesSamples  %in% normalIDs]
      etfDiseaseSamples = etFirmicutesCountriesSamples[etFirmicutesCountriesSamples  %in% caseIDs]
      etbDiseaseSamples = etBacteroidesCountriesSamples[etBacteroidesCountriesSamples  %in% caseIDs]
      
      etfNormalEnterotypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% etfNormalSamples]
      etbNormalEnterotypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% etbNormalSamples]
      etfDiseaseEnterotypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% etfDiseaseSamples]
      etbDiseaseEnterotypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% etbDiseaseSamples]
      
      etfNormalCounts = getEntCount(etfNormalEnterotypes)
      etbNormalCounts = getEntCount(etbNormalEnterotypes)
      etfDiseaseCounts = getEntCount(etfDiseaseEnterotypes)
      etbDiseaseCounts = getEntCount(etbDiseaseEnterotypes)
      
      etDiseaseTab = rbind(ETF_norm=etfNormalCounts,
                           ETF_case=etfDiseaseCounts,
                           ETB_norm=etbNormalCounts,
                           ETB_case=etbDiseaseCounts)
      
      colnames(etDiseaseTab) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
      etDiseaseTab = sweep(etDiseaseTab, MARGIN=1, rowSums(etDiseaseTab)/100, FUN="/")
      
      barplot(t(etDiseaseTab),
              col=c("#77e2c799","#bee24599","#f35b38"), las=2)
      
    }
    
    ### MAIN FIGURE FOR ET-types by richness ###
    richnessBinningWithEnterotypeMode = T
    if (richnessBinningWithEnterotypeMode) {
      
      loadMode = T
      if (loadMode) {
        diseaseMetaMap = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newDiseaseIds,]
        normalMetaMap = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNormalAllIds,]
        westernizedMetaMap = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newWesternIds,]
        indiMetaMap = basicMetaMapUpdated[basicMetaMapUpdated$sample.ID %in% newNonWesternIds,]
        
        
      }
      histMode = F
      if (histMode) {
        
        plot(density(normalMetaMap$GeneRichness), 
             ylim=c(0,2.5e-6), xaxt="n", yaxt="n", main="", xlab="")
        lines(density(diseaseMetaMap$GeneRichness), col="red") 
      }
      
      
      etTypesByRichnessMode = T
      if (etTypesByRichnessMode) {
        
        allMode = T
        diseaseMode = T
        normalMode = T
        westernizedMode = T
        nonWesternizedMode = T
        if (allMode) {
          allEntTabByRichnes = getBinTab(basicMetaMapUpdated)
          
          colnames(allEntTabByRichnes) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
          rownames(allEntTabByRichnes) = NULL
          #par(mar=c(10,10,10,10))
          par(mar=c(5,4,4,1))
          barplot(t(allEntTabByRichnes),
                  col=c(eurCol,cjpCol,tradCol), las=2)
        }
        if (diseaseMode) {
          diseaseEntTabByRichnes = getBinTab(diseaseMetaMap)
          
          colnames(diseaseEntTabByRichnes) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
          rownames(diseaseEntTabByRichnes) = NULL
          #par(mar=c(10,10,10,10))
          par(mar=c(5,4,4,1))
          barplot(t(diseaseEntTabByRichnes),
                  col=c(eurCol,cjpCol,tradCol), las=2)
          # legend("topright", 
          #        bty="n",
          #        inset = c(1.2,0),xpd=T,cex=0.7,
          #        legend=c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella"),
          #        fill=c("#77e2c799","#bee24599","#f35b38"))
          
        }
        if (normalMode) {
          normalEntTabByRichnes = getBinTab(normalMetaMap)
          
          colnames(normalEntTabByRichnes) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
          rownames(normalEntTabByRichnes) = NULL
          #par(mar=c(10,10,10,10))
          par(mar=c(5,4,4,1))
          barplot(t(normalEntTabByRichnes),
                  col=c(eurCol,cjpCol,tradCol), las=2)
          # legend("topright", 
          #        bty="n",
          #        inset = c(1.2,0),xpd=T,cex=0.7,
          #        legend=c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella"),
          #        fill=c("#77e2c799","#bee24599","#f35b38"))
        }
        if (westernizedMode) {
          westernizedEntTabByRichnes = getBinTab(westernizedMetaMap)
          
          colnames(westernizedEntTabByRichnes) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
          rownames(westernizedEntTabByRichnes) = NULL
          #par(mar=c(10,10,10,10))
          par(mar=c(5,4,4,1))
          barplot(t(westernizedEntTabByRichnes),
                  col=c(eurCol,cjpCol,tradCol), las=2)
        }
        if (nonWesternizedMode) {
          indiEntTabByRichnes = getBinTab(indiMetaMap)
          
          colnames(indiEntTabByRichnes) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
          rownames(indiEntTabByRichnes) = NULL
          #par(mar=c(10,10,10,10))
          par(mar=c(5,4,4,1))
          barplot(t(indiEntTabByRichnes),
                  col=c(eurCol,cjpCol,tradCol), las=2)
        }
        
        
      }
      
      etTypesByRichnessMode = T
      if (etTypesByRichnessMode) {
        
      }
      
      
      
      
      
      
    }
    
    checkOverallEnterotype = F
    if (checkOverallEnterotype) {
      
      ### BREIF CHECK OVERALL ENTEROTYPES ###
      enterotypeCount = do.call(rbind, lapply(diseaseSets, function(datasetId) {
        
        caseEnteroTypes = basicMetaMap[basicMetaMap$dataset.ID==datasetId & basicMetaMap$type=="Case", "enteroType"]
        contEnteroTypes = basicMetaMap[basicMetaMap$dataset.ID==datasetId & basicMetaMap$type=="Control", "enteroType"]
        
        caseEtf <- sum(caseEnteroTypes == "ET-Firmicutes")
        caseEtb <- sum(caseEnteroTypes == "ET-Bacteroides")
        caseEtp <- sum(caseEnteroTypes == "ET-Prevotella")
        
        contEtf <- sum(contEnteroTypes == "ET-Firmicutes")
        contEtb <- sum(contEnteroTypes == "ET-Bacteroides")
        contEtp <- sum(contEnteroTypes == "ET-Prevotella")
        
        resEnteroTypes = c(caseEtf, caseEtb, caseEtp, contEtf, contEtb, contEtp)
        return(resEnteroTypes)
        
      }))
      colnames(enterotypeCount) = c("CaseEtF", "CaseEtB", "CaseEtP", "ContEtF", "ContEtB", "ContEtP")
      rownames(enterotypeCount) =  diseaseSets
      
      out = apply(enterotypeCount, 1, function(currRow) {
        currTab = makeCurrTab(currRow)
        chisqOut = chisq.test(currTab)
        pvalue = chisqOut$p.value
        pvalue = ifelse(is.nan(pvalue), 1, pvalue)
        return(pvalue)
      })
      sort(out[out<0.05])
      diseaseEnteroTypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% caseIDs2]
      normalEnteroTypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% normalIDs]
      indiEnteroTypes = basicMetaMap$enteroType[basicMetaMap$sample.ID %in% indiIDs]
      
      overallEnteroTypes = list(disease=diseaseEnteroTypes, 
                                normal = normalEnteroTypes, 
                                indigenous = indiEnteroTypes)
      
      overallEnteroTypeCount = do.call(rbind, lapply(overallEnteroTypes, function(currTypes) {
        caseEtf <- sum(currTypes == "ET-Firmicutes")
        caseEtb <- sum(currTypes == "ET-Bacteroides")
        caseEtp <- sum(currTypes == "ET-Prevotella")
        
        resEnteroTypes = c(caseEtf, caseEtb, caseEtp)
        return(resEnteroTypes)
      }))
      colnames(overallEnteroTypeCount) = c("ET-Firmicutes", "ET-Bacteroides", "ET-Prevotella")
      overallEnteroTypeCount = t(overallEnteroTypeCount)
      overallEnteroTypeCount = sweep(overallEnteroTypeCount, MARGIN=2, FUN="/", STAT=colSums(overallEnteroTypeCount))
      
    }
    
    if (F) {
      diseaseDatasetIds = mgsInfo$ID[mgsInfo$Classification=="disease"]
    } 
    
    enterotypeByPerturbation = T
    if (enterotypeByPerturbation) {
      immuneDataSets = c("id26","id41")
      enterotypeCountImmune = do.call(rbind, lapply(immuneDataSets, function(datasetId) {
        
        caseEnteroTypes = basicMetaMap[basicMetaMap$dataset.ID==datasetId & tolower(basicMetaMap$subtype) == "after", "enteroType"]
        contEnteroTypes = basicMetaMap[basicMetaMap$dataset.ID==datasetId & tolower(basicMetaMap$subtype) == "before", "enteroType"]
        
        caseEtf <- sum(caseEnteroTypes == "ET-Firmicutes")
        caseEtb <- sum(caseEnteroTypes == "ET-Bacteroides")
        caseEtp <- sum(caseEnteroTypes == "ET-Prevotella")
        
        contEtf <- sum(contEnteroTypes == "ET-Firmicutes")
        contEtb <- sum(contEnteroTypes == "ET-Bacteroides")
        contEtp <- sum(contEnteroTypes == "ET-Prevotella")
        
        resEnteroTypes = c(caseEtf, caseEtb, caseEtp, contEtf, contEtb, contEtp)
        return(resEnteroTypes)
        
      }))
      colnames(enterotypeCountImmune) = c("CaseEtF", "CaseEtB", "CaseEtP", "ContEtF", "ContEtB", "ContEtP")
      rownames(enterotypeCountImmune) =  immuneDataSets
      
      outImmune = apply(enterotypeCountImmune, 1, function(currRow) {
        currTab = makeCurrTab(currRow)
        chisqOut = chisq.test(currTab)
        pvalue = chisqOut$p.value
        pvalue = ifelse(is.nan(pvalue), 1, pvalue)
        return(pvalue)
      })
      sort(outImmune[outImmune<0.05])
    }
    
  }
  
}

updatePcoaOut = T
if (updatePcoaOut) {
  
  genusMatUpdated = getTaxaSumMat(mergeMatUpdated, taxo, "genus")
  familyMatUpdated = getTaxaSumMat(mergeMatUpdated, taxo, "family")
  orderMatUpdated = getTaxaSumMat(mergeMatUpdated, taxo, "order")
  classMatUpdated = getTaxaSumMat(mergeMatUpdated, taxo, "class")
  phylumMatUpdated = getTaxaSumMat(mergeMatUpdated, taxo, "phylum")
  
  mgsCs = rowSums(mergeMatUpdated)
  genusCs = rowSums(genusMatUpdated)
  familyCs = rowSums(familyMatUpdated)
  orderCs = rowSums(orderMatUpdated)
  classCs = rowSums(classMatUpdated)
  phylumCs = rowSums(phylumMatUpdated)
  
  table(mgsCs==0)
  table(genusCs==0)
  table(familyCs==0)
  table(orderCs==0)
  table(classCs==0)
  table(phylumCs==0)
  
  vegOut=vegdist(t(mergeMatUpdated),"bray")
  vegGenusOut=vegdist(t(genusMatUpdated[genusCs!=0,]),"bray")
  vegFamilyOut=vegdist(t(familyMatUpdated[familyCs!=0,]),"bray")
  vegOrderOut=vegdist(t(orderMatUpdated),"bray")
  vegClassOut=vegdist(t(classMatUpdated),"bray")
  vegPhylumOut=vegdist(t(phylumMatUpdated),"bray")
  
  pcoaOut=pcoa(vegOut)
  pcoaGenusOut=pcoa(vegGenusOut)
  pcoaFamilyOut=pcoa(vegFamilyOut)
  pcoaOrderOut=pcoa(vegOrderOut)
  pcoaClassOut=pcoa(vegClassOut)
  pcoaPhylumOut=pcoa(vegPhylumOut) 
  
  pcoaAllUpdatedData = "C:/Data/comparative.analysis.healthy.sweden/pcoaAllData.updated.20190805.RData"
  #load(pcoaAllUpdatedData)
  
  
  save(pcoaOut, 
       pcoaGenusOut, 
       pcoaFamilyOut,
       pcoaOrderOut, 
       pcoaClassOut, 
       pcoaPhylumOut, 
       file=pcoaAllUpdatedData)
  
  pcoaGenusMat = pcoaGenusOut$vectors[,1:2]
  
  pcoaCols = rep("#c0c0c088", dim(pcoaGenusMat)[1])
  names(pcoaCols) = rownames(pcoaGenusMat)
  pcoaCols[newNormalIds] = "#00bfff88"
  pcoaCols[newNonWesternIds] = "#ff240077"
  #pcoaCols[allEtbKeySamples] = "#80008077"
  plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16) #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
  
  pcoaCols = rep("#c0c0c088", dim(pcoaGenusMat)[1])
  names(pcoaCols) = rownames(pcoaGenusMat)
  #pcoaCols[normalIDs] = "#00bfff88"
  pcoaCols[caseList$Case.id48.US] = "#ff240077"
  #pcoaCols[allEtbKeySamples] = "#80008077"
  plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16) #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
  
}

richnessBoxPlotMode = T
if (richnessBoxPlotMode) {
  
  loadSampleInfo = T
  if (loadSampleInfo) {
    
    newCaseContByGeoList = list(Sweden_D=c(newDiseaseList$Atherosclerosis_SE, 
                                           newDiseaseList$T2D_SE),
                                Sweden_N=c(newNormalList$Sweden_id28,
                                           newNormalList$Sweden_id1,
                                           newNormalList$Sweden_id14),
                                Denmark_D = c(newDiseaseList$Obesity_DK_16,
                                              newDiseaseList$Obesity_DK_25),
                                Denmark_N=c(newNormalList$Denmark_id16, 
                                            newNormalList$Denmark_id25),
                                UK_D=c(newDiseaseList$Cirrhosis_GB),
                                UK_N=c(newNormalList$UK_id29,
                                       newNormalList$UK_id39),
                                Spain_D=c(newDiseaseList$IBD_ES,
                                          newDiseaseList$T2D_ES),
                                Spain_N=c(newNormalList$Spain_id17),
                                Italy_D=c(newDiseaseList$NAFLD_IT,
                                          newDiseaseList$CRC_IT),
                                Italy_N=c(newNormalList$Italy_id37,
                                          newNormalList$Italy_id46),
                                Germany_D=c(newDiseaseList$CRC_DE,
                                            newDiseaseList$Parkinson_DE),
                                Germany_N=c(newNormalList$Germany_id44,
                                            newNormalList$Germany_id8),
                                Luxembourg_D=c(newDiseaseList$T1D_LU),
                                Luxembourg_N=c(newNormalList$Luxembourg_id35),
                                Japan_D=c(newDiseaseList$CRC_JP),
                                Japan_N=c(newNormalList$Japan_id45),
                                US_D=c(newDiseaseList$CRC_US,
                                       newDiseaseList$Melanoma_US,
                                       newDiseaseList$CFS_US,
                                       newDiseaseList$NAFLD_US),
                                US_N=c(newNormalList$US_id11,
                                       newNormalList$US_id32,
                                       newNormalList$US_id34,
                                       newNormalList$US_id36,
                                       newNormalList$US_id43),
                                China_D=c(newDiseaseList$CVD_CN,
                                          newDiseaseList$GDM_CN,
                                          newDiseaseList$T2D_CN,
                                          newDiseaseList$Cirrhosis_CN,
                                          newDiseaseList$Rheumatoid_CN,
                                          newDiseaseList$Crohn_CN,
                                          newDiseaseList$Ankylosing_CN,
                                          newDiseaseList$Behcet_CN),
                                China_N=c(newNormalList$China_id10,
                                          newNormalList$China_id12,
                                          newNormalList$China_id2,
                                          newNormalList$China_id20,
                                          newNormalList$China_id31,
                                          newNormalList$China_id52,
                                          newNormalList$China_id6,
                                          newNormalList$China_id9))
    
    newCaseContByGeoListNew = list(Sweden = list(Atherosclerosis= newDiseaseList$Atherosclerosis_SE,
                                                 T2D=newDiseaseList$T2D_SE,
                                                 IGT=diseaseUniqSubjectList$IGT_SW_id14,
                                                 Normal=c(newNormalList$Sweden_id28,
                                                          newNormalList$Sweden_id1,
                                                          newNormalList$Sweden_id14)),
                                   Denmark = list(Obesity=newDiseaseList$Obesity_DK_16,
                                                  Obesity=newDiseaseList$Obesity_DK_25,
                                                  Normal=c(newNormalList$Denmark_id16, 
                                                           newNormalList$Denmark_id25)),
                                   UK = list(Cirrhosis=newDiseaseList$Cirrhosis_GB,
                                             Normal = c(newNormalList$UK_id29,
                                                        newNormalList$UK_id39)),
                                   Spain = list(IBD=newDiseaseList$IBD_ES, 
                                                T2D=newDiseaseList$T2D_ES,
                                                Normal=newNormalList$Spain_id17),
                                   Italy = list(NAFLD=newDiseaseList$NAFLD_IT,
                                                CRC=newDiseaseList$CRC_IT,
                                                Normal=c(newNormalList$Italy_id37,
                                                         newNormalList$Italy_id46)),
                                   Germany = list(CRC=newDiseaseList$CRC_DE,
                                                  Parkinson=newDiseaseList$Parkinson_DE,
                                                  Normal=c(newNormalList$Germany_id44,
                                                           newNormalList$Germany_id8)),
                                   Finland = list(T1D = newDiseaseList$T1D_FI,
                                                  Normal = newNormalList$Finland_id21),
                                   Luxembourg = list(T1D=newDiseaseList$T1D_LU,
                                                     Normal=newNormalList$Luxembourg_id35),
                                   Japan = list(CRC=newDiseaseList$CRC_JP,
                                                Normal=newNormalList$Japan_id45),
                                   US = list(CRC=newDiseaseList$CRC_US,
                                             Melanoma=newDiseaseList$Melanoma_US,
                                             CFS=newDiseaseList$CFS_US,
                                             NAFLD=newDiseaseList$NAFLD_US,
                                             Normal=c(newNormalList$US_id11,
                                                      newNormalList$US_id32,
                                                      newNormalList$US_id34,
                                                      newNormalList$US_id36,
                                                      newNormalList$US_id43)),
                                   China = list(CVD=newDiseaseList$CVD_CN,
                                                GDM=newDiseaseList$GDM_CN,
                                                T2D=newDiseaseList$T2D_CN,
                                                Cirrhosis=newDiseaseList$Cirrhosis_CN,
                                                Rhematoid=newDiseaseList$Rheumatoid_CN,
                                                Crohn=newDiseaseList$Crohn_CN,
                                                Ankylosing=newDiseaseList$Ankylosing_CN,
                                                Behcet=newDiseaseList$Behcet_CN,
                                                Normal=c(newNormalList$China_id10,
                                                         newNormalList$China_id12,
                                                         newNormalList$China_id2,
                                                         newNormalList$China_id20,
                                                         newNormalList$China_id31,
                                                         newNormalList$China_id52,
                                                         newNormalList$China_id6,
                                                         newNormalList$China_id9)))
    
    newCaseContByGeoListNewWithFoodborne = list(Sweden = list(Atherosclerosis= newDiseaseList$Atherosclerosis_SE,
                                                              T2D=newDiseaseList$T2D_SE,
                                                              IGT=diseaseUniqSubjectList$IGT_SW_id14,
                                                              Normal=c(newNormalList$Sweden_id28,
                                                                       newNormalList$Sweden_id1,
                                                                       newNormalList$Sweden_id14)),
                                                Denmark = list(Obesity=newDiseaseList$Obesity_DK_16,
                                                               Obesity=newDiseaseList$Obesity_DK_25,
                                                               Antibiotics = newPerturbCaseList$antibiotics_DK,
                                                               Normal=c(newNormalList$Denmark_id16, 
                                                                        newNormalList$Denmark_id25)),
                                                UK = list(Cirrhosis=newDiseaseList$Cirrhosis_GB,
                                                          Normal = c(newNormalList$UK_id29,
                                                                     newNormalList$UK_id39)),
                                                Spain = list(IBD=newDiseaseList$IBD_ES, 
                                                             T2D=newDiseaseList$T2D_ES,
                                                             Normal=newNormalList$Spain_id17),
                                                Italy = list(NAFLD=newDiseaseList$NAFLD_IT,
                                                             CRC=newDiseaseList$CRC_IT,
                                                             Normal=c(newNormalList$Italy_id37,
                                                                      newNormalList$Italy_id46)),
                                                Germany = list(CRC=newDiseaseList$CRC_DE,
                                                               Parkinson=newDiseaseList$Parkinson_DE,
                                                               Normal=c(newNormalList$Germany_id44,
                                                                        newNormalList$Germany_id8)),
                                                Finland = list(T1D = newDiseaseList$T1D_FI,
                                                               Normal = newNormalList$Finland_id21),
                                                Luxembourg = list(T1D=newDiseaseList$T1D_LU,
                                                                  Normal=newNormalList$Luxembourg_id35),
                                                France = list(RCC=newDiseaseList$RCC_FR,
                                                              NSCLC=newDiseaseList$NSCLC_FR),
                                                
                                                Japan = list(CRC=newDiseaseList$CRC_JP,
                                                             Normal=newNormalList$Japan_id45),
                                                US = list(CRC=newDiseaseList$CRC_US,
                                                          Melanoma=newDiseaseList$Melanoma_US,
                                                          CFS=newDiseaseList$CFS_US,
                                                          NAFLD=newDiseaseList$NAFLD_US,
                                                          Food=newPerturbCaseList$foodborne_US,
                                                          Preterm=newPerturbCaseList$Preterm_US,
                                                          Normal=c(newNormalList$US_id11,
                                                                   newNormalList$US_id32,
                                                                   newNormalList$US_id34,
                                                                   newNormalList$US_id36,
                                                                   newNormalList$US_id43)),
                                                China = list(CVD=newDiseaseList$CVD_CN,
                                                             GDM=newDiseaseList$GDM_CN,
                                                             T2D=newDiseaseList$T2D_CN,
                                                             Cirrhosis=newDiseaseList$Cirrhosis_CN,
                                                             Rhematoid=newDiseaseList$Rheumatoid_CN,
                                                             Crohn=newDiseaseList$Crohn_CN,
                                                             Ankylosing=newDiseaseList$Ankylosing_CN,
                                                             Behcet=newDiseaseList$Behcet_CN,
                                                             Normal=c(newNormalList$China_id10,
                                                                      newNormalList$China_id12,
                                                                      newNormalList$China_id2,
                                                                      newNormalList$China_id20,
                                                                      newNormalList$China_id31,
                                                                      newNormalList$China_id52,
                                                                      newNormalList$China_id6,
                                                                      newNormalList$China_id9)))
    
    if (F) {
      newContByGeoListNew = list(Sweden = c(newNormalList$Sweden_id28,
                                            newNormalList$Sweden_id1,
                                            newNormalList$Sweden_id14),
                                 Denmark = c(newNormalList$Denmark_id16, 
                                             newNormalList$Denmark_id25),
                                 UK = c(newNormalList$UK_id29,
                                        newNormalList$UK_id39),
                                 Spain = newNormalList$Spain_id17,
                                 Italy = c(newNormalList$Italy_id37,
                                           newNormalList$Italy_id46),
                                 Germany = c(newNormalList$Germany_id44,
                                             newNormalList$Germany_id8),
                                 Finland = newNormalList$Finland_id21,
                                 Luxembourg = newNormalList$Luxembourg_id35,
                                 Japan = newNormalList$Japan_id45,
                                 US = c(newNormalList$US_id11,
                                        newNormalList$US_id32,
                                        newNormalList$US_id34,
                                        newNormalList$US_id36,
                                        newNormalList$US_id43),
                                 China = c(newNormalList$China_id10,
                                           newNormalList$China_id12,
                                           newNormalList$China_id2,
                                           newNormalList$China_id20,
                                           newNormalList$China_id31,
                                           newNormalList$China_id52,
                                           newNormalList$China_id6,
                                           newNormalList$China_id9),
                                 Peru=newNonWesternList$Peru_id32,
                                 Mongolia=newNonWesternList$Mongo_id33,
                                 Thailand = newNonWesternList$Thai_id43,
                                 Fiji = newNonWesternList$Fiji_id49,
                                 India = newNonWesternList$India_id50,
                                 Tanzania = newNonWesternList$Tanzania_id37,
                                 Madagascar = newNonWesternList$Madagascar_id42
      )
      
    }
    
    
    
  }
  
  newCirclePlotMode = T
  if (newCirclePlotMode) {
    #newCaseContRichnessNewByGeoList = grepGeneRichnesssFromNew(newCaseContByGeoListNew, basicMetaMapUpdated)
    newCaseContRichnessNewByGeoList = grepGeneRichnesssFromNew(newCaseContByGeoListNewWithFoodborne, basicMetaMapUpdated)
    #newCaseContRichnessNewByGeoList = grepGeneRichnesssFromNew(newContByGeoListNew, basicMetaMapUpdated)
    
    newCaseContMeanRichnessNewByGeoList = lapply(newCaseContRichnessNewByGeoList, function(cohortList){
      meanRichness = lapply(cohortList, function(currRichness){
        return(mean(currRichness))
      })
      return(meanRichness)
    })
    newCaseContMeanRichnessNewByGeoMelt = melt(newCaseContMeanRichnessNewByGeoList)
    
    newCaseContCvRichnessNewByGeoList = lapply(newCaseContRichnessNewByGeoList, function(cohortList){
      meanRichness = lapply(cohortList, function(currRichness){
        meanRichness = mean(currRichness)
        sdRichness = sd(currRichness)
        return(sdRichness/meanRichness)
      })
      return(meanRichness)
    })
    newCaseContCvRichnessNewByGeoMelt = melt(newCaseContCvRichnessNewByGeoList)
    
    newCaseContStatsRichnessNewByGeoTab = data.frame(Geo=newCaseContMeanRichnessNewByGeoMelt$L1, 
                                                     Cohort=newCaseContMeanRichnessNewByGeoMelt$L2,
                                                     Mean=newCaseContMeanRichnessNewByGeoMelt$value,
                                                     CV=newCaseContCvRichnessNewByGeoMelt$value)
    newCaseContStatsRichnessNewByGeoTab$Geo = factor(newCaseContStatsRichnessNewByGeoTab$Geo, 
                                                     levels=c("Denmark", 
                                                              "Sweden",
                                                              "Spain",
                                                              "Luxembourg",
                                                              "Germany",
                                                              "UK",
                                                              "Italy",
                                                              "France",
                                                              "US",
                                                              "China",
                                                              "Japan",
                                                              "Finland"))
    newCaseContStatsRichnessNewByGeoTab$Class = "Disease"
    newCaseContStatsRichnessNewByGeoTab$Class[newCaseContStatsRichnessNewByGeoTab$Cohort=="Normal"]="Normal"
    newCaseContStatsRichnessNewByGeoTab$Class = factor(newCaseContStatsRichnessNewByGeoTab$Class)
    disCol="#FFD70044"
    
    ggplot(newCaseContStatsRichnessNewByGeoTab, aes(x=Geo, y=Mean, size=CV)) +
      geom_point(aes(fill=factor(Class)), colour="black", shape=21) +
      xlab("") + ylab("") +
      scale_size(range = c(2, 10))+
      scale_fill_manual(values=c(Normal=eurCol, Disease=disCol)) +
      theme(panel.grid.major = element_line(colour="gray"), 
            axis.text.x = element_text(angle = 60, hjust=1),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
  }
  
  oldBoxplotMode = F
  if (oldBoxplotMode) {
    
    newCaseContRichnessByGeoList = grepGeneRichnesssFrom(newCaseContByGeoList, basicMetaMapUpdated)
    newCaseContRichnessByGeoTab = melt(newCaseContRichnessByGeoList)
    boxConSeqs = c("Denmark_N",
                   "Denmark_D",
                   "Sweden_N",
                   "Sweden_D",
                   "Spain_N",
                   "Spain_D",
                   "Luxembourg_N",
                   "Luxembourg_D",
                   "Germany_N",
                   "Germany_D",
                   "UK_N",
                   "UK_D",
                   "Italy_N",
                   "Italy_D",
                   "US_N",
                   "US_D",
                   "China_N",
                   "China_D",
                   "Japan_N",
                   "Japan_D")
    gsub(".*_","","Japan_N")
    
    disCol="#FFD70044"
    
    richColVec = c("Denmark_N"=eurCol,
                   "Denmark_D"=disCol,
                   "Sweden_N"=eurCol,
                   "Sweden_D"=disCol,
                   "Spain_N"=eurCol,
                   "Spain_D"=disCol,
                   "Germany_N"=eurCol,
                   "Germany_D"=disCol,
                   "Luxembourg_N"=eurCol,
                   "Luxembourg_D"=disCol,
                   "UK_N"=eurCol,
                   "UK_D"=disCol,
                   "Italy_N"=eurCol,
                   "Italy_D"=disCol,
                   "US_N"=cjpCol,
                   "US_D"=disCol,
                   "China_N"=cjpCol,
                   "China_D"=disCol,
                   "Japan_N"=cjpCol,
                   "Japan_D"=disCol)
    newCaseContRichnessByGeoTab$L1=factor(newCaseContRichnessByGeoTab$L1, levels = boxConSeqs)
    newCaseContRichnessByGeoTab$type=gsub(".*_","", newCaseContRichnessByGeoTab$L1)
    
    ggplot(newCaseContRichnessByGeoTab, aes(x=L1, y=value, fill=L1)) + geom_boxplot() + 
      scale_fill_manual(values = richColVec) +
      xlab("") + ylab("") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 60, hjust=1),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
    t.test(newCaseContRichnessByGeoList$Denmark_D, newCaseContRichnessByGeoList$Denmark_N)
    t.test(newCaseContRichnessByGeoList$Sweden_D, newCaseContRichnessByGeoList$Sweden_N)
    t.test(newCaseContRichnessByGeoList$Spain_D, newCaseContRichnessByGeoList$Spain_N)
    t.test(newCaseContRichnessByGeoList$China_D, newCaseContRichnessByGeoList$China_N)
    t.test(newCaseContRichnessByGeoList$US_D, newCaseContRichnessByGeoList$US_N)
    t.test(newCaseContRichnessByGeoList$Japan_D, newCaseContRichnessByGeoList$Japan_N)
    t.test(newCaseContRichnessByGeoList$UK_D, newCaseContRichnessByGeoList$UK_N)
    t.test(newCaseContRichnessByGeoList$Italy_D, newCaseContRichnessByGeoList$Italy_N)
    t.test(newCaseContRichnessByGeoList$Germany_D, newCaseContRichnessByGeoList$Germany_N)
    t.test(newCaseContRichnessByGeoList$Luxembourg_D, newCaseContRichnessByGeoList$Luxembourg_N)
    
    par(mfrow=c(2,1))
    par(mar=c(4,4,2,1))
    boxplot(rev(newCaseContRichnessByGeoList), las=2)
    plot.new()
    
  }
  
  oldBoxPlotMode2 = F
  if (oldBoxPlotMode2) {
    orderedCountries = c("Tanzania","Madagascar","Peru","India","Fiji","Mongo","Thai",
                         "Japan","China","US", "Finlan",
                         "Sweden","Denmark","UK","Spain","Germany","Luxembourg","Italy")
    
    europeCountries = c("Denmark", "Germany", "Italy", "Luxembourg", "Spain", "Sweden", "UK")
    tradCountries = c("Fiji", "India", "Madagascar", "Mongo", "Peru", "Tanzania", "Thai")
    cjpCountries = c("US","China","Japan")
    
    
    newNormalMetadata = basicMetaMapUpdated[match(newNormalAllIds_wo_FI, basicMetaMapUpdated$sample.ID),]
    newNormalRichnessByGeoList = split(newNormalMetadata$GeneRichness, newNormalMetadata$Geography )
    newNormalRichnessByGeoTab = melt(newNormalRichnessByGeoList)
    newNormalRichnessByGeoTab$L1[newNormalRichnessByGeoTab$L1 %in% cjpCountries ] = "UsChinaJapan"
    newNormalRichnessByGeoTab$L1[newNormalRichnessByGeoTab$L1 %in% tradCountries ] = "Traditional"
    newNormalRichnessByGeoTab$L1[newNormalRichnessByGeoTab$L1 %in% europeCountries ] = "European"
    newNormalRichnessByGeoTab$L1 = factor(newNormalRichnessByGeoTab$L1, levels = c("European","UsChinaJapan","Traditional"))
    ggplot(newNormalRichnessByGeoTab, aes(x=L1, y=value, fill=L1)) + geom_boxplot() + scale_fill_brewer("Set1") +
      theme(legend.position = "none")
    ss = split(newNormalRichnessByGeoTab$value, newNormalRichnessByGeoTab$L1)
    t.test(ss$European, ss$Traditional)
  }
  
}

richnessVolcanoModeUpdated = T
if (richnessVolcanoModeUpdated) {
  
  afterUpdate = T
  if (afterUpdate) {
    statsGenusForAllSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.genus.for.all.txt"
    statsGenusForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.genus.for.normal.txt"
    statsGenusForDiseaseSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.genus.for.disease.txt"
    statsGenusForIndustrialSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.genus.for.industrial.txt"
    statsGenusForTraditionalSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.genus.for.traditional.txt"
    
    statsMgsForAllSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.all.txt"
    statsMgsForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.normal.txt"
    statsMgsForDiseaseSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.disease.txt"
    statsMgsForIndustrialSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.industrial.txt"
    statsMgsForTraditionalSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.traditional.txt"
    
    lCut=0.25
    hCut=0.75
    
    
    richnessHist = T
    if (richnessHist) {
      plot(density(basicMetaMapUpdated$GeneRichness),  main="")
      abline(v=hAllCutNumber, col="red")
      abline(v=lAllCutNumber, col="red")
      
    }
    
    richness = basicMetaMapUpdated$GeneRichness
    names(richness) = basicMetaMapUpdated$sample.ID
    
    normalRichness = richness[names(richness) %in% newNormalAllIds]
    diseaseRichness = richness[names(richness) %in% newDiseaseIds]
    industrialRichness = richness[names(richness) %in% newWesternIds]
    
    hHalfCutNumber=as.numeric(quantile(basicMetaMapUpdated$GeneRichness,0.5))
    
    hAllCutNumber=as.numeric(quantile(basicMetaMapUpdated$GeneRichness,hCut))
    lAllCutNumber=as.numeric(quantile(basicMetaMapUpdated$GeneRichness,lCut))
    
    print(hAllCutNumber) #727622.5
    print(lAllCutNumber) #441994.5
    
    hgcSamples = names(richness[richness>=hAllCutNumber])
    lgcSamples = names(richness[richness<=lAllCutNumber])
    
    
    
    
    hgcHalfSamples = names(richness[richness>=hHalfCutNumber])
    lgcHalfSamples = names(richness[richness<=hHalfCutNumber])
    
    sort(table(basicMetaMapUpdated$Geography[basicMetaMapUpdated$sample.ID %in% hgcSamples] ))
    sort(table(basicMetaMapUpdated$Geography[basicMetaMapUpdated$sample.ID %in% lgcSamples] ))
    
    par(mar=c(5,10,4,1))
    barplot(sort(table(basicMetaMapUpdated$Geography[basicMetaMapUpdated$sample.ID %in% lgcSamples] )),
            horiz = T,
            las=1)
    barplot(sort(table(basicMetaMapUpdated$Geography[basicMetaMapUpdated$sample.ID %in% hgcSamples] )),
            horiz = T,
            las=1)
    
    caseContFraction = (cbind(hgc=c(sum(hgcSamples %in% newDiseaseIds), sum(hgcSamples %in% newNormalAllIds)),
                              lgc=c(sum(lgcSamples %in% newDiseaseIds), sum(lgcSamples %in% newNormalAllIds))))
    barplot(cbind(HGR=caseContFraction[,1]/1299, LGR=caseContFraction[,2]/1238))
    
    
    
    length(newDiseaseIds)
    
    hgcPrevMat=mergeMatUpdated[rownames(mergeMatUpdated) %in% taxo$MSP[taxo$genus=="Prevotella"],hgcHalfSamples]
    lgcPrevMat=mergeMatUpdated[rownames(mergeMatUpdated) %in% taxo$MSP[taxo$genus=="Prevotella"],lgcHalfSamples]
    hgcPrevMat=hgcPrevMat[,colnames(hgcPrevMat) %in% basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$enteroType=="ET-Prevotella"]]
    lgcPrevMat=lgcPrevMat[,colnames(lgcPrevMat) %in% basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$enteroType=="ET-Prevotella"]]
    
    hgcPrevMeans = rowMeans(hgcPrevMat)
    lgcPrevMeans = rowMeans(lgcPrevMat)
    data.frame(hgcPrevMeans, lgcPrevMeans)
    
    
    plot(hgcPrevMeans, lgcPrevMeans, log="xy")
    abline(a=0,b=1)
    length(hgcSamples)
    length(lgcSamples)
    
    normalHgcSamples = newNormalAllIds[newNormalAllIds %in% hgcSamples]
    normalLgcSamples = newNormalAllIds[newNormalAllIds %in% lgcSamples]
    
    normHgcSampleTab = data.frame(sample=normalHgcSamples, type="HGC.normal")
    normLgcSampleTab = data.frame(sample=normalLgcSamples, type="LGC.normal")
    normSampleTab = rbind(normHgcSampleTab, normLgcSampleTab)
    
    writeMode = F
    if (writeMode) write.table(normSampleTab, file="C://Data//comparative.analysis.healthy.sweden//reza.HGC.LGC.normal.samples.txt", sep="\t")
    
    checkInflowOutflowMode = T
    if (checkInflowOutflowMode) {
      volcano.stats.mgs.all = getVolcanoStat(mergeMatUpdated, hgcSamples, lgcSamples, taxo, T)
      volcano.stats.mgs.normal = getVolcanoStat(mergeMatUpdated, normalHgcSamples, normalLgcSamples, taxo, T)
      adjpCut=1e-3
      volcano.stats.mgs.normal.sel = volcano.stats.mgs.normal[volcano.stats.mgs.normal$qvalue <= adjpCut,]
      volcano.stats.mgs.normal.sel.hgc = volcano.stats.mgs.normal.sel[volcano.stats.mgs.normal.sel$lfc > 2,]
      volcano.stats.mgs.normal.sel.lgc = volcano.stats.mgs.normal.sel[volcano.stats.mgs.normal.sel$lfc < -2,]
      
      hgcSpecies = rownames(volcano.stats.mgs.normal.sel.hgc)
      lgcSpecies = rownames(volcano.stats.mgs.normal.sel.lgc)
      
      
      hgcLgcInflowList=list(HGR=inflow[hgcSpecies], LGR=inflow[lgcSpecies])
      hgcLgcOutflowList=list(HGR=outflow[hgcSpecies], LGR=outflow[lgcSpecies])
      hgcLgcInflowMelt = melt(hgcLgcInflowList)
      hgcLgcOutflowMelt = melt(hgcLgcOutflowList)
      
      t.test(hgcLgcInflowList$HGR, hgcLgcInflowList$LGR)
      t.test(hgcLgcOutflowList$HGR, hgcLgcOutflowList$LGR)
      
      wilcox.test(hgcLgcInflowList$HGR, hgcLgcInflowList$LGR, alternative = "less")
      wilcox.test(hgcLgcOutflowList$HGR, hgcLgcOutflowList$LGR, alternative = "less")
      
      
      ggplot(hgcLgcInflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
        scale_fill_manual(values = c("#ff000055","#2c528c55")) +
        xlab("") +
        ylab("") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 60, hjust=1),
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      ggplot(hgcLgcOutflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
        scale_fill_manual(values = c("#ff000055","#2c528c55")) +
        xlab("") +
        ylab("") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 60, hjust=1),
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
    }
    
    checkInflowOutflowFractionAbundance = T
    if (checkInflowOutflowFractionAbundance) {
      
      diseaseHgcSamples = newDiseaseIds[newDiseaseIds %in% hgcSamples]
      diseaseLgcSamples = newDiseaseIds[newDiseaseIds %in% lgcSamples]
      
      length(diseaseHgcSamples)
      length(diseaseLgcSamples)
      
      
      industrialHgcSamples = newWesternIds[newWesternIds %in% hgcSamples]
      industrialLgcSamples = newWesternIds[newWesternIds %in% lgcSamples]
      
      length(industrialHgcSamples)
      length(industrialLgcSamples)
      
      traditionalHgcSamples = newNonWesternIds[newNonWesternIds %in% hgcSamples]
      traditionalLgcSamples = newNonWesternIds[newNonWesternIds %in% lgcSamples]
      
      
      
      boxplot(list(hgc=richness[hgcSamples], lgc=richness[lgcSamples]))
      
      table(hgcSpecies %in% inflowSpecies)
      table(hgcSpecies %in% outflowSpecies)
      
      table(lgcSpecies %in% inflowSpecies)
      table(lgcSpecies %in% outflowSpecies)
      
      mergeMatUpdatedDigit = mergeMatUpdated != 0
      
      
      getSigFraction <- function(mgsMatDigit, targetSamples, targetSpecies,pCut=0.05) {
        targetMat = mgsMatDigit[,targetSamples]
        rowSpecies = rownames(targetMat)
        result = apply(targetMat, 2, function(currCol){
          currSpecies = rowSpecies[currCol != 0]
          overlapSpecies = currSpecies[currSpecies %in% targetSpecies]
          totalSpecies = unique(c(currSpecies, targetSpecies))
          hyperP = getHyperP(targetSpecies, currSpecies, rowSpecies,F)
          #jaccard = length(overlapSpecies)/length(currSpecies)
          return((hyperP))
        })
        
        return(mean(result < pCut))
      }
      
      inflow.outflow.sig.fraction.list = c(inflow.hgc.health = getSigFraction(mergeMatUpdatedDigit, industrialHgcSamples, inflowSpecies),
                                           outflow.hgc.health = getSigFraction(mergeMatUpdatedDigit, industrialHgcSamples, outflowSpecies),
                                           inflow.lgc.health = getSigFraction(mergeMatUpdatedDigit, industrialLgcSamples, inflowSpecies),
                                           outflow.lgc.health = getSigFraction(mergeMatUpdatedDigit, industrialLgcSamples, outflowSpecies),
                                           inflow.hgc.disease = getSigFraction(mergeMatUpdatedDigit, diseaseHgcSamples, inflowSpecies),
                                           outflow.hgc.disease = getSigFraction(mergeMatUpdatedDigit, diseaseHgcSamples, outflowSpecies),
                                           inflow.lgc.disease = getSigFraction(mergeMatUpdatedDigit, diseaseLgcSamples, inflowSpecies),
                                           outflow.lgc.disease = getSigFraction(mergeMatUpdatedDigit, diseaseLgcSamples, outflowSpecies))
      
      save(melted, file = "C:\\Data\\comparative.analysis.healthy.sweden\\melted.inflow.outflow.sig.fraction.0.05.RData")
      
      melted = melt(inflow.outflow.sig.fraction.list)
      melted$health = c(rep("health",4),rep("disease",4))
      melted$Richness = c(rep(c("HGC","HGC","LGC","LGC"),2))
      melted$in.outflow = c(rep(c("inflow","outflow"),4))
      melted$group = paste(melted$health, melted$Richness)
      melted$group = factor(melted$group)
      ggplot(melted, aes(x=group, y=value, fill=in.outflow)) + geom_bar(position="stack", stat="identity")
      
      par(mar=c(10,4,4,1))
      barplot(inflow.outflow.sig.fraction.list, ylab="fractions per inflow/outflow species", las=2)
      
      par(mar=c(7,4,4,1))
      barplot(inflow.outflow.sig.fraction.list[c(1,3)], ylab="fractions per inflow/outflow species", las=2)
      
      par(mar=c(7,4,4,1))
      barplot(inflow.outflow.sig.fraction.list[c(2,4)], ylab="fractions per inflow/outflow species", las=2)
      
      inflow.outflow.presence.list = list(inflow.hgc.health = colSums(mergeMatUpdatedDigit[inflowSpecies, industrialHgcSamples])/length(inflowSpecies),
                                          outflow.hgc.health = colSums(mergeMatUpdatedDigit[outflowSpecies, industrialHgcSamples])/length(outflowSpecies),
                                          inflow.lgc.health = colSums(mergeMatUpdatedDigit[inflowSpecies, industrialLgcSamples])/length(inflowSpecies),
                                          outflow.lgc.health = colSums(mergeMatUpdatedDigit[outflowSpecies, industrialLgcSamples])/length(outflowSpecies),
                                          inflow.hgc.disease = colSums(mergeMatUpdatedDigit[inflowSpecies, diseaseHgcSamples])/length(inflowSpecies),
                                          outflow.hgc.disease = colSums(mergeMatUpdatedDigit[outflowSpecies, diseaseHgcSamples])/length(outflowSpecies),
                                          inflow.lgc.disease = colSums(mergeMatUpdatedDigit[inflowSpecies, diseaseLgcSamples])/length(inflowSpecies),
                                          outflow.lgc.disease = colSums(mergeMatUpdatedDigit[outflowSpecies, diseaseLgcSamples])/length(outflowSpecies))
      
      par(mar=c(7,4,4,1))
      boxplot(inflow.outflow.presence.list, ylab="fractions per inflow/outflow species", las=2)
      
      mergeMatUpdatedRelAbd =  getRelAbd(mergeMatUpdated)
      
      inflow.outflow.abd.list = list(inflow.hgc.health = colSums(mergeMatUpdatedRelAbd[inflowSpecies, industrialHgcSamples]),
                                     outflow.hgc.health = colSums(mergeMatUpdatedRelAbd[outflowSpecies, industrialHgcSamples]),
                                     inflow.lgc.health = colSums(mergeMatUpdatedRelAbd[inflowSpecies, industrialLgcSamples]),
                                     outflow.lgc.health = colSums(mergeMatUpdatedRelAbd[outflowSpecies, industrialLgcSamples]),
                                     inflow.hgc.disease = colSums(mergeMatUpdatedRelAbd[inflowSpecies, diseaseHgcSamples]),
                                     outflow.hgc.disease = colSums(mergeMatUpdatedRelAbd[outflowSpecies, diseaseHgcSamples]),
                                     inflow.lgc.disease = colSums(mergeMatUpdatedRelAbd[inflowSpecies, diseaseLgcSamples]),
                                     outflow.lgc.disease = colSums(mergeMatUpdatedRelAbd[outflowSpecies, diseaseLgcSamples]))
      
      #inflow.outflow.abd.list = lapply(inflow.outflow.abd.list, function(x) log10(x + 1e-10))
      
      
      par(mar=c(10,4,4,1))
      boxplot(inflow.outflow.abd.list, ylab="log10 (cumulative sum abundance)", las=2)
      
      boxplot(inflow.outflow.abd.list[c(1,3,5,7)], ylab="log10 (cumulative sum abundance)", las=2)
      boxplot(inflow.outflow.abd.list[c(2,4,6,8)], ylab="log10 (cumulative sum abundance)", las=2)
      
      wilcox.test(inflow.outflow.abd.list$outflow.hgc.health, inflow.outflow.abd.list$outflow.lgc.health, alternative = "greater")
      wilcox.test(inflow.outflow.abd.list$outflow.hgc.health, inflow.outflow.abd.list$outflow.hgc.disease, alternative = "less")
      wilcox.test(inflow.outflow.abd.list$outflow.hgc.health, inflow.outflow.abd.list$outflow.lgc.disease, alternative = "less")
      wilcox.test(inflow.outflow.abd.list$outflow.hgc.disease, inflow.outflow.abd.list$outflow.lgc.disease, alternative = "less")
      
      wilcox.test(inflow.outflow.abd.list$inflow.hgc, inflow.outflow.abd.list$outflow.hgc)
      wilcox.test(inflow.outflow.abd.list$inflow.hgc, inflow.outflow.abd.list$inflow.lgc)
      wilcox.test(inflow.outflow.abd.list$inflow.lgc, inflow.outflow.abd.list$outflow.hgc)
      
      
      boxplot(inflow.outflow.abd.list[c(1,3)])
      boxplot(inflow.outflow.abd.list[c(2,4)])
      
      
    }
    
    length(normalHgcSamples)
    length(normalLgcSamples)
    
    diseaseHgcSamples = newDiseaseIds[newDiseaseIds %in% hgcSamples]
    diseaseLgcSamples = newDiseaseIds[newDiseaseIds %in% lgcSamples]
    
    length(diseaseHgcSamples)
    length(diseaseLgcSamples)
    
    
    industrialHgcSamples = newWesternIds[newWesternIds %in% hgcSamples]
    industrialLgcSamples = newWesternIds[newWesternIds %in% lgcSamples]
    
    length(industrialHgcSamples)
    length(industrialLgcSamples)
    
    traditionalHgcSamples = newNonWesternIds[newNonWesternIds %in% hgcSamples]
    traditionalLgcSamples = newNonWesternIds[newNonWesternIds %in% lgcSamples]
    
    length(traditionalHgcSamples)
    length(traditionalLgcSamples)
    
    
    normHgcSampleTab = data.frame(sample=normalHgcSamples, type="HGC.normal")
    normLgcSampleTab = data.frame(sample=normalLgcSamples, type="LGC.normal")
    diseaseHgcSampleTab = data.frame(sample=diseaseHgcSamples, type="HGC.disease")
    diseaseLgcSampleTab = data.frame(sample=diseaseLgcSamples, type="LGC.disease")
    industrialHgcSampleTab = data.frame(sample=industrialHgcSamples, type="HGC.industrial")
    industrialLgcSampleTab = data.frame(sample=industrialLgcSamples, type="LGC.industrial")
    traditionalHgcSampleTab = data.frame(sample=traditionalHgcSamples, type="HGC.traditional")
    traditionalLgcSampleTab = data.frame(sample=traditionalLgcSamples, type="LGC.traditional")
    allHgcSampleTab = data.frame(sample=traditionalHgcSamples, type="HGC.all")
    allLgcSampleTab = data.frame(sample=traditionalLgcSamples, type="LGC.all")
    
    
    normSampleTab = rbind(normHgcSampleTab, 
                          normLgcSampleTab,
                          diseaseHgcSampleTab,
                          diseaseLgcSampleTab,
                          industrialHgcSampleTab,
                          industrialLgcSampleTab,
                          traditionalHgcSampleTab,
                          traditionalLgcSampleTab,
                          allHgcSampleTab,
                          allLgcSampleTab)
    writeMode = F
    if (writeMode) write.table(normSampleTab, file="C://Data//comparative.analysis.healthy.sweden//reza.HGC.LGC.samples.all.types.txt", sep="\t")
    
    #newNormalAllIds
    #newDiseaseIds
    #newWesternIds = newNormalAllIds[!newNormalAllIds %in% newNonWesternIds]
    
    
    #and then volcano plots
    
    
    genusMode = F
    if (genusMode) {
      
      ### all modes ####
      
      
      volcano.stats.genus.all = getVolcanoStat(genusMatUpdated, hgcSamples, lgcSamples, taxo)
      volcano.stats.genus.normal = getVolcanoStat(genusMatUpdated, normalHgcSamples, normalLgcSamples, taxo)
      volcano.stats.genus.disease = getVolcanoStat(genusMatUpdated, diseaseHgcSamples, diseaseLgcSamples, taxo)
      volcano.stats.genus.industrial = getVolcanoStat(genusMatUpdated, industrialHgcSamples, industrialLgcSamples, taxo)
      volcano.stats.genus.traditional = getVolcanoStat(genusMatUpdated, traditionalHgcSamples, traditionalLgcSamples, taxo)
      
      write.table(volcano.stats.genus.all, statsGenusForAllSamples, sep="\t")
      write.table(volcano.stats.genus.normal, statsGenusForNormSamples, sep="\t")
      write.table(volcano.stats.genus.disease, statsGenusForDiseaseSamples, sep="\t")
      write.table(volcano.stats.genus.industrial, statsGenusForIndustrialSamples, sep="\t")
      write.table(volcano.stats.genus.traditional, statsGenusForTraditionalSamples, sep="\t")
      
      plotVolcano(volcano.stats.genus.all, F)
      
      pairs(cbind(all=volcano.stats.genus.all$lfc,
                  normal=volcano.stats.genus.normal$lfc,
                  disease=volcano.stats.genus.disease$lfc,
                  industrial=volcano.stats.genus.industrial$lfc))
      
      
    }
    
    mgsMode = T
    if (mgsMode) {
      
      generateMode = F
      if (generateMode) {
        volcano.stats.mgs.all = getVolcanoStat(mergeMatUpdated, hgcSamples, lgcSamples, taxo, T)
        volcano.stats.mgs.normal = getVolcanoStat(mergeMatUpdated, normalHgcSamples, normalLgcSamples, taxo, T)
        volcano.stats.mgs.disease = getVolcanoStat(mergeMatUpdated, diseaseHgcSamples, diseaseLgcSamples, taxo, T)
        volcano.stats.mgs.industrial = getVolcanoStat(mergeMatUpdated, industrialHgcSamples, industrialLgcSamples, taxo, T)
        volcano.stats.mgs.traditional = getVolcanoStat(mergeMatUpdated, traditionalHgcSamples, traditionalLgcSamples, taxo, T)
        
        write.table(volcano.stats.mgs.all, statsMgsForAllSamples, sep="\t")
        write.table(volcano.stats.mgs.normal, statsMgsForNormSamples, sep="\t")
        write.table(volcano.stats.mgs.disease, statsMgsForDiseaseSamples, sep="\t")
        write.table(volcano.stats.mgs.industrial, statsMgsForIndustrialSamples, sep="\t")
        write.table(volcano.stats.mgs.traditional, statsMgsForTraditionalSamples, sep="\t")
      }
      
      
      loadMode = T
      if (loadMode) {
        volcano.stats.mgs.all = read.table(statsMgsForAllSamples, sep="\t")
        volcano.stats.mgs.normal = read.table(statsMgsForNormSamples, sep="\t")
        volcano.stats.mgs.disease = read.table(statsMgsForDiseaseSamples, sep="\t")
        
      }
      
      
      
      plotVolcano(volcano.stats.mgs.all, F)
      
      pairs(cbind(all=volcano.stats.mgs.all$lfc,
                  normal=volcano.stats.mgs.normal$lfc,
                  disease=volcano.stats.mgs.disease$lfc,
                  industrial=volcano.stats.mgs.industrial$lfc))
      
      ### figure for scatter plot ### 
      
      normDiseaseScatterTab = data.frame(normal=volcano.stats.mgs.normal$lfc, 
                                         disease=volcano.stats.mgs.disease$lfc,
                                         stringsAsFactors = F)
      rownames(normDiseaseScatterTab)=rownames(volcano.stats.mgs.normal)
      normDiseaseScatterTab = normDiseaseScatterTab[bothSpecies,]
      
      normDiseaseScatterTab[normDiseaseScatterTab >= 15]=NA
      normDiseaseScatterTab[normDiseaseScatterTab <= -15]=NA
      
      ggplot(normDiseaseScatterTab, aes(x=disease, y=normal, fill=normal)) + 
        geom_point(size=2, shape=21, colour="black")  +
        xlab("") + ylab("")+
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
        geom_smooth(method='lm',formula=y~x, se=FALSE, color="red")+
        scale_fill_gradient2(midpoint = 0.2, low="#FFD70044", mid = "white", high = "purple") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      normWellnessScatterTab = data.frame(normal=volcano.stats.mgs.normal$lfc, 
                                          wellness=- wellnessMgsWilcoxStatsDecreasedPair$l2fc,
                                          stringsAsFactors = F)
      rownames(normWellnessScatterTab)=rownames(volcano.stats.mgs.normal)
      normWellnessScatterTab = normWellnessScatterTab[bothSpecies,]
      
      normWellnessScatterTab[normWellnessScatterTab >= 10]=NA
      normWellnessScatterTab[normWellnessScatterTab <= -10]=NA
      normWellnessScatterTab[normWellnessScatterTab == 0]=NA
      
      
      ggplot(normWellnessScatterTab, aes(x=wellness, y=normal, fill=normal)) + 
        geom_point(size=2, shape=21, colour="black")  +
        xlab("") + ylab("")+ xlim(c(-3.5,6))+
        geom_smooth(method='lm',formula=y~x, se=FALSE, color="red")+
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
        scale_fill_gradient2(midpoint = 0.2, low="#FFD70044", mid = "white", high = "purple") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      normAntibioticsScatterTab = data.frame(normal=volcano.stats.mgs.normal$lfc, 
                                             antibiotics=- wilcoxStatsDecrease4$l2fc,
                                             stringsAsFactors = F)
      
      
      rownames(normAntibioticsScatterTab)=rownames(volcano.stats.mgs.normal)
      normAntibioticsScatterTab = normAntibioticsScatterTab[bothSpecies,]
      
      normAntibioticsScatterTab[normAntibioticsScatterTab >= 10]=NA
      normAntibioticsScatterTab[normAntibioticsScatterTab <= -10]=NA
      normAntibioticsScatterTab[normAntibioticsScatterTab == 0]=NA
      
      
      ggplot(normAntibioticsScatterTab, aes(x=antibiotics, y=normal, fill=normal)) + 
        geom_point(size=2, shape=21, colour="black")  +
        xlab("") + ylab("")+ 
        geom_smooth(method='lm',formula=y~x, se=FALSE, color="red")+
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
        scale_fill_gradient2(midpoint = 0.2, low="#FFD70044", mid = "white", high = "purple") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      normAntibioticsShortScatterTab = data.frame(normal=volcano.stats.mgs.normal$lfc, 
                                                  antibiotics=- wilcoxStatsDecrease1$l2fc,
                                                  stringsAsFactors = F)
      
      rownames(normAntibioticsShortScatterTab)=rownames(volcano.stats.mgs.normal)
      normAntibioticsShortScatterTab = normAntibioticsShortScatterTab[bothSpecies,]
      
      normAntibioticsShortScatterTab[normAntibioticsShortScatterTab >= 10]=NA
      normAntibioticsShortScatterTab[normAntibioticsShortScatterTab <= -10]=NA
      normAntibioticsShortScatterTab[normAntibioticsShortScatterTab == 0]=NA
      
      
      ggplot(normAntibioticsShortScatterTab, aes(x=antibiotics, y=normal, fill=normal)) + 
        geom_point(size=2, shape=21, colour="black")  +
        xlab("") + ylab("")+ 
        geom_smooth(method='lm',formula=y~x, se=FALSE, color="red")+
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
        scale_fill_gradient2(midpoint = 0.2, low="#FFD70044", mid = "white", high = "purple") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      cor(normDiseaseScatterTab$normal, 
          normDiseaseScatterTab$disease, 
          method = "spearman", use="pairwise.complete.obs")
      
      cor(normWellnessScatterTab$normal, 
          normWellnessScatterTab$wellness, 
          method = "spearman", use="pairwise.complete.obs")
      
      cor(normAntibioticsScatterTab$normal, 
          normAntibioticsScatterTab$antibiotics, 
          method = "spearman", use="pairwise.complete.obs")
      
      cor(normAntibioticsShortScatterTab$normal, 
          normAntibioticsShortScatterTab$antibiotics, 
          method = "spearman", use="pairwise.complete.obs")
      
    }
    
    
  }
  
  woWellness = T
  if (woWellness) {
    
    lCut=0.25
    hCut=0.75
    
    basicMetaMapUpdatedWoWellness = basicMetaMapUpdated[!basicMetaMapUpdated$sample.ID %in% wellnessSamples, ]
    
    richness = basicMetaMapUpdatedWoWellness$GeneRichness
    names(richness) = basicMetaMapUpdatedWoWellness$sample.ID
    
    normalRichness = richness[names(richness) %in% newNormalAllIds]
    diseaseRichness = richness[names(richness) %in% newDiseaseIds]
    industrialRichness = richness[names(richness) %in% newWesternIds]
    
    hHalfCutNumber=as.numeric(quantile(basicMetaMapUpdatedWoWellness$GeneRichness,0.5))
    
    hAllCutNumber=as.numeric(quantile(basicMetaMapUpdatedWoWellness$GeneRichness,hCut))
    lAllCutNumber=as.numeric(quantile(basicMetaMapUpdatedWoWellness$GeneRichness,lCut))
    
    print(hAllCutNumber) #708843.2
    print(lAllCutNumber) #431824.5
    
    hgcSamples = names(richness[richness>=hAllCutNumber])
    lgcSamples = names(richness[richness<=lAllCutNumber])
    
    normalHgcSamples = newNormalAllIds[newNormalAllIds %in% hgcSamples]
    normalLgcSamples = newNormalAllIds[newNormalAllIds %in% lgcSamples]
    
    
    table(normalHgcSamples %in% wellnessSamples)
    table(normalLgcSamples %in% wellnessSamples)
    
    
    checkInflowOutflowMode = T
    if (checkInflowOutflowMode) {
      volcano.stats.mgs.all = getVolcanoStat(mergeMatUpdated, hgcSamples, lgcSamples, taxo, T)
      volcano.stats.mgs.normal = getVolcanoStat(mergeMatUpdated, normalHgcSamples, normalLgcSamples, taxo, T)
      adjpCut=1e-3
      volcano.stats.mgs.normal.sel = volcano.stats.mgs.normal[volcano.stats.mgs.normal$qvalue <= adjpCut,]
      volcano.stats.mgs.normal.sel.hgc = volcano.stats.mgs.normal.sel[volcano.stats.mgs.normal.sel$lfc > 2,]
      volcano.stats.mgs.normal.sel.lgc = volcano.stats.mgs.normal.sel[volcano.stats.mgs.normal.sel$lfc < -2,]
      
      hgcSpecies = rownames(volcano.stats.mgs.normal.sel.hgc)
      lgcSpecies = rownames(volcano.stats.mgs.normal.sel.lgc)
      
      
      hgcLgcInflowList=list(HGC=inflow[hgcSpecies], LGC=inflow[lgcSpecies])
      hgcLgcOutflowList=list(HGC=outflow[hgcSpecies], LGC=outflow[lgcSpecies])
      hgcLgcInflowMelt = melt(hgcLgcInflowList)
      hgcLgcOutflowMelt = melt(hgcLgcOutflowList)
      
      t.test(hgcLgcInflowList$HGC, hgcLgcInflowList$LGC)
      t.test(hgcLgcOutflowList$HGC, hgcLgcOutflowList$LGC)
      
      CPCOLS <- c("#493D8C91", "#FF735786")
      
      
      ggplot(hgcLgcInflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
        scale_fill_manual(values = c("#ff000055","#2c528c55")) +
        xlab("") +
        ylab("Inflow") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 60, hjust=1),
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
      
      ggplot(hgcLgcOutflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
        scale_fill_manual(values = c("#ff000055","#2c528c55")) +
        xlab("") +
        ylab("Outflow") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 60, hjust=1),
              panel.background = element_rect(fill = "white", 
                                              colour = "black", 
                                              size = 0.5, 
                                              linetype = "solid"))
    }
    
    
    
  }
  
  beforeUpdate = F
  if (beforeUpdate) {
    normalMat = mergeMat[,colnames(mergeMat) %in% normalIDs]
    diseaseMat = mergeMat[,colnames(mergeMat) %in% caseIDs2]
    diseaseFamilyMat = getTaxaSumMat(diseaseMat,taxo,"family",T)
    
    indiMat = mergeMat[,colnames(mergeMat) %in% indiIDs]
    indiFamilyMat = getTaxaSumMat(indiMat,taxo,"family",T)
    
    targetMat = mergeMat[,colnames(mergeMat) %in% normalIDs]
    phylaMat = getTaxaSumMat(targetMat,taxo,"phylum",T)
    familyMat = getTaxaSumMat(targetMat,taxo,"family",T)
    
    #healthy controls by geometry
    notIndiRichness = basicMetaMap$GeneRichness[match(normalIDs, basicMetaMap$sample.ID)]
    names(notIndiRichness) = basicMetaMap$sample.ID[match(normalIDs, basicMetaMap$sample.ID)]
    
    notIndiAge = basicMetaMap$Age[match(normalIDs, basicMetaMap$sample.ID)]
    notIndiBMI = basicMetaMap$BMI[match(normalIDs, basicMetaMap$sample.ID)]
    notIndiGender = basicMetaMap$Gender[match(normalIDs, basicMetaMap$sample.ID)]
    notIndiGeos = basicMetaMap$Geography[match(normalIDs, basicMetaMap$sample.ID)]
    
    notIndiRichnessByGeos = split(notIndiRichness,notIndiGeos)
    notIndiAgeByGeos = split(notIndiAge,notIndiGeos)
    notIndiBmiByGeos = split(notIndiBMI,notIndiGeos)
    notIndiGenderByGeos = split(notIndiGender,notIndiGeos)
    
    notIndiRichGeoMedian = do.call(c,lapply(notIndiRichnessByGeos,median))
    notIndiRichOrdered = names(sort(notIndiRichGeoMedian))
    
    diseaseRichness = basicMetaMap$GeneRichness[match(caseIDs2, basicMetaMap$sample.ID)]
    names(diseaseRichness) = basicMetaMap$sample.ID[match(caseIDs2, basicMetaMap$sample.ID)]
    
    diseaseAge = basicMetaMap$Age[match(caseIDs2, basicMetaMap$sample.ID)]
    diseaseBMI = basicMetaMap$BMI[match(caseIDs2, basicMetaMap$sample.ID)]
    diseaseGender = basicMetaMap$Gender[match(caseIDs2, basicMetaMap$sample.ID)]
    diseaseGeos = basicMetaMap$Geography[match(caseIDs2, basicMetaMap$sample.ID)]
    diseaseDsIds = basicMetaMap$dataset.ID[match(caseIDs2, basicMetaMap$sample.ID)]
    diseaseDsNames = caseName2Tab$Case[match(diseaseDsIds, caseName2Tab$ID)]
    
    diseaseRichnessByDs = split(diseaseRichness,diseaseDsNames)
    diseaseAgeByDs = split(diseaseAge,diseaseDsNames)
    diseaseBmiByDs = split(diseaseBMI,diseaseDsNames)
    diseaseGenderByDs = split(diseaseGender,diseaseDsNames)
    
    diseaseRichDsMedian = do.call(c,lapply(diseaseRichnessByDs,median))
    diseaseRichOrdered = names(sort(diseaseRichDsMedian))
    
    indiRichness = basicMetaMap$GeneRichness[match(indiIDs, basicMetaMap$sample.ID)]
    names(indiRichness) = basicMetaMap$sample.ID[match(indiIDs, basicMetaMap$sample.ID)]
    
    indiAge = basicMetaMap$Age[match(indiIDs, basicMetaMap$sample.ID)]
    indiBMI = basicMetaMap$BMI[match(indiIDs, basicMetaMap$sample.ID)]
    indiGender = basicMetaMap$Gender[match(indiIDs, basicMetaMap$sample.ID)]
    indiGeos = basicMetaMap$Geography[match(indiIDs, basicMetaMap$sample.ID)]
    
    indiRichnessByGeos = split(indiRichness,indiGeos)
    indiAgeByGeos = split(indiAge,indiGeos)
    indiBmiByGeos = split(indiBMI,indiGeos)
    indiGenderByGeos = split(indiGender,indiGeos)
    
    indiRichGeoMedian = do.call(c,lapply(indiRichnessByGeos,median))
    indiRichOrdered = names(sort(indiRichGeoMedian))
    
    lCut=0.25
    hCut=0.75
    
    recapMode = T
    if (recapMode) {
      load("healthy.HGC.LGC.RData")
      load("diseased.HGC.LGC.RData")
      load("indi.HGC.LGC.RData")
      
      # statSpeciesTab
      # statDiseaseSpeciesTab
      # statIndiSpeciesTab
      # 
      statSpeciesTab$lfc[statSpeciesTab$lfc==15]=16
      statSpeciesTab$lfc[statSpeciesTab$lfc==-15]=-16
      
      statIndiSpeciesTab$lfc[statIndiSpeciesTab$lfc==15]=16
      statIndiSpeciesTab$lfc[statIndiSpeciesTab$lfc==-15]=-16
      
      lfcMat = cbind(healthy=statSpeciesTab$lfc,
                     diseased=statDiseaseSpeciesTab$lfc,
                     indigenous=statIndiSpeciesTab$lfc)
      pairs(lfcMat)
      
      
      cor(statSpeciesTab$lfc, statDiseaseSpeciesTab$lfc)
      cor(statSpeciesTab$lfc, statIndiSpeciesTab$lfc)
      
      plot(statSpeciesTab$lfc, statDiseaseSpeciesTab$lfc, 
           pch=16, col="#80808088",las=1, tck = 0.02,yaxt="n",xaxt="n",
           xlab="LFC, healthy", ylab="LFC, diseased", xlim=c(-13,13),ylim=c(-13,13))
      abline(a=0,b=1,col="red")
      grid(lty=3,col="#80808077")
      
      plot(statSpeciesTab$lfc, statIndiSpeciesTab$lfc, 
           pch=16, col="#80808088",las=1,tck = 0.02,yaxt="n",xaxt="n",
           xlab="LFC, healthy", ylab="LFC, indigenous", xlim=c(-13,13),ylim=c(-13,13))
      abline(a=0,b=1,col="red")
      grid(lty=3,col="#80808077")
      
    }
    
    healthyMode = T
    if (healthyMode) {
      
      
      healthyHcutNumber=as.numeric(quantile(notIndiRichness,hCut))
      healthyLcutNumber=as.numeric(quantile(notIndiRichness,lCut))
      
      speciesMode = T
      if (speciesMode) {
        specMean = rowMeans(normalMat)
        sPart = specMean/sum(specMean)
        
        hgcSamples = notIndiRichness[notIndiRichness > as.numeric(quantile(notIndiRichness,hCut))]
        lgcSamples = notIndiRichness[notIndiRichness < as.numeric(quantile(notIndiRichness,lCut))]
        
        hgcMat = normalMat[,names(hgcSamples)]
        lgcMat = normalMat[,names(lgcSamples)]
        
        lenSpecies = dim(hgcMat)[1]
        pValues = sapply(1:lenSpecies, function(ind) {
          mergedVec = c(hgcMat[ind,],lgcMat[ind,])
          mergedFac = c(rep("HGC",length(hgcMat[ind,])),
                        rep("LGC",length(lgcMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenSpecies,function(ind){
          hmeans = mean(hgcMat[ind,])
          lmeans = mean(lgcMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(normalMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-15
        lfcs[is.infinite(lfcs)]=2^15
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statSpeciesTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=sPart, label=getSpeciesName(names(lfcs),taxo), stringsAsFactors = F)
        
        lfcInd = abs(statSpeciesTab$lfc)>2
        sigInd = statSpeciesTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statSpeciesTabModLabel  = statSpeciesTab
        statSpeciesTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        #statTabModLabel$label[!sigInd | !lfcInd] = ""
        
        plotVolcano(statSpeciesTabModLabel, F)
        
        
      }
      
      familyMode = F
      if (familyMode) {
        
        familyMean = rowMeans(familyMat)
        fPart = familyMean/sum(familyMean)
        
        as.numeric(quantile(notIndiRichness,hCut))
        as.numeric(quantile(notIndiRichness,lCut))
        
        hgcSamples = notIndiRichness[notIndiRichness > as.numeric(quantile(notIndiRichness,hCut))]
        lgcSamples = notIndiRichness[notIndiRichness < as.numeric(quantile(notIndiRichness,lCut))]
        
        hgcFamilyMat = familyMat[,names(hgcSamples)]
        lgcFamilyMat = familyMat[,names(lgcSamples)]
        
        lenFamilies = dim(hgcFamilyMat)[1]
        pValues = sapply(1:lenFamilies, function(ind) {
          mergedVec = c(hgcFamilyMat[ind,],lgcFamilyMat[ind,])
          mergedFac = c(rep("HGC",length(hgcFamilyMat[ind,])),
                        rep("LGC",length(lgcFamilyMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenFamilies,function(ind){
          hmeans = mean(hgcFamilyMat[ind,])
          lmeans = mean(lgcFamilyMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(familyMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-15
        lfcs[is.infinite(lfcs)]=2^15
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=fPart, label=names(lfcs), stringsAsFactors = F)
        
        lfcInd = abs(statTab$lfc)>2
        sigInd = statTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statTabModLabel  = statTab
        statTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        #statTabModLabel$label[!sigInd | !lfcInd] = ""
        
        
        
        
        plotVolcano(statTabModLabel, F)
        
      }
      
      save(statTab, statSpeciesTab, file="healthy.HGC.LGC.RData")
      
      
    }
    
    diseaseMode = T
    if (diseaseMode) {
      diseaseHcutNumber=as.numeric(quantile(diseaseRichness,hCut))
      diseaseLcutNumber=as.numeric(quantile(diseaseRichness,lCut))
      
      targetHCut = healthyHcutNumber
      targetLCut = healthyLcutNumber
      
      hgcDiseaseSamples = diseaseRichness[diseaseRichness > targetHCut]
      lgcDiseaseSamples = diseaseRichness[diseaseRichness < targetLCut]
      
      
      speciesMode = T
      if (speciesMode){
        specDiseaseMean = rowMeans(diseaseMat)
        sDiseasePart = specDiseaseMean/sum(specDiseaseMean)
        
        ################
        hgcDiseaseMat = diseaseMat[,names(hgcDiseaseSamples)]
        lgcDiseaseMat = diseaseMat[,names(lgcDiseaseSamples)]
        
        
        lenSpecies = dim(hgcDiseaseMat)[1]
        pValues = sapply(1:lenSpecies, function(ind) {
          mergedVec = c(hgcDiseaseMat[ind,],lgcDiseaseMat[ind,])
          mergedFac = c(rep("HGC",length(hgcDiseaseMat[ind,])),
                        rep("LGC",length(lgcDiseaseMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenSpecies,function(ind){
          hmeans = mean(hgcDiseaseMat[ind,])
          lmeans = mean(lgcDiseaseMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(diseaseMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-16
        lfcs[is.infinite(lfcs)]=2^16
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statDiseaseSpeciesTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=sPart, label=getSpeciesName(names(lfcs),taxo), stringsAsFactors = F)
        
        lfcInd = abs(statDiseaseSpeciesTab$lfc)>2
        sigInd = statDiseaseSpeciesTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statDiseaseSpeciesTabModLabel  = statDiseaseSpeciesTab
        statDiseaseSpeciesTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        #statTabModLabel$label[!sigInd | !lfcInd] = ""
        
        
        
        
        plotVolcano(statDiseaseSpeciesTabModLabel, F)
        
        ###############
      }
      
      familyMode = T
      if (familyMode) {
        diseaseFamilyMean = rowMeans(diseaseFamilyMat)
        fDiseasePart = diseaseFamilyMean/sum(diseaseFamilyMean)
        
        hgcDiseaseFamilyMat = diseaseFamilyMat[,names(hgcDiseaseSamples)]
        lgcDiseaseFamilyMat = diseaseFamilyMat[,names(lgcDiseaseSamples)]
        
        lenFamilies = dim(diseaseFamilyMat)[1]
        
        pValues = sapply(1:lenFamilies, function(ind) {
          mergedVec = c(hgcDiseaseFamilyMat[ind,],lgcDiseaseFamilyMat[ind,])
          mergedFac = c(rep("HGC",length(hgcDiseaseFamilyMat[ind,])),
                        rep("LGC",length(lgcDiseaseFamilyMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenFamilies,function(ind){
          hmeans = mean(hgcDiseaseFamilyMat[ind,])
          lmeans = mean(lgcDiseaseFamilyMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(diseaseFamilyMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-15
        lfcs[is.infinite(lfcs)]=2^15
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statDiseaseTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=fPart, label=names(lfcs), stringsAsFactors = F)
        
        lfcInd = abs(statDiseaseTab$lfc)>2
        sigInd = statDiseaseTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statDiseaseTabModLabel  = statDiseaseTab
        statDiseaseTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        #statTabModLabel$label[!sigInd | !lfcInd] = ""
        
        
        plotVolcano(statDiseaseTabModLabel, F)
        
      }
      
      save(statDiseaseSpeciesTab, statDiseaseSpeciesTab, file="diseased.HGC.LGC.RData")
      
    }
    
    indiMode = T
    if (indiMode) {
      indiHcutNumber=as.numeric(quantile(indiRichness,hCut))
      indiLcutNumber=as.numeric(quantile(indiRichness,lCut))
      
      targetHCut = healthyHcutNumber
      targetLCut = healthyLcutNumber
      
      hgcIndiSamples = indiRichness[indiRichness > targetHCut]
      lgcIndiSamples = indiRichness[indiRichness < targetLCut]
      
      speciesMode = T
      if (speciesMode) {
        specIndiMean = rowMeans(indiMat)
        sIndiPart = specIndiMean/sum(specIndiMean)
        
        ################
        hgcIndiMat = indiMat[,names(hgcIndiSamples)]
        lgcIndiMat = indiMat[,names(lgcIndiSamples)]
        
        
        lenSpecies = dim(hgcIndiMat)[1]
        pValues = sapply(1:lenSpecies, function(ind) {
          mergedVec = c(hgcIndiMat[ind,],lgcIndiMat[ind,])
          mergedFac = c(rep("HGC",length(hgcIndiMat[ind,])),
                        rep("LGC",length(lgcIndiMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenSpecies,function(ind){
          hmeans = mean(hgcIndiMat[ind,])
          lmeans = mean(lgcIndiMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(indiMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-16
        lfcs[is.infinite(lfcs)]=2^16
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statIndiSpeciesTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=sPart, label=getSpeciesName(names(lfcs),taxo), stringsAsFactors = F)
        
        lfcInd = abs(statIndiSpeciesTab$lfc)>2
        sigInd = statIndiSpeciesTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statIndiSpeciesTabModLabel  = statIndiSpeciesTab
        statIndiSpeciesTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        #statTabModLabel$label[!sigInd | !lfcInd] = ""
        
        plotVolcano(statIndiSpeciesTabModLabel, F)
        
        
      }
      
      familyMode = T
      if (familyMode) {
        indiFamilyMean = rowMeans(indiFamilyMat)
        fIndiPart = indiFamilyMean/sum(indiFamilyMean)
        
        hgcIndiFamilyMat = indiFamilyMat[,names(hgcIndiSamples)]
        lgcIndiFamilyMat = indiFamilyMat[,names(lgcIndiSamples)]
        
        lenFamilies = dim(indiFamilyMat)[1]
        
        pValues = sapply(1:lenFamilies, function(ind) {
          mergedVec = c(hgcIndiFamilyMat[ind,],lgcIndiFamilyMat[ind,])
          mergedFac = c(rep("HGC",length(hgcIndiFamilyMat[ind,])),
                        rep("LGC",length(lgcIndiFamilyMat[ind,])))
          
          wilcoxOut = wilcox.test(mergedVec ~ mergedFac)
          return(wilcoxOut$p.value)
        })
        pValues[is.nan(pValues)]=1
        adjPs = p.adjust(pValues,method="BH")
        logAdjPs = -log10(adjPs)
        
        lfcs = sapply(1:lenFamilies,function(ind){
          hmeans = mean(hgcIndiFamilyMat[ind,])
          lmeans = mean(lgcIndiFamilyMat[ind,]) 
          return(hmeans/lmeans)
        })
        names(lfcs)=rownames(indiFamilyMat)
        lfcs[is.nan(lfcs)]=1
        lfcs[lfcs==0]=2^-15
        lfcs[is.infinite(lfcs)]=2^15
        lfcs = log2(lfcs)
        summary(lfcs)
        
        plot(lfcs, -log10(adjPs))
        abline(h=-log10(1e-10))
        abline(v=3)
        abline(v=-3)
        
        statIndiTab = data.frame(lfc=lfcs, sig=logAdjPs, relAbd=fPart, label=names(lfcs), stringsAsFactors = F)
        
        lfcInd = abs(statIndiTab$lfc)>2
        sigInd = statIndiTab$sig>-log10(1e-3)
        unclassifiedInd = !grepl("unclassified", names(lfcs))
        statIndiTabModLabel  = statIndiTab
        statIndiTabModLabel$label[!sigInd | !lfcInd | !unclassifiedInd] = ""
        
        plotVolcano(statIndiTabModLabel, F)
        
        
      }
      save(statIndiTab, statIndiSpeciesTab, file="indi.HGC.LGC.RData")
      
      
    }
    #480735.3
    #767550.3
    #save.image("today.2019.03.11.RData")
    
    hgcDSamples = diseaseRichness[diseaseRichness > as.numeric(quantile(diseaseRichness,hCut))]
    lgcDSamples = diseaseRichness[diseaseRichness < as.numeric(quantile(diseaseRichness,lCut))]
    
  }
  
}

richnessFuncCheckMode = T
if (richnessFuncCheckMode) {
  statsMgsForAllSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.all.txt"
  statsMgsForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.normal.txt"
  statsMgsForDiseaseSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.disease.txt"
  statsMgsForIndustrialSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.industrial.txt"
  statsMgsForTraditionalSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.traditional.txt"
  funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
  
  
  newCutMode = T
  if (newCutMode) {
    adjpCut = 1e-3
    statsMgsForNorm = read.table(statsMgsForNormSamples, sep="\t", header=T, stringsAsFactors = F)
    statsMgsForNormSel = statsMgsForNorm[statsMgsForNorm$qvalue <= adjpCut,]
    # statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 0,]
    # statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -0,]
    statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 2,]
    statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -2,]
    
    hgcSpecies = rownames(statsMgsForNormSelHgc)
    lgcSpecies = rownames(statsMgsForNormSelLgc)
    
    statsMgsForNormSelHgcNew = statsMgsForNormSelHgc[order(statsMgsForNormSelHgc$qvalue),]
    statsMgsForNormSelLgcNew = statsMgsForNormSelLgc[order(statsMgsForNormSelLgc$qvalue),]
    hgcTopSpecies = rownames(statsMgsForNormSelHgcNew[1:25,])
    lgcTopSpecies = rownames(statsMgsForNormSelLgcNew[1:25,])
    
    
  }
  
  loadFuncAnnot = T
  if (loadFuncAnnot) {
    ### virulence factors ###
    vfMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/vfBestMat.RData"
    load(vfMatRData) #vfMatDigit
    vfTerms = colnames(vfMatDigit)
    patricVfDescTabFile = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/PATRIC_sp_gene.txt"
    patricVfDescTab = read.table(patricVfDescTabFile, header=T, sep="\t", stringsAsFactors = F)
    
    ### JGI phenotype ###
    jgiMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/jgiMat.RData"
    load(jgiMatRData)
    jgiTerms = colnames(jgiMat)
    
    ### antibiotic resistance ###
    mustardMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/mustard.mat.RData"
    load(mustardMatRData)
    arTerms = colnames(mustardMat)
    
    ### antismash - secondary metabolism ###
    antismashMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/msp1992.igc2.antismashMat.RData"
    load(antismashMatRData)
    antismashMelt = melt(antismashMat)
    antismashMelt = antismashMelt[antismashMelt$value!=0,]
    antismashTerms = colnames(antismashMat)
    
    ### Cazy  ### 
    newGutCazyMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/CAZY.INRA/igc2.new.cazy.mat.RData"
    load(newGutCazyMatRData)
    colnames(cazyMat) = paste("cazy", colnames(cazyMat), sep=".")
    cazyTerms = colnames(cazyMat)
    
  }
  
  loadFuncMode = T
  if (loadFuncMode) {
    load(funcMatRData)
    
    funcMatTranspose = t(funcMat)
    # hgcFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% hgcSpecies]
    # lgcFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% lgcSpecies]
    hgcFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% hgcSpecies]
    lgcFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% lgcSpecies]
    
    funcs = rownames(funcMatTranspose)
    pvalues = sapply(funcs, function(currFunc) {
      hgcCurrFunc = hgcFuncMat[currFunc,]
      lgcCurrFunc = lgcFuncMat[currFunc,]
      
      hgcAbsent = sum(hgcCurrFunc==0)
      hgcPresent = sum(hgcCurrFunc==1)
      
      lgcAbsent = sum(lgcCurrFunc==0)
      lgcPresent = sum(lgcCurrFunc==1)
      
      contigTab = cbind(HGC=c(present=hgcPresent, absent=hgcAbsent),
                        LGC=c(present=lgcPresent, absent=lgcAbsent))
      pvalue=chisq.test(contigTab)$p.value
      return(pvalue)
    }, simplify = F)
    pvalues = unlist(pvalues)
    
    oddsRatios = sapply(funcs, function(currFunc) {
      hgcCurrFunc = hgcFuncMat[currFunc,]
      lgcCurrFunc = lgcFuncMat[currFunc,]
      
      hgcAbsent = sum(hgcCurrFunc==0)
      hgcPresent = sum(hgcCurrFunc==1)
      
      lgcAbsent = sum(lgcCurrFunc==0)
      lgcPresent = sum(lgcCurrFunc==1)
      
      oddsRatio = (hgcPresent*lgcAbsent)/(hgcAbsent*lgcPresent)
      if (lgcPresent*hgcAbsent == 0) return(Inf) 
      return(oddsRatio)
    }, simplify = F)
    oddsRatios = unlist(oddsRatios)
    
    chisqStatTab = data.frame(term=names(pvalues),
                              pvalue=pvalues,
                              oddsRatio=oddsRatios,
                              logP = -log10(pvalues),
                              logOdds = log(oddsRatios))
    
    # View(chisqStatTab[grepl("cazy",rownames(chisqStatTab)),])
    # View(chisqStatTab[chisqStatTab$pvalue<1e-5,])
    
    View(chisqStatTab[mgePfamTerms,])
    View(chisqStatTab[arTerms,])
    
  }
  
  #inflow
  hgcLgcInflowList=list(HGC=inflow[hgcSpecies], LGC=inflow[lgcSpecies])
  hgcLgcOutflowList=list(HGC=outflow[hgcSpecies], LGC=outflow[lgcSpecies])
  hgcLgcInflowMelt = melt(hgcLgcInflowList)
  hgcLgcOutflowMelt = melt(hgcLgcOutflowList)
  
  t.test(hgcLgcInflowList$HGC, hgcLgcInflowList$LGC)
  t.test(hgcLgcOutflowList$HGC, hgcLgcOutflowList$LGC)
  
  CPCOLS <- c("#493D8C91", "#FF735786")
  
  
  p = ggplot(hgcLgcInflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
    scale_fill_manual(values = c("#ff000055","#2c528c55")) +
    xlab("") +
    ylab("Inflow") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust=1),
          panel.background = element_rect(fill = "white", 
                                          colour = "black", 
                                          size = 0.5, 
                                          linetype = "solid"))
  ggsave("hgc.lgc.inflow.boxplot.pdf", p, width=1.56, height=2.92, units="in")
  
  p = ggplot(hgcLgcOutflowMelt) + geom_boxplot(aes(x=L1, y=value, fill=L1)) +
    scale_fill_manual(values = c("#ff000055","#2c528c55")) +
    xlab("") +
    ylab("Outflow") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust=1),
          panel.background = element_rect(fill = "white", 
                                          colour = "black", 
                                          size = 0.5, 
                                          linetype = "solid"))
  ggsave("hgc.lgc.outflow.boxplot.pdf", p, width=1.56, height=2.92, units="in")
  
  
  boxplot(hgcLgcInflowList)
  boxplot(hgcLgcOutflowList)
  
  t.test(hgcLgcInflowList$hgc, hgcLgcInflowList$lgc)
  t.test(hgcLgcOutflowList$hgc, hgcLgcOutflowList$lgc)
  
  plot(inflow, statsMgsForNorm$lfc[match(names(inflow), rownames(statsMgsForNorm))] )
  plot(outflow, statsMgsForNorm$lfc[match(names(outflow), rownames(statsMgsForNorm))]  )
  
  cor.test(inflow, statsMgsForNorm$lfc[match(names(inflow), rownames(statsMgsForNorm))] , method="spearman" )
  cor.test(outflow, statsMgsForNorm$lfc[match(names(outflow), rownames(statsMgsForNorm))]  , method="spearman")
  
  
}

richnessSpeciesPCoaMode = F
if (richnessSpeciesPCoaMode) {
  statsMgsForAllSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.all.txt"
  statsMgsForNormSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.normal.txt"
  statsMgsForDiseaseSamples = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.disease.txt"
  statsMgsForIndustrialSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.industrial.txt"
  statsMgsForTraditionalSamples  = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.volcanot.out\\richness.volcano.stats.mgs.for.traditional.txt"
  
  prevCutMode = F
  if (prevCutMode) {
    adjpCut = 1e-3
    statsMgsForNorm = read.table(statsMgsForNormSamples, sep="\t", header=T, stringsAsFactors = F)
    statsMgsForNormSel = statsMgsForNorm[statsMgsForNorm$qvalue <= adjpCut,]
    statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 0,]
    statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < 0,]
    # statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 2,]
    # statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -2,]
    
    if (F) {
      hgcCut = quantile(statsMgsForNormSelHgc$relAbd_HGC, 0.75)
      lgcCut = quantile(statsMgsForNormSelLgc$relAbd_LGC, 0.75)
      
      hgcSpecies = rownames(statsMgsForNormSelHgc[statsMgsForNormSelHgc$relAbd_HGC >= hgcCut,])
      lgcSpecies = rownames(statsMgsForNormSelLgc[statsMgsForNormSelLgc$relAbd_LGC >= lgcCut,])
      
      bothSpecies = c(hgcSpecies, lgcSpecies)
    }
    
    
    if (F) {
      hgcCut = quantile(statsMgsForNormSelHgc$qvalue, 0.1)
      lgcCut = quantile(statsMgsForNormSelLgc$qvalue, 0.1)
      
      hgcSpecies = rownames(statsMgsForNormSelHgc[statsMgsForNormSelHgc$qvalue <= hgcCut,])
      lgcSpecies = rownames(statsMgsForNormSelLgc[statsMgsForNormSelLgc$qvalue <= lgcCut,])
      
      bothSpecies = c(hgcSpecies, lgcSpecies)
    }
    
    
    sort(getSpeciesName(hgcSpecies, taxo))
    sort(getSpeciesName(lgcSpecies, taxo))
    
    ###
    statsMgsForNormSelHgcNew = statsMgsForNormSelHgc[order(statsMgsForNormSelHgc$qvalue),]
    statsMgsForNormSelLgcNew = statsMgsForNormSelLgc[order(statsMgsForNormSelLgc$qvalue),]
    hgcTopSpecies = rownames(statsMgsForNormSelHgcNew[1:25,])
    lgcTopSpecies = rownames(statsMgsForNormSelLgcNew[1:25,])
    
    
    getSpeciesName(hgcTopSpecies, taxo)
    getSpeciesName(lgcTopSpecies, taxo)
  }
  
  newCutMode = T
  if (newCutMode) {
    adjpCut = 1e-3
    statsMgsForNorm = read.table(statsMgsForNormSamples, sep="\t", header=T, stringsAsFactors = F)
    statsMgsForNormSel = statsMgsForNorm[statsMgsForNorm$qvalue <= adjpCut,]
    statsMgsForNormSelHgc = statsMgsForNormSel[statsMgsForNormSel$lfc > 2,]
    statsMgsForNormSelLgc = statsMgsForNormSel[statsMgsForNormSel$lfc < -2,]
    
    hgcSpecies = rownames(statsMgsForNormSelHgc)
    lgcSpecies = rownames(statsMgsForNormSelLgc)
    
  }
  
  
  
  diseaseUpSpecies = unique(allPairMeltSel$msp)
  
  allPairMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairMelt.disease.signatures.RData"
  load(allPairMeltRData)
  allPairMeltSel = allPairMelt[allPairMelt$effect_size>=0.5,]
  allPairMeltSel = allPairMeltSel[!is.na(allPairMeltSel$effect_size),]
  diseaseUpSpecies = unique(allPairMeltSel$msp)
  
  
  allPairDownMeltRData = "C://Data//comparative.analysis.healthy.sweden//allPairDownMelt.disease.signatures.RData"
  load(allPairDownMeltRData)
  allPairDownMeltSel = allPairDownMelt[allPairDownMelt$effect_size>=0.5,]
  allPairDownMeltSel = allPairDownMeltSel[!is.na(allPairDownMeltSel$effect_size),]
  diseaseDownSpecies = unique(allPairDownMeltSel$msp)
  
  richSpecMgsMat = mergeMatUpdated[diseaseUpSpecies,]
  richSpecMgsMat = richSpecMgsMat[rowMeans(richSpecMgsMat)!=0,]
  richSpecMgsMat = richSpecMgsMat[,colMeans(richSpecMgsMat)!=0]
  
  pcoaMgsMat = pcoaMat[,1:2]
  
  out = plotContourBoth(pcoaMgsMat, 
                        newCaseByEnteroList$ETF,
                        newContByEnteroList$ETF)
  
  
  par(mfrow=c(1,2))
  zD2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in% newCaseByEnteroList$ETF,])
  zC2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in%newContByEnteroList$ETF,])
  
  zAll = list(x=c(zD2$x, zC2$x),
              y=c(zD2$y, zC2$y),
              z=c(zD2$z, zC2$z))
  
  contourMat = pcoaMgsMat[rownames(pcoaMgsMat) %in% newCaseByEnteroList$ETF,]
  plot(contourMat, xlab="X label", ylab="Y label", col="#80808033", pch=19, cex=.4, xlim=c(-1,1),ylim=c(-1,1))
  
  contour(zC2, drawlabels=FALSE, nlevels=11, col=my.cols, add=TRUE)
  
  
  par(mfrow=c(1,2))
  zD2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in% newCaseByEnteroList$ETB,])
  zC2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in%newContByEnteroList$ETB,])
  
  par(mfrow=c(1,2))
  zD2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in% newCaseByEnteroList$ETF,])
  zC2 = plotContour(pcoaMgsMat[rownames(pcoaMgsMat) %in%newContByEnteroList$ETF,])
  
  vegOut=vegdist(t(richSpecMgsMat),"bray")
  pcoaRichSpecMgsOut=pcoa(vegOut)
  pcoaRichSpecMgsMat= pcoaRichSpecMgsOut$vectors[,1:2]
  
  checkAnosim = F
  if (checkAnosim) {
    ### ETB-country
    etbSampels = unique(c(newCaseByEnteroList$ETB,
                          newContByEnteroList$ETB))
    
    etbClasses = rep("case",length(etbSampels))
    etbClasses[colnames(richSpecMgsMat[,etbSampels]) %in% newContByEnteroList$ETB] = "control"
    
    vegEtbOut=vegdist(t(richSpecMgsMat[,etbSampels]),"bray")
    pcoaEtbOut=pcoa(vegOut)
    anosimEtbOut = anosim(vegEtbOut, etbClasses) 
    # ANOSIM statistic R: 0.05184 
    # Significance: 0.001 
    
    
    ### ETF-country
    
  }
  
  out = plotContourBoth(pcoaRichSpecMgsMat, 
                        newCaseByEnteroList$ETF,
                        newContByEnteroList$ETF)
  out2 = plotContourBoth(pcoaRichSpecMgsMat, 
                         newCaseByEnteroList$ETB,
                         newContByEnteroList$ETB)
  
  
  
  
  
  
  
  
  
  par(mfrow=c(1,2))
  zD2 = plotContour(pcoaRichSpecMgsMat[rownames(pcoaRichSpecMgsMat) %in% newCaseByEnteroList$ETF,])
  zC2 = plotContour(pcoaRichSpecMgsMat[rownames(pcoaRichSpecMgsMat) %in%newContByEnteroList$ETF,])
  
  par(mfrow=c(1,2))
  zD = plotContour(pcoaRichSpecMgsMat[rownames(pcoaRichSpecMgsMat) %in% newCaseByEnteroList$ETB,])
  zC = plotContour(pcoaRichSpecMgsMat[rownames(pcoaRichSpecMgsMat) %in%newContByEnteroList$ETB,])
  
  par(mfrow=c(1,1))
  plotContour(pcoaRichSpecMgsMat)
  
  
  
  
}

checkIndJaccMode = T
if (checkIndJaccMode) {
  
  #same individual from wellness
  #T1D same individual
  #metformin
  #immunotherapy
  #antibiotics
  
  mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
  load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
  dim(mergeMatUpdatedJacc)
  
  
  loadSampleInfo = T
  if (loadSampleInfo) {
    
    
    getJaccPairList <- function(samples, basicMetadata, speciesJaccMat) {
      
      currBasicMetadata = basicMetadata[basicMetadata$sample.ID %in% samples, ]
      samplesByIndividualList = split(currBasicMetadata$sample.ID, currBasicMetadata$metadata.ID)
      
      jaccVecs = lapply(samplesByIndividualList, function(currSamples){
        if (length(currSamples)<=1) return(NA)
        currMat = speciesJaccMat[currSamples, currSamples]
        return(currMat[lower.tri(currMat)])
      })
      return(jaccVecs)
    }
    intraMat=mergeMatUpdatedJacc[!duplicated(basicMetaMapUpdated$metadata.ID),!duplicated(basicMetaMapUpdated$metadata.ID)]
    intraMat = intraMat[rownames(intraMat)%in% newNormalAllIds, colnames(intraMat)%in% newNormalAllIds ]
    allJaccs = intraMat[lower.tri(intraMat)]
    hist(allJaccs)
    
    wellnessSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id28"]
    wellnessJaccs = unlist(getJaccPairList(wellnessSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(wellnessJaccs)
    
    antibioticsSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51"]
    if (F) {
      id51MetaIds = basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id51"]
      id51MetaIds = gsub("_Dag0$","",id51MetaIds)
      id51MetaIds = gsub("_Dag42$","",id51MetaIds)
      id51MetaIds = gsub("_Dag4opt$","",id51MetaIds)
      id51MetaIds = gsub("_Dag180$","",id51MetaIds)
      id51MetaIds = gsub("_Dag4$","",id51MetaIds)
      id51MetaIds = gsub("_Dag8$","",id51MetaIds)
      id51MetaIds = gsub("_Dag8opt$","",id51MetaIds)
      basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id51"] = id51MetaIds
      
    }
    
    antibioticsJaccs = unlist(getJaccPairList(antibioticsSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(antibioticsJaccs)
    
    t1dLuSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id35" 
                                                 & basicMetaMapUpdated$type == "Case"]
    t1dLuJaccs = unlist(getJaccPairList(t1dLuSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(t1dLuJaccs)
    
    contLuSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id35" 
                                                  & basicMetaMapUpdated$type == "Control"]
    contLuJaccs = unlist(getJaccPairList(contLuSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(contLuJaccs)
    
    ##update basicmetadata for metformin data
    id53MetaIds = basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id53"]
    id53MetaIds = gsub("V[0-9]+$","",id53MetaIds)
    basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id53"] = id53MetaIds
    
    exMetfominSamples = c("SRR5580332", "SRR5580333", "SRR5580334", "SRR5580335", "SRR5579968", "SRR5579973", "SRR5579974", "SRR5580118", "SRR5580113", "SRR5580120", "SRR5580121", "SRR5580272", "SRR5580392", "SRR5580393", "SRR5580259", "SRR5580260", "SRR5580394", "SRR5580395", "SRR5579978", "SRR5579979", "SRR5579980", "SRR5579983", "SRR5580316", "SRR5580317", "SRR5580318", "SRR5580320", "SRR5558211", "SRR5558212", "SRR5558213", "SRR5558214", "SRR5558297", "SRR5558298", "SRR5558304", "SRR5580103", "SRR5580313", "SRR5580314", "SRR5580327", "SRR5580039", "SRR5580040", "SRR5580045", "SRR5580046", "SRR5579956", "SRR5580397", "SRR5580399", "SRR5558182", "SRR5558183", "SRR5558191", "SRR5558192", "SRR5558092", "SRR5558093", "SRR5558094", "SRR5558157", "SRR5558158", "SRR5558219", "SRR5558220", "SRR5558169", "SRR5558178", "SRR5558179", "SRR5558181", "SRR5558190", "SRR5558309", "SRR5558313", "SRR5558330", "SRR5558331", "SRR5558061", "SRR5558062", "SRR5558063", "SRR5558068", "SRR5558069", "SRR5558194", "SRR5558195", "SRR5558196", "SRR5558206", "SRR5558102", "SRR5558107", "SRR5558108", "SRR5558232", "SRR5558071", "SRR5558072", "SRR5558073", "SRR5558154", "SRR5558284", "SRR5558386", "SRR5558387", "SRR5558388", "SRR5580390", "SRR5580391", "SRR5580392", "SRR5580393")
    metforminT2dSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" 
                                                        & basicMetaMapUpdated$subtype %in% c("M0","M2","M4")]
    metforminT2dSamples = metforminT2dSamples[!metforminT2dSamples %in% exMetfominSamples]
    
    metforminT2dJaccs = unlist(getJaccPairList(metforminT2dSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(metforminT2dJaccs)
    
    contT2dSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" 
                                                   & basicMetaMapUpdated$subtype %in% c("P0","P2","P4")]
    contT2dSamples = contT2dSamples[!contT2dSamples %in% exMetfominSamples]
    
    contT2dJaccs = unlist(getJaccPairList(contT2dSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(contT2dJaccs)
    
    exRccSamples = c("ERR2235250", "ERR2235262", "ERR2235267", "ERR2235270", "ERR2235273", "ERR2235286", "ERR2235292", "ERR2235295", "ERR2235299", "ERR2235308")
    
    rccSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id41"]
    rccSamples = rccSamples[!rccSamples %in% exRccSamples]
    if (F) {
      rccMetaIds = basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id41"]
      rccMetaIds = gsub("_T[0-9]$","",rccMetaIds)
      
    }
    basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id41"] = rccMetaIds
    
    rccJaccs = unlist(getJaccPairList(rccSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(rccJaccs)
    exNsclcSamples = c("ERR2213662", "ERR2213668", "ERR2213684", "ERR2213687", "ERR2213690", "ERR2213693", "ERR2213698", "ERR2213703", "ERR2213709", "ERR2213730", "ERR2213739", "ERR2213743")
    
    nsclcSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id26"]
    nsclcSamples = nsclcSamples[!nsclcSamples %in% exNsclcSamples]
    if (F) {
      nsclcMetaIds = basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id26"]
      nsclcMetaIds = gsub("_T[0-9]$","",nsclcMetaIds)
    }
    basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id26"] = nsclcMetaIds
    
    nsclcJaccs = unlist(getJaccPairList(nsclcSamples, basicMetaMapUpdated, mergeMatUpdatedJacc))
    hist(nsclcJaccs)
    
    jaccList = list(all=allJaccs[!is.na(allJaccs)], 
                    antibiotics = antibioticsJaccs[!is.na(antibioticsJaccs)],
                    wellness=wellnessJaccs[!is.na(wellnessJaccs)],
                    t1d_LU=t1dLuJaccs[!is.na(t1dLuJaccs)],
                    cont_LU=contLuJaccs[!is.na(contLuJaccs)],
                    metformin_T2D=metforminT2dJaccs[!is.na(metforminT2dJaccs)],
                    cont_T2D=contT2dJaccs[!is.na(contT2dJaccs)],
                    rcc_immunotherapy=rccJaccs[!is.na(rccJaccs)],
                    nsclc_immunotherapy=nsclcJaccs[!is.na(nsclcJaccs)])
    par(mar=c(10,4,4,1))
    boxplot(jaccList,las=2, ylab="Jaccard")
    jaccMeans = unlist(lapply(jaccList, mean))
    jaccSds = unlist(lapply(jaccList, sd))
    
    newNames = c("Inter-individual",
                 "Antibiotics (DK)",
                 "Intra-individual (ES)",
                 "Intra-individual (T1D, LU)",
                 "Intra-individual (LU)",
                 "Metformin (T2D, ES)",
                 "Intra-individual (T2D, ES)",
                 "Immunotherapy (RCC, FR)",
                 "Immunotherapy (NSCLC, FR)")
    
    orderedNames = c("Inter-individual",
                     "Intra-individual (LU)",
                     "Intra-individual (ES)",
                     "Intra-individual (T2D, ES)",
                     "Intra-individual (T1D, LU)",
                     "Antibiotics (DK)",
                     "Metformin (T2D, ES)",
                     "Immunotherapy (NSCLC, FR)",
                     "Immunotherapy (RCC, FR)")
    
    classes = c("Inter",
                "Perturb",
                "Intra.normal",
                "Intra.disease",
                "Intra.normal",
                "Perturb",
                "Intra.disease",
                "Perturb",
                "Perturb")
    jaccStatTab = data.frame(name=newNames,
                             mean=jaccMeans,
                             sd=jaccSds,
                             class=classes)
    
    jaccStatTab$name = factor(jaccStatTab$name, levels=orderedNames)
    
    ggplot(jaccStatTab, aes(x=name, y=mean, size=sd)) +
      geom_point(aes(fill=factor(class)), colour="black", shape=21) +
      geom_vline(xintercept = c(1.5, 3.5, 5.5))+
      xlab("") + ylab("Jaccard") +
      scale_size(range = c(2, 10))+
      ylim(c(0.1,0.8))+
      theme(panel.grid.major = element_line(colour="gray"), 
            axis.text.x = element_text(angle = 60, hjust=1),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background = element_rect(fill = "white", 
                                            colour = "black", 
                                            size = 0.5, 
                                            linetype = "solid")) 
    
  }
}


####### UNUSED SOURCE CODES #######
checkPlotWorldMap  = F
if (checkPlotWorldMap) {
  library(rgdal)
  #download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="DATA/world_shape_file.zip")
  
  
  my_spdf <- readOGR( 
    dsn= path.expand("J://Deposit//Project//2018_microbiome_atlas//functional.annotation"), 
    layer="TM_WORLD_BORDERS_SIMPL-0.3",
    verbose=FALSE
  )
  
  library(broom)
  spdf_fortified <- tidy(my_spdf, region = "NAME")
  
  # Plot it
  library(ggplot2)
  ggplot() +
    geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white") +
    theme_void() 
  
  par(mar=c(0,0,0,0))
  plot(my_spdf, col="#f2f2f2", bg="skyblue", lwd=0.25, border=0 )
  
}

checkDatasetInfoNew2 = F
if (checkDatasetInfoNew2) {
  
  datasets = unique(basicMetaMapUpdated$dataset.ID)
  
  checkDuplicates = F
  if (checkDuplicates) {
    for( ds in datasets ) {
      isDuplicated = duplicated(basicMetaMapUpdated[basicMetaMapUpdated$dataset.ID==ds,"metadata.ID"])
      print(ds)
      print(table(isDuplicated))
    }
    
  }
  
  checkScatterPlotMode = T 
  if (checkScatterPlotMode) {
    
    diseaseMgsMat = mergeMatUpdated[,diseaseUniqAll] 
    industrialMgsMat = mergeMatUpdated[,industrializedSamples]
    traditionalMgsMat = mergeMatUpdated[,traditionalSamples]
    
    diseaseMGsRowMs = rowMeans(diseaseMgsMat)
    industrialMGsRowMs = rowMeans(industrialMgsMat)
    traditionalMGsRowMs = rowMeans(traditionalMgsMat)
    
    industrialSpecies = names(industrialMGsRowMs[industrialMGsRowMs!=0])
    traditionalSpecies = names(traditionalMGsRowMs[traditionalMGsRowMs!=0])
    
    diseaseTop10 =  names(sort(diseaseMGsRowMs, decreasing = T)[1:10])
    industrialTop10 =  names(sort(industrialMGsRowMs, decreasing = T)[1:10])
    traditionalTop10 =  names(sort(traditionalMGsRowMs, decreasing = T)[1:10])
    
    tradOnlySpecies = traditionalTop10[!traditionalTop10 %in% industrialTop10]
    industOnlySpecies = industrialTop10[!industrialTop10 %in% traditionalTop10]
    commSpecies = traditionalTop10[traditionalTop10 %in% industrialTop10]
    
    
    freqs = rowSums(mergeMatUpdated > 0)
    avgs = rowMeans(mergeMatUpdated)
    maxs = rowMaxs(mergeMatUpdated)
    
    cols = rep("#A9A9A944", length(freqs))
    names(cols) = names(freqs)
    cols[tradOnlySpecies] = "#FF000088"
    cols[industOnlySpecies] = "#0000FF55"
    cols[commSpecies] = "#800080"
    
    
    plot(avgs, freqs, col=cols, pch=16, log='x', las=1, xlab="",ylab="")
    abline(v=quantile(avgs, 0.5))
    abline(v=quantile(avgs, 0.99))
    
    ###
    plot(maxs, avgs, log="xy")
    
    
    pairs(cbind(diseaseMGsRowMs,
                industrialMGsRowMs,
                traditionalMGsRowMs), log="xy")
    
    diseaseSpecies = rowSums(mergeMatUpdated[,diseaseUniqAll] >0)
    diseaseSpecies = names(diseaseSpecies[diseaseSpecies>0])
    
    normSpecies = rowSums(mergeMatUpdated[,normUniqAll] >0)
    normSpecies = names(normSpecies[normSpecies>0])
    
    
    industSpecies = rowSums(mergeMatUpdated[,industrializedSamples] >0)
    industSpecies = names(industSpecies[industSpecies>0])
    
    tradSpecies = rowSums(mergeMatUpdated[,traditionalSamples] >0)
    tradSpecies = names(tradSpecies[tradSpecies>0])
    
    tradOnlySpecies = tradSpecies[!tradSpecies %in% industSpecies]
    industOnlySpecies = industSpecies[!industSpecies %in% tradSpecies]
    diseaseOnlySpecies = diseaseSpecies[!diseaseSpecies %in% normSpecies]
    
    
  }
  
  getDiseaseSigByPlsDaMode = T
  if (getDiseaseSigByPlsDaMode) {
    
    genusMode = T
    if (genusMode) {
      
      diseaseGenusMat = genusMatUpdated[,diseaseUniqAll] 
      normalGenusMat = genusMatUpdated[,normUniqAll]
      industrialGenusMat = genusMatUpdated[,industrializedSamples]
      traditionalGenusMat = genusMatUpdated[,traditionalSamples]
      
      diseaseVipCollectList = list()
      for( ind in seq_along(cohortIdByDiseases)) {
        currDiseaseCohorts = cohortIdByDiseases[[ind]]
        print(names(cohortIdByDiseases)[ind])
        currDiseaseSamples = unique(do.call(c,diseaseUniqSubjectList[currDiseaseCohorts]))
        currDiseaseGenusMat = genusMatUpdated[,currDiseaseSamples]
        currDataMat = cbind(industrialGenusMat, currDiseaseGenusMat)
        
        currPlsClasses = c(rep("health", dim(industrialGenusMat)[2]),
                           rep("disease", dim(currDiseaseGenusMat)[2]))
        currPlsCols = c(rep("#A9A9A988",dim(industrialGenusMat)[2]), #health: gray
                        rep("#FF000088",dim(currDiseaseGenusMat)[2])) #disease: red
        
        currPlsDaOut = opls(t(currDataMat), currPlsClasses, plotL=F, predI=2)
        currVips  = getVipVn(currPlsDaOut)
        currVipsAboveTwo = names(currVips[currVips>=2])
        currVipList = getVipList(currVips, currDataMat, currPlsClasses)
        
        diseaseVipCollectList[[ind]] = currVipsAboveTwo
        plotMode = F
        if (plotMode) {
          fileName = paste(names(cohortIdByDiseases)[ind], "pls_da" , "pdf", sep=".")
          pdf(fileName, width = 2, height = 2)
          #plot(currPlsDaOut@scoreMN[,1:2], col=currPlsCols, pch=16,cex=1,las=1,xlab="",ylab="", xaxt="n",yaxt="n")
          s.class(currPlsDaOut@scoreMN[,1:2], fac = factor(currPlsCols), sub = names(cohortIdByDiseases)[ind],
                  cell = 2, axesell = FALSE, csta = 0, grid=F,addaxes = F,
                  col = unique(currPlsCols), clabel = FALSE)
          dev.off()
        }
      }
      names(diseaseVipCollectList) = names(cohortIdByDiseases)
      
      diseaseVipTab = melt(diseaseVipCollectList)
      diseaseVipTab$value = as.character(diseaseVipTab$value)
      diseaseVipTab$L1 = as.character(diseaseVipTab$L1)
      diseaseVipTab$count = 1
      diseaseVipMat = acast(diseaseVipTab, value ~ L1, value.var = "count")
      diseaseVipMat[is.na(diseaseVipMat)]=0
      colSums(diseaseVipMat)
      dim(diseaseVipMat)
      heatmap.2(diseaseVipMat, trace="none",
                col = colorRampPalette(c("white","red"))(256),
                margins = c(5,12),
                key=F,lhei=c(1,20),
                cexCol = 0.5,cexRow = 0.5,
                sepcolor = "gray",
                colsep=seq(0, dim(diseaseVipMat)[2]),
                rowsep=seq(0, dim(diseaseVipMat)[1]),
                sepwidth = c(0.01,0.01))
      
      
      allInputList=c(newDiseaseList, newNormalAllList)
      allInputList = newNormalAllList[c("US_id11",
                                        "US_id36",
                                        "US_id34",
                                        "US_id43")] 
      
      allInputList = newNormalAllList[c("UK_id29",
                                        "UK_id39")] 
      
      allInputList = newNormalAllList[c("Italy_id46",
                                        "Italy_id37")] 
      
      allInputList = newNormalAllList[c("Sweden_id28",
                                        "Sweden_id1",
                                        "Sweden_id14")] 
      
      allInputList = newNormalAllList[c("China_id10",
                                        "China_id12",
                                        "China_id2",
                                        "China_id20",
                                        "China_id31",
                                        "China_id6",
                                        "China_id9")] 
      
      allInputEnteroList = getEntGeoPlot(grepEnteroTypesFrom(allInputList, basicMetaMapUpdated))
      
      entGeoTab = getEntGeoPlot(allInputEnteroList)
      
      avgGenusDiseaseMat = getAvgAbdMatBy(genusMatUpdated, c(diseaseUniqList) )
      avgGenusDiseaseMeans = rowMeans(avgGenusDiseaseMat)
      topGenera = names(sort(avgGenusDiseaseMeans,decreasing = T)[1:108])
      
      avgGenusDiseaseMatSel = avgGenusDiseaseMat[rownames(diseaseVipMat),]
      avgGenusDiseaseMatSel2 = avgGenusDiseaseMat[topGenera,]
      
      cor(avgGenusDiseaseMatSel)
      cc = melt( cor(avgGenusDiseaseMatSel))
      cc = cc[cc$value>=0.99,]
      cc = cc[cc$Var1 != cc$Var2,]
      View(cc)
      
      heatmap.2(cor(avgGenusDiseaseMatSel),trace="none",key=F,
                col = colorRampPalette(c("blue","white","red"))(256))
      heatmap.2(cor(avgGenusDiseaseMat),trace="none",key=F,
                col = colorRampPalette(c("blue","white","red"))(256))
      
      
      set.seed(1)
      targetMat = avgGenusDiseaseMatSel2
      diseaseTsneOut = Rtsne(t(targetMat), dims=2, perplexity = 5)
      rownames(diseaseTsneOut$Y) = colnames(targetMat)
      plot(diseaseTsneOut$Y)
      
      tsneTab = data.frame(diseaseTsneOut$Y, type=rownames(diseaseTsneOut$Y))
      ggplot(tsneTab, aes(x=X1 , y=X2, label=type)) + 
        geom_point(aes(colour = type, alpha=.02)) + 
        theme(panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text.x = element_blank(),
              axis.text.y = element_blank(),
              strip.text.y = element_blank(),
              legend.position="none", 
              panel.border = element_rect(colour = "black", fill=NA, size=.5))  + 
        geom_text_repel(size=3, force=2) + xlab("") + ylab("")
      
    }
    
    mgsMode = T
    if (mgsMode) {
      
      countryMode = T
      if (countryMode) {
        
        
        diseaseMgsCountryVipCollectList = list()
        for (ind in 1:30) {
          currDiseaseCohortId = diseaseWithCountryPairTab2[ind,1]
          currCountryId = diseaseWithCountryPairTab2[ind,2]
          currDiseaseSamples = diseaseUniqSubjectList[[currDiseaseCohortId]]
          currContSamples = normUniqCountryList[[currCountryId]] 
          currDiseaseMat = mergeMatUpdated[,currDiseaseSamples]
          currContMat = mergeMatUpdated[,currContSamples]
          currDataMat = cbind(currContMat, currDiseaseMat)
          
          currPlsClasses = c(rep("health", dim(currContMat)[2]),
                             rep("disease", dim(currDiseaseMat)[2]))
          currPlsCols = c(rep("#A9A9A988",dim(currContMat)[2]), #health: gray
                          rep("#FF000088",dim(currDiseaseMat)[2])) #disease: red
          
          currPlsDaOut = opls(t(currDataMat), currPlsClasses, plotL=F, predI=2)
          currVips  = getVipVn(currPlsDaOut)
          currVipsAboveOne = (currVips[currVips>=1])
          currVipList = getVipList(currVips, currDataMat, currPlsClasses)
          diseaseMgsCountryVipCollectList[[ind]] = currVipsAboveOne
          
          
          plotMode = T
          if (plotMode) {
            currTopVipList = takeTopOfList(currVipList, 10)
            currTopVipList = labelVipList(currTopVipList, taxo)
            
            fileName = paste(currDiseaseCohortId, "country.MGS.VIP.top10" , "pdf", sep=".")
            pdf(fileName, width = 7, height = 5)
            plotVipList(currTopVipList, marVec = c(5,25,4,1))
            dev.off()
            
            
            fileName = paste(currDiseaseCohortId, "country.MGS.pls_da" , "pdf", sep=".")
            pdf(fileName, width = 2, height = 2)
            #plot(currPlsDaOut@scoreMN[,1:2], col=currPlsCols, pch=16,cex=1,las=1,xlab="",ylab="", xaxt="n",yaxt="n")
            s.class(currPlsDaOut@scoreMN[,1:2], fac = factor(currPlsCols), sub = currDiseaseCohortId,
                    cell = 2, axesell = FALSE, csta = 0, grid=F,addaxes = F,
                    col = unique(currPlsCols), clabel = FALSE)
            dev.off()
          }
          
          
          
        }
        names(diseaseMgsCountryVipCollectList) = diseaseWithCountryPairTab2[,1]
        
        diseaseMgsCountryVipCollectList2 = lapply(diseaseMgsCountryVipCollectList, names)
        diseaseMgsVipTab = melt(diseaseMgsCountryVipCollectList2)
        diseaseMgsVipTab$value = as.character(diseaseMgsVipTab$value)
        diseaseMgsVipTab$L1 = as.character(diseaseMgsVipTab$L1)
        diseaseMgsVipTab$count = 1
        diseaseMgsVipTab$species = sapply(diseaseMgsVipTab$value, function(x) getSpeciesName(x, taxo))
        
        
        
        
      }
      allCompMode = T
      if (allCompMode) {
        diseaseMat = mergeMatUpdated[,diseaseUniqAll] 
        normalMat = mergeMatUpdated[,normUniqAll]
        industrialMat = mergeMatUpdated[,industrializedSamples]
        traditionalMat = mergeMatUpdated[,traditionalSamples]
        
        diseaseMgsVipCollectList = list()
        for( ind in seq_along(cohortIdByDiseases)) {
          currDiseaseCohorts = cohortIdByDiseases[[ind]]
          print(names(cohortIdByDiseases)[ind])
          currDiseaseSamples = unique(do.call(c,diseaseUniqSubjectList[currDiseaseCohorts]))
          currDiseaseMat = mergeMatUpdated[,currDiseaseSamples]
          currDataMat = cbind(industrialMat, currDiseaseMat)
          
          currPlsClasses = c(rep("health", dim(industrialMat)[2]),
                             rep("disease", dim(currDiseaseMat)[2]))
          currPlsCols = c(rep("#A9A9A988",dim(industrialMat)[2]), #health: gray
                          rep("#FF000088",dim(currDiseaseMat)[2])) #disease: red
          
          currPlsDaOut = opls(t(currDataMat), currPlsClasses, plotL=F, predI=2)
          currVips  = getVipVn(currPlsDaOut)
          currVipsAboveTwo = names(currVips[currVips>=2])
          currVipList = getVipList(currVips, currDataMat, currPlsClasses)
          
          diseaseMgsVipCollectList[[ind]] = currVipsAboveTwo
          
          plotMode = T
          if (plotMode) {
            fileName = paste(names(cohortIdByDiseases)[ind], "MGS.pls_da" , "pdf", sep=".")
            pdf(fileName, width = 2, height = 2)
            #plot(currPlsDaOut@scoreMN[,1:2], col=currPlsCols, pch=16,cex=1,las=1,xlab="",ylab="", xaxt="n",yaxt="n")
            s.class(currPlsDaOut@scoreMN[,1:2], fac = factor(currPlsCols), sub = names(cohortIdByDiseases)[ind],
                    cell = 2, axesell = FALSE, csta = 0, grid=F,addaxes = F,
                    col = unique(currPlsCols), clabel = FALSE)
            dev.off()
          }
        }
        names(diseaseMgsVipCollectList) = names(cohortIdByDiseases)
        
        diseaseVipTab = melt(diseaseMgsVipCollectList)
        diseaseVipTab$value = as.character(diseaseVipTab$value)
        diseaseVipTab$L1 = as.character(diseaseVipTab$L1)
        diseaseVipTab$count = 1
        diseaseVipTab$species = sapply(diseaseVipTab$value, function(x) getSpeciesName(x, taxo))
        
        diseaseVipMat = acast(diseaseVipTab, value ~ L1, value.var = "count")
        diseaseVipMat[is.na(diseaseVipMat)]=0
        colSums(diseaseVipMat)
        dim(diseaseVipMat)
        
        
        heatmap.2(diseaseVipMat, trace="none",
                  col = colorRampPalette(c("white","red"))(256),
                  margins = c(5,12),
                  key=F,lhei=c(1,20),
                  cexCol = 0.5,cexRow = 0.5,
                  sepcolor = "gray",
                  colsep=seq(0, dim(diseaseVipMat)[2]),
                  rowsep=seq(0, dim(diseaseVipMat)[1]),
                  sepwidth = c(0.01,0.01))
        
      }
      
      
      allInputList=c(newDiseaseList, newNormalAllList)
      allInputList = newNormalAllList[c("US_id11",
                                        "US_id36",
                                        "US_id34",
                                        "US_id43")] 
      
      allInputList = newNormalAllList[c("UK_id29",
                                        "UK_id39")] 
      
      allInputList = newNormalAllList[c("Italy_id46",
                                        "Italy_id37")] 
      
      allInputList = newNormalAllList[c("Sweden_id28",
                                        "Sweden_id1",
                                        "Sweden_id14")] 
      
      allInputList = newNormalAllList[c("China_id10",
                                        "China_id12",
                                        "China_id2",
                                        "China_id20",
                                        "China_id31",
                                        "China_id6",
                                        "China_id9")] 
      
      allInputEnteroList = getEntGeoPlot(grepEnteroTypesFrom(allInputList, basicMetaMapUpdated))
      
      entGeoTab = getEntGeoPlot(allInputEnteroList)
      
      avgGenusDiseaseMat = getAvgAbdMatBy(genusMatUpdated, c(diseaseUniqList) )
      avgGenusDiseaseMeans = rowMeans(avgGenusDiseaseMat)
      topGenera = names(sort(avgGenusDiseaseMeans,decreasing = T)[1:108])
      
      avgGenusDiseaseMatSel = avgGenusDiseaseMat[rownames(diseaseVipMat),]
      avgGenusDiseaseMatSel2 = avgGenusDiseaseMat[topGenera,]
      
      cor(avgGenusDiseaseMatSel)
      cc = melt( cor(avgGenusDiseaseMatSel))
      cc = cc[cc$value>=0.99,]
      cc = cc[cc$Var1 != cc$Var2,]
      View(cc)
      
      heatmap.2(cor(avgGenusDiseaseMatSel),trace="none",key=F,
                col = colorRampPalette(c("blue","white","red"))(256))
      heatmap.2(cor(avgGenusDiseaseMat),trace="none",key=F,
                col = colorRampPalette(c("blue","white","red"))(256))
      
      
      set.seed(1)
      targetMat = avgGenusDiseaseMatSel2
      diseaseTsneOut = Rtsne(t(targetMat), dims=2, perplexity = 5)
      rownames(diseaseTsneOut$Y) = colnames(targetMat)
      plot(diseaseTsneOut$Y)
      
      tsneTab = data.frame(diseaseTsneOut$Y, type=rownames(diseaseTsneOut$Y))
      ggplot(tsneTab, aes(x=X1 , y=X2, label=type)) + 
        geom_point(aes(colour = type, alpha=.02)) + 
        theme(panel.background = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              strip.text.x = element_blank(),
              axis.text.y = element_blank(),
              strip.text.y = element_blank(),
              legend.position="none", 
              panel.border = element_rect(colour = "black", fill=NA, size=.5))  + 
        geom_text_repel(size=3, force=2) + xlab("") + ylab("")
      
    }
  }
  
  getDiseaseSigMode = F
  if (getDiseaseSigMode) {
    
    sampleTab = melt(mixedUniqSubjectList)
    colnames(sampleTab) = c("sample.ID", "cohort.desc")
    sampleTab$cohort.ID =  gsub(".*_id","id",sampleTab$cohort.desc)
    
    sampleTab$type = NA
    sampleTab$type[sampleTab$cohort.desc %in% names(diseaseUniqSubjectList)] = "disease"
    sampleTab$type[sampleTab$cohort.desc %in% names(normUniqSubjectList)] = "health"
    sampleTab$type = as.factor(sampleTab$type)
    
    sampleTab$geography = NA
    sapply(seq_along(cohortIdByGeos), function(id) {
      currCountry = names(cohortIdByGeos)[id]
      currCohorts = cohortIdByGeos[[id]]
      sampleTab$geography[sampleTab$cohort.desc %in% currCohorts] <<- currCountry
    })
    sampleTab$geography = as.factor(sampleTab$geography)
    
    sampleTab$sequencer = NA
    seqInfoTab = unique(basicMetaMapUpdated[,c("dataset.ID","Sequencer")])
    apply(seqInfoTab, 1, function(currRow){
      currId = currRow[1]
      currSeq = currRow[2]
      sampleTab$sequencer[sampleTab$cohort.ID == currId] <<- currSeq
    })
    
    sampleTab$Imputation =  basicMetaMapUpdated$Imputation[match(sampleTab$sample.ID, basicMetaMapUpdated$sample.ID)]
    sampleTab$WhatImputed =  basicMetaMapUpdated$WhatImputed[match(sampleTab$sample.ID, basicMetaMapUpdated$sample.ID)]
    sampleTab$Age = basicMetaMapUpdated$Age[match(sampleTab$sample.ID, basicMetaMapUpdated$sample.ID)]
    sampleTab$BMI = basicMetaMapUpdated$BMI[match(sampleTab$sample.ID, basicMetaMapUpdated$sample.ID)]
    
    matchedMgsTab = t(mergeMatUpdated[,match(sampleTab$sample.ID, colnames(mergeMatUpdated))])
    matchedGenusTab = t(genusMatUpdated[,match(sampleTab$sample.ID, colnames(genusMatUpdated))])
    matchedFamilyTab = t(familyMatUpdated[,match(sampleTab$sample.ID, colnames(familyMatUpdated))])
    
    matchedMgsTab = matchedMgsTab*10^9
    matchedGenusTab = matchedGenusTab*10^9
    matchedFamilyTab = matchedFamilyTab*10^9
    
    matchedMgsTab = matchedMgsTab+1
    matchedGenusTab = matchedGenusTab+1
    matchedFamilyTab = matchedFamilyTab+1
    
    matchedMgsTab = log10(matchedMgsTab)
    matchedGenusTab = log10(matchedGenusTab)
    matchedFamilyTab = log10(matchedFamilyTab)
    
    # 
    # inputMgsTab = cbind(sampleTab, matchedMgsTab)
    # inputGenusTab = cbind(sampleTab, matchedGenusTab)
    # inputFamilyTab = cbind(sampleTab, matchedFamilyTab)
    featureTab = matchedMgsTab
    
    lmerSignatureOut = lmerDiseaseSignature(sampleTab, featureTab)
    geoMat = lmerSignatureOut$geography
    diseaseMat = lmerSignatureOut$disease
    diseaseMat[is.na(diseaseMat[,1]),1] = 0
    diseaseMat[is.na(diseaseMat[,2]),2] = 1
    
    featureTab = matchedGenusTab
    lmerGenusSignatureOut = lmerDiseaseSignature(sampleTab, featureTab)
    geoGenusMat = lmerGenusSignatureOut$geography
    diseaseGenusMat = lmerGenusSignatureOut$disease
    
    genusCorrected = getCorrectedFeatureTab(lmerGenusSignatureOut, sampleTab, featureTab)
    
    getCorrectedFeatureTab <- function(signatureOut, sampleTab, featureTab) {
      featureMat = as.matrix(featureTab)
      geoMat = signatureOut$geography
      sampleTab = sampleTab[match(rownames(featureTab), sampleTab$sample.ID),]
      geoFeatureMat = featureMat
      geoFeatureMat = geoFeatureMat[,match(colnames(geoMat), colnames(geoFeatureMat))]
      rownames(geoFeatureMat) = sampleTab$geography
      
      
      invisible(sapply(rownames(geoMat), function(currGeo) {
        currGeoVec = geoMat[currGeo,]
        geoFeatureMat[rownames(geoFeatureMat) == currGeo] <<- geoFeatureMat[rownames(geoFeatureMat) == currGeo] - currGeoVec
      }))
      rownames(geoFeatureMat) = sampleTab$sample.ID
      
      set.seed(1)
      tsneOut = Rtsne(geoFeatureMat, dims=2, perplexity = 100)
      set.seed(1)
      tsneUncorrectedOut = Rtsne(featureMat, dims=2, perplexity = 100 )
      
      
      geoCols = rep("gray",3781)
      geoCols[sampleTab$geography=="France"]="red"
      geoCols[sampleTab$geography=="China"]="blue"
      
      plot(tsneOut$Y, col=geoCols)
      plot(tsneUncorrectedOut$Y, col=geoCols)
      
      avgMgsCorrectMat = getAvgAbdMatBy(t(geoFeatureMat), mixedUniqSubjectList )
      avgMgsUnCorrectMat = getAvgAbdMatBy(t(featureMat), mixedUniqSubjectList )
      
      
      heatmap.2(cor(avgMgsCorrectMat), 
                trace="none", 
                margins = c(7,7), 
                col = colorRampPalette(c("blue","white","red"))(256))
      
      
      heatmap.2(cor(avgMgsUnCorrectMat), 
                trace="none", 
                margins = c(7,7), 
                col = colorRampPalette(c("blue","white","red"))(256))
      
      
      set.seed(1)
      tsneCorrectedAvgOut = Rtsne(t(avgMgsCorrectMat), dims=2, perplexity = 10)
      set.seed(1)
      tsneUncorrectedAvgOut = Rtsne(t(avgMgsUnCorrectMat), dims=2, perplexity = 10 )
      
      dataCols = rep("#80808077", 65)
      dataCols[colnames(avgMgsCorrectMat)%in%cohortIdByGeos$Sweden]="blue" 
      
      plot(tsneCorrectedAvgOut$Y, col=dataCols)
      plot(tsneUncorrectedAvgOut$Y, col=dataCols)
      
      
      return(list(corrected=geoFeatureMat))
      
    }
    
    datasetID = "Ankylosing_CN_id9"
    
    ##### contructing input matrix ####
    # rows: 
    # all possible samples --> DONE
    
    # columns: 
    # sample ID
    # cohort ID
    # healthy/disease --> DONE
    # subtypes
    # cohorts
    # geography
    # sequencer
    # all MGSs
    # all genus
    # all family
    
    # additional variables
    # index indicating MGSs
    # index indicating genus
    # index indicating family
    
    
  }
  
  jaccNetMode = T
  if (jaccNetMode) {
    
    cohortWiseJaccNetMode = T
    if (cohortWiseJaccNetMode) {
      #wellness id28
      #kcl twin?? id39
      #luxembourg id35?? 
      #metformin treatment??? id53??
      
      ###
      
      mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
      load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
      
      wellnessSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id28"]
      wellnessSubjs = basicMetaMapUpdated$metadata.ID[basicMetaMapUpdated$dataset.ID=="id28"]
      wellnessVisits = basicMetaMapUpdated$subtype[basicMetaMapUpdated$dataset.ID=="id28"]
      
      kclTwinSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id39"] 
      # intereting six subjects were connected each other
      
      us43Samples = normUniqSubjectList$US_id43
      hmpSamples = normUniqSubjectList$US_id36
      
      # intereting six subjects were connected each other
      
      t1dLuSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id35" 
                                                   & basicMetaMapUpdated$type == "Case"]
      
      contLuSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id35" 
                                                    & basicMetaMapUpdated$type == "Control"]
      
      t1dContLuSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id35"]
      
      metforminT2dSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" 
                                                          & basicMetaMapUpdated$subtype %in% c("M0","M2","M4")]
      
      contT2dSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" 
                                                     & basicMetaMapUpdated$subtype %in% c("P0","P2","P4")]
      
      metforminContT2dSamples = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" 
                                                              & basicMetaMapUpdated$subtype %in% c("M0","M2","M4","P0","P2","P4")]
      
      
      
      wellnessNodeFile = "wellness.jacc.net.nodeTable.txt"
      wellnessEdgeFile = "wellness.jacc.net.edgeTable.txt"
      wellnessOut = writeJaccNetToCytoscape(wellnessSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                            wellnessNodeFile, wellnessEdgeFile)
      
      kclTwinNodeFile = "kcl.twin.jacc.net.nodeTable.txt"
      kclTwinEdgeFile = "kcl.twin.jacc.net.edgeTable.txt"
      kclTwinOut = writeJaccNetToCytoscape(kclTwinSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                           kclTwinNodeFile, kclTwinEdgeFile)
      
      t1dLuNodeFile = "t1dLu.jacc.net.nodeTable.txt"
      t1dLuEdgeFile = "t1dLu.jacc.net.edgeTable.txt"
      t1dLuOut = writeJaccNetToCytoscape(t1dLuSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                         t1dLuNodeFile, t1dLuEdgeFile)
      
      contLuNodeFile = "contLu.jacc.net.nodeTable.txt"
      contLuEdgeFile = "contLu.jacc.net.edgeTable.txt"
      contLuOut = writeJaccNetToCytoscape(contLuSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                          contLuNodeFile, contLuEdgeFile)
      
      us43NodeFile = "us43.jacc.net.nodeTable.txt"
      us43EdgeFile = "us43.jacc.net.edgeTable.txt"
      us43Out = writeJaccNetToCytoscape(us43Samples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                        us43NodeFile, us43EdgeFile)
      
      t1dContLuNodeFile = "t1dContLu.jacc.net.nodeTable.txt"
      t1dContLuEdgeFile = "t1dContLu.jacc.net.edgeTable.txt"
      t1dContLuOut = writeJaccNetToCytoscape(t1dContLuSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                             t1dContLuNodeFile, t1dContLuEdgeFile)
      
      metContT2dNodeFile = "metContT2d.jacc.net.nodeTable.txt"
      metContT2dEdgeFile = "metContT2d.jacc.net.edgeTable.txt"
      metContT2dOut = writeJaccNetToCytoscape(metforminContT2dSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                              metContT2dNodeFile, metContT2dEdgeFile)
      
      
      hmpNodeFile = "hmp.jacc.net.0.5.nodeTable.txt"
      hmpEdgeFile = "hmp.jacc.net.0.5.edgeTable.txt"
      hmpOut = writeJaccNetToCytoscape(hmpSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                       hmpNodeFile, hmpEdgeFile)
      
      hmpNodeFile = "hmp.jacc.net.0.33.nodeTable.txt"
      hmpEdgeFile = "hmp.jacc.net.0.33.edgeTable.txt"
      hmpOut = writeJaccNetToCytoscape(hmpSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                       hmpNodeFile, hmpEdgeFile, 0.33)
      
      hmpNodeFile = "hmp.jacc.net.0.4.nodeTable.txt"
      hmpEdgeFile = "hmp.jacc.net.0.4.edgeTable.txt"
      hmpOut = writeJaccNetToCytoscape(hmpSamples, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                       hmpNodeFile, hmpEdgeFile, 0.4)
      
      japanNodeFile = "japan.jacc.net.0.33.nodeTable.txt"
      japanEdgeFile = "japan.jacc.net.0.33.edgeTable.txt"
      japanOut = writeJaccNetToCytoscape(normUniqCountryList$Japan, mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                         japanNodeFile, japanEdgeFile, 0.33)
      
      antibioticsNodeFile = "antibiotics.jacc.net.0.5.nodeTable.txt"
      antibioticsEdgeFile = "antibiotics.jacc.net.0.5.edgeTable.txt"
      antibioticsOut = writeJaccNetToCytoscape(basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51"], mergeMatUpdatedJacc, basicMetaMapUpdated, 
                                               antibioticsNodeFile, antibioticsEdgeFile, 0.5)
      
    }
    
    
    geoJaccNetMode = T
    if (geoJaccNetMode) {
      mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
      load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
      mergeMatUpdatedJaccSelected = mergeMatUpdatedJacc[normUniqAll,normUniqAll]
      
      jacc05Mode = T
      if (jacc05Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.5)
        
        geoList = normUniqSubjectList
        geoClusterScores = getClusterScores(jaccNet, geoList)
        
        jaccClusterInfo = annotateModulesByCC(jaccNet, geoList)
        writeModuleNetworkAnnot(jaccClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGeo.pruned.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGeo.pruned.txt")
        
        
        geoMergedList = normUniqCountryList
        geoMergedClusterScores = getClusterScores(jaccNet, geoMergedList)
        
        jaccMergedClusterInfo = annotateModulesByCC(jaccNet, geoMergedList)
        writeModuleNetworkAnnot(jaccMergedClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGeoMerged.pruned.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGeoMerged.pruned.txt")
        
      }
      
      jacc033Mode = F
      if (jacc033Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.33)
        
        geoList = normUniqSubjectList
        geoClusterScores = getClusterScores(jaccNet, geoList)
        
        jaccClusterInfo = annotateModulesByCC(jaccNet, geoList)
        writeModuleNetworkAnnot(jaccClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGeo.pruned.0.33.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGeo.pruned.0.33.txt")
        
        
        geoMergedList = normUniqCountryList
        geoMergedClusterScores = getClusterScores(jaccNet, geoMergedList)
        
        jaccMergedClusterInfo = annotateModulesByCC(jaccNet, geoMergedList)
        writeModuleNetworkAnnot(jaccMergedClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGeoMerged.pruned.0.33.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGeoMerged.pruned.0.33.txt")
        
      }
      
    }
    
    geoGenusJaccNetMode = T
    if (geoGenusJaccNetMode) {
      genusMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//genusMatUpdatedJaccard.RData"
      load(genusMatUpdatedJaccRData)
      genusMatUpdatedJaccSelected = genusMatUpdatedJacc[normUniqAll, normUniqAll]
      
      jacc05Mode = T
      if (jacc05Mode) {
        
        genusCut = 0.33
        jaccNet = makeJaccNet(genusMatUpdatedJaccSelected, cutJacc = genusCut)
        
        geoList = normUniqSubjectList
        geoClusterScores = getClusterScores(jaccNet, geoList)
        
        jaccClusterInfo = annotateModulesByCC(jaccNet, geoList)
        writeModuleNetworkAnnot(jaccClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGenusGeo.pruned.0.33.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGenusGeo.pruned.0.33.txt")
        
        
        geoMergedList = normUniqCountryList
        geoMergedClusterScores = getClusterScores(jaccNet, geoMergedList)
        
        jaccMergedClusterInfo = annotateModulesByCC(jaccNet, geoMergedList)
        writeModuleNetworkAnnot(jaccMergedClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGenusGeoMerged.pruned.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGenusGeoMerged.pruned.txt")
        
        
        
      }
      
    }
    
    diseaseJaccNetMode = T
    if (diseaseJaccNetMode) {
      mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
      load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
      mergeMatUpdatedJaccSelected = mergeMatUpdatedJacc[diseaseUniqAll, diseaseUniqAll]
      
      jacc05Mode = T
      if (jacc05Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.5)
        
        diseaseList = diseaseUniqSubjectList
        diseaseClusterScores = getClusterScores(jaccNet, diseaseList)
        
        jaccDiseaseClusterInfo = annotateModulesByCC(jaccNet, diseaseList)
        writeModuleNetworkAnnot(jaccDiseaseClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccDisease.pruned2.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccDisease.pruned2.txt")
        
      }
      jacc033Mode = F
      if (jacc033Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.33)
        
        diseaseList = diseaseUniqSubjectList
        diseaseClusterScores = getClusterScores(jaccNet, diseaseList)
        
        jaccDiseaseClusterInfo = annotateModulesByCC(jaccNet, diseaseList)
        writeModuleNetworkAnnot(jaccDiseaseClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccDisease.pruned.0.33.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccDisease.pruned.0.33.txt")
        
      }
      
      
    }
    
    mixedJaccNetMode = F
    if (mixedJaccNetMode) {
      
      mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
      load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
      mergeMatUpdatedJaccSelected = mergeMatUpdatedJacc[mixedUniqAll, mixedUniqAll]
      
      jacc05Mode = T
      if (jacc05Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.5)
        
        mixedList = mixedUniqSubjectList
        mixedClusterScores = getClusterScores(jaccNet, mixedList)
        
        jaccDiseaseClusterInfo = annotateModulesByCC(jaccNet, mixedList)
        writeModuleNetworkAnnot(jaccDiseaseClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccMixed.pruned.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccMixed.pruned.txt")
        
      }
      
      jacc033Mode = T
      if (jacc033Mode) {
        jaccNet = makeJaccNet(mergeMatUpdatedJaccSelected, cutJacc = 0.4)
        
        mixedList = mixedUniqSubjectList
        mixedClusterScores = getClusterScores(jaccNet, mixedList)
        
        jaccDiseaseClusterInfo = annotateModulesByCC(jaccNet, mixedList)
        writeModuleNetworkAnnot(jaccDiseaseClusterInfo, 
                                "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccMixed.pruned.0.4.txt",
                                "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccMixed.pruned.0.4.txt")
        
        
      }
      
    }
  }
  
  
}

findFunctionCluster = F
if (findFunctionCluster) {
  
  loadFuncAnnot = T
  if (loadFuncAnnot) {
    
    gemFile = "J://Deposit/Project/2018_microbiome_atlas//atlas.GEM.model.related/GEM.from.Reza.amazing/absentPresentInModels.txt"
    gemRData = "J://Deposit/Project/2018_microbiome_atlas//atlas.GEM.model.related/GEM.from.Reza.amazing/absentPresentInModels.RData"
    
    if (F) {
      gemTab = read.delim(gemFile, sep="\t", stringsAsFactors = F, header = T)
      gemMat = as.matrix(gemTab[,-1])
      rownames(gemMat) = gemTab[,1]
      gemMat = t(gemMat)
      save(gemMat, file=gemRData)
    }
    
    antismashMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/msp1992.igc2.antismashMat.RData"
    vfMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/vfBestMat.RData"
    mustardMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/mustard.mat.RData"
    igc2PfamMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PFAM/igc2PfamMat.RData"
    newGutCazyMatRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/CAZY.INRA/igc2.new.cazy.mat.RData"
    gutKoBestMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/gutKoBestMat.1990.RData"
    jgiMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/jgiMat.RData"
    funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
    funcGemMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcGemMat.20190820.RData"
    
    igc2PfamDescRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PFAM/igc2PfamDesc.RData"
    
    patricVfDescTabFile = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/PATRIC/PATRIC_sp_gene.txt"
    patricVfDescTab = read.table(patricVfDescTabFile, header=T, sep="\t", stringsAsFactors = F)
    gut1992MspVfBestDescMapRData = "j://Deposit/Project/2018_microbiome_atlas/igc2/gut1992MspVfBestDescMap.RData"
    
    keggModuleDescMatFile = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/KEGG module/moduleDescTab.RData"
    koModuleListRData = "J://Deposit/Project/2018_microbiome_atlas/functional.annotation/KEGG module/koModuleList.20190219.RData"
    
    koJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//koJaccMat.RData"
    cazyJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//cazyJaccMat.RData"
    pfamJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//pfamJaccMat.RData"
    funcJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//funcJaccMat.20191227.RData"
    funcGemJaccMatRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//funcGemJaccMat.20190820.RData"
    
    load(keggModuleDescMatFile)
    microbeModuleTab = moduleDescTab[moduleDescTab$isBacteriaModule|moduleDescTab$isArchaeaModule,]
    nonMicrobeModuleTab = moduleDescTab[!(moduleDescTab$isBacteriaModule|moduleDescTab$isArchaeaModule),]
    
    microbeModuleNames =  paste(microbeModuleTab$moduleId, microbeModuleTab$moduleName, sep=":")
    nonMicrobeModuleNames =  paste(nonMicrobeModuleTab$moduleId, nonMicrobeModuleTab$moduleName, sep=":")
    
    load(koModuleListRData)
    load(igc2PfamDescRData)
    uniqPfamNames = unique(igc2PfamDesc$pfam_name)
    uniqPfamDesc = unique(igc2PfamDesc$pfam_desc)
    compDescs = c(uniqPfamDesc[grep("competence",uniqPfamDesc)], uniqPfamDesc[grep("Competence",uniqPfamDesc)])
    compPfamNames = unique(igc2PfamDesc$pfam_name[igc2PfamDesc$pfam_desc %in% compDescs])
    
    
    load(gut1992MspVfBestDescMapRData)
    
    load(igc2PfamMatRData) #pfamMat
    load(newGutCazyMatRData) #cazyMat
    load(gutKoBestMatRData) #gutKoBestMat
    
    forNeeluMode = T
    if (forNeeluMode) {
      pfamMat = pfamMat[!rownames(pfamMat) %in% suppressMsps2,]
      pfamMelt = melt(pfamMat)
      pfamMelt = pfamMelt[pfamMelt$value==1,]
      pfamMelt = pfamMelt[,1:2]
      pfamMelt$Var1  = as.character(pfamMelt$Var1)
      pfamMelt$Var2  = as.character(pfamMelt$Var2)
      pfamMelt$species = getSpeciesName(pfamMelt$Var1, taxo)
      gutPfamMeltFile="C:\\Data\\comparative.analysis.healthy.sweden\\gutPfamMelt.for.neelu.txt"
      write.table(pfamMelt, gutPfamMeltFile, sep="\t", quote=F, row.names = F)
      
      igc2PfamDesc$species = getSpeciesName(igc2PfamDesc$msp, taxo)
      gutPfamDescFile = "C:\\Data\\comparative.analysis.healthy.sweden\\gutPfamDesc.for.neelu.txt"
      write.table(igc2PfamDesc, gutPfamDescFile, sep="\t", quote=F, row.names = F)
      
    }
    
    if (F) {
      
      cazyMelt = melt(cazyMat)
      cazyMelt = cazyMelt[cazyMelt$value!=0,]
      cazyMelt$species = getSpeciesName(as.character(cazyMelt$Var1), taxo)
      colnames(cazyMelt)[1:2] = c("MSP_ID","CAZY")
      cazyMelt = cazyMelt[!cazyMelt$MSP_ID %in% suppressMsps2, ]
      write.table(cazyMelt[,-3], 
                  "C://Data//comparative.analysis.healthy.sweden//gut.cazy.tab.for.neelu.txt",
                  row.names = F, quote = F,
                  sep="\t")
    }
    
    load(antismashMatRData) #antismashMat
    load(vfMatRData) #vfMatDigit
    load(mustardMatRData) #mustardMat --> antibiotic resistance
    load(jgiMatRData) #jgiMat
    
    checkPathMode = F
    if (checkPathMode) {
      load(koJaccMatRData)
      allKos = rownames(koJaccMat)
      
      quantile(koJaccMat[lower.tri(koJaccMat)], 0.5) #0
      quantile(koJaccMat[lower.tri(koJaccMat)], 0.90) #0.1111111
      quantile(koJaccMat[lower.tri(koJaccMat)], 0.95) #0.2377622
      quantile(koJaccMat[lower.tri(koJaccMat)], 0.99) #0.5909091 
      
      moduleLevelMode = T
      if (moduleLevelMode) {
        
        countKoModuleList <- do.call(c,lapply(koModuleList, function(kos){
          kos = kos[kos %in% allKos]
          return(length(kos))
        }))
        
        fractionKoModuleList <- do.call(c,lapply(koModuleList, function(kos){
          total = length(kos)
          kos = kos[kos %in% allKos]
          remained = length(kos)
          return(remained/total)
        }))
        
        koModuleCatalogStat = data.frame(total=do.call(c,lapply(koModuleList,length)),
                                         number=countKoModuleList, 
                                         fraction=fractionKoModuleList)
        
        nonSmallModules = rownames(koModuleCatalogStat[koModuleCatalogStat$total>=5,])
        
        coveredKoModules = rownames(koModuleCatalogStat[koModuleCatalogStat$total>=5&
                                                          koModuleCatalogStat$fraction>=0.5,])
        
        microbeModuleCoveredNames = microbeModuleNames[microbeModuleNames %in% nonSmallModules]
        nonMicrobeModuleCoveredNames = nonMicrobeModuleNames[nonMicrobeModuleNames %in% nonSmallModules]
        
        jaccKoModuleList <- lapply(koModuleList, function(kos){
          kos = kos[kos %in% allKos]
          if (length(kos)==0) return(NA)
          if (length(kos)==1) return(NA)
          currJaccMat = koJaccMat[kos,kos]
          return(currJaccMat[lower.tri(currJaccMat)])
        })
        
        jaccMeanKoModules = do.call(c,lapply(jaccKoModuleList, mean))
        koBoxList = list(microbe= jaccMeanKoModules[names(jaccMeanKoModules)%in% microbeModuleNames],
                         nonMicrobe= jaccMeanKoModules[names(jaccMeanKoModules) %in% nonMicrobeModuleNames])
        boxplot(rev(koBoxList),horizontal = T,las=2)
        
        koBoxList = list(microbe= jaccMeanKoModules[names(jaccMeanKoModules)%in% microbeModuleCoveredNames],
                         nonMicrobe= jaccMeanKoModules[names(jaccMeanKoModules) %in% nonMicrobeModuleCoveredNames])
        boxplot(rev(koBoxList),horizontal = T,las=2)
        
        jaccKoStats = 3
        
      }
      
      pathLevelMode = T
      if (pathLevelMode) {
        
      }
      
    }
    
    generateMode = F
    if (generateMode) {
      koJaccMat =  getJaccMat(gutKoBestMat)
      save(koJaccMat, file = koJaccMatRData)
      
      cazyJaccMat =  getJaccMat(cazyMat)
      save(cazyJaccMat, file = cazyJaccMatRData)
      
      pfamJaccMat =  getJaccMat(pfamMat)
      pfamJaccMat = pfamJaccMat[!rownames(pfamJaccMat) %in% suppressMsps,]
      save(pfamJaccMat, file = pfamJaccMatRData)
      
      colnames(cazyMat) = paste("cazy", colnames(cazyMat), sep=".")
      
      funcMat = mustardMat
      funcMat = mergeMatrix(funcMat, cazyMat)
      funcMat = mergeMatrix(funcMat, antismashMat)
      funcMat = mergeMatrix(funcMat, vfMatDigit)
      funcMat = mergeMatrix(funcMat, gutKoBestMat)
      funcMat = mergeMatrix(funcMat, pfamMat)
      funcMat = mergeMatrix(funcMat, jgiMat)
      funcMat = funcMat[!rownames(funcMat) %in% suppressMsps,]
      
      save(funcMat, file=funcMatRData)
      
      
      gemSpecies = rownames(gemMat)
      table(gemSpecies %in% rownames(funcMat))
      
      funcGemMat = funcMat[rownames(funcMat) %in% gemSpecies,]
      funcGemMat = mergeMatrix(funcGemMat, gemMat)
      save(funcGemMat, file=funcGemMatRData)
      
      funcJaccMat =  getJaccMat(funcMat)
      save(funcJaccMat, file = funcJaccMatRData)
      
    }
    
  }
  
  loadFuncJaccData = F
  if (loadFuncJaccData) {
    load(funcJaccMatRData)
    funcJaccTabSelRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//funcJaccTabSel.20190815.RData"
    load(funcJaccTabSelRData)
    View(funcJaccTabSel[is.na(funcJaccTabSel$Cat2),])
  }
  
  checkFuncJaccMat = F
  if (checkFuncJaccMat) {
    load(funcJaccMatRData)
    funcJaccNet = makeJaccNet(funcJaccMat, cutJacc = 0.9)
    funcJaccModList = makeModuleList(funcJaccNet)
    #annotateModulesByCC
    writeModuleNetworkAnnot(jaccClusterInfo, 
                            "C://Data//comparative.analysis.healthy.sweden//nodeTableJaccGeo.pruned.txt",
                            "C://Data//comparative.analysis.healthy.sweden//edgeTableJaccGeo.pruned.txt")
    
    
    
    makeTabMode = F
    if (makeTabMode) {
      require(reshape2)
      funcJaccTab = melt(funcJaccMat)
      jaccCut = 0.5
      funcJaccTabSel = funcJaccTab[funcJaccTab$value>=jaccCut,]
      funcJaccTabSel = funcJaccTabSel[funcJaccTabSel$Var1 != funcJaccTabSel$Var2,]
      funcJaccTabSel$Cat1 = NA
      funcJaccTabSel$Cat2 = NA
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(pfamMat)] = "Pfam"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(pfamMat)] = "Pfam"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(cazyMat)] = "Cazy"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(cazyMat)] = "Cazy"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(gutKoBestMat)] = "KO"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(gutKoBestMat)] = "KO"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(antismashMat)] = "antiSMASH"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(antismashMat)] = "antiSMASH"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(vfMatDigit)] = "virulence"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(vfMatDigit)] = "virulence"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(mustardMat)] = "AR"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(mustardMat)] = "AR"
      funcJaccTabSel$Cat1[funcJaccTabSel$Var1 %in% colnames(jgiMat)] = "Phenotype"
      funcJaccTabSel$Cat2[funcJaccTabSel$Var2 %in% colnames(jgiMat)] = "Phenotype"
      funcJaccTabSel = funcJaccTabSel[!is.na(funcJaccTabSel$value),]
      
      funcJaccTabSelRData = "J://Deposit/Project/2018_microbiome_atlas//functional.annotation//funcJaccTabSel.20190815.RData"
      save(funcJaccTabSel, file=funcJaccTabSelRData)
      
    }
    
  }
  
}

makeFunctionCluster = F ### performed on UPPMAX rackham cluster
if (makeFunctionCluster) {
  
  source("~/rscripts/2019_funcType/funcType.analysis.functions.r")
  source("~/rscripts/2019_funcType/funcType.data.files.r")
  
  funcGemJaccNetMode = T
  if (funcGemJaccNetMode) {
    funcGemModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.nodeTable.20190828.txt"
    funcGemModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.edgeTable.20190828.txt"
    
    load(funcGemJaccMatRData)
    funcGemJaccNet = makeJaccNet(funcGemJaccMat)
    save(funcGemJaccNet, file = funcGemJaccNetRData)
    
    funcGemModules = makeModuleList(funcGemJaccNet)
    save(funcGemModules, file = funcGemModulesRData)
    
    funcGemModuleAnnots = annotateModulesByCC(funcGemJaccNet, funcGemModules, 3)
    writeModuleNetworkAnnot(funcGemModuleAnnots, funcGemModNetNodeTableFile, funcGemModNetEdgeTableFile)
    save(funcGemModuleAnnots, file=funcGemModuleAnnotsRData)
    
  }
  
  rm(list=ls())
  
  source("~/rscripts/2019_funcType/funcType.analysis.functions.r")
  source("~/rscripts/2019_funcType/funcType.data.files.r")
  
  funcJaccNetMode = T
  if (funcJaccNetMode) {
    funcModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.nodeTable.20190828.txt"
    funcModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.edgeTable.20190828.txt"
    
    load(funcJaccMatRData)
    funcJaccNet = makeJaccNet(funcJaccMat)
    save(funcJaccNet, file = funcJaccNetRData)
    
    funcModules = makeModuleList(funcJaccNet)
    save(funcModules, file = funcModulesRData)
    
    funcModuleAnnots = annotateModulesByCC(funcJaccNet, funcModules, 3)
    writeModuleNetworkAnnot(funcModuleAnnots, funcModNetNodeTableFile, funcModNetEdgeTableFile)
    save(funcModuleAnnots, file=funcModuleAnnotsRData)
    
  }
  
}

mainFig1PiePlotMode = F
if (mainFig1PiePlotMode) {
  
  diseaseMgsMat = mergeMatUpdated[,diseaseUniqAll] 
  normalMgsMat = mergeMatUpdated[,normUniqAll]
  industrialMgsMat = mergeMatUpdated[,industrializedSamples]
  traditionalMgsMat = mergeMatUpdated[,traditionalSamples]
  
  familyMode = T
  if (familyMode) {
    normalFamilyRelAbd = getRelAbdOfTaxa(industrialMgsMat, "family", taxo)
    diseaseFamilyRelAbd = getRelAbdOfTaxa(diseaseMgsMat, "family", taxo)
    indiFamilyRelAbd = getRelAbdOfTaxa(traditionalMgsMat, "family", taxo)
    
    normalFamilySortedNames = names(sort(normalFamilyRelAbd,decreasing = T))
    normalFamilySortedNames = normalFamilySortedNames[!grepl("unclassified",normalFamilySortedNames)]
    diseaseFamilySortedNames = names(sort(diseaseFamilyRelAbd,decreasing = T))
    diseaseFamilySortedNames = diseaseFamilySortedNames[!grepl("unclassified",diseaseFamilySortedNames)]
    indiFamilySortedNames = names(sort(indiFamilyRelAbd,decreasing = T))
    indiFamilySortedNames = indiFamilySortedNames[!grepl("unclassified",indiFamilySortedNames)]
    
    normalTop20 = normalFamilySortedNames[1:20]
    diseaseTop20 = diseaseFamilySortedNames[1:20]
    indiTop20 = indiFamilySortedNames[1:20]
    
    mergedTop20 = unique(c(normalTop20, diseaseTop20, indiTop20))
    
    normalFamilyTop20Imputed = getRelAbdWithOthersImputed(normalFamilyRelAbd, mergedTop20)
    diseaseFamilyTop20Imputed = getRelAbdWithOthersImputed(diseaseFamilyRelAbd, mergedTop20)
    indiFamilyTop20Imputed = getRelAbdWithOthersImputed(indiFamilyRelAbd, mergedTop20)
    
    mergedList = names(sort(normalFamilyTop20Imputed,decreasing = T))
    lenMerged = length(mergedList)
    
    normalFamilyTop20Imputed = normalFamilyTop20Imputed[match(mergedList, names(normalFamilyTop20Imputed))]
    diseaseFamilyTop20Imputed = diseaseFamilyTop20Imputed[match(mergedList, names(diseaseFamilyTop20Imputed))]
    indiFamilyTop20Imputed = indiFamilyTop20Imputed[match(mergedList, names(indiFamilyTop20Imputed))]
    
    normalFamilyTop20Imputed = c(normalFamilyTop20Imputed[names(normalFamilyTop20Imputed)!="others"],
                                 normalFamilyTop20Imputed[names(normalFamilyTop20Imputed)=="others"])
    diseaseFamilyTop20Imputed = c(diseaseFamilyTop20Imputed[names(diseaseFamilyTop20Imputed)!="others"],
                                  diseaseFamilyTop20Imputed[names(diseaseFamilyTop20Imputed)=="others"])
    indiFamilyTop20Imputed = c(indiFamilyTop20Imputed[names(indiFamilyTop20Imputed)!="others"],
                               indiFamilyTop20Imputed[names(indiFamilyTop20Imputed)=="others"])
    
    piePlotMode = T
    if (piePlotMode) {
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      familyCols = c(col_vector[1:(lenMerged-1)],"lightgray")
      familyCols = col_vector[1:(lenMerged)] 
      
      diseaseTarget = diseaseFamilyTop20Imputed
      normalTarget = normalFamilyTop20Imputed
      indiTarget = indiFamilyTop20Imputed
      
      iniR = 0.2
      pie(1, radius=0.01, init.angle=90, col=c('white'), border = NA, labels='')
      floating.pie(0,0,diseaseTarget ,
                   radius=5*iniR, 
                   startpos=pi/2, 
                   col = familyCols,
                   border="black")
      
      floating.pie(0,0,normalTarget ,  
                   radius=4*iniR, 
                   startpos=pi/2, 
                   col = familyCols,
                   border="black")
      
      floating.pie(0,0,indiTarget,  
                   radius=3*iniR, 
                   startpos=pi/2, 
                   col = familyCols,
                   border="black")
      floating.pie(0,0,1,  
                   radius=2*iniR, 
                   startpos=pi/2, 
                   col="white",
                   border="black")
      floating.pie(0,0,1,  
                   radius=1.99*iniR, 
                   startpos=pi/2, 
                   col="white",
                   border=NA)
      
      
    }
  }
}

pcoaCheckMode = F
if (pcoaCheckMode) {
  
  pcoaEnteroCheck = T
  if (pcoaEnteroCheck) {
    
  }
  
  simpleCheck = F
  if (simpleCheck) {
    pcoaGenusMat= pcoaGenusOut$vectors[,1:2]
    
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    barplot(c(1,2),col = c("#00bfff","#ff2400"))
    
    pcoaCols[normUniqCountryList$Thailand] = "blue"
    pcoaCols[thaiShort] = "#00bfff"
    pcoaCols[thaiLong] = "#ff2400"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[normUniqCountryList$Thailand] = "blue"
    pcoaCols[thaiShort] = "#00bfff"
    pcoaCols[thaiLong] = "#ff2400"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    beforeAntibiotics=basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51" & basicMetaMapUpdated$subtype=="Dag0"]
    afterAntibiotics=basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51" & basicMetaMapUpdated$subtype=="Dag4"]
    afterAntibiotics2=basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51" & basicMetaMapUpdated$subtype=="Dag8"]
    afterAntibiotics3=basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id51" & basicMetaMapUpdated$subtype=="Dag180"]
    
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[afterAntibiotics] = "blue"
    pcoaCols[beforeAntibiotics] = "#00bfff"
    pcoaCols[afterAntibiotics2] = "purple"
    pcoaCols[afterAntibiotics3] = "black"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    poisonedSamples =  basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id48" ]
    
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[poisonedSamples] = "blue"
    pcoaCols[normUniqCountryList$US] = "#00bfff"
    
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    
    beforeMetformin= basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" & basicMetaMapUpdated$subtype=="M0"]
    afterMetformin= basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id53" & basicMetaMapUpdated$subtype=="M4"]
    
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[afterMetformin] = "blue"
    pcoaCols[beforeMetformin] = "#00bfff"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    renalBefore = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id41" & 
                                                  basicMetaMapUpdated$type=="Case" &
                                                  basicMetaMapUpdated$subtype=="before"]
    renalAfter = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id41" & 
                                                 basicMetaMapUpdated$type=="Case" &
                                                 basicMetaMapUpdated$subtype=="after"]
    pcoaCols = rep("#c0c0c033", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[renalAfter] = "blue"
    pcoaCols[renalBefore] = "#00bfff"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
    lungBefore = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id26" & 
                                                 basicMetaMapUpdated$type=="Case" &
                                                 basicMetaMapUpdated$subtype=="Before"]
    lungAfter = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID=="id26" & 
                                                basicMetaMapUpdated$type=="Case" &
                                                basicMetaMapUpdated$subtype=="after"]
    pcoaCols = rep("#c0c0c022", dim(pcoaGenusMat)[1])
    names(pcoaCols) = rownames(pcoaGenusMat)
    
    pcoaCols[lungAfter] = "blue"
    pcoaCols[lungBefore] = "#00bfff"
    plot(pcoaGenusMat[,1], pcoaGenusMat[,2], col=pcoaCols, pch=16, xaxt="n", yaxt="n", xlab="", ylab="") #, xlim=c(-0.55,0.4), ylim=c(-0.3,0.7)
    
  }
  
}

updateResidualRichness = F
if (updateResidualRichness) {
  
  noImputeMetaMapUpdated = basicMetaMapUpdated[!is.na(basicMetaMapUpdated$Age) & !is.na(basicMetaMapUpdated$BMI),]
  fitGeneRichnessUpdated = lmer(GeneRichness ~ Age + BMI + (1|Sequencer), data=noImputeMetaMapUpdated, na.action = na.exclude)
  noImputeMetaMapUpdated$ResidualRichness = noImputeMetaMapUpdated$GeneRichness - predict(fitGeneRichnessUpdated)
  
  basicMetaMapUpdated$ResidualRichness = NA
  basicMetaMapUpdated$ResidualRichness[match(noImputeMetaMapUpdated$sample.ID, basicMetaMapUpdated$sample.ID)] = noImputeMetaMapUpdated$ResidualRichness
  
  basicMetaCleanRDataUpdated2 = "C:\\Data/comparative.analysis.healthy.sweden/all.basic.clean.metadata.behcet.20190805.RData"
  save(basicMetaMapUpdated, file=basicMetaCleanRDataUpdated2)
  
}

checkMgsStat = F
if (checkMgsStat) {
  
  out = apply(mergeMatUpdated, 2, function(currCol){
    names(currCol) = rownames(mergeMatUpdated)
    currCol = sort(currCol, decreasing = T)
    return(names(currCol[1:10]))
  })
  outVec = as.vector(out)
  uniqTop10Species = unique(outVec)
  
  topSpeciesMergeMat = mergeMatUpdated[uniqTop10Species,]
  #corMat = cor(topSpeciesMergeMat)
  avgMixedTopMat = getAvgAbdMatBy(topSpeciesMergeMat, mixedUniqSubjectList)
  corMxiedTopMat = cor(avgMixedTopMat)
  heatmap.2(corMxiedTopMat, trace="none",
            col = colorRampPalette(c("white","red"))(256),
            margins = c(7,7),key=F)
  
  set.seed(1)
  targetMat = avgMixedTopMat
  targetTsneOut = Rtsne(t(targetMat), dims=2, perplexity = 5)
  rownames(targetTsneOut$Y) = colnames(targetMat)
  plot(targetTsneOut$Y)
  
  tsneTab = data.frame(targetTsneOut$Y, type=rownames(targetTsneOut$Y))
  ggplot(tsneTab, aes(x=X1 , y=X2, label=type)) + 
    geom_point(aes(colour = type, alpha=.02)) + 
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text.y = element_blank(),
          legend.position="none", 
          panel.border = element_rect(colour = "black", fill=NA, size=.5))  + 
    geom_text_repel(size=3, force=2) + xlab("") + ylab("")
  
  # ,cexCol = 0.5,cexRow = 0.5,
  # sepcolor = "gray",
  # colsep=seq(0, dim(diseaseVipMat)[2]),
  # rowsep=seq(0, dim(diseaseVipMat)[1]),
  # sepwidth = c(0.01,0.01)
}

checkJaccForCooccurence = F
if (checkJaccForCooccurence) {
  
  etfCountriesUpdated = c("UK","Germany","Italy","Spain","Sweden","Denmark","Finland","Luxembourg","France")
  etpCountriesUpdated = c("Peru","Mongo","Thai","Fiji","Madagascar","Tanzania","India") #,"Fiji","India"
  etbCountriesUpdated = c("China","Japan","US")
  
  etfSamples = do.call(c,normalSamplesByCountries[etfCountriesUpdated])
  etpSamples = do.call(c,normalSamplesByCountries[etpCountriesUpdated])
  etbSamples = do.call(c,normalSamplesByCountries[etbCountriesUpdated])
  
  etfMat = mergeMatUpdated[,etfSamples]
  etpMat = mergeMatUpdated[,etpSamples]
  etbMat = mergeMatUpdated[,etbSamples]
  etfMat = etfMat[rowMeans(etfMat)!=0,]
  etpMat = etpMat[rowMeans(etpMat)!=0,]
  etbMat = etbMat[rowMeans(etbMat)!=0,]
  
  industrialMat = mergeMatUpdated[,newWesternIds]
  industrialMat = industrialMat[rowMeans(industrialMat)!=0,]
  traditionalMat = mergeMatUpdated[,newNonWesternIds]
  traditionalMat = traditionalMat[rowMeans(traditionalMat)!=0,]
  diseaseMat = mergeMatUpdated[,newDiseaseIds]
  diseaseMat = diseaseMat[rowMeans(diseaseMat)!=0,]
  
  mgsJaccMat = 1-jaccard(t(mergeMatUpdated>0))
  mgsJaccIndustrialMat = 1-jaccard(t(industrialMat>0))
  mgsJaccTraditionalMat = 1-jaccard(t(traditionalMat>0))
  mgsJaccDiseaseMat = 1-jaccard(t(diseaseMat>0))
  mgsJaccEtFMat = 1-jaccard(t(etfMat>0))
  mgsJaccEtBMat = 1-jaccard(t(etbMat>0))
  mgsJaccEtPMat = 1-jaccard(t(etpMat>0))
  
  mgsJaccNet = makeJaccNet(mgsJaccMat, cutJacc = 0.75)
  mgsJaccIndustrialNet = makeJaccNet(mgsJaccIndustrialMat, cutJacc = 0.75)
  mgsJaccTraditionalNet = makeJaccNet(mgsJaccTraditionalMat, cutJacc = 0.75)
  mgsJaccDiseaseNet = makeJaccNet(mgsJaccDiseaseMat, cutJacc = 0.75)
  
  
  #### all samples ####
  mgsJaccModule = makeModuleList(mgsJaccNet)
  mgsJaccAnnot = annotateModulesByCC(mgsJaccNet, mgsJaccModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccAnnot, 
                          "mgs.jacc.0.75.nodeTable.txt",
                          "mgs.jacc.0.75.edgeTable.txt")
  
  #### industrial samples ####
  mgsJaccIndustrialModule = makeModuleList(mgsJaccIndustrialNet)
  mgsJaccIndustrialAnnot = annotateModulesByCC(mgsJaccIndustrialNet, mgsJaccIndustrialModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccIndustrialAnnot, 
                          "mgs.jacc.0.75.industrial.nodeTable.txt",
                          "mgs.jacc.0.75.industrial.edgeTable.txt")
  
  #### traditional samples ####
  mgsJaccTraditionalModule = makeModuleList(mgsJaccTraditionalNet)
  mgsJaccTraditionalAnnot = annotateModulesByCC(mgsJaccTraditionalNet, mgsJaccTraditionalModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccTraditionalAnnot, 
                          "mgs.jacc.0.75.traditional.nodeTable.txt",
                          "mgs.jacc.0.75.traditional.edgeTable.txt")
  
  #### diseased samples ####
  mgsJaccDiseaseModule = makeModuleList(mgsJaccDiseaseNet)
  mgsJaccDiseaseAnnot = annotateModulesByCC(mgsJaccDiseaseNet, mgsJaccDiseaseModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccDiseaseAnnot, 
                          "mgs.jacc.0.75.disease.nodeTable.txt",
                          "mgs.jacc.0.75.disease.edgeTable.txt")
  
  #### ET-firmiuctes samples ####
  mgsJaccEtFNet = makeJaccNet(mgsJaccEtFMat, cutJacc = 0.5)
  mgsJaccEtFModule = makeModuleList(mgsJaccEtFNet)
  mgsJaccEtFAnnot = annotateModulesByCC(mgsJaccEtFNet, mgsJaccEtFModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccEtFAnnot, 
                          "mgs.jacc.0.5.EtF.nodeTable.txt",
                          "mgs.jacc.0.5.EtF.edgeTable.txt")
  
  #### ET-bacteroides samples ####
  mgsJaccEtBNet = makeJaccNet(mgsJaccEtBMat, cutJacc = 0.5)
  mgsJaccEtBModule = makeModuleList(mgsJaccEtBNet)
  mgsJaccEtBAnnot = annotateModulesByCC(mgsJaccEtBNet, mgsJaccEtBModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccEtBAnnot, 
                          "mgs.jacc.0.5.EtB.nodeTable.txt",
                          "mgs.jacc.0.5.EtB.edgeTable.txt")
  
  #### ET-prevotella samples ####
  mgsJaccEtPNet = makeJaccNet(mgsJaccEtPMat, cutJacc = 0.5)
  mgsJaccEtPModule = makeModuleList(mgsJaccEtPNet)
  mgsJaccEtPAnnot = annotateModulesByCC(mgsJaccEtPNet, mgsJaccEtPModule, cutCluster = 3)
  
  writeModuleNetworkAnnot(mgsJaccEtPAnnot, 
                          "mgs.jacc.0.5.EtP.nodeTable.txt",
                          "mgs.jacc.0.5.EtP.edgeTable.txt")
  
  #### check ####
  targetList = mgsJaccModule$`2`
  targetList = mgsJaccIndustrialModule$`9`
  targetList = mgsJaccTraditionalModule$`29`
  sort(taxo$species[taxo$MSP %in% targetList])
  
  lapply(mgsJaccIndustrialModule[mgsJaccIndustrialAnnot$nodeTable$node], function(x) sort(taxo$species[taxo$MSP %in% x]) )
  lapply(mgsJaccTraditionalModule[mgsJaccTraditionalAnnot$nodeTable$node], function(x) sort(taxo$species[taxo$MSP %in% x]) )
  lapply(mgsJaccDiseaseModule[mgsJaccDiseaseAnnot$nodeTable$node], function(x) sort(taxo$species[taxo$MSP %in% x]) )
  
  table(mgsJaccDiseaseModule$`1` %in% mgsJaccIndustrialModule$`1`)
  table(mgsJaccDiseaseModule$`1` %in% mgsJaccIndustrialModule$`2`)
  
  statCheck = F
  if (statCheck) {
    mm = mgsJaccMat[lower.tri(mgsJaccMat)]
    mgsJaccMaxs = sapply(1:1977, function(ind){
      currVec = mgsJaccMat[ind,]
      currVec = currVec[-ind]
      return(max(currVec))
    })
    plot(density(mgsJaccMaxs))
    median(mgsJaccMaxs)
    quantile(mgsJaccMaxs,0.999)
    quantile(mgsJaccMaxs,0.001)
    
    quantile(mgsJaccMaxs,0.99)
    quantile(mgsJaccMaxs,0.01)
    
    plot(density(mm[mm!=0]), log="x")
  }
  
}

loadCoSymptomNetwork = F
if (loadCoSymptomNetwork) {
  coSymptomNetFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.disease.network.data\\co-symptom-network.txt"
  coSymtpomTab = read.delim(coSymptomNetFile, stringsAsFactors = F, sep="\t")
  
  targetDiseaseFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.disease.network.data\\co-symptom-network-diseaseNameMapping.txt"
  targetDiseaseTab = read.delim(targetDiseaseFile, stringsAsFactors = F, sep="\t")
  
  cut05Mode = F
  if (cut05Mode) {
    coSymtpomTabSel = coSymtpomTab[coSymtpomTab$MeSH.Disease.Term %in% targetDiseaseTab$Diseases.in.symptom.disease.network &
                                     coSymtpomTab$MeSH.Disease.Term.1 %in% targetDiseaseTab$Diseases.in.symptom.disease.network &
                                     coSymtpomTab$symptom.similarity.score >= 0.5, ] 
    
    dim(coSymtpomTabSel)  
    
    coSymptomNetSelFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.disease.network.data\\co-symptom-network.selected.atlas.txt"
    write.table(coSymtpomTabSel, file=coSymptomNetSelFile, sep="\t", quote = F, row.names = F)
    
  }
  
  cut033Mode = T
  if (cut033Mode) {
    coSymtpomTabSel = coSymtpomTab[coSymtpomTab$MeSH.Disease.Term %in% targetDiseaseTab$Diseases.in.symptom.disease.network &
                                     coSymtpomTab$MeSH.Disease.Term.1 %in% targetDiseaseTab$Diseases.in.symptom.disease.network &
                                     coSymtpomTab$symptom.similarity.score >= 1/3, ] 
    
    dim(coSymtpomTabSel)  
    
    coSymptomNetSelFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.disease.network.data\\co-symptom-network.selected.atlas.033.txt"
    write.table(coSymtpomTabSel, file=coSymptomNetSelFile, sep="\t", quote = F, row.names = F)
    
  }
}

loadCoMorbidGeneticNetwork = F
if (loadCoMorbidGeneticNetwork) {
  coNetworkFile = "J:\\Deposit\\Project\\2018_microbiome_atlas\\atlas.disease.network.data\\co-host-genetic-network.txt"
  coNetworkTab = read.delim(coNetworkFile, sep="\t", header=T, stringsAsFactors = F)
  coNetworkNameMapFile = "J://Deposit//Project//2018_microbiome_atlas//atlas.disease.network.data/co-host-genetic-network-diseaseNameMapping.txt"
  coNetworkNameTab = read.delim(coNetworkNameMapFile, sep="\t", header=T, stringsAsFactors = F)
  
  coNetworkTab$X..ICD1. = gsub("\\[","",coNetworkTab$X..ICD1.)
  coNetworkTab$X..ICD1. = gsub("\\]","",coNetworkTab$X..ICD1.)
  coNetworkTab$X.Name1. = gsub("\\[","",coNetworkTab$X.Name1.)
  coNetworkTab$X.Name1. = gsub("\\]","",coNetworkTab$X.Name1.)
  coNetworkTab$X.ICD2.  = gsub("\\[","",coNetworkTab$X.ICD2.)
  coNetworkTab$X.ICD2.  = gsub("\\]","",coNetworkTab$X.ICD2.)
  coNetworkTab$X.Name2. = gsub("\\[","",coNetworkTab$X.Name2.)
  coNetworkTab$X.Name2. = gsub("\\]","",coNetworkTab$X.Name2.)
  
  targetDiseaseNames = coNetworkNameTab$Diseases.in.host.genetic.disease.network
  
  table(coNetworkTab$X.Name1. %in% coNetworkNameTab$Diseases.in.host.genetic.disease.network)
  coNetworkNameTab$Diseases.in.host.genetic.disease.network[!coNetworkNameTab$Diseases.in.host.genetic.disease.network  %in%  coNetworkNameTab$X.Name1.]
  
  coNetworkSel =  coNetworkTab[coNetworkTab$X.Name1. %in% targetDiseaseNames
                               & coNetworkTab$X.Name2. %in% targetDiseaseNames,]
  
  coMorbidNetSel  = coNetworkSel[coNetworkSel$RHO>0,]
  coGeneticNetSel = coNetworkSel[coNetworkSel$NG>0,]
  
}

### pouyan data ###
compModeForPouyan = F
if (compModeForPouyan) {
  dsIds = c("id9", "id1","id10","id34",
            "id11","id44","id46","id45",
            "id6","id12","id17","id20",
            "id16","id25","id8","id31",
            "id21","id35","id2","id14", 
            "id52")
  
  statOut = list()
  for (dsId in dsIds) {
    
    currCaseIds = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID==dsId & basicMetaMapUpdated$type == "Case"]
    currContIds = basicMetaMapUpdated$sample.ID[basicMetaMapUpdated$dataset.ID==dsId & basicMetaMapUpdated$type == "Control"]
    currSamples = c(currCaseIds, currContIds)
    print(table(currCaseIds %in% currContIds))
    print(paste(dsId, length(currCaseIds), length(currContIds)))
    
    currCaseMat = mergeMatUpdated[,colnames(mergeMatUpdated)%in%currCaseIds]
    currContMat = mergeMatUpdated[,colnames(mergeMatUpdated)%in%currContIds]
    
    currDataMat = cbind(currCaseMat, currContMat)
    currDataLab = c(rep("case",dim(currCaseMat)[2]),
                    rep("cont",dim(currContMat)[2]))
    
    currCaseMeans = rowMeans(currCaseMat)
    currContMeans = rowMeans(currContMat)
    currFoldChange = currCaseMeans / currContMeans
    
    currStats = testRelations(currDataMat, currDataLab, "wilcoxon")
    currStats$dataset = dsId
    currStats$msp = rownames(currStats)
    currStats$foldChange = currFoldChange
    currStats$Case = currCaseMeans
    currStats$Cont = currContMeans
    
    
    currStats = currStats[,3:10]
    currStats$p[is.nan(currStats$p)]=1
    currStats$q[is.nan(currStats$q)]=1
    
    statOut[[dsId]] = currStats
    print(head(currStats))
    
  }
  statFinal = do.call(rbind, statOut)
  
  gssOssTabNewFile = "J://Deposit//Project//2018_oral_catalog_paper//catalogs (gut, oral), 20190523//MSP_set_class_tropism_20190523_short.txt"
  gssOssTabNew = read.table(gssOssTabNewFile, header=T, sep="\t", stringsAsFactors = F)
  rownames(gssOssTabNew) = gssOssTabNew$MSP
  
  statFinal$species = gssOssTabNew$species[match(statFinal$msp, gssOssTabNew$MSP)]
  
  write.table(statFinal, 
              #"C://Data//comparative.analysis.healthy.sweden//atlas.case.cont.comparison.pouyan.20190626.txt",
              "C://Data//comparative.analysis.healthy.sweden//atlas.case.cont.comparison.pouyan.20190827.txt",
              sep="\t")
  
  
  
}

atlasCompModeNew = F
if (atlasCompModeNew) {
  
  normalMat = mergeMatUpdated[,colnames(mergeMatUpdated) %in% newNormalIds]
  
  targetMat = mergeMatUpdated
  phylaMat = getTaxaSumMat(targetMat,taxo,"phylum",T)
  classMat = getTaxaSumMat(targetMat,taxo,"class",T)
  orderMat = getTaxaSumMat(targetMat,taxo,"order",T)
  familyMat = getTaxaSumMat(targetMat,taxo,"family",T)
  genusMat = getTaxaSumMat(targetMat,taxo,"genus",T)
  speciesMat = getTaxaSumMat(targetMat,taxo,"species",T)
  
  
  checkEnrichedByGeo <- function(newNormalAllList, mgsMat) {
    
    numDatasets = length(newNormalAllList)
    nameDatasets = names(newNormalAllList)
    
    
    
    for (targetInd in seq_along(newNormalAllList)) {
      targetSamples = newNormalAllList[[targetInd]]
      targetName = nameDatasets[targetInd]
      
      targetMat = mgsMat[, colnames(mgsMat) %in% targetSamples]
      targetMeans = rowMeans(targetMat)
      
      foldChangeList = list()
      effectSizeList = list()
      
      for (refInd in seq_along(newNormalAllList)) {
        refSamples = newNormalAllList[[refInd]]
        refName = nameDatasets[refInd]
        
        refMat = mgsMat[, colnames(mgsMat) %in% refSamples]
        refMeans = rowMeans(refMat)
        currFoldChange = refMeans/targetMeans
        currEffectSize = getEffectSizeMat(targetMat, refMat)
        
        foldChangeList[[refName]] = currFoldChange
        effectSizeList[[refName]] = currEffectSize
        
        
      }
      
      ##effectSizeMat = do.call(rbind, effectSizeList)
      ##effectSizeMatDigit = effectSizeMat>
      
    }
    
    
  } 
  
  oldMode = F
  if (oldMode) {
    getSampleListByDisease <- function(basicMetaMap, mgsInfo) {
      
      caseDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="disease"])
      caseDatasets2 = caseDatasets[caseDatasets!="id38"]
      caseIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% caseDatasets]
      caseIDs2 = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% caseDatasets2]
      diseaseOfCaseIDs = sapply(caseIDs, function(currID){
        currDatasetId = basicMetaMap$dataset.ID[basicMetaMap$sample.ID==currID]
        currDisease = unique(mgsInfo$Case[mgsInfo$ID==currDatasetId])
        return(currDisease)
      })
      caseIDsByDisease= split(names(diseaseOfCaseIDs), diseaseOfCaseIDs)
      lenDisease = length(caseIDsByDisease)
      basicMetaMap$disease = ""
      for (ind in 1:lenDisease){
        basicMetaMap$disease[basicMetaMap$sample.ID%in%caseIDsByDisease[[ind]]] = names(caseIDsByDisease)[ind]
      }
      
      
      selectedMetaMap = basicMetaMap[basicMetaMap$sample.ID %in% caseIDs,]
      samplesByDisease = split(selectedMetaMap$sample.ID, selectedMetaMap$disease)
      
      exDisease = c("atherosclerosis")# those samples having less than 5 samples
      samplesByDisease = samplesByDisease[!names(samplesByDisease)%in% exDisease]
      
      return(samplesByDisease)
      
    }
    getSampleListByCountry <- function(basicMetaMap, mgsInfo) {
      mContDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="disease"&!is.na(mgsInfo$Control)])
      hContDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="healthy"&!is.na(mgsInfo$Control)])
      
      healthyDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="healthy"])
      indiDatasets = c("id32","id37","id33","id42")
      notIndiDatasets = healthyDatasets[!healthyDatasets %in% indiDatasets]
      
      healthyIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% healthyDatasets]
      indiIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% indiDatasets]
      NotIndiHealtyIDs = healthyIDs[!healthyIDs %in% indiIDs]
      
      mContIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Control" & basicMetaMap$dataset.ID %in% mContDatasets]
      hContIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Control" & basicMetaMap$dataset.ID %in% hContDatasets]
      mhContIDs = unique(c(mContIDs, hContIDs))
      normalIDs = unique(c(NotIndiHealtyIDs, mhContIDs))
      
      
      allHealthyIds = unique(c(mhContIDs, normalIDs, indiIDs))
      selectedMetaMap = basicMetaMap[basicMetaMap$sample.ID %in% allHealthyIds,]
      
      samplesByCountry = split(selectedMetaMap$sample.ID, selectedMetaMap$Geography)
      exCountries = c("Finland","Italy")# those samples having less than 5 samples
      samplesByCountry = samplesByCountry[!names(samplesByCountry)%in% exCountries]
      return(samplesByCountry)
    }
    getSampleListByIndiCountry <- function(basicMetaMap, mgsInfo) {
      indiDatasets = c("id32","id37","id33","id42")
      indiIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% indiDatasets]
      selectedMetaMap = basicMetaMap[basicMetaMap$sample.ID %in% indiIDs,]
      
      samplesByCountry = split(selectedMetaMap$sample.ID, selectedMetaMap$Geography)
      return(samplesByCountry)
    }
    getSampleListByNotIndiCountry <- function(basicMetaMap, mgsInfo) {
      mContDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="disease"&!is.na(mgsInfo$Control)])
      hContDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="healthy"&!is.na(mgsInfo$Control)])
      
      healthyDatasets = unique(mgsInfo$ID[mgsInfo$Classification=="healthy"])
      notIndiDatasets = healthyDatasets[!healthyDatasets %in% indiDatasets]
      
      healthyIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Case" & basicMetaMap$dataset.ID %in% healthyDatasets]
      NotIndiHealtyIDs = healthyIDs[!healthyIDs %in% indiIDs]
      
      mContIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Control" & basicMetaMap$dataset.ID %in% mContDatasets]
      hContIDs = basicMetaMap$sample.ID[basicMetaMap$type=="Control" & basicMetaMap$dataset.ID %in% hContDatasets]
      mhContIDs = unique(c(mContIDs, hContIDs))
      normalIDs = unique(c(NotIndiHealtyIDs, mhContIDs))
      
      
      allHealthyIds = unique(c(mhContIDs, normalIDs))
      selectedMetaMap = basicMetaMap[basicMetaMap$sample.ID %in% allHealthyIds,]
      
      samplesByCountry = split(selectedMetaMap$sample.ID, selectedMetaMap$Geography)
      exCountries = c("Finland","Italy")# those samples having less than 5 samples
      samplesByCountry = samplesByCountry[!names(samplesByCountry)%in% exCountries]
      return(samplesByCountry)
    }
    
    samplesByDisease = getSampleListByDisease(basicMetaMap, mgsInfo)
    samplesByAllCountry = getSampleListByCountry(basicMetaMap, mgsInfo)
    sampleListByIndiCountry = getSampleListByIndiCountry(basicMetaMap, mgsInfo)
    sampleListByNotIndiCountry = getSampleListByNotIndiCountry(basicMetaMap, mgsInfo)
    
    showFourClassPerBarPlotByCountry <- function(samplesByCountry, target, mergeMat, cutOffs=c(1e-8,  1e-6)) {
      
      perList <- lapply(samplesByCountry, function(currSamples){
        currVec = mergeMat[target,currSamples]
        absent = sum(currVec < cutOffs[1])/length(currVec)
        medium = sum(currVec < cutOffs[2] & currVec >= cutOffs[1] )/length(currVec)
        high = sum(currVec >= cutOffs[2] )/length(currVec)
        #present = sum(currVec >= cutOff)/length(currVec)
        return(c(high, medium, absent))
      })
      perMat <- do.call(cbind, perList)
      rownames(perMat) = c("high","medium","absent")
      cols = c("red","gold","gray")
      cols =  adjustcolor( cols, alpha.f = 0.8)
      par(mar=c(5,10,4,2))
      barplot(perMat, horiz = T, las=1, col=cols)
      legend("topleft",legend=c("high","medium","lowly"),
             fill=cols)
      
      par(mar=c(5.1,4.1,4.1,2.1))
      return(perMat)
    }
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
    showPerBarPlotByDisease <- function(samplesByDisease, target, normalMat, mergeMat, cutOff=1e-7) {
      perList <- lapply(samplesByDisease, function(currSamples){
        currVec = mergeMat[target,currSamples]
        present = sum(currVec >= cutOff)/length(currVec)
        absent = sum(currVec < cutOff)/length(currVec)
        return(c(present, absent))
      })
      normMode = T
      if (normMode) {
        currVec = normalMat[target,]
        present = sum(currVec >= cutOff)/length(currVec)
        absent = sum(currVec < cutOff)/length(currVec)
        perList$normal = c(present, absent)
      }
      perMat <- do.call(cbind, perList)
      
      par(mar=c(5,10,4,2))
      barplot(perMat, horiz = T, las=1)
      par(mar=c(5.1,4.1,4.1,2.1))
      return(perMat)
    }
    showBoxPlotsByDisease <- function(samplesByDisease, target, normalMat, mergeMat, ylims = c(0,1e-5)) {
      abdList <- lapply(samplesByDisease, function(currSamples){
        return(mergeMat[target,currSamples])
      })
      normMode = T
      if (normMode) {
        currVec = normalMat[target,]
        abdList$normal = currVec
      }
      par(mar=c(5,10,4,2))
      boxplot(abdList,horizontal = T,las=1, ylim=ylims)
      par(mar=c(5.1,4.1,4.1,2.1))
      return(abdList)
    }
    
    targetSpecies = "msp_0025"
    # perList=showPerBarPlotByDisease(samplesByDisease, targetSpecies, normalMat, mergeMat )
    getSpeciesName(targetSpecies, taxo)
    ylims= c(0,1e-5)
    outList=showBoxPlotsByDisease(samplesByDisease, targetSpecies, normalMat, mergeMat,ylims )
    
    do.call(c,lapply(outList, mean))
    do.call(c,lapply(outList, median))
    
    
    targetSpecies = "msp_0025"
    perList=showFourClassPerBarPlotByCountry(samplesByAllCountry, targetSpecies, mergeMat, c(1e-8, 1e-6) )
    
    getSpeciesName(targetSpecies, taxo)
    
    outList=showBoxPlotsByCountry(samplesByAllCountry, targetSpecies, familyMat )
    do.call(c,lapply(outList, mean))
    do.call(c,lapply(outList, median))
    
    speciesDiseaseList = statisticsByDisease(samplesByDisease, sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.2, numCut=17, allCut=8)
    speciesDiseaseGeneralList = statisticsByDiseaseGeneral(samplesByDisease, sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.3, numCut=17)
    
    speciesSpecificList = statisticsByCountrySpecific(samplesByAllCountry, mergeMat, esCut=0.2, numCut=11)
    
    ####
    speciesUrbanSpecificList = statisticsByUrbanSpecific(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    speciesIndiSpecificList = statisticsByIndiSpecific(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    speciesDisSpecificList = statisticsByDiseaseSpecific(sampleListByNotIndiCountry, samplesByDisease,mergeMat, esCut=0.2)
    
    speciesUrbanSpecificTab = melt(speciesUrbanSpecificList)
    speciesIndiSpecificTab = melt(speciesIndiSpecificList)
    speciesDisSpecificTab = melt(speciesDisSpecificList)
    
    speciesUrbanSpecificTab$type="healthy"
    speciesIndiSpecificTab$type="indi"
    speciesDisSpecificTab$type="disease"
    
    speciesAllSpecificTab = rbind(speciesUrbanSpecificTab,
                                  speciesIndiSpecificTab,
                                  speciesDisSpecificTab)
    speciesAllSpecificTab$name = getSpeciesName(speciesAllSpecificTab$value,taxo)
    write.table(speciesAllSpecificTab,"speciesAllSpecificTab.txt",sep='\t')
    
    
    save(speciesUrbanSpecificList,
         speciesIndiSpecificList,
         speciesDisSpecificList,
         file="specific.species.by.category.RData")
    
    load("specific.species.by.category.RData")
    
    
    speciesUrbanSpecificList_es_02 = unique(do.call(c,speciesUrbanSpecificList))
    #177 species
    speciesIndiSpecificList_es_02 = unique(do.call(c,speciesIndiSpecificList))
    #187 species
    speciesDisSpecificList_es_02 = unique(do.call(c,speciesDisSpecificList))
    #72 species
    
    table(speciesUrbanSpecificList_es_02 %in% speciesIndiSpecificList_es_02)
    table(speciesUrbanSpecificList_es_02 %in% speciesDisSpecificList_es_02)
    table(speciesIndiSpecificList_es_02 %in% speciesDisSpecificList_es_02)
    
    
    urbanSpecies = unique(do.call(c,speciesUrbanSpecificList))
    indiSpecies = unique(do.call(c,speciesIndiSpecificList))
    
    getSpeciesName(urbanSpecies, taxo)
    getSpeciesName(indiSpecies, taxo)
    
    
    speciesUrbanList = statisticsByUrbanCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    speciesIndiList = statisticsByIndiCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    
    familyUrbanList = statisticsByUrbanCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, familyMat, esCut=0.2)
    familyIndiList = statisticsByIndiCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, familyMat, esCut=0.2)
    
    genusUrbanList = statisticsByUrbanCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, genusMat, esCut=0.2)
    genusIndiList = statisticsByIndiCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, genusMat, esCut=0.2)
    
    
    speciesDiseaseDownList = statisticsDownByDisease(samplesByDisease[names(samplesByDisease)!="preterm"], sampleListByNotIndiCountry, normalMat, mergeMat, esCut=0.2)
    speciesUrbanDownList = statisticsDownByUrbanCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    speciesIndiDownList = statisticsDownByIndiCountry(sampleListByNotIndiCountry, sampleListByIndiCountry, mergeMat, esCut=0.2)
    
    lapply(speciesUrbanList, function(x)getSpeciesName(x,taxo))
    
    sankeyMode = T
    targetList = speciesIndiSpecificList
    if (sankeyMode) {
      library(networkD3)
      library(tidyverse)
      require(RColorBrewer)
      
      taretMelt = melt(targetList)
      taretMelt$name = sapply(taretMelt$value, function(x) getSpeciesName(x,taxo))
      
      taretMelt$newName = paste(taretMelt$name, taretMelt$value, sep=",")
      taretMelt$value=1
      
      plotMode=T
      if (plotMode) {
        fsSize=12
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        
        targetOut = taretMelt
        elements = c(unique(targetOut$newName),unique(targetOut$L1))
        elementCols = col_vector[1:length(elements)]#
        
        interactions = targetOut[,c("newName","L1","value")]
        numInteractions = dim(interactions)[1]
        
        sourceInds = match(interactions$newName, elements)
        targetInds = match(interactions$L1, elements)
        interactionValues = interactions$value
        
        links = interactions
        colnames(links) = c("source","target","value")
        
        
        # From these flows we need to create a node data frame: it lists every entities involved in the flow
        nodes=data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique())
        links$IDsource=match(links$source, nodes$name)-1 
        links$IDtarget=match(links$target, nodes$name)-1
        
        
        color_scale <- data.frame(
          range = elementCols,
          domain = elements, 
          nodes = as.character(nodes$name),
          stringsAsFactors = FALSE
        )
        
        # Make the Network. I call my colour scale with the colourScale argument
        sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "name", fontFamily = "arial", fontSize = fsSize, colourScale = JS(
          sprintf(
            'd3.scaleOrdinal()  
            .domain(%s)
            .range(%s)
            ',
            jsonlite::toJSON(color_scale$domain),
            jsonlite::toJSON(color_scale$range)
          )
        ))
      }
      
    }
    
    
    getSampleListByDisease <- function(basicMetaMap) {
      
    }
    
    getSampleListByAge <- function(basicMetaMap) {
      
    }
    
    getSampleListByBmi <- function(basicMetaMap) {
      
    }
    
    
    
    ### country-wise comparison
    getStats <- function()
      
      ### 
  }
  
  
  }

checkLucasKoProfile = F
if (checkLucasKoProfile) {
  gutKoBestMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/gutKoBestMat.1990.RData"
  load(gutKoBestMatRData) #gutKoBestMat
  gutKoBestSums = colSums(gutKoBestMat)
  gutKoBestSums = sort(gutKoBestSums, decreasing = T)
  gutKoCountTab = data.frame(KO=names(gutKoBestSums),
                             MSP_COUNT = gutKoBestSums,
                             KO_gene = koDescMap$gene[match(names(gutKoBestSums), koDescMap$KO)],
                             DESCRIPTION=koDescMap$desc[match(names(gutKoBestSums), koDescMap$KO)])
  
  gutKoCountTabFile = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/gutKoCountTab.for.lucas.txt"
  write.table(gutKoCountTab, file=gutKoCountTabFile, row.names = F, quote = F, sep="\t")
}

diseaseNetworkMode = F
if (diseaseNetworkMode) {
  
  mergeMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//mergeMatUpdatedJaccard.RData"
  genusMatUpdatedJaccRData = "J:\\Deposit\\Project\\2018_microbiome_atlas\\microbiome.atlas.Rproj//genusMatUpdatedJaccard.RData"
  
  dim(genusMatUpdated)
  
  generateJacc = F
  if (generateJacc) {
    mergeMatUpdatedDigit = (mergeMatUpdated > 0)*1
    mergeMatUpdatedJacc = getJaccMat(mergeMatUpdatedDigit)
    save(mergeMatUpdatedJacc, 
         file=mergeMatUpdatedJaccRData)
    
    genusMatUpdatedDigit = (genusMatUpdated > 0)*1
    genusMatUpdatedJacc = getJaccMat(genusMatUpdatedDigit)
    save(genusMatUpdatedJacc, 
         file=genusMatUpdatedJaccRData)
    
  }
  
  loadJacc = T
  if (loadJacc) {
    load(mergeMatUpdatedJaccRData) #mergeMatUpdatedJacc
    jaccNet = makeJaccNet(mergeMatUpdatedJacc, cutJacc = 0.5)
    
    ### geography network ###
    geoList = newNormalAllList
    geoClusterScores = getClusterScores(jaccNet, geoList)
    
    jaccClusterInfo = annotateModulesByCC(jaccNet, geoList)
    writeModuleNetworkAnnot(jaccClusterInfo, 
                            "nodeTableJaccGeo.txt",
                            "edgeTableJaccGeo.txt")
    
    ### disease network ###
    diseaseList = newDiseaseList
    diseaseClusterScores = getClusterScores(jaccNet, diseaseList)
    
    jaccClusterInfo = annotateModulesByCC(jaccNet, diseaseList)
    writeModuleNetworkAnnot(jaccClusterInfo, 
                            "nodeTableJaccDisease.txt",
                            "edgeTableJaccDisease.txt")
    
    
  }
  
  
}

transmissionMap = F
if (transmissionMap) {
  
  generateMode = F
  if (generateMode) {
    
    freqs = rowSums(mergeMatUpdated>0)
    median(freqs) 
    
    freqSpecies = names(freqs[freqs > median(freqs)])
    rareSpecies = names(freqs[freqs < median(freqs)])
    
    digitMat = (mergeMatUpdated>0)*1
    digitFreqMat = (mergeMatUpdated[freqSpecies,]>0)*1
    digitRareMat = (mergeMatUpdated[rareSpecies,]>0)*1
    digitFreqJaccMat = 1-jaccard(digitFreqMat)
    digitRareJaccMat = 1-jaccard(digitRareMat)
    digitRareJaccMat[is.nan(digitRareJaccMat)]=0
    
    #digitJaccMat = getJaccMat(digitMat)
    
    #require(prabclus)
    digitJaccMat = 1-jaccard(digitMat)
    
    digitJaccMatRData = "C:\\Data\\comparative.analysis.healthy.sweden\\digitJaccMat.RData"
    
    save(digitJaccMat, file=digitJaccMatRData)
    #save(digitJaccMat, digitFreqJaccMat, digitRareJaccMat, file="C:\\Data\\comparative.analysis.healthy.sweden\\digitJaccMat.all.RData")
    #load(digitJaccMatRData)
    
    
    #genusup = getTaxaSumMat(mergeMat, taxo, "genus")
    #digitGenusMat = (genusMat>0)*1
    digitGenusMat = (genusMatUpdated>0)*1
    digitGenusJaccMat = 1-jaccard(digitGenusMat)
    save(digitGenusJaccMat, file="C:\\Data\\comparative.analysis.healthy.sweden\\digitGenusJaccMat.RData")
    
  }
  
  loadMode = T
  if (loadMode) {
    mgsJaccRData = "C:/Data/comparative.analysis.healthy.sweden/digitJaccMat.RData"
    genusJaccRData = "C:/Data/comparative.analysis.healthy.sweden/digitGenusJaccMat.RData"
    load(mgsJaccRData) #digitJaccMat
    load(genusJaccRData)   #digitGenusJaccMat 
    
    plot(density(digitJaccMat[lower.tri(digitJaccMat)]), ylim=c(0,5), main="")
    lines(density(digitGenusJaccMat[lower.tri(digitGenusJaccMat)]), col="red")
    
    samplesBySubjects = split(basicMetaMapUpdated$sample.ID,
                              basicMetaMapUpdated$metadata.ID)
    lenSamplesBySubjects = do.call(c,lapply(samplesBySubjects,length))
    samplesBySubjects = samplesBySubjects[lenSamplesBySubjects>1]
    
    
    allMgsJaccs = digitJaccMat[lower.tri(digitJaccMat)]
    allGenusJaccs = digitGenusJaccMat[lower.tri(digitGenusJaccMat)]
    
    subjectMgsJaccs = do.call(c,lapply(samplesBySubjects, function(currSamples){
      currMat = digitJaccMat[rownames(digitJaccMat)%in%currSamples, colnames(digitJaccMat)%in%currSamples]
      return(currMat[lower.tri(currMat)])
    }))
    plot(density(digitJaccMat[lower.tri(digitJaccMat)]), ylim=c(0,10), main="")
    lines(density(subjectMgsJaccs), col="red")
    
    
    subjectGenusJaccs = do.call(c,lapply(samplesBySubjects, function(currSamples){
      currMat = digitGenusJaccMat[rownames(digitGenusJaccMat)%in%currSamples, colnames(digitGenusJaccMat)%in%currSamples]
      return(currMat[lower.tri(currMat)])
    }))
    plot(density(digitGenusJaccMat[lower.tri(digitGenusJaccMat)]), ylim=c(0,10), main="")
    lines(density(subjectGenusJaccs), col="red")
    
  }
  
  mgsLevel = T
  if (mgsLevel) {
    load("digitJaccMat.RData")
    names(datasetList)
    
    #digitJaccMatCut = (digitJaccMat > 0.5)*1
    digitJaccMatCut = (digitRareJaccMat > 0.5)*1
    
    countsMat = lapply(datasetList, function(cohortI) {
      countsRow = lapply(datasetList,  function(cohortJ) {
        counts = sum(digitJaccMatCut[cohortI, cohortJ])
        return(counts)
        
      })
      
      return(countsRow)
    })
    countsMelt = melt(countsMat)
    countsTransmissionMat = acast(countsMelt, L1~ L2, value.var = "value")
    heatmap.2((countsTransmissionMat>0)*1, trace="none", key=F,
              col = colorRampPalette(c("white","gold"))(256),
              sepcolor = "gray",
              colsep=0:32,rowsep=0:32,
              sepwidth = c(0.01,0.01)
    )
    
    comms = c("Firmicutes","Bacteroidetes","unclassified","Proteobacteria",
              "Actinobacteria","Fusobacteria","Verrucomicrobia","Tenericutes",
              "Euryarchaeota","Synergistetes","unclassified Eukaryota",
              "Elusimicrobia","Candidatus Melainabacteria","unclassified Bacteria")
    
    freqVec = table(getGenusName(freqSpecies,taxo))
    rareVec = table(getGenusName(rareSpecies,taxo))
    comms = names(freqVec)[names(freqVec) %in% names(rareVec)]
    
    freqVec = table(getPhylaName(freqSpecies,taxo))
    rareVec = table(getPhylaName(rareSpecies,taxo))
    
    plot(as.numeric(freqVec[comms]), as.numeric(rareVec[comms]),log="xy")
    
    freqVec[names(rareVec)]
    
    wellnessIDs = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id28"]
    wellnessBySubjs = split(wellnessIDs, getSubectIds(wellnessIDs)  )
    
    wellnessJaccBySubs = lapply(wellnessBySubjs, function(x){
      currMat  = digitJaccMat[x,x]
      return(currMat[lower.tri(currMat)])
    })
    
    twinIDs = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id39"]
    madaIds = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id42"]
    t1dIds = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id35"]
    hmpIds = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id36"]
    rccIds = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id41"]
    rheIds = basicMetaMap$sample.ID[basicMetaMap$dataset.ID=="id31"]
    
    digitTwinMat = digitMat[,twinIDs]
    digitHmpMat = digitMat[,hmpIds]
    digitT1dMat = digitMat[,t1dIds]
    digitRccMat = digitMat[,rccIds]
    digitRheMat = digitMat[,rheIds]
    digitMadaMat = digitMat[,madaIds]
    digitWellMat = digitMat[,wellnessIDs]
    
    
    plot(rowSums(digitTwinMat), rowSums(digitHmpMat))
    plot(rowSums(digitRccMat), rowSums(digitHmpMat))
    plot(rowSums(digitRheMat), rowSums(digitHmpMat))
    plot(rowSums(digitT1dMat), rowSums(digitHmpMat))
    plot(rowSums(digitWellMat), rowSums(digitHmpMat))
    
    
    plot(rowSums(digitTwinMat), rowSums(digitMat))
    plot(rowSums(digitRccMat), rowSums(digitMat))
    plot(rowSums(digitRheMat), rowSums(digitMat))
    plot(rowSums(digitT1dMat), rowSums(digitMat))
    plot(rowSums(digitMadaMat), rowSums(digitMat))
    plot(rowSums(digitWellMat), rowSums(digitMat))
    plot(rowSums(digitHmpMat), rowSums(digitMat))
    
    plot(rowSums(digitMadaMat), rowSums(digitHmpMat))
    
    
    
    
    getSpeciesName(names(sort(rowSums(digitTwinMat),decreasing = T)[1:10]),taxo)
    getSpeciesName(names(sort(rowSums(digitHmpMat),decreasing = T)[1:10]),taxo)
    getSpeciesName(names(sort(rowSums(digitWellMat),decreasing = T)[1:10]),taxo)
    getSpeciesName(names(sort(rowSums(digitWellMat),decreasing = F)[1:10]),taxo)
    
    twinMat = digitJaccMat[twinIDs,twinIDs]
    wellMat = digitJaccMat[wellnessIDs,wellnessIDs]
    madaMat =  digitJaccMat[madaIds,madaIds]
    t1dMat = digitJaccMat[t1dIds,t1dIds]
    hmpMat = digitJaccMat[hmpIds,hmpIds]
    rccMat = digitJaccMat[rccIds,rccIds]
    rheMat = digitJaccMat[rheIds,rheIds]
    
    
    jaccRegList = list(twin=twinMat[lower.tri(twinMat)],
                       wellness=wellMat[lower.tri(wellMat)],
                       madagascar=madaMat[lower.tri(madaMat)],
                       hmp = hmpMat[lower.tri(hmpMat)],
                       t1d=t1dMat[lower.tri(t1dMat)],
                       rcc=rccMat[lower.tri(rccMat)],
                       rhe=rheMat[lower.tri(rheMat)])
    boxplot(jaccRegList,horizontal = T,las=1)
    
    
    boxplot(twinMat[lower.tri(twinMat)])
    boxplot(wellMat[lower.tri(wellMat)])
  }
  
  genusLevel = T
  if (genusLevel) {
    names(datasetList)
    
    digitGenusJaccMatCut = (digitGenusJaccMat > 0.9)*1
    
    countsMat = lapply(datasetList, function(cohortI) {
      countsRow = lapply(datasetList,  function(cohortJ) {
        counts = sum(digitGenusJaccMatCut[cohortI, cohortJ])
        return(counts)
        
      })
      
      return(countsRow)
    })
    countsMelt = melt(countsMat)
    countsTransmissionMat = acast(countsMelt, L1~ L2, value.var = "value")
    heatmap.2((countsTransmissionMat>0)*1, trace="none", key=F,
              col = colorRampPalette(c("white","gold"))(256),
              sepcolor = "gray",
              colsep=0:32,rowsep=0:32,
              sepwidth = c(0.01,0.01)
    )
    
  }
  
}

tsneAfterPcaReductionMode = F
if (tsneAfterPcaReductionMode) {
  
  getTsneAfterPca <- function(targetMat, plotMode = F, numComps = 5, perplexity=50) {
    pcaOut = prcomp(t(targetMat))
    pcaMat = pcaOut$x
    eigs <- pcaOut$sdev^2
    pcaSummary = rbind(
      SD = sqrt(eigs),
      Proportion = eigs/sum(eigs),
      Cumulative = cumsum(eigs)/sum(eigs))
    print(pcaSummary[,1:numComps])
    print("")
    print("... pca done")
    
    set.seed(1)
    tsneOut = Rtsne(pcaMat[,1:numComps], dims = 2, perplexity = perplexity)
    tsneMat = tsneOut$Y
    print("... tSNE done")
    
    if (plotMode) {
      tsneCols = rep("#80808088",dim(tsneMat)[1])
      names(tsneCols) = colnames(targetMat)
      
      
    }
    return(tsneMat)
  }
  
  normMat = mergeMatUpdated[,normUniqAll]
  normTsneMat = getTsneAfterPca(normMat)
  
  tsneCols = rep("#80808088",dim(tsneMat)[1])
  names(tsneCols) = colnames(targetMat)
  tsneCols[names(tsneCols) %in% newNonWesternIds] = "red"
  
  plot(tsneMat, col=tsneCols)
  
}

updateEnteroType = F
if (updateEnteroType) {
  
  genusMode = T
  if (genusMode) {
    
    genusMat = getTaxaSumMat(mergeMatUpdated, taxoNew, "genus")
    genusMatT = round(t(genusMat)*10^9)
    
    generateMode = T
    if (generateMode) {
      set.seed(1)
      genusFit <- mclapply(1:3, dmn, count=genusMatT, verbose=TRUE) 
      save(genusFit, 
           file="J:\\Deposit\\Project\\2018_microbiome_atlas\\enterotypeData\\dirichlet.multinomial.genusFit3.2019.0726.RData")
      
      processMode = T
      if (processMode) {
        
        lplc <- sapply(genusFit, laplace)
        plot(lplc, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")
        best <- genusFit[[which.min(lplc)]]
        
        p0 <- fitted(genusFit[[1]], scale=TRUE) 
        p3 <- fitted(best, scale=TRUE)
        colnames(p3) <- paste("m", 1:3, sep="")
        (meandiff <- colSums(abs(p3 - as.vector(p0))))
        diff <- rowSums(abs(p3 - as.vector(p0)))
        o <- order(diff, decreasing=TRUE)
        cdiff <- cumsum(diff[o]) / sum(diff)
        df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 10)
        
        heatmapdmn(genusMatT, genusFit[[1]], best, 10, lblwidth = 4000)
        
        clusterAssigned = apply(genusFit[[3]]@group, 1, function(x) which.max(x))
        
      }
      
      
    }
    
    loadMode = F
    if (loadMode) {
      dmnNewRData ="J:\\Deposit\\Project\\2018_microbiome_atlas\\enterotypeData\\dirichlet.multinomial.genusFit3.2019.0726.RData"
      load(dmnNewRData)
      
      processMode = T
      if (processMode) {
        lplc <- sapply(genusFit, laplace)
        plot(lplc, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")
        best <- genusFit[[which.min(lplc)]]
        
        p0 <- fitted(genusFit[[1]], scale=TRUE) 
        p3 <- fitted(best, scale=TRUE)
        colnames(p3) <- paste("m", 1:3, sep="")
        (meandiff <- colSums(abs(p3 - as.vector(p0))))
        diff <- rowSums(abs(p3 - as.vector(p0)))
        o <- order(diff, decreasing=TRUE)
        cdiff <- cumsum(diff[o]) / sum(diff)
        df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 10)
        #m1 --> firmicutes
        #m2 --> bacteroides
        #m3 --> prevotella
        
        
        heatmapdmn(genusMatT, genusFit[[1]], best, 10, lblwidth = 4000)
        
        clusterAssigned = apply(genusFit[[3]]@group, 1, function(x) which.max(x))
        
        clusterAssignedList = split(names(clusterAssigned), clusterAssigned)
        
        names(clusterAssignedList) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
        enteroTypeTab = melt(clusterAssignedList)
        colnames(enteroTypeTab) = c("sample.ID","enteroType")
        basicMetaMapUpdated$enteroType = enteroTypeTab$enteroType[match(basicMetaMapUpdated$sample.ID, enteroTypeTab$sample.ID)]
        save(basicMetaMapUpdated, file=basicMetaCleanRDataUpdated2)
        
        View(basicMetaMapUpdated)
        
      }
    }
    
    
  }
  
  mgsMode = T
  if (mgsMode) {
    mgsMatT = round(t(mergeMatUpdated)*10^9)
    
    set.seed(1)
    mgsFit <- mclapply(1:3, dmn, count=mgsMatT, verbose=TRUE) 
    save(mgsFit, 
         file="J:\\Deposit\\Project\\2018_microbiome_atlas\\enterotypeData\\dirichlet.multinomial.mgsFit3.2019.0726.RData")
    processMode = T
    if (processMode) {
      lplc <- sapply(mgsFit, laplace)
      plot(lplc, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")
      best <- mgsFit[[which.min(lplc)]]
      
      p0 <- fitted(mgsFit[[1]], scale=TRUE) 
      p3 <- fitted(best, scale=TRUE)
      colnames(p3) <- paste("m", 1:3, sep="")
      (meandiff <- colSums(abs(p3 - as.vector(p0))))
      diff <- rowSums(abs(p3 - as.vector(p0)))
      o <- order(diff, decreasing=TRUE)
      cdiff <- cumsum(diff[o]) / sum(diff)
      df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 10)
      #m1 --> firmicutes
      #m2 --> bacteroides
      #m3 --> prevotella
      
      
      heatmapdmn(mgsMatT, mgsFit[[1]], best, 10, lblwidth = 4000)
      
      clusterAssigned = apply(mgsFit[[3]]@group, 1, function(x) which.max(x))
      
      clusterAssignedList = split(names(clusterAssigned), clusterAssigned)
      
      names(clusterAssignedList) = c("ET-Bacteroides","ET-Firmicutes","ET-Prevotella")
      enteroTypeTab = melt(clusterAssignedList)
      colnames(enteroTypeTab) = c("sample.ID","enteroType")
      basicMetaMapUpdated$enteroTypeM = enteroTypeTab$enteroType[match(basicMetaMapUpdated$sample.ID, enteroTypeTab$sample.ID)]
      save(basicMetaMapUpdated, file=basicMetaCleanRDataUpdated2)
      
      View(basicMetaMapUpdated)
      
    }
  }
  
}

