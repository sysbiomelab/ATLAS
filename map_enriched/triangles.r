    antismashMode = T
    if(antismashMode) {
      
      
      pCut = 0.05
      varMat = getVarScoreMat(funcMat, antismashTerms, etfSpecies, etbSpecies, etpSpecies)
      sigMat = selectVarScoreMat(varMat, pCut)
      
      colnames(sigMat)[1:3] = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% mucinCazymes] = "mucin"
      # sigMat$class[sigMat$lab %in% storageCazymes] = "storage"
      # sigMat$class[sigMat$lab %in% pectinCazymes] = "storage"
      # sigMat$class = factor(sigMat$class, levels=c("others","mucin","storage"))
      View(sigMat)
      
      
      
      
      #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(fill="#80808033",#aes(fill=class),
                   size = 4, 
                   shape = 21, 
                   color = "#808080") +
        #scale_fill_manual(values=c("#80808033","#0000ff88","#ff000088","#d4af3788")) +
        theme_ggtern(12, "") %+replace% ggplot2::theme_void(12, "") %+replace% theme(tern.axis.text = element_blank(), tern.axis.title = element_blank(), legend.position = "none" ) + 
        labs( x       = "",
              xarrow  = "European",
              y       = "",
              yarrow  = "Westernized",
              z       = "",
              zarrow  = "Non-westernized") 
      
      
      if (F) {
        funcMatTranspose = t(funcMat)
        funcMatTranspose = funcMatTranspose[antismashTerms,]
        etfFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etfOnlySpecies]
        etbFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etbOnlySpecies]
        etpFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etpOnlySpecies]
        
        chisq_FB = getChisqStat(etfFuncMat, etbFuncMat, antismashTerms)
        chisq_BP = getChisqStat(etbFuncMat, etpFuncMat, antismashTerms)
        chisq_FP = getChisqStat(etfFuncMat, etpFuncMat, antismashTerms)
        
        pCut=0.05
        sig_FB = rownames(chisq_FB[!is.nan(chisq_FB$pvalue) & chisq_FB$pvalue<pCut,])
        sig_BP = rownames(chisq_BP[!is.nan(chisq_BP$pvalue) & chisq_BP$pvalue<pCut,])
        sig_FP = rownames(chisq_FP[!is.nan(chisq_FP$pvalue) & chisq_FP$pvalue<pCut,])
        sig_union = unique(c(sig_FB, sig_BP, sig_FP))
        
        sigMat = getFractionMat(etfFuncMat, etbFuncMat, etpFuncMat, sig_union)
        sigMat = as.data.frame(sigMat)
        colnames(sigMat) = c("European","Pacific","NonWesternized")
        sigMat$lab = rownames(sigMat)
        
        require(ggtern)
        #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
        ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
          geom_point(aes(fill=class),
                     size = 2, 
                     shape = 21, 
                     color = "black") +
          scale_fill_manual(values=c("#80808055","#0000ff88","#ff000088","#d4af3788")) +
          theme_rgbw() 
      }
      
    }
    
    vfMode = T
    if(vfMode) {
      #invasionProteins = c("APECO1_2093","t3801", "t0781", "MMAR_2285", "lmo1745", "STM2511")
      invasionProteins = c("MMAR_2285", "APECO1_2093", "Btr_1145", "lmo1847", "t0781", "t3801", "Lmon1_020100004517", "lmo0691")
      intracellularProteins = c("Rv2987c", "lmo0055", "lmo0136", "lmo0137", "lmo0810", "lmo2825", "Lmon1_020100002767", "SL1344_1828", "MMAR_2285")
      
      pCut = 1e-2
      varMat = getVarScoreMat(funcMat, vfTerms, etfSpecies, etbSpecies, etpSpecies)
      sigMat = selectVarScoreMat(varMat, pCut)
      
      colnames(sigMat)[1:3] = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      sigMat$Prod = patricVfDescTab$Product[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      sigMat$Category = patricVfDescTab$Classification[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      sigMat$class = "others"
      sigMat$class[sigMat$lab %in% invasionProteins] = "invasion"
      sigMat$class[sigMat$lab %in% intracellularProteins] = "intracellular"
      sigMat$class = factor(sigMat$class, levels=c("others","invasion", "intracellular"))
      
      
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% mucinCazymes] = "mucin"
      # sigMat$class[sigMat$lab %in% storageCazymes] = "storage"
      # sigMat$class[sigMat$lab %in% pectinCazymes] = "storage"
      # sigMat$class = factor(sigMat$class, levels=c("others","mucin","storage"))
      View(sigMat)
      
      
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(aes(fill=class),
                   size = 4, 
                   shape = 21, 
                   color = "#80808011") +
        scale_fill_manual(values=c("#80808033","#0000ff88","#ff000088","#d4af3788")) +
        theme_ggtern(12, "") %+replace% ggplot2::theme_void(12, "") %+replace% theme(tern.axis.text = element_blank(), tern.axis.title = element_blank(), legend.position = "none" ) + 
        labs( x       = "",
              xarrow  = "European",
              y       = "",
              yarrow  = "Westernized",
              z       = "",
              zarrow  = "Non-westernized") 
      
      if (F) {
        funcMatTranspose = t(funcMat)
        funcMatTranspose = funcMatTranspose[vfTerms,]
        etfFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etfOnlySpecies]
        etbFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etbOnlySpecies]
        etpFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etpOnlySpecies]
        
        chisq_FB = getChisqStat(etfFuncMat, etbFuncMat, vfTerms)
        chisq_BP = getChisqStat(etbFuncMat, etpFuncMat, vfTerms)
        chisq_FP = getChisqStat(etfFuncMat, etpFuncMat, vfTerms)
        
        pCut=1e-2
        sig_FB = rownames(chisq_FB[!is.nan(chisq_FB$pvalue) & chisq_FB$pvalue<pCut,])
        sig_BP = rownames(chisq_BP[!is.nan(chisq_BP$pvalue) & chisq_BP$pvalue<pCut,])
        sig_FP = rownames(chisq_FP[!is.nan(chisq_FP$pvalue) & chisq_FP$pvalue<pCut,])
        sig_union = unique(c(sig_FB, sig_BP, sig_FP))
        
        sigMat = getFractionMat(etfFuncMat, etbFuncMat, etpFuncMat, sig_union)
        sigMat = as.data.frame(sigMat)
        colnames(sigMat) = c("European","Pacific","NonWesternized")
        sigMat$lab = rownames(sigMat)
        sigMat$Prod = patricVfDescTab$Product[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
        sigMat$Category = patricVfDescTab$Classification[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
        sigMat$class = "others"
        sigMat$class[sigMat$lab %in% invasionProteins] = "invasion"
        sigMat$class = factor(sigMat$class, levels=c("others","invasion"))
        
        require(ggtern)
        #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
        ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
          geom_point(aes(fill=class),
                     size = 2, 
                     shape = 21, 
                     color = "black") +
          scale_fill_manual(values=c("#80808055","#0000ff88","#ff000088","#d4af3788")) +
          theme_rgbw() 
        
        View(sigMat)
      }
      
    }
    
    jgiMode = F
    if(jgiMode) {
      
      funcMatTranspose = t(funcMat)
      funcMatTranspose = funcMatTranspose[jgiTerms,]
      etfFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etfOnlySpecies]
      etbFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etbOnlySpecies]
      etpFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etpOnlySpecies]
      
      chisq_FB = getChisqStat(etfFuncMat, etbFuncMat, jgiTerms)
      chisq_BP = getChisqStat(etbFuncMat, etpFuncMat, jgiTerms)
      chisq_FP = getChisqStat(etfFuncMat, etpFuncMat, jgiTerms)
      
      pCut=1e-5
      sig_FB = rownames(chisq_FB[!is.nan(chisq_FB$pvalue) & chisq_FB$pvalue<pCut,])
      sig_BP = rownames(chisq_BP[!is.nan(chisq_BP$pvalue) & chisq_BP$pvalue<pCut,])
      sig_FP = rownames(chisq_FP[!is.nan(chisq_FP$pvalue) & chisq_FP$pvalue<pCut,])
      sig_union = unique(c(sig_FB, sig_BP, sig_FP))
      
      sigMat = getFractionMat(etfFuncMat, etbFuncMat, etpFuncMat, sig_union)
      sigMat = as.data.frame(sigMat)
      colnames(sigMat) = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      
      
      # sigMat$Prod = patricVfDescTab$Product[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$Category = patricVfDescTab$Classification[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% invasionProteins] = "invasion"
      # sigMat$class = factor(sigMat$class, levels=c("others","invasion"))
      # 
      
      
      
      require(ggtern)
      #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(#aes(fill=class),
          size = 2, 
          shape = 21, 
          color = "black") +
        scale_fill_manual(values=c("#80808055","#0000ff88","#ff000088","#d4af3788")) +
        theme_void() 
      
      View(sigMat)
    }
    
    arMode = F
    if(arMode) {
      
      
      pCut = 1e-2
      varMat = getVarScoreMat(funcMat, arTerms, etfSpecies, etbSpecies, etpSpecies)
      sigMat = selectVarScoreMat(varMat, pCut)
      
      colnames(sigMat)[1:3] = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      
      
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% mucinCazymes] = "mucin"
      # sigMat$class[sigMat$lab %in% storageCazymes] = "storage"
      # sigMat$class[sigMat$lab %in% pectinCazymes] = "storage"
      # sigMat$class = factor(sigMat$class, levels=c("others","mucin","storage"))
      View(sigMat)
      
      
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(fill="#80808033",#aes(fill=class),
                   size = 4, 
                   shape = 21, 
                   color = "#808080") +
        #scale_fill_manual(values=c("#80808033","#0000ff88","#ff000088","#d4af3788")) +
        theme_void() +
        labs( x       = "",
              xarrow  = "European",
              y       = "",
              yarrow  = "Westernized",
              z       = "",
              zarrow  = "Non-westernized") 
      
      if (F) {
        
      }
      funcMatTranspose = t(funcMat)
      funcMatTranspose = funcMatTranspose[arTerms,]
      etfFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etfOnlySpecies]
      etbFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etbOnlySpecies]
      etpFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etpOnlySpecies]
      
      chisq_FB = getChisqStat(etfFuncMat, etbFuncMat, arTerms)
      chisq_BP = getChisqStat(etbFuncMat, etpFuncMat, arTerms)
      chisq_FP = getChisqStat(etfFuncMat, etpFuncMat, arTerms)
      
      pCut=1e-5
      sig_FB = rownames(chisq_FB[!is.nan(chisq_FB$pvalue) & chisq_FB$pvalue<pCut,])
      sig_BP = rownames(chisq_BP[!is.nan(chisq_BP$pvalue) & chisq_BP$pvalue<pCut,])
      sig_FP = rownames(chisq_FP[!is.nan(chisq_FP$pvalue) & chisq_FP$pvalue<pCut,])
      sig_union = unique(c(sig_FB, sig_BP, sig_FP))
      
      sigMat = getFractionMat(etfFuncMat, etbFuncMat, etpFuncMat, sig_union)
      sigMat = as.data.frame(sigMat)
      colnames(sigMat) = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      
      
      # sigMat$Prod = patricVfDescTab$Product[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$Category = patricVfDescTab$Classification[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% invasionProteins] = "invasion"
      # sigMat$class = factor(sigMat$class, levels=c("others","invasion"))
      # 
      
      
      
      require(ggtern)
      #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(#aes(fill=class),
          size = 2, 
          shape = 21, 
          color = "black") +
        scale_fill_manual(values=c("#80808055","#0000ff88","#ff000088","#d4af3788")) +
        theme_rgbw() 
      
      View(sigMat)
    }
    
    arMode = F
    if(arMode) {
      
      funcMatTranspose = t(funcMat)
      funcMatTranspose = funcMatTranspose[arTerms,]
      etfFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etfOnlySpecies]
      etbFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etbOnlySpecies]
      etpFuncMat = funcMatTranspose[,colnames(funcMatTranspose) %in% etpOnlySpecies]
      
      chisq_FB = getChisqStat(etfFuncMat, etbFuncMat, arTerms)
      chisq_BP = getChisqStat(etbFuncMat, etpFuncMat, arTerms)
      chisq_FP = getChisqStat(etfFuncMat, etpFuncMat, arTerms)
      
      pCut=1e-5
      sig_FB = rownames(chisq_FB[!is.nan(chisq_FB$pvalue) & chisq_FB$pvalue<pCut,])
      sig_BP = rownames(chisq_BP[!is.nan(chisq_BP$pvalue) & chisq_BP$pvalue<pCut,])
      sig_FP = rownames(chisq_FP[!is.nan(chisq_FP$pvalue) & chisq_FP$pvalue<pCut,])
      sig_union = unique(c(sig_FB, sig_BP, sig_FP))
      
      sigMat = getFractionMat(etfFuncMat, etbFuncMat, etpFuncMat, sig_union)
      sigMat = as.data.frame(sigMat)
      colnames(sigMat) = c("European","Pacific","NonWesternized")
      sigMat$lab = rownames(sigMat)
      
      
      # sigMat$Prod = patricVfDescTab$Product[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$Category = patricVfDescTab$Classification[match(sigMat$lab, patricVfDescTab$RefSeq.Locus.Tag)]
      # sigMat$class = "others"
      # sigMat$class[sigMat$lab %in% invasionProteins] = "invasion"
      # sigMat$class = factor(sigMat$class, levels=c("others","invasion"))
      # 
      
      
      
      require(ggtern)
      #ggtern(data=sigCazyMat, aes(x=F,y=B, z=P, label=lab)) + 
      ggtern(data=sigMat, aes(x=European,y=Pacific, z=NonWesternized, label=lab)) + 
        geom_point(#aes(fill=class),
          size = 2, 
          shape = 21, 
          color = "black") +
        scale_fill_manual(values=c("#80808055","#0000ff88","#ff000088","#d4af3788")) +
        theme_rgbw() 
      
      View(sigMat)
    }
    
  }
  
}

funcCompareByInflowOutflow = T
if (funcCompareByInflowOutflow) {
  
  loadFuncAnnot = T
  if (loadFuncAnnot) {
    
    ### PFAM ###
    igc2PfamMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/PFAM/igc2PfamMat.RData"
    load(igc2PfamMatRData)
    pfamTerms = colnames(pfamMat)
    
    ### KEGG ###
    gutKoBestMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/gutKoBestMat.1990.RData"
    load(gutKoBestMatRData)
    koTerms = colnames(gutKoBestMat)
    
    
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
    funcMatRData = "j://Deposit/Project/2018_microbiome_atlas/functional.annotation/funcMat.20190808.RData"
    load(funcMatRData)
  }
  
  loadInOutFlowStat = T
  if (loadInOutFlowStat) {
    markovStatTabFinalFile = "C://Data//comparative.analysis.healthy.sweden//markovStatTab.final.3.txt"
    markovStatTabFinal = read.table(markovStatTabFinalFile, sep="\t", stringsAsFactors = F)
    markovStatTabFinal$GssOss = NA
    markovStatTabFinal$GssOss[rownames(markovStatTabFinal) %in% gssSpecies] = "GSS"
    markovStatTabFinal$GssOss[rownames(markovStatTabFinal) %in% ossSpecies] = "OSS"
    #write.table(markovStatTabFinal, markovStatTabFinalFile, sep="\t")
    
    gutTaxoSel = gutTaxo[rownames(markovStatTabFinal),]
    gutGenusMat = acast(data=gutTaxoSel, MSP ~ genus )
    gutFamilyMat = acast(data=gutTaxoSel, MSP ~ family )
    gutClassMat = acast(data=gutTaxoSel, MSP ~ class )
    gutPhylumMat = acast(data=gutTaxoSel, MSP ~ phylum )
    gutGenusMat[is.na(gutGenusMat)] = 0
    gutFamilyMat[is.na(gutFamilyMat)] = 0
    gutClassMat[is.na(gutClassMat)] = 0
    gutPhylumMat[is.na(gutPhylumMat)] = 0
    gutGenusMat = gutGenusMat[,!grepl("unclassified",colnames(gutGenusMat))]
    gutFamilyMat = gutFamilyMat[,!grepl("unclassified",colnames(gutFamilyMat))]
    gutClassMat = gutClassMat[,!grepl("unclassified",colnames(gutClassMat))]
    gutPhylumMat = gutPhylumMat[,!grepl("unclassified",colnames(gutPhylumMat))]
    
    getSlideSums <- function(input, wd=10) {
      output = rep(0, length(input))
      slideSums = rollapply(input, width=wd, by=1, FUN=mean, align="left", partial=T)
      output[seq(slideSums)] = slideSums
      return(output)
    }
    getSignatures <- function(targetFlows, funcMatNew, currFcs, wd=10) {
      fcTab = data.frame(targetflow=targetFlows,
                         funcMatNew[,currFcs])
      set.seed(1)
      fcTab = fcTab[sample(nrow(fcTab)),]
      fcTab = fcTab[order(fcTab[,"targetflow"], decreasing = T),]
      #fcTab$slideSum = getSlideSums(fcTab[,"func"], wd)
      fcTab[,-1] = apply(fcTab[,-1], 2, function(x) getSlideSums(x, wd))
      fcTab = fcTab[order(fcTab[,"targetflow"], decreasing = F),]
      #print(head(fcTab))
      fcLevelInput = t(as.matrix(fcTab))
      #rownames(fcLevelInput) = c("Flow","Function", "slideSum")
      colnames(fcLevelInput) = NULL
      #outMat = fcLevelInput[c("Flow","slideSum"),]
      return(fcLevelInput)
      
    }
    getSignature <- function(targetFlows, funcMatNew, currFc, wd=10) {
      fcTab = data.frame(targetflow=targetFlows,
                         func = funcMatNew[,currFc])
      set.seed(1)
      fcTab = fcTab[sample(nrow(fcTab)),]
      fcTab = fcTab[order(fcTab[,"targetflow"], decreasing = T),]
      fcTab$slideSum = getSlideSums(fcTab[,"func"], wd)
      fcTab = fcTab[order(fcTab[,"targetflow"], decreasing = F),]
      #print(head(fcTab))
      fcLevelInput = t(as.matrix(fcTab))
      rownames(fcLevelInput) = c("Flow","Function", "slideSum")
      colnames(fcLevelInput) = NULL
      outMat = fcLevelInput[c("Flow","slideSum"),]
      return(outMat)
      
    }
    plotSignature <- function(outMat, as = 5) {
      out = levelplot(outMat, aspect = as, colorkey=F,
                      xlab="",ylab="",axes=F,
                      scales = list(tck = c(0,0), x=list(rot=90), draw=F),
                      col.regions = colorRampPalette(c("white","black"))(512))
      return(out)
    }
    
    
    
    #funcMatNew = funcMat[match(rownames(markovStatTabFinal), rownames(funcMat)),]
    #table(rownames(markovStatTabFinal) %in% rownames(funcMat))
    #all 1411 TRUE
    
    mucinCazymes = c("cazy.GT27", "cazy.GT14", "cazy.GT11", "cazy.GT10", "cazy.GH89", "cazy.GH84", "cazy.GH33", "cazy.GH20", "cazy.GH29", "cazy.GH18", "cazy.GH123", "cazy.GH109", "cazy.CBM51", "cazy.CBM50", "cazy.CBM32")
    storageCazymes = c("cazy.CBM34", "cazy.CBM48", "cazy.CBM20", "cazy.GH13", "cazy.GH133", "cazy.GH32", "cazy.GH36", "cazy.GH37", "cazy.GH77", "cazy.GH57", "cazy.GT101", "cazy.GT26", "cazy.GT3", "cazy.GT35", "cazy.GT5")
    pectinCazymes = c("cazy.CBM77","cazy.CE8","cazy.GH105","cazy.GH28", "cazy.PL1")
    
    if (F) {
      targetTerms = vfTerms
      targetTerms = arTerms
      targetTerms = antismashTerms
      targetTerms = cazyTerms
      targetTerms = pectinCazymes
      targetTerms = jgiTerms
      targetTerms = mgePfamTerms
      targetTerms = pfamTerms
      targetTerms = koTerms
      
      table(lmStatInMat[,"AdjP"]< adjCut) # 2086 terms were significant
      table(lmStatOutMat[,"AdjP"]< adjCut) # 294 terms were significant
      
      adjCut=1e-5
      table(rownames(lmStatInMat[lmStatInMat[,"AdjP"]< adjCut,]) %in% targetTerms)
      table(rownames(lmStatOutMat[lmStatOutMat[,"AdjP"]< adjCut,]) %in% targetTerms)
      
      
    }
    
    inflowAssocFunctionMat <- function(markovStatTabFinal, funcMat, targetFuncs = NULL) {
      funcMatNew = funcMat[match(rownames(markovStatTabFinal), rownames(funcMat)),]
      if (is.null(targetFuncs)) funcs = colnames(funcMatNew) 
      if (!is.null(targetFuncs)) funcs = targetFuncs
      
      inflows = markovStatTabFinal$inflow
      outflows = markovStatTabFinal$outflow
      names(inflows) = rownames(markovStatTabFinal)
      names(outflows) = rownames(markovStatTabFinal)
      
      
      lmStatInMat = do.call(rbind, sapply(funcs, function(fc) {
        funcProfiles = funcMatNew[,fc]
        lmOut  = lm(funcProfiles~inflows)
        pvalue = summary(lmOut)$coefficient["inflows","Pr(>|t|)"]
        estimate = summary(lmOut)$coefficient["inflows","Estimate"]
        return(c(estimate, pvalue))
        
      }, simplify = F))
      lmStatInMat = cbind(lmStatInMat, p.adjust(lmStatInMat[,2], "BH"))
      colnames(lmStatInMat) = c("Estimate", "Pvalue", "AdjP")
      print(dim(lmStatInMat))
      lmStatInMat = lmStatInMat[lmStatInMat[,"Estimate"]>0,]
      
      print(dim(lmStatInMat))
      return(lmStatInMat)
      #head(lmStatInMat[order(lmStatInMat[,"Pvalue"]),],n=20)
      View(lmStatInMat[lmStatInMat[,"AdjP"]< 1e-10,])
      tt = rownames(lmStatInMat[lmStatInMat[,"AdjP"]< 1e-10,])
      table(tt %in% vfTerms)
      
      View(data.frame(tt,getDescFromKoTerms(tt, koDescMap)))
      
      if (F) {
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% arTerms,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% vfTerms,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% mgePfamTerms,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% cazyTerms,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% storageCazymes,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% mucinCazymes,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% jgiTerms,]
        lmStatInMat2 = lmStatInMat[rownames(lmStatInMat) %in% antismashTerms,]
        
        
        head(lmStatInMat2[order(lmStatInMat2[,"AdjP"]),],n=20)
      }
      if (F) {
        
        currFc = "K18220"
        currFc = "K09759"
        currFc = "K06404"
        
        currFc = "Integrase_DNA"
        currFc = "MobC"
        currFc = "cazy.GT32"
        #currFc = "tetM"
        currFc = "Nonsporulating"
        currFcs = c("Anaerobe","Integrase_DNA", "MobC", "bacteriocin", "van", "cazy.GT5")
        currFcs = c("Integrase_DNA","MobC", "cazy.GT5")
        currFcs = c("Integrase_DNA","MobC")
        
        plotSignature(getSignatures(inflows, funcMatNew, currFcs, 100), 2)#, 1.2
        
        levelplot(getSignatures(outflows, funcMatNew, currFcs, 100), aspect = 1.5, 
                  xlab="",ylab="",axes=F,
                  scales = list(tck = c(0,0), x=list(rot=90), draw=F),
                  col.regions = colorRampPalette(c("white","black"))(512))
        
        
      }
      if (F) {
        
        
        plotSignature(getSignature(inflows, funcMatNew, currFc, 100))
        
        checkBoxPlot = F
        if (checkBoxPlot) {
          integraseBySplit = split(inflows, funcMatNew[,"Integrase_DNA"])
          boxplot(integraseBySplit)
        }
        
      }
      
      return(lmStatInMat)
    }
    outflowAssocFunctionMat <- function(markovStatTabFinal, funcMat, targetFuncs = NULL) {
      funcMatNew = funcMat[match(rownames(markovStatTabFinal), rownames(funcMat)),]
      if (is.null(targetFuncs)) funcs = colnames(funcMatNew) 
      if (!is.null(targetFuncs)) funcs = targetFuncs
      
      outflows = markovStatTabFinal$outflow
      inflows = markovStatTabFinal$inflow
      
      lmStatOutMat = do.call(rbind, sapply(funcs, function(fc) {
        funcProfiles = funcMatNew[,fc]
        lmOut  = lm(funcProfiles~outflows)
        pvalue = summary(lmOut)$coefficient["outflows","Pr(>|t|)"]
        estimate = summary(lmOut)$coefficient["outflows","Estimate"]
        return(c(estimate, pvalue))
        
      }, simplify = F))
      lmStatOutMat = cbind(lmStatOutMat, p.adjust(lmStatOutMat[,2], "BH"))
      colnames(lmStatOutMat) = c("Estimate", "Pvalue", "AdjP")
      lmStatOutMat = lmStatOutMat[lmStatOutMat[,"Estimate"]>0,]
      #head(lmStatOutMat[order(lmStatOutMat[,"AdjP"]),],n=20)
      return(lmStatOutMat)
      
      if (F) {
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% arTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% cazyTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% vfTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% jgiTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% mgePfamTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% antismashTerms,]
        lmStatOutMat2 = lmStatOutMat[rownames(lmStatOutMat) %in% koTerms,]
        
        head(lmStatOutMat2[order(lmStatOutMat2[,"AdjP"]),],n=20)
        
        
      }
      
      if (F) {
        
        
        
        currFc = "siderophore"
        currFc = "t2pks-ladderane"
        currFc = "K00867"
        currFc = "lmo1003"
        #lmo1267 --> cell division trigger factor
        currFcs = c("Facultative", "Microaerophilic", "lmo1267", "K02248")
        
        plotSignature(getSignatures(outflows, funcMatNew, currFcs, 100), 1.5)
        
      }
      
      checkBoxPlot = F
      if (checkBoxPlot) {
        integraseBySplit = split(inflows, funcMatNew[,"Integrase_DNA"])
        boxplot(integraseBySplit)
      }
      
      return(lmStatOutMat)
    } 
