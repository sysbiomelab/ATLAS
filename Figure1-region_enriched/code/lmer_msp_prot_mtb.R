#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(dplyr)
require(lme4)
require(ade4)
require(pls)
require(car)
require(lmerTest)
require(tidyverse)
require(MuMIn)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("At least supply 2 matrices, row are samples with name individual_visit .n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  frm=args[1]
  tom=args[2]
  outfile= paste(args[1],args[2])
} else if (length(args)==3) {
  frm=args[1]
  tom=args[2]
  outfile=args[3]
} else if (length(args)>3) {
  stop("too many arguments, matrix1 matrix2 and outfilename")
}


###########
### load datafiles
#load("sunjae/current.example.MGS.other.omics.RData")
msp_genus<-read.csv(frm, row.names=1)
#prot<-read.csv(tom, row.names=1)
mtb<-read.csv(tom, row.name=1)
#extract_features
mgsFeatures<-colnames(msp_genus)
#protFeatures<-colnames(prot)
metaboFeatures<-colnames(mtb)

msp_genus$subject =  gsub("_[1-4]", "", rownames(msp_genus))
#prot$subject =  gsub("_[1-4]", "", rownames(msp_genus))
mtb$subject =  gsub("_[1-4]", "", rownames(msp_genus))
#msp_genus$visit = gsub("^[0-9]+_", "", rownames(mtb))

#metaboFeatures
#build subject id and visit columns

### run lmer on possible pairs of variables
mixedModel= do.call(rbind,sapply(metaboFeatures, function(currMetabo) {
    print(paste(which(metaboFeatures == currMetabo), date(), sep="|||||"))
    ## RUN MIXED MODEL LMER
    mixedModelResult = do.call(rbind.data.frame, sapply(mgsFeatures, function(currMgs) {
    #print(which(mgsFeatures == currMgs))
        currData = data.frame(metabolite = mtb[,currMetabo],
                              mgs = msp_genus[,currMgs],
                              subj = msp_genus$subject,
                              stringsAsFactors = F)
        lmerMetaboMgs = lmer(metabolite ~ mgs + (1|subj), data = currData)
        if (dim(summary(lmerMetaboMgs)$coefficients)[1]<2) {
            return(list(pMgs=NA,  tMgs = NA, exp.var = NA, from=currMgs, to=currMetabo))
            return(result)
        }
        #print(summary(lmerMetaboMgs))
        #print(which(mgsFeatures == currMgs))
        varNames = rownames(summary(lmerMetaboMgs)$coefficients)
        exp.var = r.squaredGLMM(lmerMetaboMgs)[1]
        tMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Estimate"], NA)
        pMgs = ifelse( any(varNames == "mgs"), summary(lmerMetaboMgs)$coefficients["mgs","Pr(>|t|)"], NA)
        return(list(pMgs=pMgs,   tMgs = tMgs, exp.var = exp.var, from = currMgs, to = currMetabo))
    }, simplify=F))
    return(mixedModelResult)
}, simplify = F))

#mixModel_msp_prot$species = getSpeciesName(mixedModelMetaboMgsTab$mgs, taxo)
#mixedMdel_RData = "mixedModel.Rdata"
mixedModel=mixedModel[complete.cases(mixedModel),]
save(mixedModel, file=outfile)


#selectedTab = mixedModel[mixedModel$exp.var > 0.25, ]
#selectedTabPos = selectedTab[selectedTab$tMgs > 0,]
#selectedTabNeg = selectedTab[selectedTab$tMgs < 0,]

# negative relations
