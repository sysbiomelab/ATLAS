#!/bin/Rscript
require(DirichletMultinomial)
require(dplyr)
require(parallel)

#### converting MGS matrix to genus matrix and scaling data ###
#genusMat <- read.table("../results/genusAbund.csv", header=T, row.names=1, dec=".", sep=",")
genusMat <- read.table("../data/genusrelabund.csv", header=T, row.names=1, dec=".", sep=",")
#genusMat <- sample_n(genusMat, 400)
genusMatT <- round(genusMat*10^9)
#genusMatT <- genusMatT[,0:(ncol(genusMatT)-1)]

#### DMN clustering ###
set.seed(1)
genusFit <- mclapply(3, dmn, count=as.matrix(genusMatT), verbose=TRUE)
lplc <- sapply(genusFit, laplace)
#plot(lplc, type="b", xlab="Number of Dirichlet Components" ,ylab="Model Fit")
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

### plotting heatmap ###
#heatmapdmn(genusMatT, genusFit[[1]], best, 10, lblwidth = 4000)

### assigning cluster names ####
#clusterAssigned = apply(genusFit[[3]]@group, 1, function(x) which.max(x))
clusterAssigned = apply(genusFit[[1]]@group, 1, function(x) which.max(x))
clusterAssignedList = split(names(clusterAssigned), clusterAssigned)
names(clusterAssignedList) = c("ET-Firmicutes","ET-Bacteroides","ET-Prevotella")
write.csv(stack(clusterAssignedList), '../results/enterotypes.csv', quote=F)
