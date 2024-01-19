require(igraph)
require(monocle)
require(Biobase)
vect_atlas <- read.csv('../results/healthymsp.csv',
		       row.names=1,
		       header=T,
		       stringsAsFactors=F)

meta <- read.csv('../results/monoclemeta.csv',
		       row.names=1,
		       header=T,
		       stringsAsFactors=F)

normAllMat = vect_atlas*10^9

trajectorySampleTab = data.frame(Library=meta$sample.ID,
				 Well=meta$project.ID,
				 Hours=0,
				 Media="NONE")
rownames(trajectorySampleTab) = meta$sample.ID

trajectoryTaxaTab = data.frame(gene_short_name=rownames(normAllMat), 
			       biotype="protein_coding",
			       num_cells_expressed=1,
			       use_for_ordering=FALSE)
rownames(trajectoryTaxaTab) = rownames(normAllMat)

pd <- new("AnnotatedDataFrame", data=trajectorySampleTab)
fd <- new("AnnotatedDataFrame", data=trajectoryTaxaTab)
cds <- newCellDataSet(as(as.matrix(normAllMat), 'dgTMatrix'), 
		      phenoData=pd, 
		      featureData=fd,
		      lowerDetectionLimit=0.1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds, 
		       max_components=2,
		       method='DDRTree')
cds <- orderCells(cds)
save(cds, file=normDataTrajectory2RData)
