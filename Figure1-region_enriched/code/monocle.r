require(igraph)
require(monocle)
require(dplyr)
require(Biobase)

# Load data
vect_atlas <- read.csv('../../data/msp.csv', row.names=1, header=T, stringsAsFactors=F)
meta <- read.csv('../../data/sampleID.csv', row.names=1, header=T, stringsAsFactors=F)
taxo <- read.csv('../../data/gutTaxo.csv',
		       row.names=1,
		       header=T,
		       stringsAsFactors=F)

# Filter healthy metadata
meta <- meta %>% filter(Disease == 'Healthy')

# Transform msp to pseudocounts
normAllMat = vect_atlas*10^9

# Prepare monocle
trajectorySampleTab = data.frame(Library=meta$sample.ID,
				 Well=meta$project.ID,
				 Hours=0,
				 Media="NONE")

rownames(trajectorySampleTab) = meta$sample.ID

trajectoryTaxaTab = data.frame(gene_short_name=taxo$species, 
			       biotype="protein_coding",
			       num_cells_expressed=1,
			       use_for_ordering=FALSE)

rownames(trajectoryTaxaTab) = rownames(taxo)
trajectoryTaxaTab = trajectoryTaxaTab[match(rownames(normAllMat), rownames(taxo)),]

pd <- new("AnnotatedDataFrame", data=trajectorySampleTab)
fd <- new("AnnotatedDataFrame", data=trajectoryTaxaTab)
cds <- newCellDataSet(as(as.matrix(normAllMat), 'dgTMatrix'), 
		      phenoData=pd, 
		      featureData=fd,
		      lowerDetectionLimit=0.1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Run monocle
cds <- reduceDimension(cds, 
		       max_components=2,
		       method='DDRTree')
#cds <- orderCells(cds)
output <- t(data.frame(cds@reducedDimS))
write.csv(output, '../results/monocleout.csv')
