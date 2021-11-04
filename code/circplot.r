require(RColorBrewer)
require(circlize)

# new plan, just load the proteomics and do the correlations between them
#msps are vect atlas
load('../data/current.example.MGS.other.omics.RData')
load('../data/vect_atlas.RData')
circosInput <- read.csv('../newcircos.csv')
circosInput <- read.csv('../newgo.csv')

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#selectedTab = mixedModelMetaboFamilyTab[mixedModelMetaboFamilyTab$exp.var > 0.1, ]
#selectedTab = mixedModelMetaboMgsTab[mixedModelMetaboMgsTab$exp.var > 0.25, ]
selectedTab = mixedModelMetaboMgsTab[mixedModelMetaboMgsTab$exp.var > 0.5, ]
circosInput = circosInput[circosInput$value > 0.99, ]

# explained variance > 10%
selectedTabPos = selectedTab[selectedTab$tMgs > 0,]
# positive relations
selectedTabNeg = selectedTab[selectedTab$tMgs < 0,]
# negative relations

circosInput = data.frame(from=getGenusName(selectedTabPos$mgs, taxo),
                         to=selectedTabPos$metabo,
                         value=(abs(selectedTabPos$exp.var)),
                         stringsAsFactors = F)


#circosInput = circosInput[!grepl("unclassified",circosInput$from),]

elements = (unique(c(circosInput$from, circosInput$to)))
lenElements = length(elements)

circosGridCols = (col_vector)[1:length(elements)]
names(circosGridCols) = elements

#circos.par(track.margin=c(0,0)) 
set.seed(1)
#circos.initialize(elements, xlim=c(0,1))
circos.initialize(elements, xlim=c(0,1))
chordDiagram(circosInput)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
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
