#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(RColorBrewer)
require(circlize)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("Supply edge file, input for circos plot.n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  f=args[1]
}

circosInput <- read.csv(f)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#circosInput = data.frame(from=getGenusName(selectedTabPos$mgs, taxo),
#                         to=selectedTabPos$metabo,
#                         value=(abs(selectedTabPos$exp.var)),
#                         stringsAsFactors = F)


#circosInput = circosInput[!grepl("unclassified",circosInput$from),]

elements = (unique(c(circosInput$from, circosInput$to)))
lenElements = length(elements)

circosGridCols = (col_vector)[1:length(elements)]
names(circosGridCols) = elements

#circos.par(track.margin=c(0,0))
set.seed(1)
#circos.initialize(elements, xlim=c(0,1))

#pdf('pdfMeta.pdf')
#pdf('pdfMspd.pdf')
pdf(args[2])
circos.initialize(elements, xlim=c(0,1))
#chordDiagram(circosInput)
#circos.par(cell.padding = c(0.02, 0, 0.02, 0))
chordDiagram(circosInput, annotationTrack = "grid",
             preAllocateTracks = 1, grid.col = circosGridCols)

circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, cex=0.7,
        CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA)
#circos.trackPlotRegion(track.index = 1,
#  panel.fun = function(x, y) {
#    xlim = get.cell.meta.data("xlim")
#    ylim = get.cell.meta.data("ylim")
#    sector.name = get.cell.meta.data("sector.index")
#    circos.text(mean(xlim), ylim[1] , cex = 0.6,
#                sector.name, facing = "clockwise",
#                niceFacing = T, adj = c(0, 0.5))
#    circos.axis(h = "top", labels = F,
#                major.tick = F,
#                sector.index = sector.name,
#                track.index = 2)
#    }, bg.border = NA)
circos.clear()
dev.off()
