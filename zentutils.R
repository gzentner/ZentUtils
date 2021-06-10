library(ZentUtils)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

zent <- zentutils(system.file("extdata", "homer_reb1_badis_motif_scan.bed",
                              package = "ZentUtils"), genome = "sacCer3")

zent <- sort_regions(zent, by = "score", decreasing = TRUE)

zent <- expand_regions(zent, length = 10)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3

zent <- get_seqs(zent, region_type = "expanded", genome = genome)

pdf(file = "colormap.pdf")
color_map(zent)
dev.off()
