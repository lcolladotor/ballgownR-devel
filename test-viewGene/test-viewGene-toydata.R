### Check that it works

suppressMessages(library(ballgownR))
dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
samplePattern <- "sample"
gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)

### Test the visualization
info <- system.file("extdata", "sample_info.txt", package="ballgownR")
info <- read.table(info, header=TRUE)
match <- sapply(names(gown$dirs), function(x) { which(info$dir == x)})
geneInfo <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match])

## Visualize
viewGene(geneInfo=geneInfo, html="toy-teset.html", spacing=0.1)
