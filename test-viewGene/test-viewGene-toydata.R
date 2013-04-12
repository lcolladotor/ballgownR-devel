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
viewGene(geneInfo=geneInfo, html="toy-test.html", spacing=0.1)

## Coverage
geneInfo2 <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match], countIntron=FALSE, whichCount="cov")
viewGene(geneInfo=geneInfo2, html="toy-test-2.html", spacing=0.1)

## No intron
geneInfo3 <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match], countIntron=FALSE)
viewGene(geneInfo=geneInfo3, html="toy-test-3.html", spacing=0.1)

## top
viewGene(geneInfo=geneInfo, html="toy-test-4.html", spacing=0.1, location="top")

###


### Checking infoGene
geneID <- "gene_1"
whichCount <- "mrcount"
countIntron <- TRUE
group <- info$outcome[match]