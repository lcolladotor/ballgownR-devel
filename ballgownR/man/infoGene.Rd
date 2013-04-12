\name{infoGene}
\alias{infoGene}
\title{Data extraction of the ballgown results to visualize with \code{viewGene}}
\usage{
  infoGene(geneID, gown, group, countIntron = TRUE,
    whichCount = "mrcount")
}
\arguments{
  \item{geneID}{specifies the geneID to look at}

  \item{gown}{specifies the output from \code{readGown}}

  \item{group}{specifies the sample groups and will be used
  to determine the colors.}

  \item{countIntron}{specifies whether to use the intron
  data or not. If set to TRUE, \code{whichCount} has to be
  mrcount, rcount or ucount.}

  \item{whichCount}{has to be mrcount, rcount or ucount if
  \code{countIntron=TRUE}. Otherwise it can also be cov and
  mcov.}
}
\value{
  A list with the information needed to make the plot with
  \code{viewGene}
}
\description{
  \code{infoGene} extracts the information according to a
  geneID. This can be the rcounts or mrcounts per exon (and
  introns optionall). Alternatively, it can extract the raw
  coverage.
}
\examples{
dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
samplePattern <- "sample"
gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)
info <- system.file("extdata", "sample_info.txt", package="ballgownR")
info <- read.table(info, header=TRUE)
match <- sapply(names(gown$dirs), function(x) { which(info$dir == x)})
geneInfo <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match])
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}
\seealso{
  viewGene
}

