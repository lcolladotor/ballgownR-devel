\name{viewGene}
\alias{viewGene}
\title{Interactive visualization of the ballgown results}
\usage{
  viewGene(geneInfo, exon.color = "#000000",
    location = "bottom", spacing = 0.02, html = NULL,
    wdir = NULL, browse = TRUE, title = NULL, main = NULL,
    height = 500, width = 800)
}
\arguments{
  \item{geneInfo}{is the output of \code{infoGene}}

  \item{geneID}{specifies the geneID to look at}

  \item{exon.color}{is the color used for the exon lines}

  \item{location}{has to be 'bottom' or 'top' and it
  controls where the exon tracks are showed}

  \item{spacing}{is the factor by which the exons are
  separated}

  \item{html}{specifies the html output file}

  \item{wdir}{specifies the directory where the plot will
  be made. You must have writing permission.}

  \item{browse=TRUE}{if you want to open a browser window
  to view the graphic}

  \item{title}{specifies the title in the HTML file}

  \item{main}{specifies the title of the plot in the HTML
  file}

  \item{height}{specifies in pixels the height of the
  plot.}

  \item{width}{specifies in pixels the width of the plot.}
}
\value{
  The html file name of the resulting plot.
}
\description{
  \code{viewGene} plots the information extracted with
  \code{infoGene} in an interactive html
}
\examples{
dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
samplePattern <- "sample"
gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)
info <- system.file("extdata", "sample_info.txt", package="ballgownR")
info <- read.table(info, header=TRUE)
match <- sapply(names(gown$dirs), function(x) { which(info$dir == x)})
geneInfo <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match])
viewGene(geneInfo=geneInfo, html="toy-test.html", spacing=0.1)
viewGene(geneInfo=geneInfo, html="toy-test-2.html", spacing=0.1, location="top")
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}
\seealso{
  infoGene
}

