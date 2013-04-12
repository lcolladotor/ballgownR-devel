\name{viewGene}
\alias{viewGene}
\title{Interactive visualization of the ballgown results}
\usage{
  viewGene(geneInfo, exon.color = "#000000",
    location = "bottom", spacing = 0.02, html = NULL,
    wdir = NULL, browse = TRUE)
}
\arguments{
  \item{geneID}{specifies the geneID to look at}

  \item{exon.color}{is the color used for the exon lines}

  \item{location}{has to be 'bottom' or 'top' and it
  controls where the exon tracks are showed}

  \item{spacing}{is the factor by which the exons are
  separated}

  \item{html}{specifies the html output file}

  \item{wdir}{specifies the directory where the plot will
  be made. You must have writing permission.}
}
\value{
  the html file name with the resulting plot
}
\description{
  \code{viewGene} plots the information extracted with
  \code{infoGene} in an interactive html
}
\examples{
?viewGene # Read the help. Example to do!
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}

