\name{viewGene}
\alias{viewGene}
\title{Interactive visualization of the ballgown results: raw coverage or rcounts}
\usage{
  viewGene(geneID, gown, group, coverage = FALSE,
    tophatDir, exon.color = "#000000", location = "bottom",
    spacing = 1, html = NULL, wdir = NULL)
}
\arguments{
  \item{geneID}{specifies the geneID to look at}

  \item{gown}{specifies the output from \code{readGown}}

  \item{coverage=FALSE}{determines whether to make a raw
  coverage plot or show the rcounts}

  \item{tophatDir}{specifies the TopHat directory. Used
  when \code{coverage} is set to TRUE.}

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
  \code{viewGene} extracts the information according to a
  geneID and plots the coverage or the rcounts (intron and
  exon tables) in an interactive html
}
\examples{
?viewGene # Read the help. Example to do!
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}

