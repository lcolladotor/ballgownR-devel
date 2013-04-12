\name{infoGene}
\alias{infoGene}
\title{Data extraction of the ballgown results to visualize with \code{viewGene}}
\usage{
  infoGene(geneID, gown, group, coverage = FALSE,
    tophatDir = NULL)
}
\arguments{
  \item{geneID}{specifies the geneID to look at}

  \item{gown}{specifies the output from \code{readGown}}

  \item{coverage=FALSE}{determines whether to make a raw
  coverage plot or show the rcounts}

  \item{tophatDir}{specifies the TopHat directory. It is
  required when \code{coverage} is set to TRUE.}
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
?viewGene # Read the help. Example to do!
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}

