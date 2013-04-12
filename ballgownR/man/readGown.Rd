\name{readGown}
\alias{readGown}
\title{Read ballgown output *.ctab files from a data directory with multiple samples}
\usage{
  readGown(dataDir, samplePattern, verbose = TRUE)
}
\arguments{
  \item{dataDir}{must be a character that specifies the
  data directory path. Inside of it there must be a
  directory for each sample.}

  \item{samplePattern}{must be a character that specifies
  the pattern to use that identifies the sample
  directories}

  \item{verbose}{is a logical. If TRUE updates will be
  printed via message.}
}
\value{
  a list with e2t: mappting between exons and transcripts,
  i2t: mapping between introns and transcripts, intron:
  i_data.ctab merged in wide format, exon: e_data.ctab
  merged in mide format, trans: t_data.ctab merged in wide
  format, dirs: sample directories (specifies the order in
  which they are merged), mergedDate: date this function
  was run
}
\description{
  \code{readGown} reads the *.ctab files that ballgown
  creates. It then uses plyr to merge the information from
  multiple samples.
}
\examples{
#Read toy data
dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
samplePattern <- "sample"
gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)
}
\author{
  Leonardo Collado-Torres \email{lcollado@jhsph.edu}
}

