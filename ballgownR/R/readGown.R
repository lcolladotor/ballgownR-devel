#' Read ballgown output *.ctab files from a data directory with multiple samples
#'
#' \code{readGown} reads the *.ctab files that ballgown creates. It then uses plyr to merge the information from multiple samples.
#'
#' @param dataDir must be a character that specifies the data directory path. Inside of it there must be a directory for each sample.
#' @param samplePattern must be a character that specifies the pattern to use that identifies the sample directories
#' @param verbose is a logical. If TRUE updates will be printed via message.
#' @return a list with e2t: mappting between exons and transcripts, i2t: mapping between introns and transcripts, intron: i_data.ctab merged in wide format, exon: e_data.ctab merged in mide format, trans: t_data.ctab merged in wide format, dirs: sample directories (specifies the order in which they are merged), mergedDate: date this function was run
#' @export
#' @author Leonardo Collado-Torres \email{lcollado@@jhsph.edu}
#' @examples
#' #Read toy data
#' dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
#' samplePattern <- "sample"
#' gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)

#### Main function
readGown <- function(dataDir, samplePattern, verbose=TRUE) {
	if(verbose) message(date())
	
	## Load required pkgs
	suppressMessages(library(plyr))
	
	if(FALSE){
		# for testing
		dataDir <- "/amber2/scratch/lcollado/brain_rna/ballgown"
		samplePattern <- "orb"
	}
	
	## Identify the sample directories
	dirs <- list.files(path=dataDir, pattern=samplePattern, full.names=TRUE)
	names(dirs) <- list.files(path=dataDir, pattern=samplePattern)
	n <- length(dirs)
	
	if(FALSE){
		# for testing
		n <- 3 
	}
	
	## Read linking tables
	if(verbose) message("Reading linking tables")
	e2t <- read.table(list.files(dirs[1], "e2t.ctab", full.names=TRUE), header=TRUE, sep="\t", colClasses=c("integer", "integer"))
	i2t <- read.table(list.files(dirs[1], "i2t.ctab", full.names=TRUE), header=TRUE, sep="\t", colClasses=c("integer", "integer"))
	
	## Order by transcript id
	e2t <- e2t[order(e2t$t_id), ]
	i2t <- i2t[order(i2t$t_id), ]
	rownames(e2t) <- 1:nrow(e2t)
	rownames(i2t) <- 1:nrow(i2t)
	
	## Read counts for all introns in <reference_transcripts>
	if(verbose) message("Reading intron data files")
	intronFiles <- sapply(dirs, list.files, pattern="i_data.ctab", full.names=TRUE)
	intronAll <- lapply(intronFiles, .readIntron)
	
	## Merge the results
	if(verbose) message("Merging intron data")
	intron <- join_all(intronAll, by=c("i_id", "chr", "strand", "start", "end"), type="left")
	colnames(intron)  <- c("i_id", "chr", "strand", "start", "end", paste(c("rcount", "ucount", "mrcount"), rep(names(dirs), each=3), sep="."))
	

	## Read read counts and raw coverage info for all exons in <reference_transcripts>
	if(verbose) message("Reading exon data files")
	exonFiles <- sapply(dirs, list.files, pattern="e_data.ctab", full.names=TRUE)
	exonAll <- lapply(exonFiles, .readExon)
	
	## Read exon data
	if(verbose) message("Merging exon data")
	exon <- join_all(exonAll, by=c("e_id", "chr", "strand", "start", "end"), type="left")
	colnames(exon) <- c("e_id", "chr", "strand", "start", "end", paste(c("rcount", "ucount", "mrcount", "cov", "cov_sd", "mcov", "mcov_sd"), rep(names(dirs), each=7), sep="."))
	
	
	## Read transcript data
	if(verbose) message("Reading transcription data files")
	transFiles <- sapply(dirs, list.files, pattern="t_data.ctab", full.names=TRUE)
	transAll <- lapply(transFiles, .readTrans)
	
	## Merge the results
	if(verbose) message("Merging transcription data")
	trans <- join_all(transAll, by=c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name"), type="left")
	colnames(trans) <- c("t_id", "chr", "strand", "start", "end", "t_name", "num_exons", "length", "gene_id", "gene_name", paste(c("cov", "FPKM"), rep(names(dirs), each=2), sep="."))
	
	if(verbose) message("Wrapping up the results")
	result <- list(e2t=e2t, i2t=i2t, intron=intron, exon=exon, trans=trans, dirs=dirs, mergedDate=date())
	
	if(verbose) message(date())
	# Done!
	return(result)
}

#### Auxiliary functions

## Read intron files
.readIntron <- function(file){
	intron <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric"))
	intron <- intron[order(intron$i_id), ]
	rownames(intron) <- 1:nrow(intron)
	return(intron)
}

## Read counts and raw coverage
.readExon <- function(file) {
	exon <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "integer", "integer", "numeric", "numeric", "numeric", "numeric", "numeric"))
	exon <- exon[order(exon$e_id), ]
	rownames(exon) <- 1:nrow(exon)
	return(exon)
}

## Read transcript data files
.readTrans <- function(file) {
	trans <- read.table(file, header=TRUE, sep="\t", colClasses=c("integer", "character", "factor", "integer", "integer", "character", "integer", "integer", "character", "character", "numeric", "numeric"))
	trans <- trans[order(trans$t_id), ]
	rownames(trans) <- 1:nrow(trans)
	return(trans)
}
