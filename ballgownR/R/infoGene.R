#' Data extraction of the ballgown results to visualize with \code{viewGene}
#'
#' \code{infoGene} extracts the information according to a geneID. This can be the rcounts or mrcounts per exon (and introns optionall). Alternatively, it can extract the raw coverage.
#'
#' @param geneID specifies the geneID to look at
#' @param gown specifies the output from \code{readGown}
#' @param coverage=FALSE determines whether to make a raw coverage plot or show the rcounts
#' @param tophatDir specifies the TopHat directory. Used when \code{coverage} is set to TRUE.
#' @return A list with the information needed to make the plot with \code{viewGene}
#' @export
#' @author Leonardo Collado-Torres \email{lcollado@@jhsph.edu}
#' @examples
#' ?viewGene # Read the help. Example to do!

viewGene <- function(geneID, gown, group, coverage=FALSE, tophatDir=NULL, exon.color="#000000", location="bottom", spacing=0.02, html=NULL, wdir=NULL) {
	## Load required libraries
	require(colorspace)
	suppressMessages(require(Rsamtools))
	
	## Find region of interest
	idx.t <- which(gown$trans$gene_id == geneID)
	
	## Find start and end
	start <- min(gown$trans$start[idx.t])
	end <- max(gown$trans$end[idx.t])	
	nBases <- end - start + 1
	
	## Subset trans, exon and intron
	t.id <- gown$trans$t_id[idx.t]
	names(t.id) <- gown$trans$t_name[idx.t]
	e.id <- lapply(t.id, function(x) { gown$e2t$e_id[ gown$e2t$t_id == x ] })
	i.id <- lapply(t.id, function(x) { gown$i2t$i_id[ gown$i2t$t_id == x ] })
	
	## Get sample names and the exon/intron rcount columns
	rcount.e <- which(gsub("\\..*", "", names(gown$exon[e.id[[1]][1], ])) == "rcount")
	rcount.i <- which(gsub("\\..*", "", names(gown$intron[i.id[[1]][1], ])) == "rcount")
	samples <- gsub("rcount\\.", "", names(gown$exon[e.id[[1]][1], ]))[rcount.e]
	
	## Get exon information by transcript
	exons.df <- lapply(1:length(e.id), function(e) {	
		## Get transcript name
		line <- names(e.id)[[e]]		
		res <- lapply(e.id[[e]], function(exon) {	
			x <- c(gown$exon$start[exon], gown$exon$end[exon])			
			y <- rep(e, length(x))
			data.frame(line=rep(line, length(x)), x=x, y=y)				 
		})
		do.call(rbind, res)
	})
	exons.df <- do.call(rbind, exons.df)
	
	## Build data: either summarized data or coverage data
	if(coverage == FALSE) {
				
		## Initialize rcount data
		rcount.df <- data.frame(line=rep(samples, each=nBases), x=start:end, y=rep(0, each=nBases))
		
		## Count exon information
		for(e in 1:length(e.id)) {
			for(exon in e.id[[e]]) {
				x <- seq(gown$exon$start[exon], gown$exon$end[exon])			
				vals <- as.vector(as.matrix(gown$exon[exon, rcount.e]))
				y <- rep(vals, each=length(x))
				rcount.df$y[ rcount.df$x %in% x ] <- rcount.df$y[ rcount.df$x %in% x ] + y
			}
		}
		
		## Count intron information
		for(i in 1:length(i.id)) {
			for(intron in i.id[[i]]) {
				x <- seq(gown$intron$start[intron], gown$intron$end[intron])			
				vals <- as.vector(as.matrix(gown$intron[intron, rcount.i]))
				y <- rep(vals, each=length(x))
				rcount.df$y[ rcount.df$x %in% x ] <- rcount.df$y[ rcount.df$x %in% x ] + y
			}
		}
		
		toAdd <- rcount.df
		
	} else if(coverage == TRUE & !is.null(tophatDir)) {
		
		## Initialize coverage data
		coverage.df <- data.frame(line=rep(samples, each=nBases), x=start:end, y=rep(0, each=nBases))
		
		## Define region to load
		which <- RangesList(IRanges(start=start, end=end))
		chr <- gown$trans$chr[idx.t[1]]
		names(which) <- chr
		param <- ScanBamParam(which=which)
		
		for(s in names(gown$dirs)) {
			bamfile <- BamFile(paste(tophatDir, s, "accepted_hits.bam", sep="/"))
			# Read it and get the coverage. Extract only the one for the chr in question
			cov <- coverage(readBamGappedAlignments(bamfile, param=param))[[chr]]
			coverage.df$y[ coverage.df$line == s] <- as.vector(cov[start:end])
		}
		
		toAdd <- coverage.df
		
	}

	## Purge useless info
	exons.filt <- filterInfoGene(exons.df)
	toAdd.filt <- filterInfoGene(toAdd)
	
	## Merge data
	data <- list(exons=exons.filt, transInfo=toAdd.filt, start=start, end=end, group=group, geneID=geneID)
	
	## Done
	return(data)
	
}

## This function removes uninformative points to avoid overloading the html
.filterInfoGene <- function(data) {
	## The filtering is done by unique elements of data$line (for example, by sample)
	result <- lapply(unique(data$line), function(x) {
		sub <- data[ data$line == x, ]
		which.diff <- which(diff(sub$y) != 0)
		if(length(which.diff) == 0) {
			res <- rbind(sub[1, ], sub[nrow(sub), ])
		} else {
			res <- rbind(sub[1, ], sub[as.vector(sapply(which.diff, function(z) { c(z, z+1)})), ])
			res <- rbind(res, sub[nrow(sub), ])
		}
		return(res)
	}) 
	result <- do.call(rbind, result)
	rownames(result) <- 1:nrow(result)
	return(result)
}