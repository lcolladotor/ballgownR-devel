#' Interactive visualization of the ballgown results: raw coverage or rcounts
#'
#' \code{viewGene} extracts the information according to a geneID and plots the coverage or the rcounts (intron and exon tables) in an interactive html
#'
#' @param geneID specifies the geneID to look at
#' @param gown specifies the output from \code{readGown}
#' @param coverage=FALSE determines whether to make a raw coverage plot or show the rcounts
#' @param tophatDir specifies the TopHat directory. Used when \code{coverage} is set to TRUE.
#' @param exon.color is the color used for the exon lines
#' @param location has to be 'bottom' or 'top' and it controls where the exon tracks are showed
#' @param spacing is the factor by which the exons are separated
#' @param html specifies the html output file
#' @param wdir specifies the directory where the plot will be made. You must have writing permission.
#' @return the html file name with the resulting plot
#' @export
#' @author Leonardo Collado-Torres \email{lcollado@@jhsph.edu}
#' @examples
#' ?viewGene # Read the help. Example to do!

viewGene <- function(geneID, gown, group, coverage=FALSE, tophatDir, exon.color="#000000", location="bottom", spacing=1, html=NULL, wdir=NULL) {
	## Load required libraries
	require(clickme)
	require(colorspace)
	suppressMessages(require(Rsamtools))
	
	## Set working directory
	if(is.null(wdir)) wdir <- getwd()
	
	## Set default html file name
	if(is.null(html)) html <- paste0(runif(1000000), "-", geneID, "-.html")

	## Copy layout to working directory
	if(!"genome_info" %in% dir()) {
		layout.dir <- system.file(package="ballgownR", "extdata", "genome_info")
		system(paste0("cp -r ", layout.dir, " ", wdir, "/"))
	}
	
	## Change clickme's path
	set_root_path(wdir)
	
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
	
	## Get exon information
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
	exons.col <- rep(exon.color, length(unique(exons.df$line)))
	
	## Define sample colors
	group.col <- rainbow_hcl(length(unique(group)), start = 30, end = 300)
	sample.col <- rep(NA, length(group))
	for(j in 1:length(unique(group))){
		k <- unique(group)[j]
		sample.col[group %in% k] <- group.col[k]
	}
	
	
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
		
	} else {
		
		## Initialize coverage data
		coverage.df <- data.frame(line=rep(samples, each=nBases), x=start:end, y=rep(0, each=nBases))
		
		## Define region to load
		which <- RangesList(IRanges(start=start, end=end))
		chr <- gown$trans$chr[idx.t[1]]
		names(which) <- chr
		param <- ScanBamParam(which=which)
		
		for(s in gown$dirs) {
			bamfile <- BamFile(paste(tophatDir, s, "accepted_hits.bam", sep="/"))
			# Read it and get the coverage. Extract only the one for the chr in question
			coverage <- coverage(readBamGappedAlignments(bamfile, param=param))[[chr]]
			coverage.df$y[ coverage.df$line == s] <- coverage
		}
		
		toAdd <- coverage.df
		
	}
	
	## Adjust where to show the exons and it's spacing
	exons.df$y <- exons.df$y * spacing
	if(location == "bottom") {
		exons.df$y <- (-1) * exons.df$y - diff(range(toAdd$y)) * 0.1
	} else {
		exons.df$y <- exons.df$y + diff(range(toAdd$y)) * 0.1
	}
	
	## Merge data
	data <- rbind(exons.df, toAdd)
	
	## Purge useless info
	filter <- function(data) {
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
	data.filt <- filter(data)
	
	
	## Format the colors
	colors <- c(exons.col, sample.col)
	## Format for js
	colors <- paste0('["', paste0(colors, collapse='","'), '"]')
	
		
	## Visualize
	clickme(data.filt, "genome_info", params=list(color=colors), html=html)
	
	## Done
	return(list(data=data.filt, color=colors))
	
}