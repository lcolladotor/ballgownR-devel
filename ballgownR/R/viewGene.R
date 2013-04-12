#' Interactive visualization of the ballgown results
#'
#' \code{viewGene} plots the information extracted with \code{infoGene} in an interactive html
#'
#' @param geneID specifies the geneID to look at
#' @param exon.color is the color used for the exon lines
#' @param location has to be 'bottom' or 'top' and it controls where the exon tracks are showed
#' @param spacing is the factor by which the exons are separated
#' @param html specifies the html output file
#' @param wdir specifies the directory where the plot will be made. You must have writing permission.
#' @param browse=TRUE if you want to open a browser window to view the graphic
#' @return the html file name with the resulting plot
#' @export
#' @author Leonardo Collado-Torres \email{lcollado@@jhsph.edu}
#' @examples
#' ?viewGene # Read the help. Example to do!

viewGene <- function(geneInfo, exon.color="#000000", location="bottom", spacing=0.02, html=NULL, wdir=NULL, browse=TRUE) {
	## Load required libraries
	require(clickme)
	require(colorspace)
	
	## for testing
	if(FALSE) {
		exon.color="#000000"
		location="bottom"
		spacing=0.02
		html=NULL
		wdir=NULL
		browse=TRUE
	}
	
	## Assign information form geneInfo
	geneID <- geneInfo$geneID
	group <- geneInfo$group
	exons.df <- geneInfo$exons
	toAdd <- geneInfo$transInfo
	start <- geneInfo$start
	end <- geneInfo$end
	
	## Set working directory
	if(is.null(wdir)) wdir <- getwd()
	
	## Set default html file name
	if(is.null(html)) html <- paste0(runif(1, max=1000000), "-", geneID, ".html")

	## Copy layout to working directory
	if(!"genome_info" %in% dir()) {
		layout.dir <- system.file(package="ballgownR", "extdata", "genome_info")
		system(paste0("cp -r ", layout.dir, " ", wdir, "/"))
	}
	
	## Change clickme's path
	set_root_path(wdir)
	
	## Contruct the colors for the exons
	exons.col <- rep(exon.color, length(unique(exons.df$line)))
	
	## Define sample colors
	group.col <- rainbow_hcl(length(unique(group)), start = 30, end = 300)
	sample.col <- rep(NA, length(group))
	for(j in 1:length(unique(group))){
		k <- unique(group)[j]
		sample.col[group %in% k] <- group.col[k]
	}
	
	## Adjust where to show the exons and it's spacing
	exons.df$y <- exons.df$y * diff(range(toAdd$y)) * spacing
	if(location == "bottom") {
		exons.df$y <- (-1) * exons.df$y - diff(range(toAdd$y)) * spacing
		border.y <- min(c(exons.df$y, toAdd$y)) * 1.1
	} else {
		exons.df$y <- exons.df$y + diff(range(toAdd$y)) * spacing
		border.y <- max(c(exons.df$y, toAdd$y)) * 1.1
	}
	
	## Set a border	to help visualize the exons
	border <- data.frame(line=c("border", "border"), x=c(start, end), y=c(border.y, border.y))
	
	## Merge data
	data <- rbind(exons.df, toAdd, border)
	
	## Format the colors
	colors <- c(exons.col, sample.col, exons.col[1])
	## Format for js
	colors <- paste0('["', paste0(colors, collapse='","'), '"]')
	
		
	## Visualize
	result <- clickme(data, "genome_info", params=list(color=colors), html_file_name=html, browse=browse)
	## TO FIX: Error: Input to str_c should be atomic vectors
	
	## Done
	return(invisible(result))
	
}