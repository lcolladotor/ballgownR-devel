#' Interactive visualization of the ballgown results
#'
#' \code{viewGene} plots the information extracted with \code{infoGene} in an interactive html
#'
#' @param geneInfo is the output of \code{infoGene}
#' @param geneID specifies the geneID to look at
#' @param exon.color is the color used for the exon lines
#' @param location has to be 'bottom' or 'top' and it controls where the exon tracks are showed
#' @param spacing is the factor by which the exons are separated
#' @param html specifies the html output file
#' @param wdir specifies the directory where the plot will be made. You must have writing permission.
#' @param browse=TRUE if you want to open a browser window to view the graphic
#' @param title specifies the title in the HTML file
#' @param main specifies the title of the plot in the HTML file
#' @return the html file name with the resulting plot
#' @export
#' @author Leonardo Collado-Torres \email{lcollado@@jhsph.edu}
#' @examples
#' dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
#' samplePattern <- "sample"
#' gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)
#' info <- system.file("extdata", "sample_info.txt", package="ballgownR")
#' info <- read.table(info, header=TRUE)
#' match <- sapply(names(gown$dirs), function(x) { which(info$dir == x)})
#' geneInfo <- infoGene(geneID="gene_1", gown=gown, group=info$outcome[match])
#' viewGene(geneInfo=geneInfo, html="toy-test.html", spacing=0.1)
#' viewGene(geneInfo=geneInfo, html="toy-test-2.html", spacing=0.1, location="top")
#' @seealso infoGene

viewGene <- function(geneInfo, exon.color="#000000", location="bottom", spacing=0.02, html=NULL, wdir=NULL, browse=TRUE, title=NULL, main=NULL) {
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
		title <- NULL
		main <- NULL
	}
	
	## Assign information form geneInfo
	geneID <- geneInfo$geneID
	group <- geneInfo$group
	exons.df <- geneInfo$exons
	toAdd <- geneInfo$transInfo
	start <- geneInfo$start
	end <- geneInfo$end
	ylab <- geneInfo$ylab
	
	## Set working directory
	if(is.null(wdir)) wdir <- getwd()
	
	## Set default html file name
	if(is.null(html)) html <- paste0(runif(1, max=1000000), "-", geneID, ".html")
		
	## Set default title
	if(is.null(title)) title <- paste("Viewing gene", geneID, "with ballgownR")
	
	## Set default main
	if(is.null(main)) main <- paste("Viewing", ylab, "for gene", geneID, "via <a href='https://github.com/lcolladotor/ballgownR-devel/tree/master/ballgownR'>ballgownR</a>")

	## Copy layout to working directory
	if(!"genome_info" %in% dir()) {
		layout.dir <- system.file(package="ballgownR", "extdata", "genome_info")
		## TODO: make this work on Windows
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
		exons.df$y <- (-1) * exons.df$y * diff(range(toAdd$y)) * spacing
		border.y <- min(c(exons.df$y, toAdd$y)) * 1.1
	} else {
		exons.df$y <- max(toAdd$y) * 1.1 + exons.df$y * diff(range(toAdd$y)) * spacing
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
	
	## Color labels
	group.col.label <- sapply(1:length(unique(group)), function(x) {
		paste0("<font color='", group.col[x], "'>", unique(group)[x], "</font>")
	})
	group.col.label <- paste0("<p>Color code: ", paste(group.col.label, collapse=", "), ".</p>")
	
		
	## Visualize
	result <- clickme(data, "genome_info", params=list(color=colors, title=title, main=main, groupLabel=group.col.label), html_file_name=html, browse=browse)
	## TO FIX: Error: Input to str_c should be atomic vectors
	## Might be due to a newer clickme version
	
	## Done
	return(invisible(result))
	
}