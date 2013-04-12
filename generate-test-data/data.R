### Ballgown data

setwd("../ballgownR/inst/extdata")

maindir <- "ballgownData"
sampledirs <- paste0("sample", 1:4)

e2t <- data.frame("e_id"=c(1:3, 1, 3, 1, 2, 2, 3, 4, 5), "t_id"=c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5))
i2t <- data.frame("i_id"=c(1:2, 3, 1, 2, 4), "t_id"=c(1, 1, 2, 3, 4, 5))

edata <- data.frame("e_id"=1:5, "chr"=rep(1, 5), "strand"=rep(c("+", "-"), c(3, 2)), "start"=c(1, 4000, 8000, 3000, 6000), "end"=c(2000, 6000, 12000, 5000, 8000), "rcount"=rep(0, 5), "ucount"=rep(0, 5), "mrcount"=rep(0, 5), "cov"=rep(0, 5), "cov_sd"=rep(0, 5), "mcov"=rep(0, 5), "mcov_sd"=rep(0, 5))

idata <- data.frame("i_id"=1:4, "chr"=rep(1, 4), "strand"=rep(c("+", "-"), c(3, 1)), "start"=c(2001, 6001, 2001, 5001), "end"=c(3999, 7999, 7999, 5999), "rcount"=rep(0, 4), ucount=rep(0, 4), mrcount=rep(0, 4))

tdata <- data.frame("t_id"=1:5, "chr"=rep(1, 5), "strand"=rep(c("+", "-"), c(4, 1)), "start"=rep(c(1, 4000, 3000), c(3, 1, 1)), "end"=rep(c(12000, 8000, 12000, 8000), c(2, 1, 1, 1)), "t_name"=paste0("trans_", 1:5), "num_exons"=rep(c(3, 2), c(1, 4)), "length"=c(8002, 6001, 4001, 6002, 5001), "gene_id"=paste0("gene_", rep(c(1, 2), c(4, 1))), "gene_name"=paste0("name_", rep(c(1, 2), c(4, 1))), "cov"=rep(0, 5), "FPKM"=rep(0, 5))

system(paste("mkdir -p", maindir))

## For now, just easy uniform data
populateData <- function(s) {
	switch(s, 
		"sample1" = 10,
		"sample2" = 20,
		"sample3" = 15,
		"sample4" = 25
	) 
}

for(s in sampledirs) {
	dir <- paste(maindir, s, sep="/")
	system(paste("mkdir -p", dir))
	
	## Write regular files
	write.table(e2t, file=paste(dir, "e2t.ctab", sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
	write.table(i2t, file=paste(dir, "i2t.ctab", sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
		
	## Populate data and write it
	enew <- edata
	inew <- idata
	tnew <- tdata
	
	enew$rcount <- enew$ucount <- enew$mrcount <- rep(populateData(s), nrow(enew))
	inew$rcount <- inew$ucount <- inew$mrcount <- rep(populateData(s), nrow(inew))
	tnew$cov <- tnew$FPKM <- rep(populateData(s), nrow(tnew))
	
	write.table(enew, file=paste(dir, "e_data.ctab", sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
	write.table(inew, file=paste(dir, "i_data.ctab", sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
	write.table(tnew, file=paste(dir, "t_data.ctab", sep="/"), quote=FALSE, sep="\t", row.names=FALSE)
}


### Sample info
info <- data.frame("sample"=paste0("sampleName", 1:4), dir=paste0("sample", 1:4), outcome=rep(c("control", "case"), each=2))
write.table(info, file="sample_info.txt", quote=FALSE, sep="\t", row.names=FALSE)



### Check that it works

dataDir <- system.file("extdata", "ballgownData", package="ballgownR")
samplePattern <- "sample"
gown <- readGown(dataDir=dataDir, samplePattern=samplePattern)