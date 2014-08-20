# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Feb 17, 2014


mzErrCounts <- read.delim("maize_errors.txt",header=FALSE,row.names=1)
fgErrCounts <- read.delim("fungi_errors.txt",header=FALSE,row.names=1)

barplot(colSums(fnErrCounts), xlab="Library name", ylab="Read count",col="yellow",main="Number of maize reads mapped to St28A transcripts per sample")
barplot(colSums(fgErrCounts), names.arg=sort(fgErrCountGroups),xlab="Library name", ylab="Read count",col="yellow",main="Number of maize reads mapped to St28A transcripts per sample",las=2)

rawCounts = read.columns(file="raw_counts.txt",sep="\t")


#TODO fungi/maize counts per library or condition