# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Feb 26, 2014

#11702 fungal transcripts in 4 timepoints x 3 reps and 1 axenic x 3 reps

#shell
#head -n 1 ../counts.txt > fungi_counts_header.txt
#grep '^jgi' ../counts.txt >> fungi_counts_header.txt
#cut --complement -f 2,3,8-11,14-16,25-29 fungi_counts_header.txt > fungi_counts_header_nocontrol.txt

#data
fnCounts <- read.delim("fungi_counts_header_nocontrol.txt",header=TRUE,row.names=1)
design = data.frame(row.names = colnames (fnCounts), 
                    condition = c("3d", "3d", "7d", "7d", "3d", "10d", "10d", "10d", "axenic", "axenic","7d", "axenic", "5d", "5d", "5d"),
                    libType = c("paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired"))

#filtering: DESeq suggests using a quantile method. Instead I'll use the method used with edgeR for comparison purposes.
#removing genes with < 1 read/million reads in every library, 15 libs
library("edgeR")
keep = rowSums ( cpm(fnCounts)>1) >= 15
dim(fnCounts[keep, ]) 
#[1] 3597   15
fnCountsFl = fnCounts[keep, ]

#quantile counts just to compare, look
rs <- rowSums(fnCounts)
theta <- 0.5
use <- rs >quantile(rs, probs = theta)
table(use)
#use
#FALSE  TRUE 
#5852  5850  

#DEseq
library("DESeq")
cds = newCountDataSet(fnCountsFl,design)
cds = estimateSizeFactors (cds) #est normalization factors
sizeFactors(cds) #print norm factor for each lib, large ones are axenic
#OYNG        OYNH        OYNN        OYNO        OYNW        OYNX        OYOB 
#0.15877283  0.04864184  0.25830195  0.38721655  0.08752413  3.44974220  1.76835773 
#OYOC        OYOG        OYOH        OYON        OYOO        OYOP        OYOS 
#1.57892713 38.23951728 68.24167231  1.90249824 47.37053757  0.18604037  0.46115122 
#OYOZ 
#0.11882388

#For PCA plot
cdsFullBlind = estimateDispersions(cds, method = "blind")
#Error in parametricDispersionFit(means, disps) : 
#  Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')
cdsFullBlindLocal = estimateDispersions(cds, method = "blind",fitType="local") #Used the locfit package to fit a dispersion-mean relation due to error
#https://stat.ethz.ch/pipermail/bioconductor/2012-March/044230.html
cdsFullBlindLocalfitonly = estimateDispersions(cds, method = "blind",sharingMode="fit-only",fitType="local")

vsdFullLocal = varianceStabilizingTransformation (cdsFullBlindLocal)
vsdFullLocalfitonly = varianceStabilizingTransformation (cdsFullBlindLocalfitonly)

plotPCA(vsdFullLocal,intgroup=c("condition","libType"))
#same plot, no benefit using sharingMode param
plotPCA(vsdFullLocalfitonly,intgroup=c("condition","libType"))

#For Est Disps
cds = estimateDispersions(cds)
#OLD ERROR due to missing meta-data in design
#Error in if (nr < 2) stop("nrow(modelMatrix) must be >=2.") : 
#  argument is of length zero
plotDispEsts(cds)

#heatmap of indiv genes
library("RColorBrewer")
library ("gplots")
select = order(rowMeans(counts(cds)), decreasing = TRUE)[1:3000]  
hmcol = colorRampPalette(brewer.pal(8, "GnBu"))(100)
heatmap.2(exprs(vsdFullLocal)[select,], col = hmcol, trace="none", margin=c(10,6),labRow=NA)
heatmap.2(exprs(vsdFullLocal)[select,], col = redblue(75), trace="none", margin=c(10,6),labRow=NA,labCol=design$condition)

#heatmap of dataset to dataset
dists = dist( t( exprs(vsdFullLocal) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cdsFullBlindLocal), paste(rownames(design), design$condition, sep=" : "))
heatmap.2(mat, trace="none", col = rev(redblue(75)), margin=c(13, 13))

#generate DE tables (comparison is expt to axenic)
res10toax = nbinomTest (cds, "10d", "axenic")
plotMA(res10toax)
#Error in MA[, array] - x : non-numeric argument to binary operator
#In addition: There were 50 or more warnings (use warnings() to see the first 50)
#> warnings()
#Warning messages:
#  1: In mean.default(sort(x, partial = half + 0L:1L)[half +  ... :
#         argument is not numeric or logical: returning NA
#  2: In mean.default(sort(x, partial = half + 0L:1L)[half +  ... :
#         argument is not numeric or logical: returning NA
res7toax = nbinomTest (cds, "7d", "axenic")
res5toax = nbinomTest (cds, "5d", "axenic")
res3toax = nbinomTest (cds, "3d", "axenic")
write.csv(res10toax,"fungi_10d_alltags_deseq.csv")
write.csv(res7toax,"fungi_7d_alltags_deseq.csv")
write.csv(res5toax,"fungi_5d_alltags_deseq.csv")
write.csv(res3toax,"fungi_3d_alltags_deseq.csv")

#shell
#awk -F"," '{if($4 < .05) print $1}' fungi_5d_alltags_edgeR.csv | sed 's,\",,g'>edgeR
#awk -F"," '{if($9 < .05) print $2}' fungi_5d_alltags_deseq.csv| sed 's,\",,g' > adj
#awk -F"," '{if($8 < .05) print $2}' fungi_5d_alltags_deseq.csv | sed 's,\",,g' > raw
