# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Feb 15, 2014


#63540 maize transcripts in 24 diff conditions ( 4 time points x 3 reps x ( 1 expt + 1 control )).

#shell
#grep -v '^jgi' counts.txt > maize_counts.txt
#cut --complement -f 2,3,19,20,22  maize_counts.txt > maize_counts_noaxenic.txt


mzCounts <- read.delim("maize_counts_noaxenic_noheader.txt",header=FALSE,row.names=1)
hist(floor(colSums(mzCounts)), labels=TRUE,xlab="Reads in sample",
     main="Histogram of reads mapped per sample")
coltf = floor(colSums(mzCounts) /1e06) #mill reads per lib, could have used cpm()
hist(coltf, labels=TRUE,xlab="Million reads in sample",
     main="Histogram of reads mapped per sample")

#boxplot of all counts, las just print X axis labels vertically
boxplot(mzCounts,main="Spread of all counts per library - maize", las=2)

collabels = c("3d1", "3d2", "7d1", "7d2","3dctrl1", "10dctrl1", "3dctrl2", 
              "10dctrl2", "3d3", "10d1", "5dctrl1", "3dctrl3", "10dctrl3", 
              "10d2", "10d3", "7d3", "5d1", "5d2", "5dctrl2", "5dctrl3", 
              "7dctrl1", "7dctrl2", "7dctrl3", "5d3")
colnames(mzCounts) = collabels
collabels_ord = c("3d1", "3d2", "3d3", "3dctrl1", "3dctrl2" ,"3dctrl3", 
                  "5d1", "5d2", "5d3", "5dctrl1", "5dctrl2" ,"5dctrl3", 
                  "7d1", "7d2", "7d3", "7dctrl1", "7dctrl2" ,"7dctrl3", 
                  "10d1", "10d2", "10d3", "10dctrl1", "10dctrl2" ,"10dctrl3")


#filter out low counts
#Initially removed genes with < 1 read/million in ANY library
# dim(mzCounts[rowSums(cpm(mzCounts) >1) ==24,])
# [1] 17357    24
# mzCountsCln = subset(mzCounts,V2>coltf[1] & V3>coltf[2] & V4>coltf[3] & V5>coltf[4] 
#                      & V6>coltf[5] & V7>coltf[6] & V8>coltf[7] & V9>coltf[8]& V10>coltf[9]& 
#                        V11>coltf[10]& V12>coltf[11]& V13>coltf[12]& V14>coltf[13]& 
#                        V15>coltf[14]& V16>coltf[15]& V17>coltf[16]& V18>coltf[17]& 
#                        V19>coltf[18]& V20>coltf[19]& V21>coltf[20]& V22>coltf[21]& 
#                        V23>coltf[22]& V24>coltf[23]& V25>coltf[24])
# colnames(mzCountsCln) = collabels

#tutorial method removes only genes with < 1 read/million in 3 or less libs. 
#3 being the (min) number of replicates
# > dim(mzCounts[rowSums(cpm(mzCounts) >1) >=3,])
# [1] 27200    24
mzCpms = cpm(mzCounts)
mzCountsCln = mzCounts[rowSums(cpm(mzCounts) >1) >=3,] #get subset
boxplot(mzCountsCln,main="Spread of filtered counts per library - maize", las=2)

#correlation matrix
library("lattice")
levelplot(cor(mzCounts[,collabels_ord]), aspect="iso", scales=list(x=list(rot=90)),main="Correlation Matrix", cuts=50)
library("corrplot")
corrplot(cor(mzCounts[,collabels_ord]), method="square", tl.col="black", 
         addgrid.col="black", is.corr=FALSE, main="Maize raw counts")
corrplot(cor(mzCountsCln[,collabels_ord]), method="square", tl.col="black", 
         addgrid.col="black", is.corr=FALSE, main="Maize filtered counts")


groups<- c("3d", "3d", "7d", "7d", "3dctrl", "10dctrl", "3dctrl", "10dctrl", "3d", "10d", "5dctrl", "3dctrl", "10dctrl", "10d", "10d", "7d", "5d", "5d", "5dctrl", "5dctrl", "7dctr", "7dctrl", "7dctrl", "5d")
hist(colSums(mzCountsCln), labels=TRUE,xlab="Million reads in sample",
     main="Histogram of reads mapped per sample - cleaned")
barplot(colSums(mzCountsCln),names.arg=groups, xlab="Library name", ylab="Read count",
        las=2,col="yellow",main="Number of million reads mapped per sample - cleaned")

library("edgeR")
dge <- DGEList(count=mzCountsCln,group=groups)
dge <- calcNormFactors(dge) # normalize libs to prevent over expressed genes from blanking out rest
barplot(dge$counts,las=2,main="Normalized total counts per library - maize")

#http://cgrlucb.wikispaces.com/edgeR+spring2013
library(RColorBrewer)
colors <- brewer.pal(9, "Set1")
scale = dge$samples$lib.size * dge$samples$norm.factors
normCounts = round(t(t(dge$counts)/scale) * mean(scale))
boxplot(log2(normCounts+1), las=2, col=colors[dge$samples$group],
        main="Spread of normalized counts per library - maize")

plotMDS(dge,labels=groups,col=c("darkblue","lightblue","darkgreen","lightgreen","darkred","red","black","gray")[factor(groups)], main="PCA of maize RNAseq") #diff colors for each group

dge <- estimateCommonDisp(dge) # Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
dge <- estimateTagwiseDisp(dge) # Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
plotMeanVar(dge,show.tagwise.vars=TRUE,NBline=TRUE,
            main="Variance vs Expr for maize transcripts") #each dot represents the estimated mean and variance for each gene, with binned variances as well as the trended common dispersion overlaid. Explore the mean-variance relationship for DGE data
plotBCV(dge,main="maize BCV vs Counts per Million") #Plots genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).

dge_3 = exactTest(dge,pair=c('3d','3dctrl'))# Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
dge_5 = exactTest(dge,pair=c('5d','5dctrl'))
dge_7 = exactTest(dge,pair=c('7d','7dctrl'))
dge_10 = exactTest(dge,pair=c('10d','10dctrl'))

dge_3_tt = topTags(dge_3,n=nrow(dge_3)) #Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change
# > dim(dge_3_tt[dge_3_tt$table$FDR < .05,])
# [1] 12530     4

dge_5_tt = topTags(dge_5,n=nrow(dge_5))
dge_7_tt = topTags(dge_7,n=nrow(dge_7))
dge_10_tt = topTags(dge_10,n=nrow(dge_10))

dge_nc <- cpm(dge,normalized.lib.sizes=TRUE) #Computes counts per million (CPM)
head(dge_nc[rownames(dge_3_tt$table),order(groups)],5) # depth-adjusted reads per million for some of the top 5 differentially expressed genes

plotSmear(dge,de.tags=rownames(dge_3_tt$table$FDR < .05),main="3days, log-Fold Change vs log-Conc, FDR < .05")
plotSmear(dge,de.tags=rownames(dge_5_tt$table$FDR < .05),main="5days, log-Fold Change vs log-Conc, FDR < .05")
plotSmear(dge,de.tags=rownames(dge_7_tt$table$FDR < .05),main="7days, log-Fold Change vs log-Conc, FDR < .05")
plotSmear(dge,de.tags=rownames(dge_10_tt$table$FDR < .05),main="10days, log-Fold Change vs log-Conc, FDR < .05")

write.csv(dge_3_tt$table, file="edgeR_maize_3d_alltags.csv")
write.csv(dge_5_tt$table, file="edgeR_maize_5d_alltags.csv")
write.csv(dge_7_tt$table, file="edgeR_maize_7d_alltags.csv")
write.csv(dge_10_tt$table, file="edgeR_maize_10d_alltags.csv")

#TODO for prez
hist(dge_3_tt$table$FDR)
hist(dge_3_tt$table$FDR[dge_3_tt$table$FDR < .05])#only , .05
dge_3_tt_matrix = as.matrix(dge_3_tt$table)
rc <- rainbow(nrow(dge_3_tt_matrix), start = 0, end = .3)
cc <- rainbow(ncol(dge_3_tt_matrix), start = 0, end = .3)
#poor viz, first try
hv <- heatmap(dge_3_tt_matrix, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = " variables", ylab =  "Transcripts",
              main = "heatmap of all transcripts")
heatmap(dge_3_tt[dge_3_tt$table$FDR < .05,], col = cm.colors(256), scale = "column",
        RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
        xlab = " variables", ylab =  "Transcripts",
        main = "heatmap of all transcripts")

library("gplots")
heatmap.2(dge_3_tt[dge_3_tt$table$FDR < .05,]) #better, should try for pval < .05 only

