# PPath@Cornell
# Surya Saha ss2489@cornell.edu
# Feb 17, 2014

#11702 fungal transcripts in 4 timepoints x 3 reps and 1 axenic x 3 reps

#shell
#cut --complement -f 2,3,8-11,14-16,25-29 fungi_counts.txt > fungi_counts_nocontrol_noheader.txt 

fnCounts <- read.delim("fungi_counts_nocontrol_noheader.txt",header=FALSE,row.names=1)
hist(floor(colSums(fnCounts)), labels=TRUE,xlab="Reads in sample",main="Histogram of reads mapped per sample")
coltf = floor(colSums(fnCounts) /1e06) #mill reads per lib
hist(coltf, labels=TRUE,xlab="Million reads in sample",main="Histogram of reads mapped per sample")

#filter out low counts
#removing reads with < 1 read/million reads in library
fnCountsCln = subset(fnCounts,V2>coltf[1] & V3>coltf[2] & V4>coltf[3] & V5>coltf[4] & V6>coltf[5] & V7>coltf[6] & V8>coltf[7] & V9>coltf[8]& V10>coltf[9]& V11>coltf[10]& V12>coltf[11]& V13>coltf[12]& V14>coltf[13]& V15>coltf[14]& V16>coltf[15])
coltfcln = floor(colSums(fnCountsCln) /1e06) #mill reads per lib
hist(coltfcln, labels=TRUE,xlab="Million reads in sample",main="Histogram of reads mapped per sample after removing low counts")
sort(floor(colSums(fnCountsCln)))

#labels
groups<- c("3d", "3d", "7d", "7d", "3d", "10d", "10d", "10d", "axenic", "axenic","7d", "axenic", "5d", "5d", "5d")
barplot(coltfcln,names.arg=groups, xlab="Library name", ylab="Read count",las=2,col="yellow",main="Number of million reads mapped per sample")

#edgeR
dge <- DGEList(count=fnCountsCln,group=groups)
dge <- calcNormFactors(dge) # normalize libs to prevent over expressed genes from blanking out rest
barplot(dge$counts,names.arg=groups,las=2,main="Normalized fungal counts per library")
plotMDS(dge,labels=groups,col=c("darkblue","darkgreen","darkred","black","orange")[factor(groups)]) #diff colors for each group

dge <- estimateCommonDisp(dge) # Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags
dge <- estimateTagwiseDisp(dge) # Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
plotMeanVar(dge,show.tagwise.vars=TRUE,NBline=TRUE) #each dot represents the estimated mean and variance for each gene, with binned variances as well as the trended common dispersion overlaid. Explore the mean-variance relationship for DGE data
plotBCV(dge) #Biological Coefficient of Variation

dge_3 = exactTest(dge,pair=c('3d','axenic'))# Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts
#Warning message:
#In mglmOneGroup(y, dispersion = dispersion, offset = offset) :
#  max iteractions exceeded for 1 tags
dge_5 = exactTest(dge,pair=c('5d','axenic'))
dge_7 = exactTest(dge,pair=c('7d','axenic'))
dge_10 = exactTest(dge,pair=c('10d','axenic'))

dge_3_tt = topTags(dge_3,n=nrow(dge_3)) #Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change
dge_5_tt = topTags(dge_5,n=nrow(dge_5))
dge_7_tt = topTags(dge_7,n=nrow(dge_7))
dge_10_tt = topTags(dge_10,n=nrow(dge_10))

dge_nc <- cpm(dge,normalized.lib.sizes=TRUE) #Computes counts per million (CPM)
head(dge_nc[rownames(dge_3_tt$table),order(groups)],5) # depth-adjusted reads per million for some of the top 5 differentially expressed genes on day 3
head(dge_nc[rownames(dge_5_tt$table),order(groups)],5) # depth-adjusted reads per million for some of the top 5 differentially expressed genes on day 5

plotSmear(dge,de.tags=rownames(dge_3_tt$table$FDR < .05),main="Fungi 3days, log-Fold Change versus log-Concentration for FDR < .05")
plotSmear(dge,de.tags=rownames(dge_5_tt$table$FDR < .05),main="Fungi 5days, log-Fold Change versus log-Concentration for FDR < .05")
plotSmear(dge,de.tags=rownames(dge_7_tt$table$FDR < .05),main="Fungi 7days, log-Fold Change versus log-Concentration for FDR < .05")
plotSmear(dge,de.tags=rownames(dge_10_tt$table$FDR < .05),main="Fungi 10days, log-Fold Change versus log-Concentration for FDR < .05")

write.csv(dge_3_tt$table, file="fungi_3d_alltags_edgeR.csv")
write.csv(dge_5_tt$table, file="fungi_5d_alltags_edgeR.csv")
write.csv(dge_7_tt$table, file="fungi_7d_alltags_edgeR.csv")
write.csv(dge_10_tt$table, file="fungi_10d_alltags_edgeR.csv")

#shell
# awk -F"," '{if($4 < .05) print }' fungi_3d_alltags_edgeR.csv | wc -l
# awk -F"," '{if($4 < .05) print }' fungi_5d_alltags_edgeR.csv | wc -l
# awk -F"," '{if($4 < .05) print }' fungi_7d_alltags_edgeR.csv | wc -l
# awk -F"," '{if($4 < .05) print }' fungi_10d_alltags_edgeR.csv | wc -l
