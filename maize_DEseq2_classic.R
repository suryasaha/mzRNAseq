# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Adapted from Tyr's script https://github.com/nelsonlab/rnaseq/blob/master/deseq2_code_tyr.R
# Aug 19, 2014

#63540 maize transcripts in 24 diff conditions ( 4 time points x 3 reps x ( 1 expt + 1 control )).

#shell
#grep -v '^jgi' counts.txt > maize_counts.txt
#cut --complement -f 2,3,19,20,22  maize_counts.txt > maize_counts_noaxenic.txt

#################################################################
#data
#################################################################
mzCounts <- read.delim("maize_counts_noaxenic.txt",header=TRUE,row.names=1)

#read in metadata about conditions
metadata <- read.delim("metadata.txt", row.names = 1)

# drop unused levels (axenic samples) in the metadata dataframe
mzMeta <- droplevels(metadata[row.names(metadata) %in% names(mzCounts),])
dim(metadata)
# [1] 27  4
dim(mzMeta)
# [1] 24  4

###### Filtering ########

#removing genes with < 1 read/million reads in every library, 24 libs
library("edgeR")
keepcpm = rowSums (cpm(mzCounts)>1) >= 24
dim(mzCounts[keepcpm, ])
#[1] 17357 24
#filter out low counts
mzCountsFlCPM = mzCounts[keepcpm, ]

#quantile counts
rs <- rowSums(mzCounts)
#remove the genes in the lowest 50% quantile (as indicated by the parameter theta)
theta <- 0.5
keeptheta <- rs >quantile(rs, probs = theta)
table(keeptheta)
#use
#FALSE TRUE
#31778 31762
mzCountsFlTheta = mzCounts[keeptheta,]

#################################################################
#DESeq2
#################################################################
library("DESeq2") #version 1.4.5

ddsmzCountFlCPM = DESeqDataSetFromMatrix(
  countData = mzCountsFlCPM,
  colData = mzMeta,
  design = ~ condition
)

ddsmzCountFlCPM = DESeq(ddsmzCountFlCPM)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# Warning message:
#   In estimateDispersionsFit(object, fitType = fitType, quiet = quiet) :
#   the parametric fit of dispersion estimates over the mean of counts
# failed, which occurs when the trend is not well captured by the
# function y = a/x + b. A local regression fit is automatically performed,
# and the analysis can continue. You can specify fitType='local' or 'mean'
# to avoid this message if re-running the same data.
# When using local regression fit, the user should examine plotDispEsts(dds)
# to make sure the fitted line is not sharply curving up or down based on
# the position of individual points.
png("/home/surya/work/RJN_GT_Setur_Maize/R/maize_seseq2_plotDispEsts_ddsmzCountFlCPM.png")
plotDispEsts(ddsmzCountFlCPM, main = "plotDispEsts_ddsmzCountFlCPM_localfit")
dev.off()
resmzCountFlCPM10 = results(ddsmzCountFlCPM, contrast = c("condition", "10day", "10daycontrol"))
nrow(resmzCountFlCPM[which(resmzCountFlCPM$padj <= 0.05),])
# [1] 9744
resmzCountFlCPM7 = results(ddsmzCountFlCPM, contrast = c("condition", "7day", "7daycontrol"))
resmzCountFlCPM5 = results(ddsmzCountFlCPM, contrast = c("condition", "5day", "5daycontrol"))
resmzCountFlCPM3 = results(ddsmzCountFlCPM, contrast = c("condition", "3day", "3daycontrol"))
nrow(resmzCountFlCPM5[which(resmzCountFlCPM5$padj <= 0.05),])
# [1] 5397
# day 5 is again low, like in deseq
nrow(resmzCountFlCPM3[which(resmzCountFlCPM3$padj <= 0.05),])
# [1] 10086
nrow(resmzCountFlCPM7[which(resmzCountFlCPM7$padj <= 0.05),])
# [1] 9483

# MA plots
png("/home/surya/work/RJN_GT_Setur_Maize/R/maize_deseq2_plotMA_ddsmzCountFlCPM3.png")
plotMA(resmzCountFlCPM3, main = "maize_deseq2_plotMA_ddsmzCountFlCPM3", ylim = c(-2,2))
dev.off()
png("/home/surya/work/RJN_GT_Setur_Maize/R/maize_deseq2_plotMA_ddsmzCountFlCPM5.png")
plotMA(resmzCountFlCPM5, main = "maize_deseq2_plotMA_ddsmzCountFlCPM5", ylim = c(-2,2))
dev.off()


ddsmzCountFlTheta = DESeqDataSetFromMatrix(
  countData = mzCountsFlTheta,
  colData = mzMeta,
  design = ~ condition
)
ddsmzCountFlTheta = DESeq(ddsmzCountFlTheta)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
png("/home/surya/work/RJN_GT_Setur_Maize/R/maize_seseq2_plotDispEsts_ddsmzCountFlTheta.png")
plotDispEsts(ddsmzCountFlCPM, main = "plotDispEsts_ddsmzCountFlTheta_Parafit")
dev.off()
resmzCountFlTheta10 = results(ddsmzCountFlTheta, contrast = c("condition", "10day", "10daycontrol"))
nrow(resmzCountFlTheta10[which(resmzCountFlTheta10$padj <= 0.05),])
# [1] 14582
nrow(resmzCountFlTheta7[which(resmzCountFlTheta7$padj <= 0.05),])
resmzCountFlTheta7 = results(ddsmzCountFlTheta, contrast = c("condition", "7day", "7daycontrol"))
nrow(resmzCountFlTheta7[which(resmzCountFlTheta7$padj <= 0.05),])
# [1] 14258
resmzCountFlTheta5 = results(ddsmzCountFlTheta, contrast = c("condition", "5day", "5daycontrol"))
nrow(resmzCountFlTheta5[which(resmzCountFlTheta5$padj <= 0.05),])
# [1] 7097
# low count again!
resmzCountFlTheta3 = results(ddsmzCountFlTheta, contrast = c("condition", "3day", "3daycontrol"))
nrow(resmzCountFlTheta3[which(resmzCountFlTheta3$padj <= 0.05),])
# [1] 14678



