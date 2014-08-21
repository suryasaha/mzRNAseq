# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Aug 19, 2014

#63540 maize transcripts in 24 diff conditions ( 4 time points x 3 reps x ( 1 expt + 1 control )).

#shell
#grep -v '^jgi' counts.txt > maize_counts.txt
#cut --complement -f 2,3,19,20,22  maize_counts.txt > maize_counts_noaxenic.txt

#################################################################
#data
#################################################################

mzCounts <- read.delim("maize_counts_noaxenic.txt",header=TRUE,row.names=1)
# this design trips up nbiomTest, using simple condition vector
# design = data.frame(row.names = colnames (mzCounts), 
#                     condition = c("3d", "3d", "7d", "7d", "3dctrl", "10dctrl", "3dctrl", "10dctrl", "3d", "10d", "5dctrl", "3dctrl", "10dctrl", "10d", "10d", "7d", "5d", "5d", "5dctrl", "5dctrl", "7dctrl", "7dctrl", "7dctrl", "5d" ),
#                     libType = c("paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired"))
condition = c("3d", "3d", "7d", "7d", "3dctrl", "10dctrl", "3dctrl", "10dctrl", "3d", "10d", "5dctrl", "3dctrl", "10dctrl", "10d", "10d", "7d", "5d", "5d", "5dctrl", "5dctrl", "7dctrl", "7dctrl", "7dctrl", "5d" )


########### filtering ##############
#removing genes with < 1 read/million reads in every library, 24 libs
library("edgeR")
keepcpm = rowSums (cpm(mzCounts)>1) >= 24
dim(mzCounts[keepcpm, ]) 
#[1] 17357    24
#filter out low counts
mzCountsFlCPM = mzCounts[keepcpm, ]

#quantile counts 
rs <- rowSums(mzCounts)
#remove the genes in the lowest 50% quantile (as indicated by the parameter theta)
theta <- 0.5
keeptheta <- rs >quantile(rs, probs = theta)
table(keeptheta)
#use
#FALSE  TRUE 
#31778 31762 
mzCountsFlTheta = mzCounts[keeptheta,]

#################################################################
#DESeq 
#################################################################
library("DESeq") #version 1.16.0

# Trying 2 solutions
#   filtered using CPM>1, blind method, local fit
#   filtered using theta=0.5, blind method, parametric fit

cdsCPM = newCountDataSet(mzCountsFlCPM,condition)
cdsCPM = estimateSizeFactors (cdsCPM) #est normalization factors
sizeFactors(cdsCPM)
# OYNG      OYNH      OYNN      OYNO      OYNP      OYNS      OYNT      OYNU      OYNW 
# 1.4492061 0.9644800 1.2003134 1.1093406 0.6829350 0.9734716 0.8719655 1.0909045 0.8871821 
# OYNX      OYNY      OYNZ      OYOA      OYOB      OYOC      OYON      OYOP      OYOS 
# 0.6131025 0.6627994 1.0035467 0.9227379 1.0604205 0.8729568 1.3505564 1.4293639 1.0865230 
# OYOT      OYOU      OYOW      OYOX      OYOY      OYOZ 
# 1.3513753 1.5297504 0.5773062 1.0842262 0.9711836 1.2050516

##### CPM>1 filtration, blind, paramentic fit ######
#cdsCPMFullBlind = estimateDispersions(cdsCPM, method = "blind")
# Error in parametricDispersionFit(means, disps) : 
#   Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')
# In addition: There were 29 warnings (use warnings() to see them)
# warnings()
# 28: step size truncated due to divergence
# 29: glm.fit: algorithm did not converge
#No error if unfiltered counts or relaxed filtering using theta (0.4 or 0.5) are used 
#for estimating dispersion. Removes genes in the lowest 40% or 50% quantile (as indicated 
#by the parameter theta). We can use diff starting data sets for edgeR and DEseq and 
#compare the DE list. Filtering is just a part of the algorithm, not the data itself. 
#=> USING THETA 0.5 also. See Expt log in Evernote for details

##### CPM>1 filtration, blind, local fit ######

#Using local fit instead of blind. Used locfit package to fit a dispersion-mean relation
#using sharingMode = "fit-only" to avoid nbiomTest warning (see Expt log)
cdsCPMBlindLocal = estimateDispersions(cdsCPM, method = "blind",sharingMode = "fit-only", fitType="local") 
str(fitInfo(cdsCPMBlindLocal))
# List of 5
# $ perGeneDispEsts: num [1:17357] 0.349 0.1807 0.0937 0.2456 0.0959 ...
# $ dispFunc       :function (q)  
#   ..- attr(*, "fitType")= chr "local"
# $ fittedDispEsts : num [1:17357] 0.196 0.273 0.175 0.113 0.176 ...
# $ df             : num 23
# $ sharingMode    : chr "fit-only"

head(fData(cdsCPMBlindLocal))
# disp_blind
# AC147602.5_FGT004  0.1960810
# AC148152.3_FGT008  0.2734398
# AC148167.6_FGT001  0.1752220
plotDispEsts(cdsCPMBlindLocal, main='Maize cdsCPM>1 Blind method Local fit Dispersions')

#For PCA plot
vsdCPMBlindLocal = varianceStabilizingTransformation (cdsCPMBlindLocal)
plotPCA(vsdCPMBlindLocal,intgroup=c("condition"))

#generate DE tables (comparison is expt to control)
res10 = nbinomTest (cdsCPMBlindLocal, "10d", "10dctrl")
res7 = nbinomTest (cdsCPMBlindLocal, "7d", "7dctrl")
res5 = nbinomTest (cdsCPMBlindLocal, "5d", "5dctrl")
res3 = nbinomTest (cdsCPMBlindLocal, "3d", "3dctrl")

#MA-plots, i.e. a scatter plot of logarithmic fold changes (on the y-axis) 
#versus the mean of normalized counts (on the x-axis).
plotMA(res10,main='Maize day 10 CPM>1 DiffExpr vs Expr Strength')
plotMA(res7,main='Maize day 7 CPM>1 DiffExpr vs Expr Strength')
plotMA(res5,main='Maize day 5 CPM>1 DiffExpr vs Expr Strength')
plotMA(res3,main='Maize day 3 CPM>1 DiffExpr vs Expr Strength')


dim(res10[which(res10$padj < 0.1),])
#[1] 3763    8
head(res10[order(res10$log2FoldChange, decreasing = TRUE),])

write.csv(res10,"maize_10d_alltags_deseq_cpm.csv")
write.csv(res7,"maize_7d_alltags_deseq_cpm.csv")
write.csv(res5,"maize_5d_alltags_deseq_cpm.csv")
write.csv(res3,"maize_3d_alltags_deseq_cpm.csv")


##### CPM>1 filtration, per-condition, parametric - Abandoned ####

#vsdCPMcondPerCondLocal = varianceStabilizingTransformation (cdsCPMcondPerCondLocal)
# Error in getVarianceStabilizedData(cds) : 
#   Use 'estimateDispersions' with 'method="blind"' (or "pooled") before calling 'getVarianceStabilizedData'
# Abandoned per-condition method as varianceStabilizingTransformation requires a blind estimate 
# of the variance function ( i.e., one ignoring conditions)


####### Theta 0.5 filtration, blind, paramentic fit #########

cdsTheta = newCountDataSet(mzCountsFlTheta,condition)
cdsTheta = estimateSizeFactors (cdsTheta) #est normalization factors
sizeFactors(cdsTheta)
# OYNG      OYNH      OYNN      OYNO      OYNP      OYNS      OYNT      OYNU      OYNW 
# 1.4289807 0.9731867 1.2119138 1.1110971 0.6860478 0.9893979 0.8677988 1.1138187 0.9014228 
# OYNX      OYNY      OYNZ      OYOA      OYOB      OYOC      OYON      OYOP      OYOS 
# 0.6280902 0.6674862 0.9946133 0.9304562 1.0580560 0.8776981 1.3523209 1.4490049 1.0833369 
# OYOT      OYOU      OYOW      OYOX      OYOY      OYOZ 
# 1.3847217 1.5727752 0.5850125 1.0949113 0.9989436 1.2217337

cdsThetaBlind = estimateDispersions(cdsTheta, method = "blind",sharingMode = "fit-only", fitType="local") 
str(fitInfo(cdsThetaBlind))
# $ perGeneDispEsts: num [1:31762] 0.3556 0.5104 0.1812 0.0946 0.2485 ...
# $ dispFunc       :function (q)  
#   ..- attr(*, "fitType")= chr "local"
# $ fittedDispEsts : num [1:31762] 0.234 0.29 0.277 0.225 0.281 ...
# $ df             : num 23
# $ sharingMode    : chr "fit-only"

head(fData(cdsThetaBlind))
# disp_blind
# AC147602.5_FGT004  0.2337891
# AC148152.3_FGT005  0.2903896
# AC148152.3_FGT008  0.2765086
plotDispEsts(cdsThetaBlind, main='Maize cdsTheta 0.5 Blind method Local fit Dispersions')

#For PCA plot
vsdThetaBlind = varianceStabilizingTransformation (cdsThetaBlind)
plotPCA(vsdThetaBlind,intgroup=c("condition"))

#generate DE tables (comparison is expt to control)
rest10 = nbinomTest (cdsThetaBlind, "10d", "10dctrl")
rest7 = nbinomTest (cdsThetaBlind, "7d", "7dctrl")
rest5 = nbinomTest (cdsThetaBlind, "5d", "5dctrl")
rest3 = nbinomTest (cdsThetaBlind, "3d", "3dctrl")

#MA-plots, i.e. a scatter plot of logarithmic fold changes (on the y-axis) 
#versus the mean of normalized counts (on the x-axis).
plotMA(rest10,main='Maize day 10 theta 0.5 DiffExpr vs Expr Strength')
plotMA(rest7,main='Maize day 7 theta 0.5 DiffExpr vs Expr Strength')
plotMA(rest5,main='Maize day 5 theta 0.5 DiffExpr vs Expr Strength')
plotMA(rest3,main='Maize day 3 theta 0.5 DiffExpr vs Expr Strength')


