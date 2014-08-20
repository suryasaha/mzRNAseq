# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Aug 19, 2014

#63540 maize transcripts in 24 diff conditions ( 4 time points x 3 reps x ( 1 expt + 1 control )).

#shell
#grep -v '^jgi' counts.txt > maize_counts.txt
#cut --complement -f 2,3,19,20,22  maize_counts.txt > maize_counts_noaxenic.txt

#data
mzCounts <- read.delim("maize_counts_noaxenic.txt",header=TRUE,row.names=1)
design = data.frame(row.names = colnames (mzCounts), 
                    condition = c("3d", "3d", "7d", "7d", "3dctrl", "10dctrl", "3dctrl", "10dctrl", "3d", "10d", "5dctrl", "3dctrl", "10dctrl", "10d", "10d", "7d", "5d", "5d", "5dctrl", "5dctrl", "7dctrl", "7dctrl", "7dctrl", "5d" ),
                    libType = c("paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired","paired"))


#filtering: DESeq suggests using a quantile method. Instead I'll use the method used with edgeR for comparison purposes.
#removing genes with < 1 read/million reads in every library, 24 libs
library("edgeR")
keepcpm = rowSums (cpm(mzCounts)>1) >= 24
dim(mzCounts[keepcpm, ]) 
#[1] 17357    24
#filter out low counts
mzCountsFlCPM = mzCounts[keepcpm, ]

#quantile counts just to compare, look
rs <- rowSums(mzCounts)
#remove the genes in the lowest 50% quantile (as indicated by the parameter theta)
theta <- 0.5
keeptheta <- rs >quantile(rs, probs = theta)
table(keeptheta)
#use
#FALSE  TRUE 
#31778 31762 
mzCountsFlTheta = mzCounts[keeptheta,]


