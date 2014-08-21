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

#read in metadata about conditions
metadata <- read.delim("metadata.txt", row.names = 1)

# drop unused levels (axenic samples) in the metadata dataframe
mzMeta <- droplevels(metadata[row.names(metadata) %in% names(mzCounts),])
dim(metadata)
# [1] 27  4
dim(mzMeta)
# [1] 24  4
