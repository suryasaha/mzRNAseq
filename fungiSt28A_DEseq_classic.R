# PPath@Cornell
# Surya Saha ss2489@cornell.edu
# Feb 26, 2014

#11702 fungal transcripts in 4 timepoints x 3 reps and 1 axenic x 3 reps

#shell
#head -n 1 ../counts.txt > fungi_counts_header.txt
#grep '^jgi' ../counts.txt >> fungi_counts_header.txt
#cut --complement -f 2,3,8-11,14-16,25-29 fungi_counts_header.txt > fungi_counts_header_nocontrol.txt

#data
fnCounts <- read.delim("fungi_counts_header_nocontrol.txt",header=FALSE,row.names=1)
design = data.frame(row.names = colnames (fnCounts), condition = c("3d", "3d", "7d", "7d", "3d", "10d", "10d", "10d", "axenic", "axenic","7d", "axenic", "5d", "5d", "5d"))

