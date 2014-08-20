# PPath@Cornell
# Surya Saha ss2489 at cornell dot edu
# Mar 21, 2014

#11702 fungal transcripts in 4 timepoints x 3 reps and 1 axenic x 3 reps

#shell
#cut -f2-3 fungi_counts.txt --complement > fungi_counts_noheader_noLenDesc.txt

counts <- read.delim("fungi_counts_noheader_noLenDesc.txt",header=TRUE,row.names=1)
collabels = c("3d1", "3d2", "7d1", "7d2","3dctrl1", "10dctrl1", "3dctrl2", 
              "10dctrl2", "3d3", "10d1", "5dctrl1", "3dctrl3", "10dctrl3", "10d2", "10d3", "ax1","ax2",
              "7d3", "ax3","5d1", "5d2", "5dctrl2", "5dctrl3", "7dctrl1", "7dctrl2", "7dctrl3", "5d3")
colnames(counts) = collabels
collabels_ord = c("3d1", "3d2", "3d3", "3dctrl1", "3dctrl2" ,"3dctrl3", 
                  "5d1", "5d2", "5d3", "5dctrl1", "5dctrl2" ,"5dctrl3", 
                  "7d1", "7d2", "7d3", "7dctrl1", "7dctrl2" ,"7dctrl3", 
                  "10d1", "10d2", "10d3", "10dctrl1", "10dctrl2" ,"10dctrl3", 
                  "ax1", "ax2", "ax3")

library("lattice")
levelplot(cor(counts[,collabels_ord]), aspect="iso", scales=list(x=list(rot=90)),main="Correlation Matrix", cuts=50)

library("corrplot")
corrplot(cor(counts[,collabels_ord]), method="square", tl.col="black", addgrid.col="black", is.corr=FALSE, main="St28a raw counts")
corrplot(cor(counts[,c("3d1","3d2","3d3","3dctrl1","3dctrl2","3dctrl3")]), method="square", tl.col="black", addgrid.col="black", is.corr=FALSE)
