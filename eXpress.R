# ---------------------------- #
# Author: DungChi Wu
# Date: Sat Dec  9 17:56:41 2017
# Description: 
# ---------------------------- #

# Following is the exercise: 
# You would first need to set up the working directory
setwd("C:\\Users\\Administrator\\Desktop\\R_hw")

install.packages("RColorBrewer")
install.packages("survival")

source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite("pasilla")

# Load necessary library
library("DESeq2")
library("RColorBrewer")
library("ggplot2")


# And then loading the expression value data.
# value data <-read.table()
value_data <-read.table("eXpress_value.tsv")

# Also the metadata for each experiment.
# metadata <-read.table()
metadata <-read.table("eXpress_passage.tsv")


# Setup some convient color variable for future used.
colors_rev <- colorRampPalette( rev(brewer.pal(9, "Blues")))( n = 299)
colors_for <- colorRampPalette( brewer.pal(9, "Blues"))(n = 299)
colors_rgb <- colorRampPalette( c("green", "black", "red"))(n = 299)


library("pheatmap")
sampleCor <- cor(log2(1+value_data), method = "pearson")
jpeg("pheatmap.jpg",height = 1024,width= 1024)
pheatmap(sampleCor,              # correlation matrix 
         col = colors_for,       # color used in heatmap
         annotation = metadata,   # show the information
         display_numbers = TRUE  # show the correlation value
)
dev.off()

# You may pick the two group to plot again




# Exercise ----------------------------------------------------------------
# Based on the above information, the two cell lines A549 and K562 are distinguishable from each other.
# Please select these experiment data and the related
# information from rawCts and rawInfo, and do the following:

# 1. Run the DESeq pipeline on these two cell lines, A549 and K562
#    e.g. using K562 as base level

attach(metadata)

dds <- DESeqDataSetFromMatrix(countData = round(value_data),
                              colData = metadata,
                              design = ~ Passage)
dds


dds$Passage <- relevel(dds$Passage, ref = "P38")
dds <- DESeq(dds)
res <- results(dds)
res


plotMA(res, ylim=c(-5, 5) )
hist(res$padj, breaks = 100, col = "skyblue", border = "slateblue", main = "")

# 2. Report the number of differentially expressed gene using following criteria: 
#    padj < 1e-5 AND log2FoldChange > 4 or < -4 and save it into a file.


resSig <- subset(res[order(res$padj),], padj < 1e-5)
resSig 

write.csv(as.data.frame(resSig), file="eXpress_passage_results.csv")


