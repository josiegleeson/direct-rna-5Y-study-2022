# Performing differential expression analysis with DESeq2 (v1.26.0)
# Ran separately on 5Y and sequin gene counts (from featureCounts)

library(DESeq2)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(pcaExplorer)
library(ggplot2)

#Import data
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

#remove extra columns of data from featureCounts
countdata <- countdata[ ,6:ncol(countdata)]

colnames(countdata) <- c("diff2", "diff3", "diff4", "undiff3", "undiff4")
countdata <- as.matrix(countdata)

# Assign conditions
(condition <- factor(c(rep("Differentiated", 3), rep("Undifferentiated", 2))))

# Create coldata frame and make DESeq Dataset
coldata <- data.frame(row.names=colnames(countdata), condition)

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# Optional filtering step to remove very low counts (less than 5 shown here)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Run DESeq2 pipeline and get results
dds <- DESeq(dds)
res <- results(dds)

# Save normalised counts separately if desired
norm_counts <- as.data.frame(counts(dds, normalized=TRUE))
write.csv(norm_counts, file="norm_counts.txt")

## Order by adjusted p-value
table(res$padj<0.05)
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

# Write results
write.csv(resdata, file="diffexpr-results.csv")

# Regularised log transformation for visualisation
rld <- rlogTransformation(dds)

# Colours
mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]

# Heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# PCA
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomright", main="Principal Component Analysis", textcx=1, ...) {
      require(genefilter)
      require(calibrate)
      require(RColorBrewer)
      rv = rowVars(assay(rld))
      select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
      pca = prcomp(t(assay(rld)[select, ]))
      fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
    if (is.null(colors)) {
      if (nlevels(fac) >= 3) {
        colors = brewer.pal(nlevels(fac), "Paired")
      }   else {
        colors = c("black", "red")
      }
    }
    pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
    pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
    pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
    pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
    plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
    with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
    legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
png("pca-1-2.png", 1000, 1000, pointsize=30)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-10, 15), ylim=c(-10, 10))
dev.off()

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=FALSE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="blue", pch=20, cex=1))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

# Volcano Plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot-nolabel.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=2, sigthresh=0.05, textcx=0.8, xlim=c(-10, 10), legendpos="topright")
dev.off()
