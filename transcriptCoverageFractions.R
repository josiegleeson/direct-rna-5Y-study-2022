library(dplyr)
library(plotly)
library(GenomicFeatures)
library(MASS)
library(ggplot2)
library(viridis)

# Make database of known transcript lengths from gtf file
gencode <- makeTxDbFromGFF(file="gencode.v31.comprehensive.annotation.gtf", format="gtf", dataSource="gencode", organism="Homo sapiens")
gencodeLengths <- transcriptLengths(gencode, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
lengths <- data.frame(gencodeLengths$tx_name, gencodeLengths$tx_len)

# Import 5Y data from previous extractDataFromBam script 
all_5y <- read.csv(file = "all_5y.csv", header = TRUE, na.strings="#N/A")

# Merge files by transcript name for each length
merged <- merge(all_5y, lengths, by.x = "transcriptName", by.y = "gencodeLengths.tx_name", all.x = TRUE)
merged <- na.omit(merged)

# coverage = width / known length
merged <- merged %>% 
  mutate(coverages=width/gencodeLengths.tx_len)

# proportion aligned = aligned length / read length
merged <- merged %>% 
  mutate(proportionAligned=alignedLength/readLength)

# find the median coverage fraction for each unique transcript
medianmerged <- mergedex %>% 
  group_by(transcript) %>% 
  summarise_at(vars(coverages, gencodeLengths.tx_len, Counts), funs(median(., na.rm=TRUE)))

# merge again with average salmon counts per transcript
expression <- read_csv(file = "5ysalmon.csv")
expressionmerged <- merge(merged, expression, by.x = "transcript", by.y = "Name", all.x = TRUE)
expressionmerged <- expressionmerged[ !is.na(expressionmerged$Counts), ] 

# set intervals and factor
merged$gencodeLengths.tx_len <- cut(merged$gencodeLengths.tx_len, breaks=c(0,500,1000,1500,2000,200000), right=FALSE, dig.lab = 5)
merged$gencodeLengths.tx_len <- as.factor(merged$gencodeLengths.tx_len)

# violin plot


# boxplots with expression


# function for density scatterplot 
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# create df of lengths and coverage fractions
dat <- data.frame(
  x = merged$gencodeLengths.tx_len,
  y = merged$coverage
)

# call function on df
dat$density <- get_density(dat$x, dat$y, n = 300)

# plot and save
pdf("5y-density-all-reads.tiff", width = 8, height = 6)
ggplot(dat) + 
  geom_point(aes(x, y, color = density), alpha=0.2) +
  xlim(0,10000) +
  scale_color_viridis()
dev.off()
