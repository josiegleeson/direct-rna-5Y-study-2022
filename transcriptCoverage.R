library(dplyr)
library(plotly)
library(GenomicFeatures)

# Make database of known transcript lengths from gencode

gencode <- makeTxDbFromGFF(file="gencode.v31.primary_assembly.annotation.gtf", format="gtf", dataSource="gencode", organism="Homo sapiens")
gencodeLengths <- transcriptLengths(gencode, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
lengths <- data.frame(gencodeLengths$tx_name, gencodeLengths$tx_len)

# Import 5Y data from previous extractDataFromBam script 
all_5y <- read.csv(file = "5y-all.csv", header = TRUE, na.strings="#N/A")

# Merge files by transcript name for each length
merged <- merge(all_5y, lengths, by.x = "transcriptName", by.y = "gencodeLengths.tx_name", all.x = TRUE)
merged <- na.omit(merged)

# coverage = width / known length
merged <- merged %>% 
  mutate(coverages=width/gencodeLengths.tx_len)

# proportion aligned = aligned length / read length
merged <- merged %>% 
  mutate(proportionAligned=alignedLength/readLength)

# set intervals and factor
merged$gencodeLengths.tx_len <- cut(merged$gencodeLengths.tx_len, breaks=c(0,500,1000,1500,2000,200000), right=FALSE, dig.lab = 5)
merged$gencodeLengths.tx_len <- as.factor(merged$gencodeLengths.tx_len)

# violin plot of coverage of annotated transcript grouped by length of annotated transcript
q <- merged %>%
  plot_ly(
    x = ~gencodeLengths.tx_len,
    y = ~coverages,
    split = ~gencodeLengths.tx_len,
    type = 'violin',
    points=FALSE,
    box = list(visible=T)
  ) %>% 
  layout(
    title = "Human Transcripts",
    colorway='#2980B9',
    showlegend = FALSE,
    xaxis = list(
      title = "Annotated Transcript Length"
    ),
    yaxis = list(
      title = "Coverage",
      zeroline = F,
      range=c(0,1)
    )
  )
q

# violin plot of proportion aligned
p <- merged %>%
  plot_ly(
    x = ~sample,
    y = ~proportionAligned,
    split = ~sample,
    type = 'violin',
    points=FALSE,
    box = list(visible=T)
  ) %>% 
  layout(
    title = "Proportion of Read Aligned",
    colorway='#2980B9',
    showlegend = FALSE,
    xaxis = list(
      title = "Sample"
    ),
    yaxis = list(
      title = "Proportion",
      zeroline = F,
      range=c(0.85,1)
    )
  )
p

# violin plot of read length vs aligned length 
j <- merged %>%
  plot_ly(type = 'violin') %>%
  add_trace(
    x = ~sample,
    y = ~readLength,
    legendgroup = 'Read Length',
    scalegroup = 'Read Length',
    name = 'Read Length',
    side = 'negative',
    points=FALSE,
    box=list(visible=T)
  ) %>%
  add_trace(
    x = ~sample,
    y = ~alignedLength,
    legendgroup = 'Aligned Length',
    scalegroup = 'Aligned Length',
    name = 'Aligned Length',
    side = 'positive',
    points=FALSE,
    box=list(visible=T)
  ) %>% 
  layout(
    xaxis=list(title=""),
    yaxis = list(title="Length",
      zeroline = F
    ),
    violingap=0,
    violingroupgap=0,
    violinmode = 'overlay'
  )
j

# scatter plot of coverage of annotated transcript
s <- merged %>%
  plot_ly(
    x = ~gencodeLengths.tx_len,
    y = ~coverages,
    type = 'scatter',
    visible=T,
    mode='markers'
  ) %>% 
  layout(
    title = "5Y Transcripts",
    colorway='#e08d3c',
    showlegend = FALSE,
    xaxis = list(
      title = "Annotated Transcript Length" 
    ),
    yaxis = list(
      title = "Coverage",
      zeroline = F,
      range=c(0,1),
      tickvals = list(0.2,0.4,0.6,0.8,1)
    )
  )
s 
