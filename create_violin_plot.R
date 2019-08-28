library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(processx)

# Import csv from extract_data_bam script
df <- read.csv("diff2_output.csv")

# Find min and max transcript lengths
summary(df$width)

# Create column read length / region width to get proportion
df <- transform(df, fullLengthWidth = readLength / width)

# Set appropriate sequence for intervals
df$widthInterval <- cut(df$width, breaks=c(0,500,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,2000000), right=FALSE, dig.lab = 5)

# Check
summary(df$widthInterval)

# Convert intervals to factors
df$widthInterval <- as.factor(df$widthInterval)

# Produce violin plots
q <- df %>%
  plot_ly(
    x = ~widthInterval,
    y = ~fullLengthWidth,
    split = ~widthInterval,
    type = 'violin'
  ) %>% 
  layout(
    title = "Diff2",
    colorway='#FF9AA2',
    showlegend = FALSE,
    xaxis = list(
      title = "Maximum Transcript Length"
    ),
    yaxis = list(
      title = "Fraction of maximum transcript covered per read",
      zeroline = F,
      range=c(0,1)
    )
  )
q

# Save image using R Studio interface on the right