# Adapted from Soneson et al. 2019 and added accuracy, aligned fractions, alignment scores, alignment type. Included making of 
# gencode TxDb and merging this with alignemnts.
# Included extra filters to select the 'true' transcript origin of each read based on alignemnt score, aligned fraction, accuracy
# coverage fraction of transcipt.

# Usage: Rscript extractDataFromBam.R yourfile.bam gencode.gtf ouputprefix

# GenomicAlignments/Features package is from bioconductor, need to install bioconductor then run:
# BiocManager::install("GenomicFeatures")

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  bamfile <- args[1]
  gencode <- args[2]
  output <- args[3]
  suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicFeatures)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  })
  
  bam <- readGAlignments(bamfile, use.names = TRUE,
                         param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                              what = c("qname","flag", "rname", 
                                                       "pos", "mapq")))
  
  ops <- GenomicAlignments::CIGAR_OPS
  wdths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = ops)
  keep.ops <- GenomicAlignments::explodeCigarOps(cigar(bam), ops = ops)
  explodedcigars <- IRanges::CharacterList(relist(paste0(unlist(wdths), 
                                                         unlist(keep.ops)), wdths))
  for (opts in setdiff(GenomicAlignments::CIGAR_OPS, "=")) {
    mcols(bam)[[paste0("nbr", opts)]] <- 
      sapply(explodedcigars, function(cg) sum(as.numeric(gsub(paste0(opts, "$"), "", cg)), na.rm = TRUE))
  }
  mcols(bam)$readLength <- rowSums(as.matrix(mcols(bam)[, c("nbrS", "nbrH", "nbrM", "nbrI")]))
  bam
  
  tmp <- data.frame(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::rename(read = qname,
                  nbrJunctions = njunc) %>%
    dplyr::select(-cigar) %>%
    dplyr::mutate(alignedLength = nbrM + nbrI) ## equivalent to readLength-nbrS-nbrH
  
  tmp2 <- as.data.frame(table(names(subset(bam, flag %in% c(0, 16)))))
  if (nrow(tmp2) == 0) tmp2 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp2 %>% dplyr::rename(read = Var1, nbrPrimaryAlignments = Freq))
  
  tmp3 <- as.data.frame(table(names(subset(bam, flag %in% c(256, 272)))))
  if (nrow(tmp3) == 0) tmp3 <- data.frame(Var1 = tmp$read[1], Freq = 0)
  tmp <- tmp %>% 
    dplyr::left_join(tmp3 %>% dplyr::rename(read = Var1, nbrSecondaryAlignments = Freq))
  
  tmp <- tmp %>% dplyr::mutate(nbrSecondaryAlignments = replace(nbrSecondaryAlignments, 
                                                                is.na(nbrSecondaryAlignments), 0))
  tmp <- tmp %>% dplyr::mutate(alignedFraction=alignedLength/readLength)
  
  tmp <- tmp %>% 
    dplyr::mutate(MID=nbrM+nbrI+nbrD)
  
  tmp <- tmp %>% 
    dplyr::mutate(NMID=nbrM+nbrI+nbrD-NM)
  
  tmp <- tmp %>% 
    dplyr::mutate(accuracy=NMID/MID)
  
  tmp$seqnames <- as.character(tmp$seqnames)
  
  if((grepl("\\|", tmp$seqnames)) == TRUE) {
    tmp <- data.frame(tmp, do.call(rbind, strsplit(tmp$seqnames, "\\|"))[,1])
    tmp$do.call.rbind..strsplit.tmp.seqnames............1. <- as.character(tmp$do.call.rbind..strsplit.tmp.seqnames............1.)
    tmp$transcript <- tmp$do.call.rbind..strsplit.tmp.seqnames............1.
    tmp$seqnames <- NULL
    tmp$rname <- NULL
    tmp$do.call.rbind..strsplit.tmp.seqnames............1. <- NULL
  } else {
    tmp$transcript <- tmp$seqnames
    tmp$seqnames <- NULL
  }
  
  # Make gencode database
  txs <- makeTxDbFromGFF(gencode, format="gtf", dataSource="gencode", organism="Homo sapiens")
  gencodeLengths <- transcriptLengths(txs, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  lengths <- data.frame(gencodeLengths$tx_name, gencodeLengths$tx_len)
  
  merged <- merge(tmp, lengths, by.x = "transcript", by.y = "gencodeLengths.tx_name", all.x = TRUE)
  
  merged <- merged %>% 
    dplyr::mutate(coverage=width/gencodeLengths.tx_len)
  
  write.csv(merged, file = paste0(output,".csv"))
  
}

main()
