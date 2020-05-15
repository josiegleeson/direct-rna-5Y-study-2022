#!/usr/bin/env bash
# direct RNA/cDNA pipeline with sequins

# Wrap script in PBS/slurm with:
#for f in diff1 diff2tr1 diff2tr2 diff3 undiff1 undiff2
#do
#echo "Procseeing $f"
#./dRNA.sh sample_genome.fasta sample_transcriptome.fasta sample_annotation.gtf sequin_genome.fasta sequin_transcriptonme.fasta sequin_annotation.gtf $f/pass.fastq $f"_output"
#done

# Requires: python3, minimap2, samtools, bedtools, featureCounts, salmon
HGENOME=$1
HTRANSCRIPTOME=$2
HANNOTATION=$3

SGENOME=$4
STRANSCRIPTOME=$5
SANNOTATION=$6

FASTQ=$7

OUTPREF=$8

function runSample() {
mkdir "$OUTPREF"
mkdir "$OUTPREF/alignments"

echo "Performing genome alignment on $FASTQ"
# Genome
# Mapping with minimap2
minimap2 -ax splice -uf -k14 $HGENOME $FASTQ > $OUTPREF/alignments/genomic-aln.sam
# Convert sam to bam
samtools view -bS $OUTPREF/alignments/genomic-aln.sam > $OUTPREF/alignments/genomic-aln.bam
# Sort bam
samtools sort -o $OUTPREF/alignments/sorted-genomic-aln.bam $OUTPREF/alignments/genomic-aln.bam
# Index the sorted bam
samtools index $OUTPREF/alignments/sorted-genomic-aln.bam
# Create bam with primary alignment only
samtools view -b -h -F 2308 $OUTPREF/alignments/sorted-genomic-aln.bam > $OUTPREF/alignments/primary-genomic-aln.bam
# Convert to bed12 (this is the input for FLAIR)
bedtools bamtobed -bed12 -i $OUTPREF/alignments/primary-genomic-aln.bam > $OUTPREF/alignments/primary-genomic-aln.bed12

echo "Performing transcriptome alignment on $FASTQ" 
# Transcriptome
# Mapping with minimap2
minimap2 -ax map-ont $HTRANSCRIPTOME $FASTQ > $OUTPREF/alignments/transcriptomic-aln.sam
# Convert sam to bam
samtools view -bS $OUTPREF/alignments/transcriptomic-aln.sam > $OUTPREF/alignments/transcriptomic-aln.bam
# Sort bam
samtools sort -o $OUTPREF/alignments/sorted-transcriptomic-aln.bam $OUTPREF/alignments/transcriptomic-aln.bam
# Index the sorted bam
samtools index $OUTPREF/alignments/sorted-transcriptomic-aln.bam
# Create bam with primary alignment only
samtools view -b -h -F 2308 $OUTPREF/alignments/sorted-transcriptomic-aln.bam > $OUTPREF/alignments/primary-transcriptomic-aln.bam
# Conver to bed12
bedtools bamtobed -bed12 -i $OUTPREF/alignments/primary-transcriptomic-aln.bam > $OUTPREF/alignments/primary-transcriptomic-aln.bed12

mkdir "$OUTPREF/quantifications"

echo "Performing gene quantification"
# Quantifty genes with featureCounts
featureCounts -L -a $HANNOTATION -o $OUTPREF/quantifications/primary-genomic-quant $OUTPREF/alignments/genomic-aln.bam --primary
    # All alignments for differential expression analysis
featureCounts -L -a $HANNOTATION -o $OUTPREF/quantifications/genomic-quant $OUTPREF/alignments/genomic-aln.bam

echo "Performing transcript quantification"
# Quantifying transcripts with salmon
    # Primary for transcript stats such as % full-length etc
salmon quant -t $HTRANSCRIPTOME -l A -q --gencode -a $OUTPREF/alignments/primary-transcriptomic-aln.bam $OUTPREF/alignments/transcriptomic-aln.bam -o $OUTPREF/quantifications/primary-transcriptomic-quant
    # All alignments for differential expression analysis
salmon quant -t $HTRANSCRIPTOME -l A -q --gencode -a $OUTPREF/alignments/transcriptomic-aln.bam -o $OUTPREF/quantifications/transcriptomic-quant

echo "Finished for sample"

}

function runSequins() {
mkdir "$OUTPREF-sequins"
mkdir "$OUTPREF-sequins/alignments"

echo "Performing genome alignment on $FASTQ"
# Genome
# Mapping with minimap2
minimap2 -ax splice $SGENOME $FASTQ > $OUTPREF-sequins/alignments/genomic-aln.sam
# Convert sam to bam
samtools view -bS $OUTPREF-sequins/alignments/genomic-aln.sam > $OUTPREF-sequins/alignments/genomic-aln.bam
# Sort bam
samtools sort -o $OUTPREF-sequins/alignments/sorted-genomic-aln.bam $OUTPREF-sequins/alignments/genomic-aln.bam
# Index the sorted bam
samtools index $OUTPREF-sequins/alignments/sorted-genomic-aln.bam
# Create bam with primary alignment only
samtools view -b -h -F 2308 $OUTPREF-sequins/alignments/sorted-genomic-aln.bam > $OUTPREF-sequins/alignments/primary-genomic-aln.bam
# Convert to bed12 (this is the input for FLAIR)
bedtools bamtobed -bed12 -i $OUTPREF-sequins/alignments/primary-genomic-aln.bam > $OUTPREF-sequins/alignments/primary-genomic-aln.bed12

echo "Performing transcriptome alignment on $FASTQ" 
# Transcriptome
# Mapping with minimap2
minimap2 -ax map-ont $STRANSCRIPTOME $FASTQ > $OUTPREF-sequins/alignments/transcriptomic-aln.sam
# Convert sam to bam
samtools view -bS $OUTPREF-sequins/alignments/transcriptomic-aln.sam > $OUTPREF-sequins/alignments/transcriptomic-aln.bam
# Sort bam
samtools sort -o $OUTPREF-sequins/alignments/sorted-transcriptomic-aln.bam $OUTPREF-sequins/alignments/transcriptomic-aln.bam
# Index the sorted bam
samtools index $OUTPREF-sequins/alignments/sorted-transcriptomic-aln.bam
# Create bam with primary alignment only
samtools view -b -h -F 2308 $OUTPREF-sequins/alignments/sorted-transcriptomic-aln.bam > $OUTPREF-sequins/alignments/primary-transcriptomic-aln.bam
# Conver to bed12
bedtools bamtobed -bed12 -i $OUTPREF-sequins/alignments/primary-transcriptomic-aln.bam > $OUTPREF-sequins/alignments/primary-transcriptomic-aln.bed12

mkdir "$OUTPREF-sequins/quantifications"

echo "Performing gene quantification"
# All alignments for differential expression analysis
featureCounts -L -a $SANNOTATION -o $OUTPREF-sequins/quantifications/genomic-quant $OUTPREF-sequins/alignments/sorted-genomic-aln.bam

# Quantifty genes with featureCounts
featureCounts -L -a $SANNOTATION -o $OUTPREF-sequins/quantifications/primary-genomic-quant $OUTPREF-sequins/alignments/sorted-genomic-aln.bam --primary

echo "Performing transcript quantification"
# Quantifying transcripts with salmon
    # Primary for transcript stats such as % full-length etc
salmon quant -t $STRANSCRIPTOME -l A -q -a $OUTPREF-sequins/alignments/primary-transcriptomic-aln.bam -o $OUTPREF-sequins/quantifications/primary-transcriptomic-quant
    # All alignments for differential expression analysis
salmon quant -t $STRANSCRIPTOME -l A -q -a $OUTPREF-sequins/alignments/transcriptomic-aln.bam -o $OUTPREF-sequins/quantifications/transcriptomic-quant

echo "Finished for sequins"

}

# check
if [ $# == 8 ]
then
    runSample
    runSequins
elif [ $1 == "-h" ]
then
    echo "Usage: ./dRNA.sh sample_genome.fasta sample_transcriptome.fasta sample_annotation.gtf sequin_genome.fasta sequin_transcriptonme.fasta sequin_annotation.gtf all.fastq sample_name"
    echo "python3, minimap2, samtools, bedtools, featureCounts and salmon must be in your PATH"
elif [ $# != 8 ]
then
    echo "Usage: ./dRNA.sh sample_genome.fasta sample_transcriptome.fasta sample_annotation.gtf sequin_genome.fasta sequin_transcriptonme.fasta sequin_annotation.gtf all.fastq sample_name"
    echo "python3, minimap2, samtools, bedtools, featureCounts and salmon must be in your PATH"
else
    runSample
    runSequins
fi
