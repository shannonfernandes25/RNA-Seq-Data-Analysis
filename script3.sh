#!/bin/bash
'''
# HISAT2 Alignment
echo "Starting HISAT2 alignment..."
hisat2 -x ../grch38/genome -U trimmed.fastq -S aligned.sam --rna-strandness R -p 8
if [ $? -ne 0 ]; then
    echo "Error during HISAT2 alignment"
    exit 1
fi
echo "HISAT2 alignment complete"
'''
# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view -@ 8 -bS aligned.sam > aligned.bam
if [ $? -ne 0 ]; then
    echo "Error converting SAM to BAM"
    exit 1
fi
echo "SAM to BAM conversion complete"

# Sort BAM file
echo "Sorting BAM file..."
samtools sort -@ 8 -o sorted.bam aligned.bam
if [ $? -ne 0 ]; then
    echo "Error sorting BAM file"
    exit 1
fi
echo "BAM sorting complete"

# Index BAM file
echo "Indexing BAM file..."
samtools index sorted.bam
if [ $? -ne 0 ]; then
    echo "Error indexing BAM file"
    exit 1
fi
echo "BAM indexing complete"

# Assess Alignment Quality
echo "Assessing alignment quality with SAMtools flagstat..."
samtools flagstat sorted.bam
if [ $? -ne 0 ]; then
    echo "Error assessing alignment quality"
    exit 1
fi
echo "Alignment quality assessment complete"

# FeatureCounts for Gene Quantification
echo "Running featureCounts..."
featureCounts -a ../Homo_sapiens.GRCh38.110.gtf/Homo_sapiens.GRCh38.110.gtf -o gene_counts.txt -T 8 sorted.bam
if [ $? -ne 0 ]; then
    echo "Error running featureCounts"
    exit 1
fi
echo "FeatureCounts analysis complete"

