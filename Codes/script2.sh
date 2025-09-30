#!/bin/bash

# Trimming reads with Trimmomatic 

echo "Starting trimming with Trimmomatic..."

java -jar trimmomatic-0.39.jar SE -threads 6 *.fastq trimmed.fastq LEADING:3 TRAILING:3 -phred33

if [ $? -ne 0 ]; then
    echo "Error during trimming"
    exit 1
fi
echo "Trimming completed"

# Generate FastQC reports
echo "Running FastQC on trimmed reads..."
fastqc trimmed.fastq

if [ $? -ne 0 ]; then
    echo "Error running FastQC"
    exit 1
fi
echo "FastQC analysis completed"

# Generate MultiQC report
echo "Generating MultiQC report..."
multiqc .

if [ $? -ne 0 ]; then
    echo "Error generating MultiQC report"
    exit 1
fi
echo "MultiQC report generated successfully"


