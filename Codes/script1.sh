#!/bin/bash

# Download the SRA file
echo "Downloading SRA file..."
prefetch SRR25192668.sra
if [ $? -ne 0 ]; then
    echo "Error downloading SRA file"
    exit 1
fi
echo "SRA file downloaded successfully"

# Convert SRA to FASTQ files
echo "Converting SRA to FASTQ..."
fastq-dump --split-files SRR25192668/SRR25192668.sra
if [ $? -ne 0 ]; then
    echo "Error converting SRA to FASTQ"
    exit 1
fi
echo "FASTQ files generated successfully"

# Run FastQC on the FASTQ files
echo "Running FastQC..."
fastqc *.fastq
if [ $? -ne 0 ]; then
    echo "Error running FastQC"
    exit 1
fi
echo "FastQC analysis completed"

