#!/bin/bash

# Purpose:
#   To align R1 & R3 Fastqs as built and corrected by divmux.codon

# Return:
#   Aligned BAM

# Usage:
# ./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R3_correct> <threads> <path_bwa> <path_bwarefDB>
# e.g.:
# ./02_divalign.sh nanoCNT modA ScKDMA_S1 ScKDMA_S1_R1_001_correct.fastq ScKDMA_S1_R3_001_correct.fastq 20 /home/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/refBWAmem2/hg19.fa

# Install bwa-mem2 and index the ref
#   curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
#   ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index hg19.fa
# Install samtools (if not already installed)
set -e

# Variables
exp_type="${1}";
modality="${2}";
sample_name="${3}";
R1="${4}";
R3="${5}";
threads="${6}";
path_bwa="${7}";
path_bwarefDB="${8}";

if [ ${#@} -lt 8 ] ; then
    printf '\nUsage:\n';
    printf '    02_divalign.sh \\\n';
    printf '        exp_type \\\n';
    printf '        modality \\\n';
    printf '        sample_name \\\n';
    printf '        R1 \\\n';
    printf '        R3 \\\n';
    printf '        threads \\\n';
    printf '        path_bwa \\\n';
    printf '        path_bwarefDB \\\n';
    printf 'Parameters:\n';
    printf '  - exp_type: Experience type (ie. nanoCT).\n';
    printf '  - modality: Modality.\n';
    printf '  - sample_name: Sample name.\n';
    printf '  - R1:   Path to corrected FASTQ R1 (uncompressed).\n';
    printf '  - R3:   Path to corrected FASTQ R3 (uncompressed).\n';
    printf '  - threads: Number of threads to align with.\n';
    printf 'Purpose: Align R1 & R3 as built by Divmux\n';
    printf '         Output is a bam file of name <sample_name.modality.exp_type.bam> \n\n';
    exit 1
fi

RGID="${sample_name}.${modality}.${exp_type}" #sample_name.modality.exp_type
library="${sample_name}.${modality}" #sample_name.modality
platform="ILLUMINA" #technology

# Alignement
${path_bwa} mem ${path_bwarefDB} \
-t ${threads} \
-R "@RG\tID:${RGID}\tSM:${sample_name}\tLB:${library}\tPL:${platform}" \
-C \
${R1} ${R3} | samtools view -bS --threads 4 -o ${RGID}.bam -

# UNUSED CODE
# Get alignement stats
#samtools flagstat -@ ${threads} ${RGID}.bam

# Split by BC (creates too many files)
#samtools split -@ ${threads} --max-split '-1' -u NoBarcode.sam -f "%*_%\!" -d 'CB' -v ${RGID}.bam

#Get list of unique BCs in SAM file
#samtools view temp.sam | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq > uniq_BCs
