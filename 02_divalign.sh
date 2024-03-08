#!/bin/bash

# Purpose:
#   To align R1 & R2 Fastqs as built and corrected by divmux.codon

# Return:
#   Aligned BAM

# Usage:
# ./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R2_correct> <threads> <path_bwa> <path_bwarefDB> <PathGATK> <RemDups> <PathSamtools> <PathOutputBam> <PathOutputPicardDupStats> <sam_header>
# e.g.:
# ./02_divalign.sh nanoCNT modA ScKDMA_S1 R1_correct.fq R2_correct.fq 4 /home/ahrmad/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/ahrmad/refBWAmem2/hg19.fa /home/ahrmad/gatk-4.5.0.0/gatk true /home/ahrmad/micromamba/envs/ali/bin/samtools /home/ahrmad/testing/small/xs/TEST.bam /home/ahrmad/testing/small/xs/TEST_DupMetrics.txt /home/ahrmad/testing/small/xs/sam_header.txt

# Install bwa-mem2 and index the ref
#   curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
#   ./bwa-mem2-2.2.1_x64-linux/bwa-mem2 index hg19.fa
# Install samtools (if not already installed)
# Install gatk (https://github.com/broadinstitute/gatk/releases/latest)
set -e

# Variables
exp_type="${1}";
modality="${2}";
sample_name="${3}";
R1="${4}";
R2="${5}";
threads="${6}";
path_bwa="${7}";
path_bwarefDB="${8}";
PathGATK="${9}";
RemDups="${10}";
PathSamtools="${11}";
PathOutputBam="${12}";
PathOutputPicardDupStats="${13}";
PathSam_header="${14}";

if [ ${#@} -lt 14 ]; then
    printf '\nUsage:\n';
    printf '    02_divalign.sh \\\n';
    printf '        exp_type \\\n';
    printf '        modality \\\n';
    printf '        sample_name \\\n';
    printf '        R1 \\\n';
    printf '        R2 \\\n';
    printf '        threads \\\n';
    printf '        path_bwa \\\n';
    printf '        path_bwarefDB \\\n';
    printf '        PathGATK \\\n';
    printf '        RemDups \\\n';
    printf '        PathSamtools \\\n';
    printf '        PathOutputBam \\\n';
    printf '        PathOutputPicardDupStats\n';
    printf '        PathSam_header\n';
    printf 'Parameters:\n';
    printf '  - exp_type: Experience type (ie. nanoCT).\n';
    printf '  - modality: Modality.\n';
    printf '  - sample_name: Sample name.\n';
    printf '  - R1:   Path to corrected FASTQ R1 (uncompressed).\n';
    printf '  - R2:   Path to corrected FASTQ R2 (uncompressed).\n';
    printf '  - threads: Number of threads to align with.\n';
    printf '  - path_bwa: Path to BWA binary.\n';
    printf '  - path_bwarefDB: Path to BWA DB w/ extension (ex:hg19.fa).\n';
    printf '  - PathGATK: Path to gatk bin.\n';
    printf '  - RemDups: Remove Duplicates (true) or only mark them (false).\n';
    printf '  - PathSamtools: Path to samtools binary.\n';
    printf '  - PathOutputBam: Path to output BAM file.\n';
    printf '  - PathOutputPicardDupStats: Path to output Picard Duplication Stats file.\n';
    printf '  - PathSam_header: Path to SAM header file as produced by demux.codon.\n';
    printf 'Purpose: Align R1 & R2 as built by divmux codon script\n';
    printf '         Output is a bam file with marked/removed duplicates and associated duplication metrics \n\n';
    exit 1
fi

RGID="${sample_name}-${modality}-${exp_type}" #sample_name.modality.exp_type
library="${sample_name}.${modality}" #sample_name.modality
platform="ILLUMINA" #technology

# Alignement
echo "Aligning ${R1} and ${R2} with BWA-MEM2"

${path_bwa} mem ${path_bwarefDB} \
-t ${threads} \
-C \
${R1} ${R2} > ${RGID}_TEMP.sam

echo "Changing header and sorting alignment..."
# Prepare the header
${PathSamtools} view -H ${RGID}_TEMP.sam > ${RGID}_TEMPHEADER.sam
sed -e "s/.*/&\tSM:$sample_name\tLB:$library\tPL:$platform/" ${PathSam_header} > ${RGID}_RG_HEADER.sam
# Find the line number of the last occurrence of "@SQ"
last_sq_line=$(grep -n "@SQ" "${RGID}_TEMPHEADER.sam" | tail -n1 | cut -d':' -f1)
# Split the header based on the line number of the last "@SQ" occurrence
head -n "$last_sq_line" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER1.sam"
echo "" > "${RGID}_TEMPHEADER2.sam" 
tail -n +"$((last_sq_line + 1))" "${RGID}_TEMPHEADER.sam" >> "${RGID}_TEMPHEADER2.sam"

# Append the new header to the RG-replaced SAM file and convert to BAM
{ cat ${RGID}_TEMPHEADER1.sam ${RGID}_RG_HEADER.sam ${RGID}_TEMPHEADER2.sam; ${PathSamtools} view ${RGID}_TEMP.sam; } | ${PathSamtools} sort --threads ${threads} -n -o ${RGID}_TEMP.bam -

# Mark duplicates
echo "Marking duplicates..."
${PathGATK} MarkDuplicatesSpark -I ${RGID}_TEMP.bam -O ${PathOutputBam} -M ${PathOutputPicardDupStats} --remove-all-duplicates ${RemDups} --spark-master local[${threads}] 

#Remove temp files
echo "Removing temp files..."
rm ${RGID}_TEMP.sam ${RGID}_TEMPHEADER.sam ${RGID}_TEMPHEADER1.sam ${RGID}_TEMPHEADER2.sam ${RGID}_RG_HEADER.sam ${RGID}_TEMP.bam

# UNUSED CODE
# Get alignement stats
#samtools flagstat -@ ${threads} ${PathOutputBam}