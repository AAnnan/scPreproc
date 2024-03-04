#!/bin/bash

# Purpose:
#   To align R1 & R2 Fastqs as built and corrected by divmux.codon

# Return:
#   Aligned BAM

# Usage:
# ./02_divalign.sh <exp_type> <modality> <sample_name> <Path_R1_correct> <Path_R2_correct> <threads> <path_bwa> <path_bwarefDB> <PathGATK> <RemDups> <PathSamtools> <PathOutputBam> <PathOutputPicardDupStats>
# e.g.:
# ./02_divalign.sh nanoCNT modA ScKDMA_S1 R1_correct.fq R2_correct.fq 4 /home/ahrmad/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/ahrmad/refBWAmem2/hg19.fa /home/ahrmad/gatk-4.5.0.0/gatk true /home/ahrmad/micromamba/envs/ali/bin/samtools /home/ahrmad/testing/TEST.bam /home/ahrmad/testing/TEST_DupMetrics.txt

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

if [ ${#@} -lt 13 ] ; then
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
-R "@RG\tID:${RGID}\tSM:${sample_name}\tLB:${library}\tPL:${platform}" \
-C \
${R1} ${R2} > ${RGID}_TEMP.sam

# replace RG:Z: with CB:Z:
echo "Adding RG tags and header to the alignment file..."
awk '
BEGIN {
    FS="\t";    # Set field separator as tab
    OFS="\t";   # Set output field separator as tab
}
{
    # Skip lines starting with "@"
    if (substr($0, 1, 1) == "@") {
        print;  # Print the line as it is
        next;  # Skip to the next line
    }

    # Use the match function to extract the string following RG:Z: and CB:Z:
    match($0, /RG:Z:([^[:space:]]+)/, rg);
    match($0, /CB:Z:([^[:space:]]+)/, cb);

    # Replace the string following RG:Z: with the string following CB:Z:
    sub("RG:Z:" rg[1], "RG:Z:" cb[1]);
    # Print the modified line
    print
}' ${RGID}_TEMP.sam > ${RGID}_TEMPRG.sam

# Build new header with unique BCs
awk '
BEGIN {
    FS="\t";    # Set field separator as tab
    OFS="\t";   # Set output field separator as tab

    # Declare an associative array to store unique Barcodes
    # This will act as a set to keep track of seen Barcodes
    delete barcodes;
}
{
    # Skip lines starting with "@" (header lines)
    if (substr($0, 1, 1) == "@") {
        if ($1 == "@RG") {
            match($0, /SM:([^[:space:]]+)/, all_sm);
            match($0, /LB:([^[:space:]]+)/, all_lb);
            match($0, /PL:([^[:space:]]+)/, all_pl);
            next;
        } else if ($1 == "@PG") {
            last_PG_line = $0;
            next;
        } else {
            # Print other header lines as they are
            print;
            next;
        }
    }

    # extract the Barcode from CB:Z:
    match($0, /CB:Z:([^[:space:]]+)/, cb);
    if (cb[1]) {
        # Store the extracted Barcode in the associative array
        barcodes[cb[1]] = 1;
    }
}
END {
    # Output the new @RG header lines with unique Barcodes
    for (barcode in barcodes) {
        # Construct the @RG line with the unique Barcode
        rg_line = "@RG\tID:" barcode "\tSM:" all_sm[1] "\tLB:" all_lb[1] "\tPL:" all_pl[1];
        # Print the new @RG line
        print rg_line;
    }
    # Print the last @PG line
    print last_PG_line;
}' ${RGID}_TEMP.sam > ${RGID}_TEMPHEADER.sam

# Append the new header to the RG-replaced SAM file and convert to BAM
echo "Sorting alignment..."
time { cat ${RGID}_TEMPHEADER.sam; grep -v '^@' ${RGID}_TEMPRG.sam; } | ${PathSamtools} sort --threads ${threads} -n -o ${RGID}_TEMP.bam -

# Mark duplicates
echo "Marking duplicates..."
${PathGATK} MarkDuplicatesSpark -I ${RGID}_TEMP.bam -O ${PathOutputBam} -M ${PathOutputPicardDupStats} --remove-all-duplicates ${RemDups} --spark-master local[${threads}] 

#Remove temp files
echo "Removing temp files..."
rm ${RGID}_TEMP.sam ${RGID}_TEMPRG.sam ${RGID}_TEMPHEADER.sam ${RGID}_TEMP.bam

# UNUSED CODE
# Get alignement stats
#samtools flagstat -@ ${threads} ${PathOutputBam}