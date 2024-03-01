#!/bin/bash

# Purpose:
#   To create the peak/cell matrix from the aligned BAM file

# Return:
#   peak/cell matrix TSV matrix

# Usage:

set -e

### To create the peak/cell matrix
# Convert Peak BED file to SAF format
awk 'OFS="\t" {print $1"."$2"."$3"."$4, $1, $2+1, $3, "."}' ${inputBEDFile} > ${outputSAFFile}
sed -i '1s/^/GeneID\tChr\tStart\tEnd\tStrand\n/' ${outputSAFFile}

# if you're not confident in your sed version (ie using a mac), change the above line to this:
#echo -e "GeneID\tChr\tStart\tEnd\tStrand" | cat - ${inputBEDFile} > temp && mv temp ${outputSAFFile}

#Run featureCounts to count number of reads per peak
# install subread/featurecounts: "conda install bioconda::subread"
featureCounts -T ${threads} -a ${outputSAFFile} -F SAF -p --byReadGroup -O -o ${RGID}.tsv ${RGID}.bam


# UNUSED CODE
# Get alignement stats
#samtools flagstat -@ ${threads} ${RGID}.bam

# Split by BC (WILL CRASH YOUR SYSTEM creates too many files)
#samtools split -@ ${threads} --max-split '-1' -u NoBarcode.sam -f "%*_%\!" -d 'CB' -v ${RGID}.bam

#Get list of unique BCs in SAM file
#samtools view temp.sam | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort | uniq > uniq_BCs

# TESTING CODE
#macs3 callpeak --treatment TESTDONE.bam \
#--format BAMPE \
#--gsize hs \
#--nomodel \
#--nolambda \
#--keep-dup auto \
#--qvalue 0.01 \
#--call-summits \
#--outdir callpeak_TEST \
#--name merged \
#--verbose 2

#awk 'OFS="\t" {print $1"."$2"."$3"."$4, $1, $2+1, $3, "."}' TEST.bed > TEST.saf
#sed -i '1s/^/GeneID\tChr\tStart\tEnd\tStrand\n/' TEST.saf
#featureCounts -T 16 -a TEST.saf -F SAF -p --byReadGroup -O -o TESTDONE.tsv TESTDONE.bam

#PathPicard="/home/ahrmad/picard.jar"
#RemDups=false
#java -jar ${PathPicard} MarkDuplicates I=TESTDONE.bam O=TESTDONE_MarkedDup.bam M=TESTDONE_DupMetrics.txt REMOVE_DUPLICATES=${RemDups}
