# scPreproc
Fast barcode-aware demultiplexing and alignment.

## Usage
For demultiplexing (divmux.codon):
- Have [**codon**](https://github.com/exaloop/codon#install) and [**seq**](https://github.com/exaloop/seq#installation) installed.
-  Run script like this:
`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 divmux.codon R2.gz R1 R3 whitelist stats 1`

e.g.:
`codon run -plugin seq -release -D BC_LEN=16 -D BC_START=29 divmux.codon ScKDMA_S1_R2_001.fastq.gz ScKDMA_S1_R1_001.fastq ScKDMA_S1_R3_001.fastq whitelist_10xATAC_737K-cratac-v1.txt ScKDMA_S1_BC_stats.tsv 1`

Default values for BC_LEN and BC_START need to be passed when invoking codon!

For alignment (divali.sh):
- Have [**bwa**](https://github.com/lh3/bwa) or faster yet, [**bwa-mem2**](https://github.com/bwa-mem2/bwa-mem2) binaries and [**samtools**](https://github.com/samtools/samtools) installed.
- Run script like this:
`./divali.sh exp_type modality sample_name R1 R3 threads path_bwa path_refDB`

e.g.:
`./divali.sh nanoCT modA ScKDMA_S1 ScKDMA_S1_R1_001_correct.fastq ScKDMA_S1_R3_001_correct.fastq 20 /home/bwa-mem2-2.2.1_x64-linux/bwa-mem2 /home/ref-hg19-bwamem2/hg19.fa`

