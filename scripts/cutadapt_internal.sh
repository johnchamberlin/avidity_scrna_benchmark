#!/bin/bash

# we want to trim paired end reads for adapters and polyA tails
# but we also need to trim longer-than-expected polyT from read1
# this is internal sequence so we need to split the barcodes off first

ml seqtk
ml bbtools

# 1. split barcodes off of read1 (keep 59 bp on barcode segment)
# 2. trim TSO at 5' and polyA at 3' of read2
# 3. trim TSO and Truseq at 3' of read1, and polyT and 5' of read1
# 4. Paste barcodes back on, filter out 0 length reads, done?

for samp in AG1 AG2 AG3 AG4 AG5; do
    # name of parent file without path to it:
    r1=${samp}_R1.fastq.gz
    r2=${samp}_R2.fastq.gz

    # name of temporary files:
    r1_nobc_trim="${r1%.fastq.gz}.nobc_trim.fastq.gz"
    r1_trim="${r1%.fastq.gz}.trim.fastq.gz"
    r2_trim="${r2%.fastq.gz}.trim.fastq.gz"

    r1_cdna=${r1%.fastq.gz}.cdna.fastq.gz

    # trim the first 48 bp (20 bp of the primer), keep rest as 'cDNA'
    echo "extracting cDNA segment from R1"
    seqtk trimfq -b 48 ../${r1} | gzip > ${r1_cdna}

    # "On paired-end reads, --poly-a removes poly-A tails from R1 and poly-T “heads” from R2"
    # we need to reverse the read order (lame)
    echo "trimming"
    # -A: 3' end of second read input (R1 in our case - TSO and TruSeq)
    # -a: 3' end of first inputted read (R2 in our case - polyA tail)
    # -g: 5' end of first inputted read (R2 in our case - TSO)
    # --poly-a: remove poly-A tails from first read and poly-T heads from second read,
    # but we already removed polyA tails as an adapter so those shouldn't change
    cutadapt -j 24 \
         -A TSO_R1=CCCATGTACTCTGCGTTGATACCACTGCTT \
         -A TruSeq_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -g TSO_R2=AAGCAGTGGTATCAACGCAGAGTACATGGG \
         -a polya_R2=A{30} \
         --poly-a \
         -n 4 \
         -o $r2_trim \
         -p $r1_nobc_trim \
         --json="${samp}.cutadapt.json" \
         ../$r2 ${r1_cdna}
         # add the barcodes back onto r1 trimmed:

    echo "pasting barcodes back onto R1"

    paste <(seqtk trimfq -e 102 ../${r1}) <(zcat ${r1_nobc_trim}) | awk '{if (NR%4==1) {print $1 " " $2} else if (NR%4==3) {print "+"} else {print $1 $2}}' | gzip > $r1_trim

    echo "filtering short r2"
    seqtk seq -L 31 ${r2_trim} | gzip > ${r2_trim%.fastq.gz}.fil.fastq.gz

    echo "fixing read pairs"
    repair.sh in=${r1_trim} in2=${r2_trim%.fastq.gz}.fil.fastq.gz out=${r1_trim%.fastq.gz}.fix.fastq.gz out2=${r2_trim%.fastq.gz}.fix.fastq.gz changequality=f

done

rename DS DStrimint *.gz
rename trim.fix.fastq.gz fastq.gz *.gz
# remove temporary files: