
ml seqtk 

date

STAR=/scratch/ucgd/lustre/work/quinlan/software-shared/STAR-2.7.9a/source/STAR
STARREF=/scratch/ucgd/lustre-labs/quinlan/u1212967/ref/cellranger/GRCh38-2024-A_STAR_2.7.9a
WHITELIST=~/3M-february-2018.txt

name=$1 # basename of the library
dir=$2 # fastq directory
alnmode=$3 # paired single bco
tech=$4 # illumina or aviti
extra=$5

echo $name
echo $dir

BCS=${dir}/${name}_R1.fastq.gz
CDNA=${dir}/${name}_R2.fastq.gz

OUTPREFIX=/scratch/ucgd/lustre-labs/quinlan/u1212967/aviti_benchmark/star/${tech}/${extra}/${name}_${alnmode}/${name}.
echo $OUTPREFIX
file $BCS
file $CDNA

# paired data with 5' trimming of R1. Softclip first 47 bp because we have to 
# provide at least 1 alignable bp from both reads for the merging modes to work
if [ ${alnmode} == "paired47" ]; then
      extra_params="--soloBarcodeReadLength 150 \
            --readFilesCommand zcat \
            --peOverlapNbasesMin 5 \
            --alignEndsProtrude 25 ConcordantPair \
            --readFilesIn $CDNA $BCS \
            --soloBarcodeMate 2 \
            --clip5pNbases 0 47"
fi

# barcode read only alignment (R2)
if [ ${alnmode} == "bco48" ]; then
      extra_params="--readFilesIn $BCS \
            --readFilesCommand zcat \
      --soloBarcodeMate 1 \
      --soloBarcodeReadLength 150 \
      --clip5pNbases 48 \
      --soloStrand Reverse"
fi

# basic single-end mode with no trimming
if [ ${alnmode} == "single" ]; then
      extra_params="--readFilesIn $CDNA $BCS \
            --readFilesCommand zcat \
            --soloBarcodeReadLength 0"
fi

# cellranger4 mode trimming
if [ ${alnmode} == "single_null" ]; then # if we didn't trim
      extra_params="--readFilesIn $CDNA $BCS \
            --readFilesCommand zcat \
      --soloBarcodeReadLength 0 \
      --clipAdapterType CellRanger4"
fi

$STAR --soloType CB_UMI_Simple \
      --soloCBstart 1 \
      --soloCBlen 16 \
      --soloUMIlen 12 \
      --soloUMIstart 17 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIfiltering MultiGeneUMI \
      --soloUMIdedup 1MM_Directional \
      --runThreadN 12 \
      --soloCBwhitelist $WHITELIST \
      --genomeDir $STARREF \
      --outFileNamePrefix $OUTPREFIX \
      --limitBAMsortRAM 48000000000 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --outReadsUnmapped Fastx \
      --outSAMunmapped Within \
      --soloFeatures Gene GeneFull \
      --outFilterMatchNmin 31 \
      ${extra_params}
