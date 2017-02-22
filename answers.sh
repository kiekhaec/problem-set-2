
# /usr/bin/env bash

datasets="$HOME/data-sets"

#1Use BEDtools intersect to identify the size of the largest overlap
#between CTCF and H3K4me3 locations. 

#BED file containing ChIP-seq peaks for H3K4me3 in Hela cells: 
#bed/encode.h3k4me3.hela.chr22.bed.gz 

#Bed File containing peak calls for ENCODE transcription factor ChIP-seq
#experiements: bed/encode.tfbs.chr22.bed.gz

CTCF="$datasets/bed/encode.tfbs.chr22.bed.gz"
H3K4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(bedtools intersect -a $CTCF -b $H3K4me3 \
| awk 'BEGIN {OFS ="\t"} ($4 == "CTCF") {print $0, $3 - $2}'\
| sort -k5nr \
| head -n1 \
| cut -f5)

echo "answer-1: $answer_1"

#2.Use BEDtools to calculate the GC content of nucleotides 19,000,000
#to 19,000,500 on chr22 of hg19 genome build.
#Report the GC content as a fraction e.g. 0.50

#One option is to use bedtools getfasta-- this extracts sequences from a FASTA file for each of the
#intervals defined in a BED/GFF/VCF file. In this problem a different solution, bedtools nuc, will be employed to calculate GC percent.

 FASTA="$datasets/fasta/hg19.chr22.fa"
 answer_2=$(echo -e "chr22\t19000000\t19000500" > interval.bed \
 | bedtools nuc -fi $FASTA -bed interval.bed \
 | cut -f5 \
 | tail -n1)

 echo "answer-2: $answer_2"

#3.Use BEDtools to identify the length of the CTCF ChIP-seq peak (e.g.,
#interval) that has the largest mean signal in ctcf.hela.chr22.bg.

#Note:bedtools map will be used to compute the mean score of the 
#BEDGRAPH record that overlaps genes. The fourth column in BEDGRAPH is the score.

#File containing peak calls for ENCODE transcription factor ChIP-seq
#experiements: bed/encode.tfbs.chr22.bed.gz.

#Use the peaks already defined in bed/encode.tfbs.chr22.bed.gz
#to get the CTCF ChIP-Seq peaks.

CTCFPeaks="$datasets/bed/encode.tfbs.chr22.bed.gz"
Signal="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(bedtools map -c 4 -o mean -a $CTCFPeaks -b $Signal \
| sort -k5nr \
| awk 'BEGIN {OFS="\t"} ($4 == "CTCF") {print $3 - $2}' \
| head -n1)

echo "answer-3: $answer_3"

#4. Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
#of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz.
#Report the gene name (e.g., 'ABC123')

#A file containing transcription start sites (TSS) for chr22:
#bed/tss.hg19.chr22.bed.gz

#File with Bedgraph with CTCF ChIP-seq data in bedGraph format:
#bedgraph/ctcf.hela.chr22.bg.gz

CTCF="$datasets/bedtools/ctcf.hela.chr22.bg.gz"
TSS="$datasets/bed/tss.hg19.chr22.bed.gz"
hg19="$datasets/genome/hg19.genome"

answer_4=$(bedtools slop -i $TSS -g $hg19 -l 1000 -r 0 -s \
| sort \
| bedtools map -c 4 -o median -a $TSS -b $CTCF \
| sort -k7nr \
| head -n1 \
| cut -f4)

echo "answer-4: $answer_4"

#5. Use BEDtools to identify the longest interval on chr22 that is not
#covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.

#BED file with all genes in hg19: bed/genes.hg19.bed.gz.

#"genome file" with chromosome size info: genome/hg19.genome

genes="$datasets/bed/genes.hg19.bed.gz"

Genome="$datasets/genome/hg19.genome"

answer_5=$(bedtools sort -i $genes \
| bedtools complement -i - -g $Genome \
| awk 'BEGIN {OFS="\t"} ($1=="chr22") {print $1, $2, $3, $3-$2}' \
| sort -k4nr \
| head -n1 \
| awk '{print $1":"$2"-"$3}')

echo "answer-5: $answer_5"

#6.Extra credit: Use one or more bedtools that we have not covered in
#class. In Problem #1 we were asked to Use BEDtools intersect to identify
#the size of the largest overlap between CTCF and H3K4me3 locations.
#Another tool that can be used to address this question is bedtools
#subtract. 

CTCF="$datasets/bed/encode.tfbs.chr22.bed.gz"
H3K4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_6=$(bedtools subtract -a $CTCF -b $H3K4me3 \
| awk 'BEGIN {OFS ="\t"} ($4 == "CTCF") {print $0, $3 - $2}'\
| sort -k5nr \
| head -n1 \
| cut -f5)

echo "answer-6: $answer_6"
















