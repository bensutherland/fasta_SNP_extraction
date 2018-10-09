#!/bin/bash
# Use BLAST to align the marker file against the reference genome
# Note: genome must be prepared using makeblastdb

# Variables
GENOME="02_input_data/genome/oket_k101-contigs_500bp_renamed.fa"
MARKERS="02_input_data/all_input_marker_and_seq_w_pos_for_BLAST.fa"
BLAST_OUTPUT="03_blast_out/prepared_tags_v_oket-contigs_outfmt6.txt"
BAD_LOCI="03_blast_out/bad_loci.txt"
BLAST_OUTPUT_TWO_OR_LESS="03_blast_out/prepared_tags_v_oket-contigs_outfmt6_rem_multimap.txt"
BLAST_SORTED=${BLAST_OUTPUT_TWO_OR_LESS%.txt}"_sort.txt"
OUTPUT="03_blast_out/prepared_tags_v_oket-contigs_outfmt6_rem_multimap_sort_single_hit.txt"

# BLAST search against genome
blastn -db $GENOME -query $MARKERS -out $BLAST_OUTPUT -outfmt 6 -evalue 1e-20

# Identify bad loci
awk '{ print $1 }' $BLAST_OUTPUT | 
    # sort BLAST output to identify when more than two hits
    sort -n | uniq -c | sort -nk1 | 
    awk '$1 > 2 { print $2 }' - \
    > $BAD_LOCI 

# Extract bad loci from BLAST output
grep -vf $BAD_LOCI $BLAST_OUTPUT > $BLAST_OUTPUT_TWO_OR_LESS 

# Sort remaining BLAST output
sort -k1,1 -k12,12gr -k11,11g -k3,3gr $BLAST_OUTPUT_TWO_OR_LESS \
    > $BLAST_SORTED 

# Keep single record per accession
DATA=$BLAST_SORTED ; for i in $(cut -f1 $DATA | sort -u); do grep -w -m 1 "$i" $DATA; done > $OUTPUT 
