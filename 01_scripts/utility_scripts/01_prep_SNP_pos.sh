#!/bin/bash
# Take the accn,seq file w/ [A/T] alleles indicated and find the position of the SNP
# Then remove the square brackets and prepare for BLAST
# Output is "02_input_data/all_input_marker_and_seq_w_pos_for_BLAST.fa"

# Variables
RAW_DATA="02_input_data"
INPUT="all_input_marker_and_seq.csv"
SNP_POS="first_snp_position.csv"
OUTPUT="all_input_marker_and_seq_w_pos.csv"

rm $RAW_DATA/$SNP_POS 

## Make file containing position of first SNP 
awk -F"," '{ print $2 }' $RAW_DATA/$INPUT | 
    while read l
    do
        echo $l | cut -d"[" -f1 | wc -c
    done > $RAW_DATA/$SNP_POS 

## Combine position info into the file 
paste -d, $RAW_DATA/$INPUT $RAW_DATA/$SNP_POS > $RAW_DATA/$OUTPUT

## Put into fasta
awk -F"," '{ print ">" $1 "-" $3 "\n" $2 }' $RAW_DATA/$OUTPUT > $RAW_DATA/${OUTPUT%.csv}".fa"
 
## Remove characters ([]/) to prevent issues with BLAST
sed 's/\[//g; s/\///g; s/]//g' $RAW_DATA/${OUTPUT%.csv}".fa" > $RAW_DATA/${OUTPUT%.csv}"_for_BLAST.fa"

## Reporting at end
echo "In total, there are the following number of records prepared for BLAST"
grep -cE '^>' $RAW_DATA/${OUTPUT%.csv}"_for_BLAST.fa"
