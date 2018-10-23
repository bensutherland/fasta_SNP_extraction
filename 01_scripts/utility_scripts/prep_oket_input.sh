#!/bin/bash
# Prepare input files containing chum sequences for amplicon development
# Will produce "02_raw_data/all_input_marker_and_seq.csv"

# Global variables
RAW_DATA="02_input_data"

# Inputs
SEEB_INPUT="Seeb"
UW_INPUT="UW_GTseq"
WDFW_INPUT="WDFW_GTseq"
UW_SECOND_SET="UW_all" # this is to capture the 'not included' snps
MGL_INPUT="MGL_chum_amplicons_w_alleles"

# Date and extension
DATE_EXT="2018-10-05.csv" # essentially the input extension
DATE_EXT2="2018-10-23.csv" # to allow for a second date (MGL)
OUTPUT_EXT="marker_and_seq.csv" # the output extension

## Put all inputs into the same format
# SEEB input
# Remove empty lines
grep -vE '^,,$' $RAW_DATA/$SEEB_INPUT"_"$DATE_EXT | 

    # Keep only the required columns
    awk -F"," '{ print $1 "," $3 }' - |

    # Remove anything following the first space in the accession name, save out
    sed 's/\s.*,/,/' |

    # Insert the forward slash between the alleles
    sed 's/\[[A-Z]/&\//g' > $RAW_DATA/$SEEB_INPUT"_"$OUTPUT_EXT

# UW input
# Remove header
grep -vE '^Panel,' $RAW_DATA/$UW_INPUT"_"$DATE_EXT | 

    # Keep only the required columns
    awk -F"," '{ print $2","$3}' - > $RAW_DATA/$UW_INPUT"_"$OUTPUT_EXT

# UW input second set
# Remove header
grep -vE '^Panel,' $RAW_DATA/$UW_SECOND_SET"_"$DATE_EXT | 

    # Only keep lines containing the 'not included' markers
    grep 'not included' - | 

    # Keep only the required columns
    awk -F"," '{ print $2","$3 }' - > $RAW_DATA/$UW_SECOND_SET"_"$OUTPUT_EXT


# WDFW input
# Replace the hyphen with an underscore to not disrupt future steps
sed 's/\-/\_/g' $RAW_DATA/$WDFW_INPUT"_"$DATE_EXT > $RAW_DATA/$WDFW_INPUT"_"$OUTPUT_EXT

# MGL input
# Remove header
grep -vE '^mname' $RAW_DATA/$MGL_INPUT"_"$DATE_EXT2 | 
    
    # Remove fasta symbol '>'
    sed 's/>//g' - > $RAW_DATA/$MGL_INPUT"_"$OUTPUT_EXT

## Combine all data into one file
cat $RAW_DATA/$SEEB_INPUT"_"$OUTPUT_EXT \
    $RAW_DATA/$UW_INPUT"_"$OUTPUT_EXT \
    $RAW_DATA/$WDFW_INPUT"_"$OUTPUT_EXT \
    $RAW_DATA/$UW_SECOND_SET"_"$OUTPUT_EXT \
    $RAW_DATA/$MGL_INPUT"_"$OUTPUT_EXT \
    > $RAW_DATA/"all_input_"$OUTPUT_EXT

