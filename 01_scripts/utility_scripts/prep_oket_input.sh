#!/bin/bash
# Prepare three input files containing chum sequences for amplicon development
# Will produce "02_raw_data/all_input_marker_and_seq.csv"

# Variables
RAW_DATA="02_input_data"

# Inputs
SEEB_INPUT="Seeb"
UW_INPUT="UW_GTseq"
WDFW_INPUT="WDFW_GTseq"

# Date and extension
DATE_EXT="2018-10-05.csv"

OUTPUT_EXT="marker_and_seq.csv"

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

# WDFW input
# Rename only
sed 's/\-/\_/g' $RAW_DATA/$WDFW_INPUT"_"$DATE_EXT > $RAW_DATA/$WDFW_INPUT"_"$OUTPUT_EXT

## Combine all data into one file
cat $RAW_DATA/$SEEB_INPUT"_"$OUTPUT_EXT \
    $RAW_DATA/$UW_INPUT"_"$OUTPUT_EXT \
    $RAW_DATA/$WDFW_INPUT"_"$OUTPUT_EXT \
    > $RAW_DATA/"all_input_"$OUTPUT_EXT
