# fasta_SNP_extraction
Extract a section of a reference genome flanking an input locus with SNPs.

This repo is set up to identify windows within a reference genome using a marker file. This was first done for Eulachon as part of the Molecular Genetics Laboratory (MGL) at the Fisheries and Oceans Canada Pacific Biological Station (PBS).      

Note: Use at your own risk. This is provided as is, without any guarantee of usefulness.
All scripts are to be run from the existing folder.   


### 0.1 Input information for Eulachon
Original files, put in `02_input_data`:    
`tpac_adaptive_loci_2018-03-22.csv`
`tpac_neutral_loci_2018-03-22.csv`

These are in the following format:    
```
Loci_name,Loci_Sequence,SNP_location
Tpa_0002,GTTCGTTCAGTCTAAAGAGAATAATAGGACCGGGAGGTTTGTACTTATGTAATAGACACAC[A/T]CACATACATATGCGCA,62
```

Genome file:    
`tpac_assembly/tpac_assembly_v1/tpac-contigs_min200.fa`


#### 0.1a Merge and prepare input loci files
Join input files, remove header, remove everything after the occasional colon (multi-SNP location):     
`cat 02_input_data/tpac_adaptive_loci_2018-03-22.csv 02_input_data/tpac_neutral_loci_2018-03-22.csv | grep -vE 'Loci.name' - | awk -F":" '{ print $1 }' - > 02_input_data/input_loci.csv`

Re-calculate the position of the first SNP based on square bracket position:         
`awk -F, '{ print $2 }' 02_input_data/input_loci.csv | while read l ; do echo $l | cut -d[ -f1 | wc -c ; done > 02_input_data/first_snp_position.csv`

Add re-calculated position to the original file   
`paste -d, 02_input_data/input_loci.csv 02_input_data/first_snp_position.csv > 02_input_data/input_loci_w_pos.csv`   

Select columns of interest and put into FASTA:      
`awk -F"," '{ print ">"$1"-"$4"\n"$2 }' 02_input_data/input_loci_w_pos.csv > 02_input_data/input_loci_w_pos.fa`   

Remove square brackets and slash to prevent issues with BLAST:    
`sed 's/\[//g' 02_input_data/input_loci_w_pos.fa | sed 's/\///g' - | sed 's/]//g' - > 02_input_data/input_loci_w_pos_no_badchar.fa`

The output looks like this, where the accession name is '>markername-firstSNPposition':    
```
>Tpa_0012-44
TCAGTGAGGCCCACTTCCTGGGTTTCCATCGACACACACACACCACGCTAGTGGATGGAGGGAAGGACGATTCAGGGA
```

#### 0.1b Prepare genome for alignment    
Fix the names of the fasta accessions, removing spaces, commas, '+' or dashes:         
`sed 's/\ /\_/g' 02_input_data/genome/tpac-contigs_min200.fa | cut -d, -f1 | sed 's/\+//g' - | sed 's/\-//g' - > 02_input_data/genome/tpac-contigs_min200_renamed.fa`

Set up genome file as a blast database
`makeblastdb -in ./02_input_data/genome/tpac-contigs_min200_renamed.fa -parse_seqids -dbtype nucl`

Get lengths of each contig from the genome (to be used by R script below)     
`GENOME="02_input_data/genome/tpac-contigs_min200_renamed.fa" ; cat $GENOME | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n+2 > 02_input_data/genome/ref_genome_seq_lengths.txt`

This will be used by an Rscript below to ensure windows don't go beyond the contig size. 



### 1. BLAST search loci against genome 
Use BLAST to align the loci to the genome:    
`blastn -db ./02_input_data/genome/tpac-contigs_min200_renamed.fa -query ./02_input_data/input_loci_w_pos_no_badchar.fa -out 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6.txt -outfmt 6 -evalue 1e-20` 


Retain only those loci that have fewer than three significant hits:
First find the loci that map more than three times:      
`awk '{ print $1 }' 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6.txt | sort -n | uniq -c | sort -nk1 | awk '$1 > 2 { print $2 }' - > bad_loci.txt`

Then extract these from the BLAST output:   
`grep -vf bad_loci.txt 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6.txt > 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6_remove_multimap.txt`

Sort the BLAST output:
`sort -k1,1 -k12,12gr -k11,11g -k3,3gr 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6_remove_multimap.txt > 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6_remove_multimap_sorted.txt`

Keep only a single record per accession:    
`DATA="03_blast_out/input_loci_v_tpac_no_badchar_outfmt6_remove_multimap_sorted.txt" ; for i in $(cut -f1 $DATA | sort -u); do grep -w -m 1 "$i" $DATA; done > 03_blast_out/input_loci_v_tpac_no_badchar_outfmt6_single_hit.txt`



### 2. Identify the window range on the ref genome to capture the locus
Use Rscript `01_scripts/identify_req_region_w_BLAST.R`

This will output:     
`04_extraction/ranges_for_amplicon_extraction.bed`    
`04_extraction/ranges_for_amplicon_extraction.csv`

The bedfile is for extraction, the .csv gives additional useful information

### 3. Collect the amplicons
Collect the amplicon ranges using the bed file and the assembly using bedtools:   
`GENOME="02_input_data/genome/tpac-contigs_min200_renamed.fa" ; bedtools getfasta -fi $GENOME -bed 04_extraction/ranges_for_amplicon_extraction.bed -fo 04_extraction/tpac_amplicon_approx_by_BLAST.fa`

Rename the amplicons:     
`awk 'BEGIN{RS=">"}{print $1"\t"$2;}' tpac_amplicon_approx_by_BLAST.fa | tail -n+2 > 04_extraction/tpac_amplicon_approx_by_BLAST.txt`
Note: the tail -n+2 removes an empty first line
(from https://www.biostars.org/p/235052/ )


########## HERE TODAY #############






but like this:
`awk 'BEGIN{RS=">"}{print $1"\t"$2;}' tpac_amplicon_approx_by_BLAST.fa | tail -n+2 > tpac_amplicon_approx_by_BLAST.txt`


GO to Rscript and input the `tpac_amplicon_approx_by_BLAST.txt`

This exports `all_fields.txt`
which is essentially everything needed with all tab-delim











#### OLD CODE ####


Collect necessary info into csv   
`awk '{ print $1","$2","$3","$4 }' input_loci_w_pos_q30.sam > aligned_req_info.csv`

Re-collect the SNP position from the first field
`awk -F"-" '{ print $1","$2 }' aligned_req_info.csv > aligned_req_info_w_snp.csv`



#### Collect necessary info
Output a forward alignment file
`BAM="tpac_neutral_loci_w_new_snp_pos.fa.aligned.bam" ; samtools view -F 16 -q 30  $BAM > $BAM".forward.sam"`

Collect the marker name, the ref name, and the position of the start alignment into a csv
`awk '{ print $1","$3","$4 }' tpac_neutral_loci_w_new_snp_pos.fa.aligned.bam.forward.sam > tpac_neutral_req_info_forward.csv`

Gain the SNP position
`awk -F"_" '{ print $1"_"$2","$3 }' *info_forward.csv > info_forward_w_snp.csv`


### Obtain the segments needed from the fasta
Note: make sure .bed is tab-delimited
`bedtools getfasta -fi test2.fa -bed test2.bed -fo out.fa`


