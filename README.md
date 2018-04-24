# fasta_SNP_extraction
Extract a section of a reference genome assembly flanking a SNP from an input tag.   

Designed as part of a project with the Molecular Genetics Laboratory (MGL) at the Fisheries and Oceans Canada Pacific Biological Station (PBS).      

Warning: this pipeline is provided as is, without any guarantee of usefulness.    

All scripts are to be run from the main folder.   


### 0.1 Input information (here for Eulachon)    
Original files, put in `02_input_data`:    
`tpac_adaptive_loci_2018-03-22.csv`
`tpac_neutral_loci_2018-03-22.csv`

Original files are in the following format:    
```
Loci_name,Loci_Sequence,SNP_location
Tpa_0002,GTTCGTTCAGTCTAAAGAGAATAATAGGACCGGGAGGTTTGTACTTATGTAATAGACACAC[A/T]CACATACATATGCGCA,62
```

Genome contig assembly, put in `02_input_data/genome`:    
`tpac_assembly/tpac_assembly_v1/tpac-contigs_min200.fa`

#### 0.1a Combine and prepare input loci files
Concatenate input files, remove header, remove rest of line after first colon (occurs when multiple SNP locations):     
`cat 02_input_data/tpac_adaptive_loci_2018-03-22.csv 02_input_data/tpac_neutral_loci_2018-03-22.csv | grep -vE 'Loci.name' - | awk -F":" '{ print $1 }' - > 02_input_data/input_loci.csv`

Calculate the *position of the first SNP in the locus* (based on square bracket position):         
`awk -F, '{ print $2 }' 02_input_data/input_loci.csv | while read l ; do echo $l | cut -d[ -f1 | wc -c ; done > 02_input_data/first_snp_position.csv`

Add re-calculated position to the original file   
`paste -d, 02_input_data/input_loci.csv 02_input_data/first_snp_position.csv > 02_input_data/input_loci_w_pos.csv`   

Select columns of interest and put into FASTA:      
`awk -F"," '{ print ">"$1"-"$4"\n"$2 }' 02_input_data/input_loci_w_pos.csv > 02_input_data/input_loci_w_pos.fa`   

Remove characters ([]/) to prevent issues with BLAST:    
`sed 's/\[//g' 02_input_data/input_loci_w_pos.fa | sed 's/\///g' - | sed 's/]//g' - > 02_input_data/input_loci_w_pos_no_badchar.fa`

Rename and delete intermediate files
`mv 02_input_data/input_loci_w_pos_no_badchar.fa 02_input_data/prepared_tags.fa`
`mv 02_input_data/*.csv 02_input_data/input_loci_w_pos.fa 02_input_data/z-draft_input/`

The output looks like this, where the accession name is `'>markername-firstSNPposition'`:    
```
>Tpa_0012-44
TCAGTGAGGCCCACTTCCTGGGTTTCCATCGACACACACACACCACGCTAGTGGATGGAGGGAAGGACGATTCAGGGA
```

#### 0.1b Prepare genome for alignment    
Fix the names of the contig assembly accessions, removing spaces, commas, '+' or dashes:         
`sed 's/\ /\_/g' 02_input_data/genome/tpac-contigs_min200.fa | cut -d, -f1 | sed 's/\+//g' - | sed 's/\-//g' - > 02_input_data/genome/tpac-contigs_min200_renamed.fa`

Set up genome file as a blast database
`makeblastdb -in ./02_input_data/genome/tpac-contigs_min200_renamed.fa -parse_seqids -dbtype nucl`

Print name of contig, then the length of the contig (tab sep - to be used by R script below)     
`GENOME="02_input_data/genome/tpac-contigs_min200_renamed.fa" ; cat $GENOME | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | tail -n+2 > 02_input_data/genome/ref_genome_seq_lengths.txt`

This will be used by an Rscript below to ensure amplicon windows don't go beyond the contig size. 



### 1. BLAST search loci against genome 
Use BLAST to align the loci to the genome:    
`blastn -db ./02_input_data/genome/tpac-contigs_min200_renamed.fa -query ./02_input_data/prepared_tags.fa -out 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6.txt -outfmt 6 -evalue 1e-20`

Retain only those loci that have fewer than three significant hits:
```
# Find loci names that map more than two times:      
awk '{ print $1 }' 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6.txt | sort -n | uniq -c | sort -nk1 | awk '$1 > 2 { print $2 }' - > 03_blast_out/bad_loci.txt

# Then extract these from the BLAST output:   
grep -vf 03_blast_out/bad_loci.txt 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6.txt > 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6_rem_multimap.txt

# Sort the BLAST output:
sort -k1,1 -k12,12gr -k11,11g -k3,3gr 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6_rem_multimap.txt > 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6_rem_multimap_sort.txt

Keep only a single record per accession:    
DATA="03_blast_out/prepared_tags_v_tpac-contigs_outfmt6_rem_multimap_sort.txt" ; for i in $(cut -f1 $DATA | sort -u); do grep -w -m 1 "$i" $DATA; don
e > 03_blast_out/prepared_tags_v_tpac-contigs_outfmt6_rem_multimap_sort_single_hit.txt
```


### 2. Identify the window range on the ref genome to capture the locus
Use Rscript `01_scripts/identify_req_region_w_BLAST.R`

This will use as inputs the single hit per blast query, and will find the target window, making sure to not go before the start of the reference contig, nor past the end of the reference contig. In these two cases, the start of the window will be 0 or the end of the window will be the length of the contig, respectively. 

This will output:     
`04_extraction/ranges_for_amplicon_extraction.bed`    
`04_extraction/ranges_for_amplicon_extraction.csv`

The bedfile is for extraction, the .csv gives additional useful information for determining identified positions.  

### 3. Collect the amplicons
Collect the target sequence using the bed file (from above) and the genome assembly, using bedtools:      
`GENOME="02_input_data/genome/tpac-contigs_min200_renamed.fa" ; bedtools getfasta -fi $GENOME -bed 04_extraction/ranges_for_amplicon_extraction.bed -fo 04_extraction/tpac_amplicon_approx_by_BLAST.fa`


### 4. Rename amplicons
Turn fasta into tab separated file for ease of renaming amplicons:     
`awk 'BEGIN{RS=">"}{print $1"\t"$2;}' 04_extraction/tpac_amplicon_approx_by_BLAST.fa | tail -n+2 > 04_extraction/tpac_amplicon_approx_by_BLAST.txt`
Note: the tail -n+2 removes an empty first line
(from https://www.biostars.org/p/235052/ )

Use the Rscript `01_scripts/rename_amplicons.R`     


