
#### Front Matter ####
# Clean space
# rm(list=ls())

#### Read in necessary data ####
# Import blast outfmt6 results and remove duplicates from the file
sequence.filename <- paste0("04_extraction/", "amplicon_approx_by_BLAST.txt")
sequence <- read.table(sequence.filename)
colnames(sequence) <- c("match.id", "seq")
head(sequence)
sequence <- sequence[order(sequence[,'match.id']), ]
sequence <- sequence[!duplicated(sequence$match.id), ] # Remove exact duplicates of bedtools extracted segment
dim(sequence) # removes several hits
# Note that it doesn't matter here which segment is deleted, as they are exact duplicates when removed

# Import duplicates to drop
duplicates_to_drop.df <- read.csv(file = "05_amplicons/identifying_duplicate_to_drop.csv", header = T, stringsAsFactors = F)
head(duplicates_to_drop.df)
marker_to_drop <- duplicates_to_drop.df$marker_to_drop # Identify the marker to drop (contains non-priority duplicates)
length(unique(marker_to_drop))

# Import all supplemental information
data_suppl_all.df  <- read.csv(file = "05_amplicons/data_suppl_all_with_duplicates.csv", header = T, stringsAsFactors = F)
head(data_suppl_all.df)
dim(data_suppl_all.df)
# This still contains all duplicate markers

#### Dropping Duplicates ####
# Drop duplicates if they are not priority markers (as defined in previous script)
data_suppl_all_no_priority_dups.df <- data_suppl_all.df[(!data_suppl_all.df$mname %in% marker_to_drop), ]
dim(data_suppl_all.df)
dim(data_suppl_all_no_priority_dups.df)

# Remove the duplicated record from the data based on coordinates
data_suppl_all_no_priority_dups.df <- data_suppl_all_no_priority_dups.df[order(data_suppl_all_no_priority_dups.df[, 'match.id']), ]
data_suppl_all_no_dups.df <- data_suppl_all_no_priority_dups.df[!duplicated(data_suppl_all_no_priority_dups.df$match.id), ]
dim(data_suppl_all_no_dups.df)
head(data_suppl_all_no_dups.df)

#### Combine datasets ####
# Merge two dataframes
data.all <- merge(x = sequence, y = data_suppl_all_no_dups.df, by = "match.id")
dim(data.all)
head(data.all)


#### Find the snp.spot within the extracted window ####
data.all$snp.pos.in.window <- data.all$snp.spot - data.all$begin.region 
head(data.all$snp.pos.in.window)

# connect to the larger df
data.all2 <- data.all
colnames(data.all2)
head(data.all2[, 7:ncol(data.all2)], n = 20)

# Change seq to character
data.all2$seq <- as.character(data.all2$seq)


#### Correct the snp.pos.in.window when reverse complementing ####
# To match the marker file, when the reference genome has been reverse-complemented during the alignment (for.or.rev = rev),
# Calculate where the snp will be in the window after it has been reverse complemented. 

for(i in 1:nrow(data.all2)){
  
  # If alignment is reversed, to prepare for reverse complementing amplicon window, subtract from total.window
  if(data.all2$for.or.rev[i]=="rev"){
    data.all2$snp.pos.in.window[i] <- (nchar(data.all2$seq[i]) + 1) - data.all2$snp.pos.in.window[i]
  }
}

head(data.all2[,7:ncol(data.all2)], n = 20)


## How many reverse or forward orientations were there?
dim(data.all2[data.all2$for.or.rev=="rev", ])
dim(data.all2[data.all2$for.or.rev=="for", ])

# Remove soft masking for rev comp script
data.all2$seq <- toupper(x = data.all2$seq)

# Retain only necessary columns
data.export <- data.all2[,c("match.id", "mname", "snp.pos.in.window", "for.or.rev" , "seq")]
head(data.export)

write.table(x = data.export, file = "05_amplicons/all_fields.txt", sep = "\t"
            , col.names = F, row.names = F, quote = F)
