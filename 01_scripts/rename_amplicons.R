# Rename accessions in a fasta
# Specifically built for T. pacificus project by Ben Sutherland (DFO) 2018-03-26
# Also applied to other species as of 2018-10-24
# As usual, use at own risk, no guarantees on usefulness

#### Front Matter ####
# Clean space
#rm(list=ls())

# Set species shortform
species <- "oket"

# Set working directory
#setwd("~/Documents/06_tpac/fasta_SNP_extraction/")
setwd("~/Documents/07_oket/fasta_SNP_extraction/")

# Install packages
#install.packages("stringr")
library(stringr)
# install.packages("tidyr")
library(tidyr)

## Important: set window size used in 'identify_req_region_w_BLAST.R'
total.window <- as.numeric(read.table(file="04_extraction/total_window_size.txt"))
total.window

#### 0. Input data ####
# Import blast outfmt6 results
data.filename <- paste0("04_extraction/", species, "_amplicon_approx_by_BLAST.txt")
data <- read.table(data.filename)
head(data)
colnames(data) <- c("match.id", "seq")
head(data)

# Remove exact duplicates of aligned markers
dim(data)
data <- data[order(data[,'match.id']), ]
data <- data[!duplicated(data$match.id), ]
dim(data) # removes several hits

# Note that it doesn't matter which segment is deleted, as they are exact duplicates when removed

# can separate the match.id if needed
#data2 <- separate(data = data, col = "contig.name", into = c("contig", "range"), sep = "\\:")
#head(data2)

# Import other information to attach
data.suppl <- read.csv(file = "04_extraction/ranges_for_amplicon_extraction.csv")
head(data.suppl)
match.id <- paste(data.suppl$ref.name,":",data.suppl$begin.region,"-",data.suppl$end.region, sep = "")
data.suppl.all <- cbind(data.suppl, match.id)
head(data.suppl.all)

# As above, remove exact duplicates of aligned markers
data.suppl.all <- data.suppl.all[order(data.suppl.all[, 'match.id']), ]
data.suppl.all <- data.suppl.all[!duplicated(data.suppl.all$match.id), ]
dim(data.suppl.all)

# Note that this matters more about which is deleted, because certain markers are more useful than others

# Merge two dataframes
data.all <- merge(x = data, y = data.suppl.all, by = "match.id")
dim(data)
dim(data.suppl.all)
dim(data.all)
head(data.all)


#### Find the snp.spot within the extracted window ####
snp.pos.in.window <- data.all$snp.spot - data.all$begin.region 
head(snp.pos.in.window)

# connect to the larger df
data.all2 <- cbind(data.all, snp.pos.in.window)
colnames(data.all2)
head(data.all2[,7:ncol(data.all2)], n = 20)


#### Correct the snp.pos.in.window when reverse complementing ####
# To match the marker file, when the reference genome has been reverse-complemented during the alignment (for.or.rev = rev),
# Calculate where the snp will be in the window after it has been reverse complemented. 

for(i in 1:nrow(data.all2)){
  # If alignment is reversed, to prepare for reverse complementing amplicon window, subtract from total.window
  if(data.all2$for.or.rev[i]=="rev"){
    data.all2$snp.pos.in.window[i] <- total.window - data.all2$snp.pos.in.window[i]
  }
}

head(data.all2[,7:ncol(data.all2)], n = 20)


###### Correct the specific issue of reverse complement SNP being at middle position
colnames(data.all2)
data.all2[ which(data.all2$for.or.rev=="rev" & data.all2$snp.pos.in.window=="200"), "snp.pos.in.window" ] <- 201


###### Separate output for reverse complement target regions and forward regions
dim(data.all2[data.all2$for.or.rev=="rev", ])
dim(data.all2[data.all2$for.or.rev=="for", ])


# Get the needed pieces for the name
data.export <- data.all2[,c("match.id", "mname", "snp.pos.in.window", "for.or.rev" , "seq")]

write.table(x = data.export, file = "05_amplicons/all_fields.txt", sep = "\t"
            , col.names = F, row.names = F, quote = F)
