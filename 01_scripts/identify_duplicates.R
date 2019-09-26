# Identify duplicate based on genomic coordinates
# As usual, use at own risk, no guarantees on usefulness

#### Front Matter ####
# Clean space
# rm(list=ls())

# Set species shortform
#species <- "oket"
species <- "oner"

# Set working directory
#setwd("~/Documents/06_tpac/fasta_SNP_extraction/")
#setwd("~/Documents/07_oket/fasta_SNP_extraction/")
setwd("~/Documents/08_oner/fasta_SNP_extraction/")

# Install packages
#install.packages("stringr")
library(stringr)
# install.packages("tidyr")
library(tidyr)

## Important: set window size used in 'identify_req_region_w_BLAST.R'
total.window <- as.numeric(read.table(file="04_extraction/total_window_size.txt"))
total.window

#### 0. Input data ####
# Import other information about the amplicon to attach
data.suppl <- read.csv(file = "04_extraction/ranges_for_amplicon_extraction.csv")
head(data.suppl)
# Create a matching identifier to attach to the main data
data.suppl$match.id <- paste(data.suppl$ref.name,":",data.suppl$begin.region,"-",data.suppl$end.region, sep = "")
data.suppl.all <- data.suppl
head(data.suppl.all)

data.suppl.all.before_drop <- data.suppl.all

# Save this file for re-use later
write.csv(x = data.suppl.all, file = "05_amplicons/data_suppl_all_with_duplicates.csv", quote = F, row.names = F)

# Remove the duplicated record from the data
data.suppl.all <- data.suppl.all[order(data.suppl.all[, 'match.id']), ]
data.suppl.all <- data.suppl.all[!duplicated(data.suppl.all$match.id), ]
dim(data.suppl.all)
head(data.suppl.all)

## Identify the names of the amplicons that will be dropped due to duplication
dropped_duplicates.df <- data.suppl.all.before_drop[!data.suppl.all.before_drop$mname %in% data.suppl.all$mname, c("mname", "qname", "match.id")]
dim(dropped_duplicates.df)
head(dropped_duplicates.df)
write.csv(x = dropped_duplicates.df, file = "05_amplicons/dropped_duplicates.csv", row.names = F, quote = F)

# As some marker sources are more valuable than others, we need to decide which of the duplicates to delete
# use the following file for this purpose
# Saved out after the duplicates have been removed
write.csv(x = data.suppl.all, file = "05_amplicons/dropped_duplicates_suppl_info.csv", row.names = F, quote = F)

### Go to next script