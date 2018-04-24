# Renaming fasta
# Clean space
#rm(list=ls())

# Set working directory
setwd("~/Documents/06_tpac")

# Install packages
#install.packages("stringr")
library(stringr)
# install.packages("tidyr")
library(tidyr)


#### 0. Input data ####
# Import blast outfmt6 results
data <- read.table("tpac_amplicon_approx_by_BLAST.txt")
head(data)
colnames(data) <- c("match.id", "seq")
head(data)
# can separate the match.id if needed
#data2 <- separate(data = data, col = "contig.name", into = c("contig", "range"), sep = "\\:")
#head(data2)


# Import other information to attach
data.suppl <- read.csv(file = "ranges_for_amplicon_extraction.csv")
head(data.suppl)
match.id <- paste(data.suppl$ref.name,":",data.suppl$begin.region,"-",data.suppl$end.region, sep = "")
data.suppl.all <- cbind(data.suppl, match.id)
head(data.suppl.all)


# Merge two dataframes
data.all <- merge(x = data, y = data.suppl.all, by = "match.id")


# Determine exactly where the SNP is on the extracted window
snp.pos.in.window <- data.all$snp.spot - data.all$begin.region
head(snp.pos.in.window)

data.all2 <- cbind(data.all, snp.pos.in.window)

colnames(data.all2)
#


# Get the needed pieces for the name
data.export <- data.all2[,c("match.id", "mname", "snp.pos.in.window", "seq")]

write.table(x = data.export, file = "all_fields.txt", sep = "\t"
            , col.names = F, row.names = F, quote = F)
