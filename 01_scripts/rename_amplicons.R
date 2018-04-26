# Rename accessions in a fasta
# Specifically built for T. pacificus project by Ben Sutherland (DFO) 2018-03-26
# As usual, use at own risk, no guarantees on usefulness

#### Front Matter ####
# Clean space
#rm(list=ls())

# Set working directory
setwd("~/Documents/06_tpac/fasta_SNP_extraction/")

# Install packages
#install.packages("stringr")
library(stringr)
# install.packages("tidyr")
library(tidyr)




#### 0. Input data ####
# Import blast outfmt6 results
data <- read.table("04_extraction/tpac_amplicon_approx_by_BLAST.txt")
head(data)
colnames(data) <- c("match.id", "seq")
head(data)

# remove duplicates
data <- data[order(data[,'match.id']), ]
data <- data[!duplicated(data$match.id), ]
dim(data)

# can separate the match.id if needed
#data2 <- separate(data = data, col = "contig.name", into = c("contig", "range"), sep = "\\:")
#head(data2)

# Import other information to attach
data.suppl <- read.csv(file = "04_extraction/ranges_for_amplicon_extraction.csv")
head(data.suppl)
match.id <- paste(data.suppl$ref.name,":",data.suppl$begin.region,"-",data.suppl$end.region, sep = "")
data.suppl.all <- cbind(data.suppl, match.id)
head(data.suppl.all)


# remove duplicates
data.suppl.all <- data.suppl.all[order(data.suppl.all[, 'match.id']), ]
data.suppl.all <- data.suppl.all[!duplicated(data.suppl.all$match.id), ]
dim(data.suppl.all)


# Merge two dataframes
data.all <- merge(x = data, y = data.suppl.all, by = "match.id")
dim(data)
dim(data.suppl.all)
dim(data.all)
head(data.all)

# Determine exactly where the SNP is on the extracted window
snp.pos.in.window <- data.all$snp.spot - data.all$begin.region 
head(snp.pos.in.window)

data.all2 <- cbind(data.all, snp.pos.in.window)

colnames(data.all2)
head(data.all2[,7:ncol(data.all2)], n = 20)


###### Correct the snp.pos.in.window when reverse complementing
for(i in 1:nrow(data.all2)){
  # Only correct if reverse
  if(data.all2$for.or.rev[i]=="rev"){
    data.all2$snp.pos.in.window[i] <- 200 - data.all2$snp.pos.in.window[i]
  }
}

head(data.all2[,7:ncol(data.all2)], n = 20)


###### Correct the specific issue of reverse complement SNP being at position 101 when title suggests 100
colnames(data.all2)
data.all2[ which(data.all2$for.or.rev=="rev" & data.all2$snp.pos.in.window=="100"), "snp.pos.in.window" ] <- 101


###### Separate output for reverse complement target regions and forward regions
dim(data.all2[data.all2$for.or.rev=="rev", ]) # 1935
dim(data.all2[data.all2$for.or.rev=="for", ]) # 1971


# Get the needed pieces for the name
data.export <- data.all2[,c("match.id", "mname", "snp.pos.in.window", "for.or.rev" , "seq")]

write.table(x = data.export, file = "05_amplicons/all_fields.txt", sep = "\t"
            , col.names = F, row.names = F, quote = F)
