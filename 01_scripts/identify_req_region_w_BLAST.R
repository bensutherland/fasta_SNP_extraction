# Identify a window around a SNP for loci aligned to a reference genome
# Specifically built for T. pacificus project by Ben Sutherland (DFO) 2018-03-26
# Also applied to other projects starting 2018-10-24
# As usual, use at own risk, no guarantees on usefulness

#### Front Matter ####

# Clean space
#rm(list=ls())

# Todo: make possible to change species name of working directory
#species <- "oket"
species <- "oner"

# Set working directory
#setwd("~/Documents/06_tpac/fasta_SNP_extraction")
#setwd("~/Documents/07_oket/fasta_SNP_extraction")
setwd("~/Documents/08_oner/fasta_SNP_extraction/")

# Install packages
#install.packages("stringr")
library(stringr)

#### 0.a Input data ####
# Import single hit blast outfmt6 results
data.filename <- list.files(path = "03_blast_out/", pattern = "_sort_single_hit.txt", full.names = T)
data <- read.table(data.filename, stringsAsFactors = F)
head(data)
str(data)
colnames(data) <- c("qname","ref.name","ident","align.leng","mismatch","gapopen","qstart","qend","sstart","send","eval","bitscore")
head(data)

# Import contig length file
contig.length <- read.table(file = "02_input_data/genome/ref_genome_seq_lengths.txt", sep = "\t", header = F)
head(contig.length)
colnames(contig.length) <- c("ref.name", "length")
head(contig.length)


#### 0.b Format Input data ####
# Separate marker name from the position of the SNP
split.mname <- strsplit(x = data$qname, split = "-(?=[^-]+$)", perl = TRUE)
str(split.mname)

split.mname.df <- as.data.frame(matrix(unlist(split.mname), nrow=length(split.mname), byrow = T), stringsAsFactors = F)
colnames(split.mname.df) <- c("mname","snp.pos")
head(split.mname.df)

# Combine new split name with original dataframe
data.split <- cbind(data, split.mname.df)
head(data.split)
data.split$checking.col <- substr(x = data.split$qname, start = 1, stop = nchar(data.split$mname))

# Confirm that your cbind retained the correct order
dim(data.split)
table(data.split$mname==data.split$checking.col)["TRUE"] # Data checking
table(data.split$mname==data.split$checking.col)["FALSE"] #should be NA


# Make snp.pos into a numeric value
data.split[,"snp.pos"] <- as.numeric(as.character(data.split$snp.pos))
str(data.split)

# merge data
dim(data.split)
dim(contig.length)
data.collected <- merge(x = data.split, y = contig.length, by = "ref.name")
dim(data.collected)
head(data.collected)

# Rename the new merged, split data as "data"
data <- data.collected


### Characterize input data
length(unique(data$qname))
length(data$qname) # makes sure there are not any duplicates


#### 1. Identify the location of the SNP on the reference genome ####
## Important Note: 
##  alignments may be forward or reverse aligned, but in either case, 
##   the query will be in the forward direction, the subject may invert

# Identify the location of the SNP on the ref genome within the set extraction window
# note that this will make the snp.spot value, which will then be independent of whether the align is for.or.rev
for(i in 1:nrow(data)){
  
  # if subject is forward
  if(data$sstart[i] < data$send[i]){
    data$snp.spot[i] <- (data$sstart[i] + (data$snp.pos[i] - data$qstart[i]))
    data$for.or.rev[i] <- "for"
  } 
  
  # if subject is reverse 
  if(data$sstart[i] > data$send[i]){
    data$snp.spot[i] <- (data$sstart[i] - (data$snp.pos[i] - data$qstart[i]))
    data$for.or.rev[i] <- "rev"
  }
}

head(data, n = 10)


#### 2. Find total.window bp window range ####
# Identify the range, but when the start position is less than total.window/2, use 0 instead. 
total.window <- 400
half.window <- total.window/2

# For future scripts, record window size:
write.table(total.window, file="04_extraction/total_window_size.txt", quote = F)

# rename as data2 (temporary)
data2 <- data

## Avoid going passed the start of the contig:
# How many targets are more than half.window bp in from the start of contig?
table(data2$snp.spot > half.window) # FALSE means that the start of the contig is < half.window

# Identify the range using for loop
for(i in 1:nrow(data2)){
  
  # if snp.spot is greater than half.window bp into ref contig, do the full method
  if(data2$snp.spot[i] > half.window){
    data2$begin.region[i] <- data2$snp.spot[i] - half.window 
    data2$end.region[i] <- data2$snp.spot[i] + half.window 
  } else {
    
  # if snp.spot is less than half.window bp into the contig, then start extraction window at the beginning of the contig
    data2$begin.region[i] <- 1
    data2$end.region[i] <- total.window
  }
}

data3 <- data2
head(data3)

## Avoid going passed the end of the contig:
## Correct the end.regions that go beyond the length of the contig
# How many target windows are beyond the sequence length of the contig?
table(data3$end.region > data3$length)

# When position of end of target window is greater than length of the contig, use from the end to total.window bp before
for(i in 1:nrow(data3)){
  
  if(data3$end.region[i] > data3$length[i]){
    data3$end.region[i] <- data3$length[i]
    data3$begin.region[i] <- (data3$length[i] - (total.window-1))
  }
}

head(data3)

### Avoid when SNP position is outside of the contig
dim(data3)
head(data3)
table(data3$snp.spot < 1)
data3[data3$snp.spot < 1 ,]  # worth investigating this when occurs (#TODO)

# how many rows will be retained (note this shows rows and columns)
dim(data3[data3$snp.spot > 1 , ])

# Retain only the loci where the SNP is within the window
data3 <- data3[data3$snp.spot > 1 , ]
dim(data3)

## This corrects the issue, rather than just dropping it
# ### Change all windows that extend before the start of the contig to instead start at 1
# # also change snp.spot accordingly as this window shifts
# #data.backup <- data3
# tail(head(data3, n = 20), n =10)
# 
# shift.value <- NULL
# 
# for(i in 1:nrow(data3)){
#   if(data3$begin.region[i] < 0){
#     # Find the value by which to shift the begin.region and snp.spot
#     shift.value <- abs(data3$begin.region[i]) + 1
# 
#     # Change the begin.region from a negative to +1 value
#     data3$begin.region[i] <- data3$begin.region[i] + shift.value
# 
#     # In parallel, also change the snp.spot value accordingly
#     data3$snp.spot[i] <- data3$snp.spot[i] + shift.value
# 
#   }
# }
# 
# tail(head(data3, n = 20), n =10)



## Plot where on the reference contig the data is
pdf("04_extraction/start_position_extract_window.pdf", width = 8, height = 5)

hist(data3$begin.region, breaks = 200, main = "Start point of extracting region", xlab = "distance (bp)"
     , las = 1)

# The following turns off the saving out of the image
dev.off()


#### Convert to 0 based number system for bedtools extraction
head(data3)
data3$begin.region <- (data3$begin.region - 1)

# Export whole data file
write.csv(data3, file = "04_extraction/ranges_for_amplicon_extraction.csv", row.names = F, quote = F)

# Export as a bed file extractor
data3.bed <- data3[,c("ref.name","begin.region","end.region")]
head(data3.bed)
write.table(x = data3.bed, file = "04_extraction/ranges_for_amplicon_extraction.bed", sep = "\t"
            , col.names = F, row.names = F, quote = F)
