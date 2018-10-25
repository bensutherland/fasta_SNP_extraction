# Identify a window around a SNP for loci aligned to a reference genome
# Specifically built for T. pacificus project by Ben Sutherland (DFO) 2018-03-26
# Also applied to other projects starting 2018-10-24
# As usual, use at own risk, no guarantees on usefulness

#### Front Matter ####

# Clean space
#rm(list=ls())

# Todo: make possible to change species name of working directory
species <- "oket"

# Set working directory
#setwd("~/Documents/06_tpac/fasta_SNP_extraction")
setwd("~/Documents/07_oket/fasta_SNP_extraction")

# Install packages
#install.packages("stringr")
library(stringr)

#### 0.a Input data ####
# Import single hit blast outfmt6 results
data.filename <- paste0("03_blast_out/prepared_tags_v_", species, "-contigs_outfmt6_rem_multimap_sort_single_hit.txt")
data <- read.table(data.filename)
head(data)
colnames(data) <- c("qname","ref.name","ident","align.leng","mismatch","gapopen","qstart","qend","sstart","send","eval","bitscore")
head(data)

# Import contig length file
contig.length <- read.table(file = "02_input_data/genome/ref_genome_seq_lengths.txt", sep = "\t", header = F)
head(contig.length)
colnames(contig.length) <- c("ref.name", "length")
head(contig.length)


#### 0.b Format Input data ####
# Separate marker name from the position of the SNP
split.mname <- str_split_fixed(string = data$qname, pattern = "-", n = 2)
colnames(split.mname) <- c("mname","snp.pos")
head(split.mname)

# Combine new split name with original dataframe
data.split <- cbind(data, split.mname)
head(data.split)

# Make snp.pos into a numeric value
data.split[,"snp.pos"] <- as.numeric(as.character(data.split$snp.pos))
str(data.split)

# merge data
dim(data.split)
dim(contig.length)
data.collected <- merge(x = data.split, y = contig.length, by = "ref.name")
dim(data.collected) # should be 3928
head(data.collected)

# Rename the new merged, split data as "data"
data <- data.collected


### Characterize input data
length(unique(data$qname))
length(data$qname)


#### 1. Identify the location of the SNP on the reference genome ####
## Important Note: 
##  alignments may be forward or reverse aligned, but in either case, 
##   the query will be in the forward direction, the subject may invert

# Identify the location of the SNP on the ref genome within the set extraction window
# note that this will make the snp.spot value, which will then be independent of whether the align is for.or.rev
snp.spot <- rep(NA, times = nrow(data))
for.or.rev <- rep(NA, times = nrow(data))

for(i in 1:nrow(data)){
  
  # if subject is forward
  if(data$sstart[i] < data$send[i]){
    snp.spot[i] <- (data$sstart[i] + (data$snp.pos[i] - data$qstart[i]))
    for.or.rev[i] <- "for"
  } 
  
  # if subject is reverse 
  if(data$sstart[i] > data$send[i]){
    snp.spot[i] <- (data$sstart[i] - (data$snp.pos[i] - data$qstart[i]))
    for.or.rev[i] <- "rev"
  }
}

head(snp.spot)
head(for.or.rev)

# Attach to the rest of the data
head(data)
all.data <- cbind(data, snp.spot, for.or.rev)

head(all.data)


#### 2. Find X bp window range ####
# Identify the range, but when the start position is less than X/2, use 0 instead. 
total.window <- 400
half.window <- total.window/2

# For future scripts, record window size:
write.table(total.window, file="04_extraction/total_window_size.txt", quote = F)

# rename as data2 (temporary)
data2 <- all.data

## Avoid going passed the start of the contig:
# How many targets are more than half.window bp in from the start of contig?
table(data2$snp.spot > half.window)

# Set nulls
begin.region <- rep(NA, times = nrow(data2))
end.region <- rep(NA, times = nrow(data2))

# Identify the range using for loop
for(i in 1:nrow(data2)){
  
  # if snp.spot is greater than half.window bp into ref contig, do the full method
  if(data2$snp.spot[i] > half.window){
    begin.region[i] <- data2$snp.spot[i] - half.window 
    end.region[i] <- data2$snp.spot[i] + half.window 
  } else {
    
  # if snp.spot is less than half.window bp into the contig, then start extraction window at the beginning of the contig
    begin.region[i] <- 1
    end.region[i] <- total.window
  }
}

# the following can be used to confirm the math
head(begin.region, n = 20)
head(data2$snp.spot, n = 20)
head(end.region, n = 20)

# Collect into datafile
data3 <- cbind(data2, begin.region, end.region)

head(data3)

## Avoid going passed the end of the contig:
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

# how many rows will be retained (note this shows rows and columns)
dim(data3[data3$snp.spot > 1 , ])

# Retain only the loci where the SNP is within the window
data3 <- data3[data3$snp.spot > 1 , ]


### Change all windows that extend before the start of the contig to instead start at 1
# also change snp.spot accordingly as this window shifts
#data.backup <- data3
tail(head(data3, n = 20), n =10)

shift.value <- NULL

for(i in 1:nrow(data3)){
  if(data3$begin.region[i] < 0){
    # Find the value by which to shift the begin.region and snp.spot
    shift.value <- abs(data3$begin.region[i]) + 1
    
    # Change the begin.region from a negative to +1 value
    data3$begin.region[i] <- data3$begin.region[i] + shift.value
    
    # In parallel, also change the snp.spot value accordingly
    data3$snp.spot[i] <- data3$snp.spot[i] + shift.value
                               
  }
}

tail(head(data3, n = 20), n =10)



## Plot where on the reference contig the data is
pdf("04_extraction/start_position_extract_window.pdf", width = 8, height = 5)

hist(data3$begin.region, breaks = 200, main = "Start point of extracting region", xlab = "distance (bp)"
     , las = 1)

# The following turns off the saving out of the image
dev.off()


# Export whole data file
write.csv(data3, file = "04_extraction/ranges_for_amplicon_extraction.csv", row.names = F, quote = F)

# Export as a bed file extractor
data3.bed <- data3[,c("ref.name","begin.region","end.region")]
head(data3.bed)
write.table(x = data3.bed, file = "04_extraction/ranges_for_amplicon_extraction.bed", sep = "\t"
            , col.names = F, row.names = F, quote = F)

