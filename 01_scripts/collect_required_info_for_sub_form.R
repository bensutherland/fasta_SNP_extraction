# This script is used to collect various parts of the data and assemble it to 
# make the submission form requested by ThermoFisher

#rm(list=ls())

setwd("~/Documents/06_tpac/fasta_SNP_extraction_400bp/")

##### ALLELES ####
# Use the RAD tags file to get alleles
rad.tags <- read.csv(file = "06_output/tpac_amplicon_panel_rad_tags.csv", header = F
                     , col.names = c("radtag","radlocus","snp.pos"))
head(rad.tags)
rad.tags <- rad.tags[,c(1:2)] # only keep necessary columns
head(rad.tags)

# Use the allele file created from the RAD tags file via awk (first SNP of locus)
alleles <- read.csv(file = "06_output/alleles_only.txt", header = F
                    , col.names = c("ref", "var"))
head(alleles)

# Combine these two files
rad.tags.and.alleles <- cbind(rad.tags, alleles)
head(rad.tags.and.alleles)


##### AMPLICON SEQUENCE, snp.pos.in.window #####
# Use the fasta file to get sequence
fasta.txt <- read.table(file = "06_output/tpac_amplicon_panel_v0.2.txt", header = F
                    , sep = "\t", col.names = c("accn.name","seq"))
head(fasta.txt)

##
# simple.marker.name <- gsub(fasta.txt$V1, pattern = "\\__.*", replacement = "")
# head(simple.marker.name)

# Attach simple marker name to the fasta
# fasta.info <- cbind(fasta.txt, simple.marker.name)
# head(fasta.info)
##

# Extract the radtag name from accn.name the radtag (for matching purposes)
radtag <-      gsub(x = 
                               gsub(x = fasta.txt$accn.name, pattern = ".*Tpa", replacement = "Tpa")
                            , pattern = "\\__.*" , replacement = "")
head(radtag)

# Attach the radtag name to the fasta.txt
fasta.txt <- cbind(radtag, fasta.txt)
head(fasta.txt)


## Extract the position of the SNP from the accn.name
step1 <- gsub(x = fasta.txt$accn.name, pattern = ".*Tpa", replacement = "Tpa") # remove everything before Tpa
head(step1)
step2 <- gsub(x = step1, pattern = "\\__rev|\\__for", replacement = "") # remove the end 'for' or 'rev'
head(step2)
step3 <- gsub(x = step2, pattern = ".*\\__", replacement = "") # remove everything before the "__"
snp.pos.in.window <- step3

# Attach the position of the snp to the fasta.txt
fasta.txt <- cbind(fasta.txt, snp.pos.in.window)
head(fasta.txt)


##### OBTAIN OTHER INFO ####
all.fields <- read.table(file = "05_amplicons/all_fields.txt", sep = "\t", header = F
                         , col.names = c("bed.extract", "radtag", "old.pos","orient", "incor.amp"))
all.fields <- all.fields[,c("bed.extract","radtag")]
# note cannot use the sequence, as has not been reversed. Can use the marker name though

# split at the colon:
scaff <- gsub(all.fields$bed.extract, pattern = "\\:.*", replacement = "")
head(scaff)

#
# # Make the marker ID to match the fasta file
# marker.id <- gsub(all.fields$bed.extract, pattern = "\\:", replacement = "_")
# head(marker.id)
#

# combine
extra.info <- cbind(all.fields, scaff)

head(extra.info)


#### NOW COLLECT TOGETHER #####
# Sort data
rad.tags.and.alleles.sorted <- rad.tags.and.alleles[order(rad.tags.and.alleles$radtag),]
fasta.txt.sorted  <- fasta.txt[order(fasta.txt$radtag),]
extra.info.sorted <- extra.info[order(extra.info$radtag), ]

head(rad.tags.and.alleles.sorted)
head(fasta.txt.sorted)
head(extra.info.sorted)
dim(rad.tags.and.alleles.sorted)
dim(fasta.txt.sorted)
dim(extra.info.sorted)

# Merge data
submission.df <- merge(x=rad.tags.and.alleles.sorted, y = fasta.txt.sorted, by = "radtag")
colnames(submission.df)

submission.df <- merge(x=submission.df, y=extra.info.sorted, by = "radtag")
colnames(submission.df)

head(submission.df, n = 5)

## Add extra needed into
strand <- rep("NA", times = nrow(submission.df))
mtype <- rep("SNP", times = nrow(submission.df))
priority <- rep("first.priority", times = nrow(submission.df))

submission.df <- cbind(submission.df,strand,mtype,priority)

# Select columns needed
submission.df <- submission.df[,c("accn.name","scaff","snp.pos.in.window","snp.pos.in.window","ref","var","strand","mtype","seq")]

head(submission.df)

write.csv(x = submission.df, file = "06_output/submission.csv", quote = F, row.names = F)

