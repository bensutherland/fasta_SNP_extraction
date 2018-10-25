# This script is used to collect various parts of the data and assemble it to 
# make the submission form requested by ThermoFisher

#rm(list=ls())

# Set working directory
setwd("~/Documents/07_oket/fasta_SNP_extraction")
# setwd("~/Documents/06_tpac/fasta_SNP_extraction_400bp_second/")

# Set species
species <- "oket"

##### 1. Import Data ####
#### 1.A) Use the RAD tags file to get alleles ####
rad.tags.filename <- paste0("06_output/", species, "_all_amplicons_rad_tags.csv")
rad.tags <- read.csv(file = rad.tags.filename , header = F
                     , col.names = c("radtag","radlocus","snp.pos"))
head(rad.tags)
rad.tags <- rad.tags[,c("radtag","radlocus")] # only keep necessary columns
head(rad.tags)
dim(rad.tags)
length(unique(rad.tags$radtag))

# Earlier in the pipeline, an allele file was created from the first SNP in the above RAD tags file (via awk)
# Note this must still correspond to the order of the above 'tpac_all_amplicons_rad_tags.csv' file
alleles <- read.csv(file = "06_output/alleles_only.txt", header = F
                    , col.names = c("ref", "var"))
head(alleles)

# Combine these two files
rad.tags.and.alleles <- cbind(rad.tags, alleles)
head(rad.tags.and.alleles, n = 10)
# Confirm they still match (i.e. were in the same order)

#### 1. B) Use the amplicon text file to get the amplicon sequence and the position of the SNP (in amplicon window) from accn name ####
fasta.txt.filename <- paste0("06_output/", species, "_all_amplicons.txt")
fasta.txt <- read.table(file = fasta.txt.filename, header = F
                    , sep = "\t", col.names = c("accn.name","seq"))
head(fasta.txt)

# i) For matching purposes, extract the radtag name from accn.name
# The radtag name is the second field when separated by '__'
# Assumes four fields: contig name, radtag name, position of snp, original orientation
split.accn.name <- str_split_fixed(string = fasta.txt$accn.name, pattern = "__", n = 4)
radtag <- as.character(split.accn.name[, 2])
head(radtag)
length(radtag)

# OLD: Original for Tpac, not extensible. 
# radtag <-      gsub(x = 
#                       gsub(x = fasta.txt$accn.name, pattern = ".*Tpa", replacement = "Tpa")
#                             , pattern = "\\__.*" , replacement = "")
# head(radtag)

# Combine this radtag name to the original fasta.txt
fasta.txt <- cbind(radtag, fasta.txt)
head(fasta.txt)


# ii) Extract the within-window SNP position from the accession name
snp.pos.in.window <- as.numeric(as.character(split.accn.name[,3]))
head(snp.pos.in.window)

## OLD: Not very extensible
# step1 <- gsub(x = fasta.txt$accn.name, pattern = ".*Tpa", replacement = "Tpa") # remove everything before Tpa
# head(step1)
# step2 <- gsub(x = step1, pattern = "\\__rev|\\__for", replacement = "") # remove the end 'for' or 'rev'
# head(step2)
# step3 <- gsub(x = step2, pattern = ".*\\__", replacement = "") # remove everything before the "__"
# head(step3)
# snp.pos.in.window <- step3
# rm(step1, step2, step3)
# head(snp.pos.in.window)

# Combine this window position to the original fasta.txt
fasta.txt <- cbind(fasta.txt, snp.pos.in.window)
head(fasta.txt)


#### 1.C) Use the all_fields.txt file to get the scaffold identity ####
all.fields <- read.table(file = "05_amplicons/all_fields.txt", sep = "\t", header = F
                         , col.names = c("bed.extract", "radtag", "old.pos","orient", "incor.amp"))
all.fields <- all.fields[,c("bed.extract","radtag")] # retain only necessary info
# note cannot use the sequence, as has not been reversed. Can use the marker name though
head(all.fields)

# split at the colon:
scaff <- gsub(all.fields$bed.extract, pattern = "\\:.*", replacement = "")
head(scaff)

# Combine this scaff to the all.fields file
scaff.info <- cbind(all.fields, scaff)
scaff.info <- scaff.info[,c("radtag","scaff")]

head(scaff.info)


#### 2. Join All Info #####
# Sort all data by radtag
rad.tags.and.alleles.sorted <- rad.tags.and.alleles[order(rad.tags.and.alleles$radtag),]
fasta.txt.sorted  <- fasta.txt[order(fasta.txt$radtag),]
scaff.info.sorted <- scaff.info[order(scaff.info$radtag), ]

# Observe objects, composition and size
head(rad.tags.and.alleles.sorted)
head(fasta.txt.sorted)
head(scaff.info.sorted)
dim(rad.tags.and.alleles.sorted)
dim(fasta.txt.sorted)
dim(scaff.info.sorted)

# Merge data
submission.df <- merge(x=rad.tags.and.alleles.sorted, y = fasta.txt.sorted, by = "radtag")
head(submission.df)
colnames(submission.df)

submission.df <- merge(x=submission.df, y=scaff.info.sorted, by = "radtag")
colnames(submission.df)

head(submission.df, n = 5)

## Add extra needed into for the submission sheet
strand <- rep("NA", times = nrow(submission.df))
mtype <- rep("SNP", times = nrow(submission.df))
priority <- rep("first.priority", times = nrow(submission.df))

submission.df <- cbind(submission.df,strand,mtype,priority)

# Select columns needed
#submission.df <- submission.df[,c("accn.name","scaff","snp.pos.in.window","snp.pos.in.window","ref","var","strand","mtype","seq")]

head(submission.df)

save.image(file = "06_output/submission_data_part_1.RData")

# Next, go to "01_scripts/figuring_out_another_method_of_inserting_alleles.R"

#write.csv(x = submission.df, file = "06_output/submission.csv", quote = F, row.names = F)

