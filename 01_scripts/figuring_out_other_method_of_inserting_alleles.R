# Identify the location just before or just after the polymorphism in the amplicon
# Assumes that your rad locus file and amplicon file are in the same orientation

# rm(list=ls())

# Set working directory
setwd("~/Documents/07_oket/fasta_SNP_extraction")
# setwd("~/Documents/06_tpac/fasta_SNP_extraction_400bp_third_try")

# Load libraries
require(tidyr)
require(stringr)

# Import data
load(file = "06_output/submission_data_part_1.RData")

#### 1. Isolate the first polymorphism in the RAD tag ####
head(rad.tags)
# side: fix an error in oket file:
rad.tags$radlocus <- gsub(pattern = "}", replacement = "]", x = rad.tags$radlocus)
dim(rad.tags)

# Take out the section before the first square brackets
first.piece.extracted <- separate(data = rad.tags
                                  , col = "radlocus", into = c("before.snp", "after.snp"), sep = "\\[" )
# expect warnings as it will read up to the second instance then stop
head(first.piece.extracted)
colnames(first.piece.extracted)
dim(first.piece.extracted)

# Take out the section after the first square brackets
allele_isolated <- separate(data = first.piece.extracted
                            , col = "after.snp", into = c("allele","after.snp"), sep ="\\]")
head(allele_isolated, n = 2)
colnames(allele_isolated)
dim(allele_isolated)

#### 2. Determine which flanking region is longer in the rad locus (before or after) ####
# For loop to find whether before or after flanking sequence is longer
side.to.use <- NULL

for(i in 1:nrow(allele_isolated)){
  if(nchar(allele_isolated$before.snp[i]) > nchar(allele_isolated$after.snp[i])){
    side.to.use[i] <- "before"
  } else {
    side.to.use[i] <- "after"
  }
}

head(side.to.use)
length(side.to.use)

# Combine this specification to the isolated alleles
allele_isolated <- cbind(allele_isolated, side.to.use)
head(allele_isolated)
table(allele_isolated$side.to.use)

###### HAND EDIT SECTION #####
# This is a hacky solution, but will allow you to edit whether you will use the before or after manually, if it fails later #
# to.edit <- "Oki_RAD41030_31"
failed.amplicons

failed.of.interest <- NULL
for(i in 1:length(failed.amplicons)){
  
  # identify which row working on
  failed.of.interest <- failed.amplicons[i]
  print(failed.of.interest)
  
  # If statement to change before to after..
  if( allele_isolated[which(allele_isolated$radtag==failed.of.interest), "side.to.use"] == "before" ){
    allele_isolated[which(allele_isolated$radtag==failed.of.interest), "side.to.use"] <- "after"
  } else { allele_isolated[which(allele_isolated$radtag==failed.of.interest), "side.to.use"] <- "before"
  }
}


#### 3. Obtain a str.size vector for matching the region in amplicon ####
# Set size of the matching vector string
str.size <- 20

# Use for loop to obtain the matching vector from either before or after the polymorphism
extract.piece <- NULL
marker.vector <- NULL

for(i in 1:nrow(allele_isolated)){
  
  # if the 'side.to.use' for vector matching is after the polymorphism
  if(allele_isolated$side.to.use[i]=="after"){
    
    print("Using after piece")
    
    # then obtain 'str.size' characters after the allele (from 1 to str.size)
    extract.piece <- str_sub(allele_isolated$after.snp[i], start = 1, end = str.size)
    
    print(extract.piece) } else {
      
      # if before is longer, use before
      print("Using before piece")
      
      # then obtain 'str.size' characters before the allele (from -str.size to -1)
      extract.piece <- str_sub(allele_isolated$before.snp[i], start = -str.size , end = -1)
      
      print(x = extract.piece)
    }
  # put the extract.piece for this round into a vector
  marker.vector[i] <- extract.piece
  }

# Combine the data
all.info.needed <- cbind(allele_isolated, marker.vector, allele_isolated$radtag)

head(all.info.needed)


#### 4. Using the matching string, identify where the polymorphism is in the amplicon ####

## Set up ##
# head(fasta.txt) # file with the amplicon
head(submission.df) # object w/ amplicon
head(all.info.needed) # file with the matching bit

# combine these two objects
complete_info <- merge(x = submission.df, y = all.info.needed, by = "radtag")
head(complete_info)

# Check formats
str(complete_info)
# Make vector into character
complete_info$seq <- as.character(complete_info$seq)
complete_info$marker.vector <- as.character(complete_info$marker.vector)
complete_info$side.to.use <- as.character(complete_info$side.to.use)
str(complete_info)


## Loop to split up the amplicon and insert the allele ##
temp.record <- NULL
complete_info_w_amplicon <- NULL
match.success <- NULL

for(i in 1:nrow(complete_info)){
  
  # Take one record at a time
  temp.record <- complete_info[i,]
  print(paste("Working on ", temp.record$radtag))
  
  # If 'marker.vector' is expected after SNP
  if(temp.record$side.to.use=="after"){
    print("Using an 'After' marker")
    
    # Separate seq column into before and after marker
    temp.record <- separate(data = temp.record, col="seq", into = c("before_marker", "after_marker"), sep = temp.record$marker.vector, remove = F)
    print(temp.record)
    
    # Record whether it worked or not
    if(is.na(temp.record$after_marker)==TRUE){ 
      match.success[i] <- "fail"
    } else {
        match.success[i] <- "success"
      }
    
    # Paste the marker.vector back into the 'after marker' section (this is dropped during separate)
    temp.record$after_marker <- paste(temp.record$marker.vector, temp.record$after_marker, sep = "")

    # Remove the last character from 'before marker', as this is the SNP
    temp.record$before_marker <- substr(x = temp.record$before_marker, start = 1, stop = (nchar(temp.record$before_marker)-1) )
    
  } else {
  # If 'marker.vector' is expected before SNP
    
    # Separate seq column into before and after marker
    temp.record <- separate(data = temp.record, col="seq", into = c("before_marker", "after_marker"), sep = temp.record$marker.vector, remove = F)
    
    # Record whether it worked or not
    if(is.na(temp.record$after_marker)==TRUE){ 
      match.success[i] <- "fail"
    } else {
      match.success[i] <- "success"
    }
    
    # Paste the marker.vector back into the 'before marker' section (this is dropped during separate)
    temp.record$before_marker <- paste(temp.record$before_marker, temp.record$marker.vector, sep = "")
    
    # Remove the first character from 'after marker', as this is the SNP
    temp.record$after_marker <- substr(x = temp.record$after_marker, start = 2, stop = nchar(temp.record$after_marker) )
    
  }
  
  # Bring it all back together into an amplicon
  temp.record$amplicon_complete <- paste(temp.record$before_marker, "[", temp.record$allele, "]", temp.record$after_marker, sep = "")
  
  complete_info_w_amplicon <- rbind(complete_info_w_amplicon, temp.record)
  print("Good Job!")
}

# Check for matching fails, and if so, go back to just after step 2 and change 'before' to 'after' or visa versa.
table(match.success) # 19 fails, 3369 successes
# complete_info_w_amplicon[which(match.success=="fail"),]

# failed.matches <- complete_info_w_amplicon[which(match.success=="fail"), "radtag"] # these were the failing markers
# failed.amplicons <- as.character(droplevels(failed.matches))
# failed.amplicons

# Failed amplicons may manifest as NAs during the splitting, but more reliably, as a substantial difference from the 
# expected location of the snp in the window as taken from the accession name.


##### 5. Final check of the expected snp position and the actual snp position #####
colnames(complete_info_w_amplicon)

# Recalculate what position in the amplicon window has the polymorphism
final_check <- separate(data = complete_info_w_amplicon, col = "amplicon_complete", into = c("recalc_seq_before_snp", "dropseq"), sep = "\\[" , remove = F)
head(final_check)
final_check$nchar_before_snp <- nchar(final_check$recalc_seq_before_snp)
final_check$recalc_pos_snp_in_window <- final_check$nchar_before_snp + 1 
colnames(final_check)
head(final_check, n = 5)

# Keep the recal_pos_snp_in_window as an object
colnames(final_check)
necess.posit.info <- final_check[, c("radtag", "snp.pos.in.window","recalc_pos_snp_in_window")]

str(final_check)
final_check$snp.pos.in.window  <- as.numeric(as.character(final_check$snp.pos.in.window))
str(final_check)

# Calculate how far off the expected and observed snp position are
number.digits.off <- final_check$recalc_pos_snp_in_window - final_check$snp.pos.in.window
number.digits.off
failed.matches2 <- final_check[which(abs(number.digits.off) > 5 ), "radtag" ]
failed.matches2 <- as.character(droplevels(failed.matches2))
failed.matches2

#### 6. Only keep the successful insertions ####
colnames(complete_info_w_amplicon)
dim(complete_info_w_amplicon)

# Remove the failed matches from the object
j <- NULL

for(i in 1:length(failed.matches2)){
  j <- failed.matches2[i]
  print(j)

  # remove the failed match amplicon from the dataframe
  complete_info_w_amplicon <- complete_info_w_amplicon[ -c(which(complete_info_w_amplicon[,"radtag"]==j) ) , ]
}


dim(complete_info_w_amplicon)


#### 7. Choose only the top Fst markers ####
### THE FOLLOWING COMMENTED LINES WERE SPECIFIC TO TPAC ###
# adaptive.loci <- read.csv("02_input_data/adaptive_loci_fst_for_amplicon_selection.csv", col.names = c("radtag","Fst"))
# neutral.loci <- read.csv("02_input_data/neutral_loci_fst_for_amplicon_selection.csv", col.names = c("radtag","Fst"))
# head(adaptive.loci)
# head(neutral.loci)
# 
# # merge w/ existing complete_info_w_amplicons
# complete_info_w_amplicon_adaptive <- merge(x = complete_info_w_amplicon, y = adaptive.loci, by = "radtag")
# dim(adaptive.loci) # total number adaptive loci
# dim(complete_info_w_amplicon_adaptive) # total number that are found in the complete amplicons
# 
# complete_info_w_amplicon_neutral <- merge(x = complete_info_w_amplicon, y = neutral.loci, by = "radtag")
# dim(neutral.loci) # total number neutral loci
# dim(complete_info_w_amplicon_neutral)
# head(complete_info_w_amplicon_neutral) # note that Fst is no longer in order, re-order by Fst
# # Sort by Fst
# complete_info_w_amplicon_neutral <- complete_info_w_amplicon_neutral[order(complete_info_w_amplicon_neutral$Fst, decreasing = T),]
# head(complete_info_w_amplicon_neutral)
# 
# # How many adaptives do you have?
# dim(complete_info_w_amplicon_adaptive)[1]
# # How many neutrals to you need?
# number.neutral.needed <- 600 - dim(complete_info_w_amplicon_adaptive)[1]
# number.neutral.needed # 419
# 
# complete_info_w_amplicon_neutral_limited <- head(complete_info_w_amplicon_neutral, n = number.neutral.needed)
# dim(complete_info_w_amplicon_neutral_limited)

# # Combine the adaptive and neutral complete infos
# complete_info_final <- rbind(complete_info_w_amplicon_adaptive, complete_info_w_amplicon_neutral_limited)
# dim(complete_info_final)

### END SPECIFIC TO TPAC ###

complete_info_final <- complete_info_w_amplicon # just to match the steps above

# Attach the 'recalc_pos_snp_in_window'
complete_info_final <- merge(x = complete_info_final, y = necess.posit.info, by = "radtag")
head(complete_info_final)
dim(complete_info_final)

#### 8. Format the submission dataframe to only keep necessary columns ####
#NEEDED: marker ID, ref genome scaff, begin pos, end pos, ref allele, variant allele, strand, marker type, priority, seq
colnames(complete_info_final) # these are what we have in total

# Format into the submission format
final <- complete_info_final[ , c("accn.name", "scaff", "recalc_pos_snp_in_window", "recalc_pos_snp_in_window", "ref", "var", "strand", "mtype", "priority", "amplicon_complete")]
head(final)

# Shorten the amplicon name
### AGAIN, SPECIFIC TO TPAC ###
# final$accn.name <- gsub(x = 
#        gsub(x = final$accn.name, pattern = ".*Tpa", replacement = "Tpa")
#      , pattern = "\\__.*" , replacement = "")
# head(final)
### END SPECIFIC TO TPAC

# Shorten the amplicon name
head(final$accn.name)
dim(final)
split.accn.name <- str_split_fixed(string = final$accn.name, pattern = "__", n = 4)
new.accn.name <- as.character(split.accn.name[, 2])
final$accn.name <- new.accn.name

# Write out result, this will be used for the submission
amplicon.panel.version <- 0.1
filename <- paste("06_output/", species, "_amplicon_panel_v", amplicon.panel.version, ".csv", sep = "")
write.csv(x = final, file = filename, quote = F, row.names = F)


# now go to terminal and check a few to be convinced it worked..



### OLD CODE ####

## This was an attempt to re-run identifying the correct location using the other side marker, but it didn't help
# 
# allele_isolated$radtag <- as.character(allele_isolated$radtag) 
# str(allele_isolated)
# 
# for(i in 1:length(failed.matches2)){
#   j <- failed.matches2[i]
#   print(j)
#   
#   # If it was originally said to be an 'after', choose 'before'
#   if(allele_isolated[which(allele_isolated[,"radtag"]==j), "side.to.use"]=="after"){ 
#     print(c(j, " was originally an 'after' match"))
#     allele_isolated[which(allele_isolated[,"radtag"]==j), "side.to.use"] <- "before"
#     print(c(j, "is now a 'before' match"))
#   } else {
#     allele_isolated[which(allele_isolated[,"radtag"]==j), "side.to.use"] <- "after"
#     print(c(j, "was originally a 'before' match"))
#     print(c(j, "is now an 'after' match"))
#   }
# }
# 


#View radtag
#rad.tags[which(rad.tags$radtag=="Tpa_1251"),]

# #Inspection
# write.csv(x = complete_info_w_amplicon, file = "06_output/inspecting_results.csv", quote = F, row.names = F)
