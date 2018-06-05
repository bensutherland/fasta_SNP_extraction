# Identify the location just before or just after the polymorphism in the amplicon
# Assumes that your rad locus file and amplicon file are in the same orientation

# rm(list=ls())
setwd("~/Documents/06_tpac/fasta_SNP_extraction_400bp/")

# Load libraries
require(tidyr)
require(stringr)

# Import data
load(file = "06_output/submission_data_part_1.RData")

# Relevant data to use
head(rad.tags)


#### 1. Isolate the first polymorphism in the RAD tag ####
# Take out the section before the first square brackets
first.piece.extracted <- separate(data = rad.tags, col = "radlocus", into = c("before.snp", "after.snp"), sep = "\\[" )
# expect warnings as it will read up to the second instance then stop
head(first.piece.extracted)

# Take out the section after the first square brackets
allele_isolated <- separate(data = first.piece.extracted, col = "after.snp", into = c("allele","after.snp"), sep ="\\]")
head(allele_isolated)


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

# Combine this specification to the isolated alleles
allele_isolated <- cbind(allele_isolated, side.to.use)
head(allele_isolated)
table(allele_isolated$side.to.use)


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
table(match.success) # two matching fails
complete_info_w_amplicon[which(match.success=="fail"),]

failed.matches <- complete_info_w_amplicon[which(match.success=="fail"), "radtag"] # these were the failing markers
failed.amplicons <- as.character(droplevels(failed.matches))
failed.amplicons

# Failed amplicons may manifest as NAs during the splitting, but more reliably, as a substantial difference from the 
# expected location of the snp in the window as taken from the accession name.

##### 5. Final check of the expected snp position and the actual snp position #####
colnames(complete_info_w_amplicon)

# Recalculate what position in the amplicon window has the polymorphism
final_check <- separate(data = complete_info_w_amplicon, col = "amplicon_complete", into = c("recalc_seq_before_snp", "dropseq"), sep = "\\[" , remove = F)
head(final_check)
final_check$nchar_before_snp <- nchar(final_check$recalc_seq_before_snp)
final_check$recalc_pos_snp_in_window <- final_check$nchar_before_snp + 1 
head(final_check, n = 5)

# # view
# final_check[, c("snp.pos.in.window","recalc_pos_snp_in_window")]

str(final_check)
final_check$snp.pos.in.window  <- as.numeric(as.character(final_check$snp.pos.in.window))

# Calculate how far off the expected and observed snp position are
number.digits.off <- final_check$recalc_pos_snp_in_window - final_check$snp.pos.in.window
number.digits.off
failed.matches2 <- final_check[which(abs(number.digits.off) > 5 ), "radtag" ]
failed.matches2 <- as.character(droplevels(failed.matches2))
failed.matches2

#### 6. Only keep the successful insertions ####
colnames(complete_info_w_amplicon)

for(i in 1:length(failed.matches2)){
  j <- failed.matches2[i]
  print(j)
  
  # replace the poorly-located allele amplicon with an unedited amplicon
  complete_info_w_amplicon[which(complete_info_w_amplicon[,"radtag"]==j), "amplicon_complete"] <- complete_info_w_amplicon[which(complete_info_w_amplicon[,"radtag"]==j), "seq"]
}




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
rad.tags[which(rad.tags$radtag=="Tpa_1251"),]

# #Inspection
# write.csv(x = complete_info_w_amplicon, file = "06_output/inspecting_results.csv", quote = F, row.names = F)

colnames(submission.df)
colnames(complete_info_w_amplicon)
# make a couple extra columns needed for submission
begin <- rep("NA", times = nrow(complete_info_w_amplicon))
end <- rep("NA", times = nrow(complete_info_w_amplicon))

complete_info_w_amplicon <- cbind(complete_info_w_amplicon, begin, end)
colnames(complete_info_w_amplicon)

final <- complete_info_w_amplicon[ , c("accn.name", "scaff", "begin", "end", "ref", "var", "strand", "mtype", "priority", "amplicon_complete")]
head(final)

write.csv(x = final, file = "06_output/complete_info_w_amplicon.csv", quote = F, row.names = F)
