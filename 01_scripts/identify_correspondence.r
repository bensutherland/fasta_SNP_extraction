## Identify corresponding markers based on genomic coordinates
# rm(list=ls())

# Read in all info
data.suppl <- read.csv(file = "05_amplicons/dropped_duplicates_suppl_info.csv", stringsAsFactors = F)
# Read in dropped duplicates (default selection)
dropped_dups <- read.csv(file = "05_amplicons/dropped_duplicates.csv", stringsAsFactors = F)

head(dropped_dups) # contains the dropped duplicates
head(data.suppl) # does not contain the duplicates


# Find the corresponding retained marker for the dropped markers
corresp_data <- merge(x = dropped_dups, y = data.suppl, by = "match.id", all.x = TRUE)
head(corresp_data)
dim(corresp_data)
corresp_data <- corresp_data[, c("match.id", "mname.x", "qname.x", "mname.y", "qname.y")] # Reduce df

head(corresp_data)

# Rename columns
colnames(corresp_data) <- c("match.id", "mname.dropped", "qname.dropped", "mname.retained", "qname.retained")

head(corresp_data)

# These should match
dim(corresp_data)
dim(dropped_dups)

# Import priority lists
amp_mnames_adfg <- read.csv(file = "02_input_data/amp_mnames_adfg.csv", header = F, stringsAsFactors = FALSE)
head(amp_mnames_adfg)
dim(amp_mnames_adfg)
colnames(amp_mnames_adfg) <- c("mname", "source")

amp_mnames_critfc <- read.csv(file = "02_input_data/amp_mnames_critfc.csv", header = F, stringsAsFactors = FALSE)
head(amp_mnames_critfc)
dim(amp_mnames_critfc)
colnames(amp_mnames_critfc) <- c("mname", "source")

amp_mnames_other <- read.csv(file = "02_input_data/priority_source_snps_ebr_2019-09-25.csv", header = F, stringsAsFactors = FALSE)
head(amp_mnames_other)
dim(amp_mnames_other)
colnames(amp_mnames_other) <- c("mname", "source")

# Combine the priority panels
amp_mnames_priority <- rbind(amp_mnames_adfg, amp_mnames_critfc, amp_mnames_other)
dim(amp_mnames_priority)

# Identify the corresponding marker to the dropped marker, in order to keep the priority marker
# Identify the markers that were dropped but in the priority list

# fill vector with dummy character
corresp_data$duplicated_name <- rep("no", times = nrow(corresp_data))

# Loop to correct the marker that will be dropped based on the priority list
for(i in 1:nrow(corresp_data)){
  
  # If the dropped marker is in the priority list, select the other marker to drop as the duplicate
  if(corresp_data$mname.dropped[i] %in% amp_mnames_priority$mname & !corresp_data$mname.retained[i] %in% amp_mnames_priority$mname){
    
    corresp_data$marker_to_drop[i] <- corresp_data$mname.retained[i] 
    corresp_data$change_marker[i] <- "yes"
    
  # If the dropped marker is not in the priority list, keep the dropped marker as the one to drop
  } else if(!corresp_data$mname.dropped[i] %in% amp_mnames_priority$mname){
    
    corresp_data$marker_to_drop[i] <- corresp_data$mname.dropped[i]
    corresp_data$change_marker[i] <- "no"
    
  } else if(corresp_data$mname.dropped[i] %in% amp_mnames_priority$mname & corresp_data$mname.retained[i] %in% amp_mnames_priority$mname){
    
    corresp_data$marker_to_drop[i] <- corresp_data$mname.dropped[i]
    corresp_data$change_marker[i] <- "no"
    
  }
  
  # If the two names (dropped and retain) are equal, this name should not be in the drop category, otherwise both will be removed
  if(corresp_data$mname.dropped[i]==corresp_data$mname.retained[i]){
    
    corresp_data$duplicated_name[i] <- "yes"
    
  }
  
}


# If the record has an identical name problem, remove it from the drop list and let the duplicate screen remove the duplicate
dim(corresp_data)
corresp_data <- corresp_data[corresp_data$duplicated_name=="no", ]
dim(corresp_data)

head(corresp_data, n = 10)
table(corresp_data$change_marker)

# Write out results
write.csv(corresp_data, file = "05_amplicons/identifying_duplicate_to_drop.csv", quote = F, row.names = F)
