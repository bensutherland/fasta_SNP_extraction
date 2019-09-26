## Identify corresponding markers based on genomic coordinates
# rm(list=ls())

#### 01. Read in, join, and format ####
# Read in all info (needed to pair markers from the deleted list to the retained list)
data.suppl <- read.csv(file = "05_amplicons/dropped_duplicates_suppl_info.csv", stringsAsFactors = F)
# Read in dropped duplicates (default selection)
dropped_dups <- read.csv(file = "05_amplicons/dropped_duplicates.csv", stringsAsFactors = F)

head(dropped_dups) # simple mnames and coordinates of markers to be dropped due to duplicate
dim(dropped_dups)
head(data.suppl) # all_info, when duplicates have been removed
dim(data.suppl)

# To find the corresponding retained marker for the dropped markers, need to combine the two files
corresp_data <- merge(x = dropped_dups, y = data.suppl, by = "match.id", all.x = TRUE)
head(corresp_data)
dim(corresp_data)
dim(dropped_dups) # should match
corresp_data <- corresp_data[, c("match.id", "mname.x", "qname.x", "mname.y", "qname.y")] # Reduce df
colnames(corresp_data) <- c("match.id", "mname.dropped", "qname.dropped", "mname.retained", "qname.retained")

head(corresp_data)

#### 02. Import priority lists ####
# Import priority lists, priority 1
priority1.FNs <- list.files(path = "02_input_data", pattern = "_priority1_mnames.csv", full.names = T)
priority1.df <- NULL; temp <- NULL
for(i in 1:length(priority1.FNs)){
  
  temp <- read.csv(file = priority1.FNs[i], header = F, stringsAsFactors = FALSE)
  colnames(temp) <- c("mname", "source")
  priority1.df <- rbind(priority1.df, temp)
  
}

amp_mnames_priority1 <- priority1.df
head(amp_mnames_priority1)
dim(amp_mnames_priority1)
table(amp_mnames_priority1$source) # these are the requested numbers of markers from each priority

# Import priority lists, priority 2
priority2.FNs <- list.files(path = "02_input_data", pattern = "_priority2_mnames.csv", full.names = T)
priority2.df <- NULL; temp <- NULL
for(i in 1:length(priority2.FNs)){
  
  temp <- read.csv(file = priority2.FNs[i], header = F, stringsAsFactors = FALSE)
  colnames(temp) <- c("mname", "source")
  priority2.df <- rbind(priority2.df, temp)
  
}

amp_mnames_priority2 <- priority2.df
head(amp_mnames_priority2)
dim(amp_mnames_priority2)
table(amp_mnames_priority2$source) # these are the requested numbers of markers from each priority


#### Replace expected dropped markers that are from the priority lists with non-priority markers ####
head(corresp_data)

# # fill vector with dummy character
# corresp_data$duplicated_name <- rep("no", times = nrow(corresp_data))

# Loop to correct the marker that will be dropped based on the priority list
for(i in 1:nrow(corresp_data)){
  
  # Identify whether the dropped marker is in the priority1 list
  if(corresp_data$mname.dropped[i] %in% amp_mnames_priority1$mname){
    
    corresp_data$dropped.priority1[i] <- "yes"
    
  }else{
    
    corresp_data$dropped.priority1[i] <- "no"
    
  }
  
  # Identify whether the dropped marker is in the priority2 list
  if(corresp_data$mname.dropped[i] %in% amp_mnames_priority2$mname){
    
    corresp_data$dropped.priority2[i] <- "yes"
    
  }else{
    
    corresp_data$dropped.priority2[i] <- "no"
    
  }
  
  # Identify whether the retained marker is in the priority1 list
  if(corresp_data$mname.retained[i] %in% amp_mnames_priority1$mname){
    
    corresp_data$retained.priority1[i] <- "yes"
    
  }else{
    
    corresp_data$retained.priority1[i] <- "no"
    
  }
  
  # Identify whether the retained marker is in the priority2 list
  if(corresp_data$mname.retained[i] %in% amp_mnames_priority2$mname){
    
    corresp_data$retained.priority2[i] <- "yes"
    
  }else{
    
    corresp_data$retained.priority2[i] <- "no"
    
  }
  
  
  
  
}  

head(corresp_data, n = 10)

table(corresp_data$dropped.priority1)
table(corresp_data$dropped.priority2)
table(corresp_data$retained.priority1)
table(corresp_data$retained.priority2)


# Set nulls (necessary for missing conditional rows)
corresp_data$marker_to_drop <- rep("NA", times = nrow(corresp_data))
corresp_data$switched <- rep("NA", times = nrow(corresp_data))

# Switch markers away from dropping preferred markers
for(i in 1:nrow(corresp_data)){
  
  # If the marker dropped is priority one, and the retained is not, switch
  if(corresp_data$dropped.priority1[i]=="yes" & corresp_data$retained.priority1[i]=="no"){
    
    # Switch the retained marker into the drop column
    corresp_data$marker_to_drop[i] <- corresp_data$mname.retained[i]
    
    # Reporting
    corresp_data$switched[i] <- "true-keep_pri1"
  
  # If the marker dropped is priority one, and the retained marker is also priority one, do nothing  
  } else if(corresp_data$dropped.priority1[i]=="yes" & corresp_data$retained.priority1[i]=="yes"){
    
    # do nothing
    
    # Reporting
    corresp_data$switched[i] <- "cant_both_are_pri_1"
  
  # If the marker dropped is not priority one, and the marker dropped is priority two, then switch
  } else if(corresp_data$dropped.priority1[i]=="no" & corresp_data$dropped.priority2[i]=="yes" & corresp_data$retained.priority1[i]=="no"){
    
    # switch the retained marker to the drop column
    corresp_data$marker_to_drop[i] <- corresp_data$mname.retained[i]
    
    # Reporting
    corresp_data$switched[i] <- "true-keep_pri2"
    
  } else if(corresp_data$dropped.priority1[i]=="no" & corresp_data$dropped.priority2[i]=="no"){
    
    # Keep with the planned drop marker
    corresp_data$marker_to_drop[i] <- corresp_data$mname.dropped[i]
    
    # Reporting
    corresp_data$switched[i] <- "false"
    
  # 
  } else if(corresp_data$dropped.priority1[i]=="no" & corresp_data$retained.priority1[i]=="yes"){
    
    # Keep with the planned drop marker
    corresp_data$marker_to_drop[i] <- corresp_data$mname.dropped[i]
    
    # Reporting
    corresp_data$switched[i] <- "false"
    
  }
}


head(corresp_data[,-1], n = 10)
table(corresp_data$switched)


# If can't drop either as they are both priority 1 (could be because same name, or both markers are in priority one, then 
#  remove it from the drop list and let the duplicate screen remove the duplicate at random 
#  This should not happen if one of the markers is non-priority 1
dim(corresp_data)
corresp_data <- corresp_data[corresp_data$switched!="cant_both_are_pri_1", ]
dim(corresp_data)

head(corresp_data, n = 10)
table(corresp_data$switched)

# Write out results
write.csv(corresp_data, file = "05_amplicons/identifying_duplicate_to_drop.csv", quote = F, row.names = F)
