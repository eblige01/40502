# The goal of this code is to analyze paired primary breast tumors to LRR(locoregional recurrence) or met from same paitent
library(dplyr)
library(tidyr)
# Loading data
paired_metadata <- read.csv(file.choose())
pam50_metadata <- read.delim(file.choose(), header = TRUE)

# Removing NA values from the call column 
pam50_metadata <- pam50_metadata %>% filter(!is.na(Call))


# Removing any rows with exclude = Yes
paired_metadata <- subset(paired_metadata,paired_metadata$Exclude == "No")
# Identify DeIdentifiedNumbers that occur more than once
duplicate_ids <- names(which(table(paired_metadata$DeIdentifiedNumber) > 1))

# Subset the data to include only rows with duplicate DeIdentifiedNumbers
paired_metadata_sub <- subset(paired_metadata, DeIdentifiedNumber %in% duplicate_ids)

# Rearrange the table so duplicates are together
paired_metadata_sub <- paired_metadata_sub[order(paired_metadata_sub$DeIdentifiedNumber), ]

#Formatting PAM50 data before the merge 
pam50_metadata$X <- gsub("_C$", "", pam50_metadata$X)
colnames(pam50_metadata)[colnames(pam50_metadata) == "X"] <- "SlideID"
# Merging the two
merged_metadata <- merge(paired_metadata_sub , pam50_metadata, by = "SlideID")

# Subseting columns that are needed
merged_metadata <- merged_metadata[, c("SlideID", "DeIdentifiedNumber","PrimaryMet","Call")]

# Removing new uniques after combining and cleaning the data 
duplicate_ids <- names(which(table(merged_metadata$DeIdentifiedNumber) > 1))

# Subset the data to include only rows with duplicate DeIdentifiedNumbers
merged_metadata <- subset(merged_metadata, DeIdentifiedNumber %in% duplicate_ids)

# Subsetting to only include Primary/Metastatic pairs
merged_metadata <- merged_metadata %>%
  group_by(DeIdentifiedNumber) %>%
  filter(all(c("Primary","Metastatic_LRR") %in% PrimaryMet)) %>%
  ungroup()
# Removing 1017942 to make even pairs since it is the lowest confidence 
merged_metadata <- merged_metadata %>%filter(my_column != "remove_this")

pairwise_comparisons <- merged_metadata %>%
  filter(PrimaryMet %in% c("Primary", "Metastatic_LRR")) %>%  # Keep only Primary and Metastatic_LRR samples
  group_by(DeIdentifiedNumber) %>%  # Group by DeIdentifiedNumber
  spread(PrimaryMet, Call) %>%  # Reshape so Primary and Metastatic_LRR are separate columns
  ungroup() %>%  # Ungroup after spreading
  mutate(transition = paste(Primary, Metastatic_LRR, sep = " -> ")) %>%  # Create a column for transitions
  count(transition)  # Count unique transitions
