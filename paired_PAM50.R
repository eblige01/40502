# The goal of this code is to analyze paired primary breast tumors to LRR(locoregional recurrence) or met from same paitent
library(dplyr)
library(tidyr)
library("ggalluvial")
library(ggplot2)
library(reshape2)
## PAM50 anaylsis 

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
merged_metadata <- merged_metadata[, c("SlideID" ,"DeIdentifiedNumber","PrimaryMet","Call")]

# Removing new uniques after combining and cleaning the data
duplicate_ids <- names(which(table(merged_metadata$DeIdentifiedNumber) > 1))

# Subset the data to include only rows with duplicate DeIdentifiedNumbers
merged_metadata <- subset(merged_metadata, DeIdentifiedNumber %in% duplicate_ids)

# Subsetting to only include Primary/Metastatic pairs
merged_metadata <- merged_metadata %>%
  group_by(DeIdentifiedNumber) %>%
  filter(all(c("Primary","Metastatic_LRR") %in% PrimaryMet)) %>%
  ungroup()
# Removing values to make even pairs
merged_metadata <- merged_metadata %>% filter(!SlideID %in% c("1017943", "1024138", "1035479"))

# Removing SlideID for reformatting
merged_metadata <- merged_metadata %>% select(-SlideID)
rearranged_data <- merged_metadata %>%
  pivot_wider(
    names_from = PrimaryMet,
    values_from = Call,
    values_fn = list(Call = ~first(.))
  )


# Aggregate the data to get the total counts for each Primary
totals <- rearranged_data %>%
  group_by(Primary) %>%
  summarise(total = n())  # 'n()' counts the number of rows for each Primary


# Create the plot
ggplot(rearranged_data,
       aes(axis1 = Primary , axis2 = Metastatic_LRR)) +
geom_alluvium(aes(fill = Primary), width = 1/6, color = "black",aggregate.y = TRUE) +  # Change flow color based on Metastatic_LRR
geom_stratum(aes(fill =  after_stat(stratum)), width = 1/6) +
geom_label(stat = "stratum", aes(label =  after_stat(stratum)), size = 3) +  # Smaller labels
geom_text(stat = "flow", aes(label = after_stat(count)), nudge_x = 0.2, size = 3) +  # Flow counts
scale_x_discrete(limits = c("Primary", "Metastatic_LRR"), expand = c(.05, .05)) +
scale_fill_brewer(type = "qual", palette = "Set1") +
theme(legend.position = "none") +
ggtitle("PAM50 Subtype Distribution: Primary vs Metastatic_LRR")


## RNA Deconvolution analysis

rna_data <- read.csv(file.choose())
rna_metadata <- read.csv(file.choose())

# Reformating rna_metadata to match rna_data
rna_metadata$rna_decon_sampleid <- gsub("_", "", tolower(rna_metadata$rna_decon_sampleid))

# Replacing sample IDs with slide ids in rna_data
colnames(rna_data)[-1] <- rna_metadata$Slide.ID..H.E...Biobank..[match(colnames(rna_data)[-1], rna_metadata$rna_decon_sampleid)]

#

