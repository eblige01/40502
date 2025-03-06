# The goal of this code is to analyze paired primary breast tumors to LRR(locoregional recurrence) or met from same paitent
library(dplyr)
library(tidyr)
library("ggalluvial")
library(ggplot2)
library(reshape2)
library(umap)
library(writexl)
## PAM50 analysis 

# Loading data
paired_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\paired_metadata_dupadj.csv")


# Removing any rows with exclude = Yes
paired_metadata <- subset(paired_metadata,paired_metadata$Exclude == "No")
# Identify DeIdentifiedNumbers that occur more than once
duplicate_ids <- names(which(table(paired_metadata$DeIdentifiedNumber) > 1))

# Subset the data to include only rows with duplicate DeIdentifiedNumbers
paired_metadata_sub <- subset(paired_metadata, DeIdentifiedNumber %in% duplicate_ids)

# Rearrange the table so duplicates are together
paired_metadata_sub <- paired_metadata_sub[order(paired_metadata_sub$DeIdentifiedNumber), ]

### PAM50 Anaylsis
pam50_metadata <- read.delim("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE)

# Removing NA values from the call column
pam50_metadata <- pam50_metadata %>% filter(!is.na(Call))

#Formatting PAM50 data before the merge
pam50_metadata$X <- gsub("_C$", "", pam50_metadata$X)
colnames(pam50_metadata)[colnames(pam50_metadata) == "X"] <- "SlideID"
# Merging the two
pam50_pairs_merged <- merge(paired_metadata_sub , pam50_metadata, by = "SlideID")

# Subseting columns that are needed
pam50_pairs_merged  <- pam50_pairs_merged [, c("SlideID" ,"DeIdentifiedNumber","PrimaryMet","Call")]

# Identifing new duplicates after combining and cleaning the data
duplicate_ids <- names(which(table(pam50_pairs_merged $DeIdentifiedNumber) > 1))

# Subset the data to include only rows with duplicate DeIdentifiedNumbers
pam50_pairs_merged  <- subset(pam50_pairs_merged , DeIdentifiedNumber %in% duplicate_ids)

# Subsetting to only include Primary/Metastatic pairs
pam50_pairs_merged <- pam50_pairs_merged %>%
  group_by(DeIdentifiedNumber) %>%
  filter(all(c("Primary","Metastatic_LRR") %in% PrimaryMet)) %>%
  ungroup()
# Removing values to make even pairs
pam50_pairs_merged <- pam50_pairs_merged %>% filter(!SlideID %in% c("1017943", "1024138", "1035479"))

# Removing SlideID for reformatting
pam50_pairs_merged <- pam50_pairs_merged %>% dplyr::select(-SlideID)
rearranged_data <- pam50_pairs_merged %>%
  pivot_wider(
    names_from = PrimaryMet,
    values_from = Call,
    values_fn = list(Call = ~dplyr::first(.))
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
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("PAM50 Subtype Distribution: Primary vs Metastatic_LRR")



### RNA Deconvolution analysis
data <- file.choose()
# rna_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\rna_decon_matrix_40502.csv")
# rna_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\all_rna_samples_metadata.csv")
# 
# # Reformating rna_metadata to match rna_data
# rna_metadata$rna_decon_sampleid <- gsub("_", "", tolower(rna_metadata$rna_decon_sampleid))
# 
# # Replacing sample IDs with slide ids in rna_data
# colnames(rna_data)[-1] <- rna_metadata$Slide.ID..H.E...Biobank..[match(colnames(rna_data)[-1], rna_metadata$rna_decon_sampleid)]
# 
# #Transposing rna_data
# rna_data <- t(rna_data)
# 
# # Fixing formatting issues
# colnames(rna_data) <- rna_data[1,]
# rna_data <- as.data.frame(rna_data)
# rna_data <- rna_data[-1,]
# rownames(rna_data) <- gsub("X","",rownames(rna_data))
# rna_data$SlideID <- rownames(rna_data)
# #Subseting rna_data to only include pairs and removing pairs that dont have rna data
# 
# 
# # rna_data <- subset(rna_data, SlideID %in% paired_metadata_sub$SlideID)
# # paired_metadata_sub <- subset(paired_metadata_sub, SlideID %in% rna_data$SlideID)
# rna_data$SlideID <- as.integer(rna_data$SlideID)
# rna_data <- na.omit(rna_data)
# paired_rna_merged <- rna_data %>%
#   right_join(select(paired_metadata_sub, SlideID, PrimaryMet), by = "SlideID")
# 
# 
# 
# # Pairs Vs Immune
# # Initialize an empty dataframe to store module_module_results
# rna_results <- data.frame(Cell = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())
# 
# # Get the gene column names (excluding sample_id and TIL_group)
# cell_cols <- setdiff(names(paired_rna_merged), c("SlideID", "PrimaryMet"))
# 
# # Loop through each gene column
# for (cell in cell_cols) {
#   # Ensure the column is numeric
#   paired_rna_merged[[cell]] <- as.numeric(paired_rna_merged[[cell]])
#   # Perform Wilcoxon test (Mann-Whitney U test)
#   test_result <- wilcox.test(paired_rna_merged[[cell]] ~ paired_rna_merged$PrimaryMet)
# 
#   # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
#   group_means <- paired_rna_merged %>%
#     group_by(PrimaryMet) %>%
#     summarise(mean_expression = mean(!!sym(cell), na.rm = TRUE)) %>%
#     arrange(desc(mean_expression))
# 
#   # Determine which group has higher expression
#   higher_group <- group_means$PrimaryMet[1]
# 
#   # Append module_results to dataframe
#   rna_results <- rbind(rna_results, data.frame(Cell = cell, p_value = test_result$p.value,higher_expression = higher_group))
# }
# 
# rna_results$p_adj <- p.adjust(rna_results$p_value, method = "bonferroni")
# 
# # Subset significant module_results after adjustment (e.g., FDR < 0.05)
# significant_rna_results <- rna_results %>% filter(p_adj < 0.05)
# sig_cells <- significant_rna_results$Cell
# #Saving results
# write_xlsx(rna_results, "pairs_decon_no_sig_results.xlsx")




#Module analysis
M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv", dec=",")
signature_data <- read.table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\cdt.txt", header = TRUE, sep = "\t", comment.char = "", quote = "")

# Reformatting signature data for analysis
# Removing unnecessary rows and coloumns
resignature_data <- signature_data %>% select(-c("GID","CLID","GWEIGHT"))
resignature_data <- resignature_data %>% slice(-c(1,2))

# Transposing the dataframe
resignature_data <- t(resignature_data)

# Making first row colnames and making rownames a column
colnames(resignature_data) <- resignature_data[1,]
resignature_data <- resignature_data[-1,]
resignature_data <- as.data.frame(resignature_data)
resignature_data$INVESTIGATOR_SAMPLENAME <- rownames(resignature_data)

# Moving INVESTIGATOR_SAMPLENAME for to the from for convenience
resignature_data <- resignature_data[, c("INVESTIGATOR_SAMPLENAME", setdiff(names(resignature_data), "INVESTIGATOR_SAMPLENAME"))]

# Removing leading X for formatting
resignature_data <- resignature_data %>% mutate(INVESTIGATOR_SAMPLENAME = sub("^X", "", INVESTIGATOR_SAMPLENAME))

# Removing SlideIDs not found in the metadata
resignature_data <- resignature_data %>%
  filter(INVESTIGATOR_SAMPLENAME %in% M40502_joined_metadata$INVESTIGATOR_SAMPLENAME)

# Merging with metadata for the SlideID
resignature_data <- merge(resignature_data, M40502_joined_metadata[, c("INVESTIGATOR_SAMPLENAME", "Slide.ID..H.E...Biobank..")], by = "INVESTIGATOR_SAMPLENAME")
colnames(resignature_data)[colnames(resignature_data) == "Slide.ID..H.E...Biobank.."] <- "SlideID"
# Merging paired data and module data
paired_data <- merge(resignature_data, paired_metadata_sub[, c("SlideID", "PrimaryMet")], by = "SlideID")

module_results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())

# Get the gene column names (excluding sample_id and TIL_group)
gene_cols <- setdiff(names(paired_data), c("INVESTIGATOR_SAMPLENAME", "PrimaryMet"))

# Loop through each gene column
for (gene in gene_cols) {
  # Ensure the column is numeric
  paired_data[[gene]] <- as.numeric(paired_data[[gene]])
  # Perform Wilcoxon test (Mann-Whitney U test)
  test_result <- wilcox.test(paired_data[[gene]] ~ paired_data$PrimaryMet)

  # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
  group_means <- paired_data %>%
    group_by(PrimaryMet) %>%
    summarise(mean_expression = mean(!!sym(gene), na.rm = TRUE)) %>%
    arrange(desc(mean_expression))

  # Determine which group has higher expression
  higher_group <- group_means$PrimaryMet[1]

  # Append module_results to dataframe
  module_results <- rbind(module_results, data.frame(Gene = gene, p_value = test_result$p.value,higher_expression = higher_group))
}

module_results$p_adj <- p.adjust(module_results$p_value, method = "BH")

# Subset significant module_results after adjustment (e.g., FDR < 0.05)
significant_module_results <- module_results %>% filter(p_adj < 0.05)

#Reordering significant module_results
significant_module_results <- significant_module_results %>%
  mutate(PMID = ifelse(grepl("_PMID\\.\\d+$", Gene),
                       sub(".*_PMID\\.", "", Gene),
                       NA_character_)) %>%
  arrange(desc(!is.na(PMID)), PMID)

write_xlsx(significant_module_results, "Paired_gene_signatures.xlsx")