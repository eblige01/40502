library(ggplot2)
library("gplots")
library(dplyr)  # Install if not already installed
library(writexl)


# Loading in the data

M40502_joined_metadata <- read.csv("~/Desktop/Research/StoverLab_rotation/data/40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("~/Desktop/Research/StoverLab_rotation/data/D40502_data.csv", dec=",")
pam50 <-  read.table("~/Desktop/Research/StoverLab_rotation/data/PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE, sep = "\t") 
rna_seq_df <- read.csv("~/Desktop/Research/StoverLab_rotation/data/rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]


# Reformating and merging 
M40502_joined_metadata <- M40502_joined_metadata[!is.na(M40502_joined_metadata$Label.on.Curls),]
names(pam50)[1] <- "INVESTIGATOR_SAMPLENAME"
mergeddata <- merge(D40502_data, M40502_joined_metadata[, c("patid", "sTILs","rna_decon_sampleid","INVESTIGATOR_SAMPLENAME")], by = "patid")
mergeddata$rna_decon_sampleid <- gsub("Sample_", "sample",mergeddata$rna_decon_sampleid)
mergeddata <- merge(mergeddata,pam50,by="INVESTIGATOR_SAMPLENAME")
mergeddata <- mergeddata[!is.na(mergeddata$sTILs),]

# Removing duplicate paitent IDs by selecting the highest TILS
mergeddata <- mergeddata[!is.na(mergeddata$sTILs), ]
sortedData <- mergeddata[order(mergeddata$patid),]
uniqueData <- sortedData %>%
  group_by(patid) %>%
  slice_max(order_by = sTILs, n = 1) %>%
  ungroup()
head(uniqueData)
uniqueData<- unique(uniqueData)
uniqueData <- uniqueData %>%
  group_by(rna_decon_sampleid) %>%
  filter(n() == 1) %>%
  ungroup()

sampleID <- uniqueData$rna_decon_sampleid
uniqueData$rna_decon_sampleid <- gsub("Sample_", "sample",sampleID)
uniqueData$sTILs <- ifelse(uniqueData$sTILs >= 5,"High","Low")
# RNA SEQ analysis

# Subsetting RNA SEQ to only containing high TILS samples

rna_seq_df <- rna_seq_df[,colnames(rna_seq_df) %in% sampleID]
ncol(rna_seq_df)

rownames(uniqueData) <- uniqueData$rna_decon_sampleid

# Reformatting and merging
trans_rna_seq_df <- as.data.frame(t(rna_seq_df))

trans_rna_seq_df$rna_decon_sampleid <- row.names(trans_rna_seq_df)
  
trans_rna_seq_df <- merge(trans_rna_seq_df, uniqueData[, c("rna_decon_sampleid", "sTILs")], by = "rna_decon_sampleid")

# Looping and testing each cell type for each estimate

# Initialize an empty dataframe to store module_module_results
rna_results <- data.frame(Cell = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())

# Get the gene column names (excluding sample_id and TIL_group)
cell_cols <- setdiff(names(trans_rna_seq_df), c("rna_decon_sampleid", "sTILs"))

# Loop through each gene column
for (cell in cell_cols) {
  # Ensure the column is numeric
  trans_rna_seq_df[[cell]] <- as.numeric(trans_rna_seq_df[[cell]])
  # Perform Wilcoxon test (Mann-Whitney U test)
  test_result <- wilcox.test(trans_rna_seq_df[[cell]] ~ trans_rna_seq_df$sTILs)
  
  # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
  group_means <- trans_rna_seq_df %>%
    group_by(sTILs) %>%
    summarise(mean_expression = mean(!!sym(cell), na.rm = TRUE)) %>%
    arrange(desc(mean_expression))
  
  # Determine which group has higher expression
  higher_group <- group_means$sTILs[1]
  
  # Append module_results to dataframe
  rna_results <- rbind(rna_results, data.frame(Cell = cell, p_value = test_result$p.value,higher_expression = higher_group))
}

rna_results$p_adj <- p.adjust(rna_results$p_value, method = "bonferroni")

# Subset significant module_results after adjustment (e.g., FDR < 0.05)
significant_rna_results <- rna_results %>% filter(p_adj < 0.05)

# # Spearman correlation 
# ## Convert columns to numeric 
# mergeCiber$sTILs <- as.numeric(as.character(mergeCiber$sTILs))
# corrList <- list()
# for (type in columnname2) {
#   # Convert the current column to numeric
#   mergeCiber[[type]] <- as.numeric(as.character(mergeCiber[[type]]))
#   
#   # Perform Spearman correlation
#   correlation <- cor(mergeCiber$sTILs, mergeCiber[[type]], method = "spearman")
#   
#   # Store the result in the list
#   corrList[[type]] <- correlation
# }
# 
# # Convert corrList to a dataframe
# corr_df <- data.frame(
#   Cell_Type = names(corrList),
#   Correlation = unlist(corrList)
# )
# 
# # Convert corrList to a dataframe
# corr_df <- data.frame(
#   Cell_Type = names(corrList),
#   Correlation = as.numeric(unlist(corrList))
# )
# 
# # Create a bar plot using ggplot2
# ggplot(corr_df, aes(x = Cell_Type, y = Correlation)) +
#   geom_bar(stat = "identity", fill = "skyblue") +  # Create bars
#   labs(y = "Spearman Correlation") +
#   theme_minimal() +
#   geom_segment(aes(x=0,xend=0,y=0,yend=0.4),color = "black",size = .5)+
#   geom_hline(yintercept = 0,color = "black",size=.5) +
#   theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 65, hjust = 1,size =7.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())  # Rotate x-axis labels for better readability

# Making box plot with colored dots
# mergeCiber$sTILs <- as.factor(mergeCiber$sTILs)
# mergeCiber$`T cell CD8+_CIBERSORT` <- as.numeric(mergeCiber$`T cell CD8+_CIBERSORT`)
# 
# ggplot(mergeCiber, aes(x = sTILs, y = `T cell CD8+_CIBERSORT`)) +
#   geom_boxplot(fill = "black", color = "black", outlier.shape = NA) +  # Black boxes
#   geom_jitter(aes(color = factor(Call)), width = 0.2, size = 2) +   # Colored dots for subtypes
#   labs(title = "TILs Levels by Subtype",
#        x = "TILs Group",
#        y = "T cell CD8+_CIBERSORT") +
#   scale_color_manual(values = c("red", "pink", "lightblue", "blue", "green")) +  # Assign colors to subtypes
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 12),  # Customize x-axis label size
#         axis.text.y = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10)) 
# 
# # making violin plot
# mergeCiber$`Macrophage M2_CIBERSORT` <- as.numeric(mergeCiber$`Macrophage M2_CIBERSORT`)
# ggplot(mergeCiber, aes(x = sTILs, y = `Macrophage M2_CIBERSORT`)) +
#   geom_violin(fill = "black", color = "black") +  # Create black violin shapes
#   geom_jitter(aes(color = factor(Call)), width = 0.2, size = 1) +   # Colored dots for subtypes
#   labs(color = "PAM50 subtype",
#        x = "TILs Group",
#        y = "Macrophage M2_CIBERSORT") +
#   scale_color_manual(values = c("red", "pink", "lightblue", "blue", "green")) +  # Assign colors to subtypes
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 12),  # Customize x-axis label size
#         axis.text.y = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10)) +
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

# Module Analysis 

signature_data <- read.table("/Users/eblige99/Desktop/Research/StoverLab_rotation/data/cdt.txt", header = TRUE, sep = "\t", comment.char = "", quote = "")

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

# Reformatting and merging 
resignature_data <- resignature_data %>% mutate(INVESTIGATOR_SAMPLENAME = sub("^X", "", INVESTIGATOR_SAMPLENAME))

resignature_data <- merge(resignature_data,uniqueData[, c("INVESTIGATOR_SAMPLENAME", "sTILs")],by = "INVESTIGATOR_SAMPLENAME",all.y=TRUE)

# Initialize an empty dataframe to store module_module_results
module_results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())

# Get the gene column names (excluding sample_id and TIL_group)
gene_cols <- setdiff(names(resignature_data), c("INVESTIGATOR_SAMPLENAME", "sTILs"))

# Loop through each gene column
for (gene in gene_cols) {
  # Ensure the column is numeric
  resignature_data[[gene]] <- as.numeric(resignature_data[[gene]])
  # Perform Wilcoxon test (Mann-Whitney U test)
  test_result <- wilcox.test(resignature_data[[gene]] ~ resignature_data$sTILs)
  
  # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
  group_means <- resignature_data %>%
    group_by(sTILs) %>%
    summarise(mean_expression = mean(!!sym(gene), na.rm = TRUE)) %>%
    arrange(desc(mean_expression))
  
  # Determine which group has higher expression
  higher_group <- group_means$sTILs[1]
  
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

# Saving significant module_results as a excel file
write_xlsx(significant_module_results, "sTILs_gene_signatures.xlsx")
#Solo testing individual signatures
gene_to_plot <- "Fibroblasts_MCP_PMID.31942075_PMID.31942077"
test_result <- wilcox.test(resignature_data[[gene_to_plot]] ~ resignature_data$sTILs)
# Box plot
ggplot(resignature_data, aes(x = sTILs, y = !!sym(gene_to_plot), fill = sTILs)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add jittered points
  annotate("text", x = 1.5, y = max(resignature_data[[gene_to_plot]]), label = paste("p = ", format(test_result$p.value, digits = 2)), size = 5) +
  labs(title = paste("Boxplot of", gene_to_plot),
       x = "TIL Group",
       y = "Expression Level") +
  theme_minimal()
