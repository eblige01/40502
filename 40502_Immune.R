# 40502 Immune (TILSs) analysis 
library(ggplot2)
library("gplots")
library(dplyr)  
library(writexl)
library(FSA) 

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
uniqueData$sTILs_cat <- ifelse(uniqueData$sTILs >= 5,"High","Low")
### RNA SEQ analysis

# Subsetting RNA SEQ to only containing high TILS samples (Only for TILs analyses)

rna_seq_df <- rna_seq_df[,colnames(rna_seq_df) %in% sampleID]
ncol(rna_seq_df)

rownames(uniqueData) <- uniqueData$rna_decon_sampleid

# Reformatting and merging
trans_rna_seq_df <- as.data.frame(t(rna_seq_df))

trans_rna_seq_df$rna_decon_sampleid <- row.names(trans_rna_seq_df)


trans_rna_seq_df <- merge(trans_rna_seq_df, uniqueData[, c("rna_decon_sampleid","sTILs_cat","sTILs")], by = "rna_decon_sampleid")
row.names(trans_rna_seq_df) <- trans_rna_seq_df$rna_decon_sampleid
trans_rna_seq_df <- trans_rna_seq_df %>% select (-"rna_decon_sampleid")


# # sTILs Vs Immune
# # Initialize an empty dataframe to store module_module_results
# rna_results <- data.frame(Cell = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())
# 
# # Get the gene column names (excluding sample_id and TIL_group)
# cell_cols <- setdiff(names(trans_rna_seq_df), c("rna_decon_sampleid", "sTILs_cat","sTILs"))
# 
# # Loop through each gene column
# for (cell in cell_cols) {
#   # Ensure the column is numeric
#   trans_rna_seq_df[[cell]] <- as.numeric(trans_rna_seq_df[[cell]])
#   # Perform Wilcoxon test (Mann-Whitney U test)
#   test_result <- wilcox.test(trans_rna_seq_df[[cell]] ~ trans_rna_seq_df$sTILs_cat)
# 
#   # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
#   group_means <- trans_rna_seq_df %>%
#     group_by(sTILs_cat) %>%
#     summarise(mean_expression = mean(!!sym(cell), na.rm = TRUE)) %>%
#     arrange(desc(mean_expression))
# 
#   # Determine which group has higher expression
#   higher_group <- group_means$sTILs_cat[1]
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
# write_xlsx(significant_rna_results, "sTILs_decon.xlsx")
# 
# 
# ### Module Analysis 
# 
# signature_data <- read.table("/Users/eblige99/Desktop/Research/StoverLab_rotation/data/cdt.txt", header = TRUE, sep = "\t", comment.char = "", quote = "")
# 
# # Reformatting signature data for analysis 
# # Removing unnecessary rows and coloumns
# resignature_data <- signature_data %>% select(-c("GID","CLID","GWEIGHT"))
# resignature_data <- resignature_data %>% slice(-c(1,2))
# 
# # Transposing the dataframe
# resignature_data <- t(resignature_data)
# 
# # Making first row colnames and making rownames a column 
# colnames(resignature_data) <- resignature_data[1,]
# resignature_data <- resignature_data[-1,]
# resignature_data <- as.data.frame(resignature_data)
# resignature_data$INVESTIGATOR_SAMPLENAME <- rownames(resignature_data)
# 
# # Moving INVESTIGATOR_SAMPLENAME for to the from for convenience 
# resignature_data <- resignature_data[, c("INVESTIGATOR_SAMPLENAME", setdiff(names(resignature_data), "INVESTIGATOR_SAMPLENAME"))]
# 
# # Reformatting and merging 
# resignature_data <- resignature_data %>% mutate(INVESTIGATOR_SAMPLENAME = sub("^X", "", INVESTIGATOR_SAMPLENAME))
# 
# resignature_data <- merge(resignature_data,uniqueData[, c("INVESTIGATOR_SAMPLENAME", "sTILs_cat","sTILs")],by = "INVESTIGATOR_SAMPLENAME",all.y=TRUE)
# 
# # sTILs vs Modules
# # Initialize an empty dataframe to store module_results
# module_results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())
# 
# # Get the gene column names (excluding sample_id and TIL_group)
# gene_cols <- setdiff(names(resignature_data), c("INVESTIGATOR_SAMPLENAME", "sTILs_cat","sTILs"))
# 
# # Loop through each gene column
# for (gene in gene_cols) {
#   # Ensure the column is numeric
#   resignature_data[[gene]] <- as.numeric(resignature_data[[gene]])
#   # Perform Wilcoxon test (Mann-Whitney U test)
#   test_result <- wilcox.test(resignature_data[[gene]] ~ resignature_data$sTILs_cat)
# 
#   # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
#   group_means <- resignature_data %>%
#     group_by(sTILs_cat) %>%
#     summarise(mean_expression = mean(!!sym(gene), na.rm = TRUE)) %>%
#     arrange(desc(mean_expression))
# 
#   # Determine which group has higher expression
#   higher_group <- group_means$sTILs_cat[1]
# 
#   # Append module_results to dataframe
#   module_results <- rbind(module_results, data.frame(Gene = gene, p_value = test_result$p.value,higher_expression = higher_group))
# }
# 
# module_results$p_adj <- p.adjust(module_results$p_value, method = "BH")
# 
# # Subset significant module_results after adjustment (e.g., FDR < 0.05)
# significant_module_results <- module_results %>% filter(p_adj < 0.05)
# 
# #Reordering significant module_results
# significant_module_results <- significant_module_results %>%
#   mutate(PMID = ifelse(grepl("_PMID\\.\\d+$", Gene),
#                        sub(".*_PMID\\.", "", Gene),
#                        NA_character_)) %>%
#   arrange(desc(!is.na(PMID)), PMID)
# 
# write_xlsx(significant_module_results, "sTILs_module.xlsx")

# Solo testing individual signatures
# gene_to_plot <- "UNC_Scorr_LumA_Correlation_JCO.2009_PMID.19204204"
# 
# # Perform Kruskal-Wallis test
# test_result <- kruskal.test(resignature_data[[gene_to_plot]] ~ resignature_data$Call)
# 
# # Summarize data: Mean and Standard Error (SE)
# plot_data <- resignature_data %>%
#   group_by(Call) %>%
#   summarise(mean_expression = mean(!!sym(gene_to_plot), na.rm = TRUE),
#             se = sd(!!sym(gene_to_plot), na.rm = TRUE) / sqrt(n()))
# 
# # Bar plot with error bars
# ggplot(plot_data, aes(x = Call, y = mean_expression, fill = Call)) +
#   geom_bar(stat = "identity", color = "black", width = 0.6) +  # Bar plot
#   geom_errorbar(aes(ymin = mean_expression - se, ymax = mean_expression + se), 
#                 width = 0.2, color = "black") +  # Error bars
#   annotate("text", x = 1.5, y = max(plot_data$mean_expression), 
#            label = paste("p =", format(test_result$p.value, digits = 2)), size = 5) +
#   labs(title = paste("Mean Expression of", gene_to_plot),
#        x = "PAM50 Subtype",
#        y = "Mean Expression ± SE") +
#   theme_minimal()
