# 40502 Immune (TILSs) analysis 
library(ggplot2)
library("gplots")
library(dplyr)  
library(writexl)
library(FSA) 
library(DESeq2)
library(org.Hs.eg.db) 
library(enrichplot)
library("clusterProfiler")

# Loading in the data
# Replace with the latest version if needed
M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv", dec=",")
pam50 <-  read.table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE, sep = "\t") 
rna_seq_df <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]
cell_cat <- rna_seq_df$Cell_Cat

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
trans_rna_seq_df <- trans_rna_seq_df %>% dplyr::select (-"rna_decon_sampleid")


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

#Barplot for comparing significant cell types

# Making table for cell counts and significant counts.
cell_cat_table <- data.frame(
  Category = rep(unique(cell_cat), 2),
  Value = c(15,3,10,3,9,19,5,8,4,4,37,2,9,2,4,0,1,7,1,1,0,2,17,0),
  Group = rep(c("Group 1", "Group 2"), each = 12)
)

# Plot with Overlaying Bars
ggplot(cell_cat_table, aes(x = reorder(Category,-Value), y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Cell Types",    # Change X-axis label
    y = "Count",    # Change Y-axis label
  ) 
  


# Kendall's correlation Immune Decon 
## Convert columns to numeric
# cont_tils <- as.numeric(as.character(trans_rna_seq_df$sTILs))
# 
# sig_cols <- setdiff(names(trans_rna_seq_df), c("INVESTIGATOR_SAMPLENAME", "sTILs_cat","sTILs"))
# corrList <- list()
# for (type in sig_cols) {
#   # Convert the current column to numeric
#   trans_rna_seq_df[[type]] <- as.numeric(as.character(trans_rna_seq_df[[type]]))
# 
#   # Perform Spearman correlation
#   correlation <- cor(cont_tils, trans_rna_seq_df[[type]], method = "kendall")
# 
#   # Store the result in the list
#   corrList[[type]] <- correlation
# }
# 
# # Convert corrList to a dataframe
# corr_df <- data.frame(
#   Cell_Type = names(corrList),
#   Kendalls_Correlation = unlist(corrList)
# )
# 
# # Barplots for top and bottom correlations
# top_positive <- corr_df %>% arrange(desc(Kendalls_Correlation)) %>% head(15)
# top_negative <-  corr_df %>% arrange(Kendalls_Correlation) %>% head(15)
# ggplot(top_positive, aes(x = reorder(Cell_Type, Kendalls_Correlation), y = Kendalls_Correlation)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +  # Flip for better readability
#   #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   labs(title = "Top 15 Most Positively Correlated Cell Types",
#        x = "Cell Type",
#        y = "Kendall's Correlation") +
#   theme_minimal()
# # Saving results
# write_xlsx(corr_df, "kendalls_decon_sTILs_results.xlsx")

### Module Analysis
data <- file.choose()
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

# Reformatting and merging
resignature_data <- resignature_data %>% mutate(INVESTIGATOR_SAMPLENAME = sub("^X", "", INVESTIGATOR_SAMPLENAME))

resignature_data <- merge(resignature_data,uniqueData[, c("INVESTIGATOR_SAMPLENAME", "sTILs_cat","sTILs")],by = "INVESTIGATOR_SAMPLENAME",all.y=TRUE)

 
## Correlation
## Convert columns to numeric
cont_tils <- as.numeric(as.character(resignature_data$sTILs))

sig_cols <- setdiff(names(resignature_data), c("INVESTIGATOR_SAMPLENAME", "sTILs_cat","sTILs"))
corrList <- list()
for (type in sig_cols) {
  # Convert the current column to numeric
  resignature_data[[type]] <- as.numeric(as.character(resignature_data[[type]]))

  # Perform correlation
  correlation <- cor(cont_tils, resignature_data[[type]], method = "kendall")

  # Store the result in the list
  corrList[[type]] <- correlation
}

# Convert corrList to a dataframe
corr_df <- data.frame(
  Cell_Type = names(corrList),
  Kendalls_Correlation = unlist(corrList)
)
# Barplots for top and bottom correlations
top_positive <- corr_df %>% arrange(desc(Kendalls_Correlation)) %>% head(15)
top_negative <-  corr_df %>% arrange(Kendalls_Correlation) %>% head(15)
ggplot(top_positive, aes(x = reorder(Cell_Type, Kendalls_Correlation), y = Kendalls_Correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip for better readability
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Top 15 Most Positively Correlated signatures",
       x = "Cell Type",
       y = "Kendall's Correlation") +
  theme_minimal()

# Saving results
 write_xlsx(corr_df, "kendalls_sTILs_module_results.xlsx")

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
#DESEQ2 
data <- file.choose()
ge_matrix <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\ge_matrix_40502.csv")

### Removing samples that do not have sTILs data
ge_matrix <- ge_matrix %>% select(1,all_of(uniqueData$rna_decon_sampleid))
### Reformating count matrix
row.names(ge_matrix) <- ge_matrix[,1]
ge_matrix <- ge_matrix[,-1]

ge_matrix <- round(ge_matrix)
### Generating sample list for colData
samples_list <- colnames(ge_matrix)


### Metadata for DESeq2
colData <- data.frame(
  condition = uniqueData$sTILs_cat,
  row.names = samples_list
)

### Constructing DESEq DataSet
dds <- DESeqDataSetFromMatrix(countData = ge_matrix,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
### Making results table
res <- results(dds, contrast = c("condition","High","Low"))

resOrdered <- res[order(res$padj),]
res_df <- as.data.frame(resOrdered)
res_df$genes <- row.names(res_df)
### Saving results

write_xlsx(res_df, "sTILs_DESeq2.xlsx")
# Filter for significant genes (adjusted p-value < 0.05 and abs(log2FoldChange) > 1)
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

sig_upregulated_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]

sig_downregulated_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
# Identify the top 15 significant genes by adjusted p-value (or use other criteria)
top_genes <- head(sig_genes[order(sig_genes$padj), ], 15)

### Making Volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=padj < 0.05), alpha=0.5) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title="sTILs (High vs Low)", x="Log2 Fold Change", y="-Log10(p-value)") +
  theme(legend.position="none") +
  geom_text(data=top_genes, aes(x=log2FoldChange, y=-log10(pvalue), label=rownames(top_genes)), size=2.5, vjust=-1, hjust=1)

# Isolating gene names
deg_genes <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ])
upreg_genes <- rownames(sig_upregulated_genes)
downreg_genes <- rownames(sig_downregulated_genes)

go_results <- enrichGO(gene = downreg_genes,
                       OrgDb = org.Hs.eg.db,   
                       keyType = "SYMBOL",    
                       ont = "BP",            
                       pAdjustMethod = "BH",  
                       qvalueCutoff = 0.05)  
# View the results
summary(go_results)
# Saving results
write_xlsx(as.data.frame(go_results), "sTILs_GO_downregulated_results.xlsx")

# Plot the GO enrichment results
plot1 <- dotplot(go_results)  + theme(axis.text.y = element_text(angle = 0, hjust = 1))
plot2 <- barplot(go_results) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
