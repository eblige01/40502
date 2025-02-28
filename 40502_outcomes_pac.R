# 40502 Paclitaxel among HR+  analysis 
library(survival)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(writexl)
# Loading data

M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\NCTN-D3-recoded.csv", dec=",")
rna_seq_df <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]

# Reformating for the analysis
trans_rna_seq_df <- as.data.frame(t(rna_seq_df))
trans_rna_seq_df$rna_decon_sampleid <- row.names(trans_rna_seq_df)

# Reformating and merging 
M40502_joined_metadata <- M40502_joined_metadata[!is.na(M40502_joined_metadata$Label.on.Curls),]
mergeddata <- merge(D40502_data, M40502_joined_metadata[, c("patid","rna_decon_sampleid","INVESTIGATOR_SAMPLENAME")], by = "patid")
mergeddata$rna_decon_sampleid <- gsub("Sample_", "sample",mergeddata$rna_decon_sampleid)


#Subset data to only include paclitaxel treatment and patients that are HR+

pac_sub <- mergeddata %>% filter(arm == 1 , strat2_recep == 1)

# Creating response column 
pac_sub$responder_stat <- ifelse(pac_sub$bestresp <= 2,"responder","nonresponder")

# Removing patients with no tumor response data
pac_sub <- pac_sub %>% filter(!is.na(responder_stat))

pac_sub_rna <- merge(trans_rna_seq_df,pac_sub[c("responder_stat","rna_decon_sampleid")],by = "rna_decon_sampleid")
# Pac_response Vs Decon
# Initialize an empty dataframe to store module_module_results
rna_results <- data.frame(Cell = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())

# Get the gene column names (excluding sample_id and TIL_group)
cell_cols <- setdiff(names(pac_sub_rna), c("rna_decon_sampleid", "bestresp","responder_stat"))

# Loop through each gene column
for (cell in cell_cols) {
  # Ensure the column is numeric
  pac_sub_rna[[cell]] <- as.numeric(pac_sub_rna[[cell]])
  # Perform Wilcoxon test (Mann-Whitney U test)
  test_result <- wilcox.test(pac_sub_rna[[cell]] ~ pac_sub_rna$responder_stat)

  # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
  group_means <- pac_sub_rna %>%
    group_by(responder_stat) %>%
    summarise(mean_expression = mean(!!sym(cell), na.rm = TRUE)) %>%
    arrange(desc(mean_expression))

  # Determine which group has higher expression
  higher_group <- group_means$responder_stat[1]

  # Append module_results to dataframe
  rna_results <- rbind(rna_results, data.frame(Cell = cell, p_value = test_result$p.value,higher_expression = higher_group))
}

rna_results$p_adj <- p.adjust(rna_results$p_value, method = "bonferroni")

# Subset significant module_results after adjustment (e.g., FDR < 0.05)
significant_rna_results <- rna_results %>% filter(p_adj < 0.05)
sig_cells <- significant_rna_results$Cell

#Saving results
write_xlsx(rna_results, "response_pac_decon.xlsx")

# # Module analysis
# signature_data <- read.table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502\\Data\\cdt.txt", header = TRUE, sep = "\t", comment.char = "", quote = "")

# 
# # Reformatting signature data for analysis
# # Removing unnecessary rows and coloumns
# resignature_data <- signature_data %>% select(-c("GID","CLID","GWEIGHT"))
# resignature_data <- resignature_data %>% dplyr::slice(-c(1,2))
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
# resignature_data <- merge(resignature_data,pac_sub[, c("INVESTIGATOR_SAMPLENAME", "responder_stat")],by = "INVESTIGATOR_SAMPLENAME",all.y=TRUE)
# 
# # Pac_response vs Modules
# # Initialize an empty dataframe to store module_results
# module_results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())
# 
# # Get the gene column names (excluding sample_id and TIL_group)
# gene_cols <- setdiff(names(resignature_data), c("INVESTIGATOR_SAMPLENAME", "responder_stat"))
# 
# # Loop through each gene column
# for (gene in gene_cols) {
#   # Ensure the column is numeric
#   resignature_data[[gene]] <- as.numeric(resignature_data[[gene]])
#   # Perform Wilcoxon test (Mann-Whitney U test)
#   test_result <- wilcox.test(resignature_data[[gene]] ~ resignature_data$responder_stat)
# 
#   # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
#   group_means <- resignature_data %>%
#     group_by(responder_stat) %>%
#     summarise(mean_expression = mean(!!sym(gene), na.rm = TRUE)) %>%
#     arrange(desc(mean_expression))
# 
#   # Determine which group has higher expression
#   higher_group <- group_means$responder_stat[1]
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
# write_xlsx(module_results, "response_pac_modules.xlsx")

#DESEQ2 

ge_matrix <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\ge_matrix_40502.csv")

# Removing Nas from RNA IDs
pac_sub <- pac_sub %>% filter(!is.na(rna_decon_sampleid))
### Removing samples that do not have responder data
ge_matrix <- ge_matrix %>% select(1,all_of(pac_sub$rna_decon_sampleid))
### Reformating count matrix
row.names(ge_matrix) <- ge_matrix[,1]
ge_matrix <- ge_matrix[,-1]

ge_matrix <- round(ge_matrix)
samples_list <- colnames(ge_matrix)



# Removing samples that do not have any RNA data from pac_sub
pac_sub <-  pac_sub %>% filter(rna_decon_sampleid %in% samples_list)
# Removing duplicates
pac_sub <- pac_sub %>% distinct()

### Metadata for DESeq2
colData <- data.frame(
  condition = pac_sub$responder_stat,
  row.names = samples_list
)

### Constructing DESEq DataSet
dds <- DESeqDataSetFromMatrix(countData = ge_matrix,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
### Making results table
res <- results(dds, contrast = c("condition","responder","nonresponder"))

resOrdered <- res[order(res$padj),]
res_df <- as.data.frame(resOrdered)
res_df$genes <- row.names(res_df)
### Saving results

write_xlsx(res_df, "response_pac_DESeq2.xlsx")
# Filter for significant genes (adjusted p-value < 0.05 and abs(log2FoldChange) > 1)
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

# Identify the top 15 significant genes by adjusted p-value (or use other criteria)
top_genes <- head(sig_genes[order(sig_genes$padj), ], 15)

# Filter for significant upregulated genes
sig_upregulated_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]

# Identify the top 15 significant upregulated genes by adjusted p-value
top_upregulated_genes <- head(sig_upregulated_genes[order(sig_upregulated_genes$padj), ], 15)

### Making Volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=padj < 0.05), alpha=0.5) +  # Points for all genes
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title="Paclitaxel Response", x="Log2 Fold Change", y="-Log10(p-value)") +
  theme(legend.position="none") +
  geom_text(data=top_genes, aes(x=log2FoldChange, y=-log10(pvalue), label=rownames(top_genes)), size=2.5, vjust=-1, hjust=1)

deg_genes <- rownames(sig_genes)
sig_upregulated_list <- rownames(sig_upregulated_genes)
# # GO enrichment
# go_results <- enrichGO(gene = sig_upregulated_list,
#                        OrgDb = org.Hs.eg.db,   
#                        keyType = "SYMBOL",    
#                        ont = "BP",            
#                        pAdjustMethod = "bonferroni",  
#                        qvalueCutoff = 0.05)  
# # View the results
# summary(go_results)
# # Saving results
# write_xlsx(as.data.frame(go_results), "sTILs_GO_results.xlsx")
# 
# # Plot the GO enrichment results
# plot1 <- dotplot(go_results)  + theme(axis.text.y = element_text(angle = 0, hjust = 1))
# plot2 <- barplot(go_results) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


