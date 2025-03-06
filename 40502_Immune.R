# 40502 Immune (TILSs) analysis 
library(ggplot2)
library("gplots")
library(dplyr)  
library(writexl)
library(DESeq2)
library(org.Hs.eg.db)  
library("clusterProfiler")
library(extrafont)
library("corrplot")
library(ggrepel)

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
# cell_cat_table <- data.frame(
#   Category = rep(unique(cell_cat), 2),
#   Value = c(15,3,10,3,9,19,5,8,4,4,37,2,9,2,4,0,1,7,1,1,0,2,17,0),
#   Group = rep(c("Group 1", "Group 2"), each = 12)
# )
# 
# # Plot with Overlaying Bars
# ggplot(cell_cat_table, aes(x = reorder(Category,-Value), y = Value, fill = Group)) +
#   geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 11.5),
#         axis.text.y = element_text(size = 9),
#         legend.title = element_blank(),
#         panel.grid.major = element_blank(),  # Remove major gridlines
#         panel.grid.minor = element_blank(),  # Remove minor gridlines
#         axis.line = element_line(color = "black", size = 0.8), # Add axis lines
#         axis.ticks = element_line(color = "black")) +
#   scale_y_continuous(expand = c(0, 0),
#         breaks = seq(0, max(cell_cat_table$Value), by = 5)) + # Makes y = 0 touch the x axis and decreases increments
#   scale_fill_manual( values = c("Group 1" = "grey56","Group 2" = "blue"), # this line is to change the labels for legends and fix the overlapping color problem
#         labels = c("Total estimates","Different between sTILs")) +
#   labs(
#     x = "Cell Types",    # Change X-axis label
#     y = "Estimate count",    # Change Y-axis label
#   ) 
#   





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

signature_data <- read.table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\cdt.txt", header = TRUE, sep = "\t", comment.char = "", quote = "")

# Reformatting signature data for analysis
# Removing unnecessary rows and coloumns
resignature_data <- signature_data %>% select(-c("GID","CLID","GWEIGHT"))
resignature_data <- resignature_data %>% dplyr :: slice(-c(1,2))

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


# Initialize lists for correlations and p-values
corrList <- list()
pValueList <- list()

# Convert continuous sTILs to numeric
cont_tils <- as.numeric(as.character(resignature_data$sTILs))

sig_cols <- setdiff(names(resignature_data), c("INVESTIGATOR_SAMPLENAME", "sTILs_cat", "sTILs"))

# Loop through each column and compute Kendall's correlation and p-value
for (type in sig_cols) {
  resignature_data[[type]] <- as.numeric(as.character(resignature_data[[type]]))
  
  test_result <- cor.test(cont_tils, resignature_data[[type]], method = "kendall")
  
  corrList[[type]] <- test_result$estimate  # Kendall's Tau
  pValueList[[type]] <- test_result$p.value # p-value
}

# Convert to data frame
corr_df <- data.frame(
  Cell_Type = names(corrList),
  Kendalls_Correlation = unlist(corrList),
  P_Value = unlist(pValueList)
)

# Add significance stars (e.g., * p < 0.05, ** p < 0.01, *** p < 0.001)
corr_df <- corr_df %>%
  mutate(Significance = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01  ~ "**",
    P_Value < 0.05  ~ "*",
    TRUE            ~ ""
  ))



# lollipop plot for top and bottom correlations
top_positive <- corr_df %>% arrange(desc(Kendalls_Correlation)) %>% head(20)
top_negative <-  corr_df %>% arrange(Kendalls_Correlation) %>% head(20)


ggplot(top_positive, aes(x = reorder(Cell_Type, Kendalls_Correlation), y = Kendalls_Correlation)) +
  geom_segment(aes(xend = Cell_Type, y = 0, yend = Kendalls_Correlation), color = "steelblue", size = 1) +  # Lollipop stem
  geom_point(color = "red", size = 4) +  # Lollipop head
  coord_flip() +  # Flip for readability
  labs(x = "Signatures",
       y = "Kendall's Correlation") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Removes all background grid lines
    axis.line.y = element_line(color = "black"),  # Keeps only the y-axis line
    axis.line.x = element_line(color = "black"), # Removes x-axis line if unnecessary
    plot.margin = margin(20, 30, 20, 30)
  ) +
  expand_limits(y = .5) +
  scale_y_continuous(expand = c(0, 0),)
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
# Adding direction column for res_df
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$pvalue < .001] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$pvalue < .001] <- "DOWN"

# Subsetting differentially expressed genes
diff_expr <- subset(res_df, res_df$diffexpressed != "NO")
write_xlsx(res_df, "sTILs_DESeq2.xlsx")
# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
top_genes<- res_df[abs(res_df$log2FoldChange) > 1, ]  # Filter genes with |Log2FC| > 1
top_genes <- top_genes[order(top_genes$padj), ]  # Order by significance (padj)
top_genes <- head(top_genes, 30)  # Select the top 30

# Assign gene names for labeling, else NA
res_df$delabel <- ifelse(res_df$genes %in% top_genes$genes, res_df$genes, NA)

### Making Volcano plot


ggplot(res_df,aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel )) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.001), col = "gray", linetype = 'dashed') + 
  geom_point() + 
  scale_color_manual(values = c("blue","grey","red"),
                     labels = c("Downregulated","Not significant","Upregulated")) +
  coord_cartesian(ylim = c(0, 30), xlim = c(-6, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs( x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 6, 2)) + # to customise the breaks in the x axis
  geom_text_repel(max.overlaps = Inf, size = 3, force = 3 , color = "black") +
  ggtitle('TILS High(>=5%) vs Low(<5%)') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(.7), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(.7), color = 'black'),
        plot.title = element_text(hjust = 0.5,size = rel(.9),face = "bold"))
  
# Isolating gene names
upreg_genes <- (subset(diff_expr,diff_expr$diffexpressed == "UP" ))$genes
downreg_genes <- (subset(diff_expr,diff_expr$diffexpressed == "DOWN" ))$genes
go_results <- enrichGO(gene = upreg_genes,
                       OrgDb = org.Hs.eg.db,   
                       keyType = "SYMBOL",    
                       ont = "BP",            
                       pAdjustMethod = "BH",  
                       qvalueCutoff = 0.05)  
# View the results
summary(go_results)
# Saving results
write_xlsx(as.data.frame(go_results), "sTILs_GO_upregulated_results.xlsx")

# Plot the GO enrichment results
dotplot(go_results,x= "count") +  
            theme(axis.text.y = element_text(angle = 0, hjust = 1),
            plot.title = element_text(hjust = 0.5,size = rel(2),face = "bold")) + 
            ggtitle("Pathways of downregualated genes")
plot2 <- barplot(go_results) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
