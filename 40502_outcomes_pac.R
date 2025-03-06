# 40502 Paclitaxel among HR+  analysis 
library(survival)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(writexl)
# Loading data

M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv", dec=",")
rna_seq_df <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\rna_decon_matrix_40502.csv", dec=",")
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
ge_matrix <- ge_matrix %>% dplyr :: select(1,all_of(pac_sub$rna_decon_sampleid))
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
# Adding direction column for res_df
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$pvalue < .05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$pvalue < .05] <- "DOWN"

# Subsetting differentially expressed genes
diff_expr <- subset(res_df, res_df$diffexpressed != "NO")

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
top_genes<- res_df[abs(res_df$log2FoldChange) > 1, ]  # Filter genes with |Log2FC| > 1
top_genes <- top_genes[order(top_genes$padj), ]  # Order by significance (padj)
top_genes <- head(top_genes, 30)  # Select the top 30

# Assign gene names for labeling, else NA
res_df$delabel <- ifelse(res_df$genes %in% top_genes$genes, res_df$genes, NA)

### Making Volcano plot

ggplot(res_df,aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel )) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() + 
  scale_color_manual(values = c("blue","grey","red"),
                     labels = c("Downregulated","Not significant","Upregulated")) +
  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs( x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-5, 5, 2)) + # to customise the breaks in the x axis
  geom_text_repel(max.overlaps = Inf, size = 3, force = 3 , color = "black") +
  ggtitle('Paclitaxel Responders vs Nonresponders') +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(.7), color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(.7), color = 'black'),
        plot.title = element_text(hjust = 0.5,size = rel(.9),face = "bold"))

# Isolating gene names
upreg_genes <- (subset(diff_expr,diff_expr$diffexpressed == "UP" ))$genes
downreg_genes <- (subset(diff_expr,diff_expr$diffexpressed == "DOWN" ))$genes
diff_express_genes <- diff_expr$genes

