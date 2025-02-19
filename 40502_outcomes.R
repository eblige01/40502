# 40502 RNA features vs outcomes  analysis 
library(survival)
library(dplyr)

# Loading data

M40502_joined_metadata <- read.csv("~/Desktop/Research/StoverLab_rotation/data/40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("~/Desktop/Research/StoverLab_rotation/data/D40502_data.csv", dec=",")
rna_seq_df <- read.csv("~/Desktop/Research/StoverLab_rotation/data/rna_decon_matrix_40502.csv", dec=",")
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


pac_sub <- merge(trans_rna_seq_df,pac_sub[c("bestresp","rna_decon_sampleid")],by = "rna_decon_sampleid")

# Removing patients with no tumor response data
pac_sub <- na.omit(pac_sub)

# Creating response column 
pac_sub$responder_stat <- ifelse(pac_sub$bestresp <= 2,"responder","nonresponder")

# # Response Vs Decon PAc and HR+
# # Initialize an empty dataframe to store module_module_results
# rna_results <- data.frame(Cell = character(), p_value = numeric(), stringsAsFactors = FALSE,higher_expression = character())
# 
# # Get the gene column names (excluding sample_id and TIL_group)
# cell_cols <- setdiff(names(pac_sub), c("rna_decon_sampleid", "bestresp","responder_stat"))
# 
# # Loop through each gene column
# for (cell in cell_cols) {
#   # Ensure the column is numeric
#   pac_sub[[cell]] <- as.numeric(pac_sub[[cell]])
#   # Perform Wilcoxon test (Mann-Whitney U test)
#   test_result <- wilcox.test(pac_sub[[cell]] ~ pac_sub$responder_stat)
# 
#   # Calculate mean expression for each group will be used to determine the dirrection of sig module_results
#   group_means <- pac_sub %>%
#     group_by(responder_stat) %>%
#     summarise(mean_expression = mean(!!sym(cell), na.rm = TRUE)) %>%
#     arrange(desc(mean_expression))
# 
#   # Determine which group has higher expression
#   higher_group <- group_means$responder_stat[1]
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
# write_xlsx(rna_results, "response_pac_decon.xlsx")

