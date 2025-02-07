# 40502 PAM50 analysis
library(ggplot2)
library(gplots)
library(dplyr)  
library(writexl)
library(FSA) 
library(vcd)

### Loading in the data

M40502_joined_metadata <- read.csv("~/Desktop/Research/StoverLab_rotation/data/40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("~/Desktop/Research/StoverLab_rotation/data/D40502_data.csv", dec=",")
pam50 <-  read.table("~/Desktop/Research/StoverLab_rotation/data/PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE, sep = "\t") 
rna_seq_df <- read.csv("~/Desktop/Research/StoverLab_rotation/data/rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]


### Reformating and merging 
M40502_joined_metadata <- M40502_joined_metadata[!is.na(M40502_joined_metadata$Label.on.Curls),]
names(pam50)[1] <- "INVESTIGATOR_SAMPLENAME"
mergeddata <- merge(D40502_data, M40502_joined_metadata[, c("patid", "sTILs","rna_decon_sampleid","INVESTIGATOR_SAMPLENAME")], by = "patid")
mergeddata$rna_decon_sampleid <- gsub("Sample_", "sample",mergeddata$rna_decon_sampleid)
mergeddata <- merge(mergeddata,pam50,by="INVESTIGATOR_SAMPLENAME")
mergeddata_tils_sub <- mergeddata[!is.na(mergeddata$sTILs),]
mergeddata_tils_sub$sTILs_cat <- ifelse(mergeddata_tils_sub$sTILs >= 5,"High","Low")

# ## Chi-Squared analysis TILS vs PAM50
# 
# ### Creating contingency table 
# chisq_table <- structable(mergeddata_tils_sub$Call,mergeddata_tils_sub$sTILs_cat)
# 
# ### Perform chi-squared test
# chi_test <- chisq.test(chisq_table)
# 
# ### Expected Results and observed results
# chi_test$observed
# chi_test$expected
# ### Mosaic plot for visualization 
# chisq_table <- t(chisq_table)
# 
# mosaic(chisq_table, 
#        shade = TRUE, 
#        legend = TRUE,
#        labeling_args = list(
#          gp_labels = gpar(fontsize = 7, col = "black"),  # Font size and color
#          rot_labels = 30),  # Rotate labels by 45 degrees
#        gp_legend = gpar(fontsize = 6))  # Font size for the legend labels

# PAM50 vs RNA Decon estimates
### Subsetting duplicate RNA seq sample IDs 
unique_PAM <- mergeddata %>%
  distinct(rna_decon_sampleid, .keep_all = TRUE)

# Reformatting and merging
trans_rna_seq_df <- as.data.frame(t(rna_seq_df))

trans_rna_seq_df$rna_decon_sampleid <- row.names(trans_rna_seq_df)

# Change this line for other analysis TILS/PAM50
trans_rna_seq_df <- merge(trans_rna_seq_df, unique_PAM[, c("rna_decon_sampleid", "Call")], by = "rna_decon_sampleid")
row.names(trans_rna_seq_df) <- trans_rna_seq_df$rna_decon_sampleid
trans_rna_seq_df <- trans_rna_seq_df %>% select (-"rna_decon_sampleid")

# Empty data frame for results
PAM50_results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE, significant_pairs = character())

# Get the gene column names (excluding sample_id and Call)
gene_cols <- setdiff(names(trans_rna_seq_df), c("rna_decon_sampleid", "Call"))

# Loop through each gene column
for (gene in gene_cols) {
  # Ensure the column is numeric
  trans_rna_seq_df[[gene]] <- as.numeric(trans_rna_seq_df[[gene]])

  # Perform Kruskal-Wallis test
  test_result <- kruskal.test(trans_rna_seq_df[[gene]] ~ trans_rna_seq_df$Call)

  # If Kruskal-Wallis test is significant, perform Dunn's test
  significant_pairs <- NA_character_
  if (test_result$p.value < 0.05) {
    dunn_result <- dunnTest(trans_rna_seq_df[[gene]] ~ trans_rna_seq_df$Call, method = "bh")

    # Check column names and extract significant pairwise comparisons
    if ("Comparison" %in% colnames(dunn_result$res) & "P.adj" %in% colnames(dunn_result$res)) {
      sig_pairs <- dunn_result$res %>%
        filter(P.adj < 0.05) %>%
        pull(Comparison)  # Extract significant pair names

      # Store significant comparisons as a comma-separated string
      if (length(sig_pairs) > 0) {
        significant_pairs <- paste(sig_pairs, collapse = ", ")
      }
    }
  }

  # Append results to module_results dataframe
  PAM50_results <- rbind(PAM50_results,
                          data.frame(Gene = gene,
                                     p_value = test_result$p.value,
                                     significant_pairs = significant_pairs))
}

# Adjust p-values
PAM50_results$p_adj <- p.adjust(PAM50_results$p_value, method = "bonferroni")

# Subset significant module_results after adjustment (e.g., FDR < 0.05)
PAM50_results_sig <-PAM50_results %>% filter(p_adj < 0.05)
#Saving results
 write_xlsx(PAM50_results_sig, "PAM50_decon.xlsx")

