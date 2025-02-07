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