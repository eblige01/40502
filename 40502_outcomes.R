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


pac_sub <- merge(trans_rna_seq_df,pac_sub[c("survmos","pfsmos","pfsstat","survstat","rna_decon_sampleid")],by = "rna_decon_sampleid")
# Converting for cox model
pac_sub$pfsmos <- as.numeric(pac_sub$pfsmos)
# Making cell type list for loop
cell_types <-colnames(pac_sub) 
cell_types <- setdiff(cell_types,c("survmos","pfsmos","pfsstat","survstat","rna_decon_sampleid"))

# Create an empty list to store results
cox_results <- list()

# Loop through each cell type and run Cox regression
for (cell in cell_types) {
  # Create formula
  cox_formula <- as.formula(paste("Surv(pfsmos, pfsstat) ~", cell))
  
  # Fit Cox model
  cox_model <- coxph(cox_formula, data = pac_sub)
  
  # Extract summary
  cox_summary <- summary(cox_model)
  
  # Store HR, confidence interval, and p-value
  cox_results[[cell]] <- data.frame(
    CellType = cell,
    HR = exp(cox_summary$coefficients[1, 1]),
    Lower_CI = exp(cox_summary$conf.int[1, 3]),
    Upper_CI = exp(cox_summary$conf.int[1, 4]),
    p_value = cox_summary$coefficients[1, 5]
  )
}

