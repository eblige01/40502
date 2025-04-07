library("tidyverse")
library("broom")
library("pbapply")
library("future.apply")
library("bigstatsr")
library("writexl")
library("survival")
library("RColorBrewer")
library(ggrepel)
# 40502 Outcome modules
M40502_joined_metadata <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv")
D40502_data <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv")
signature_data <- read_table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\cdt.txt")


# Creating Bc class subsets
hr_pos_subset <- D40502_data |> 
  filter(patid %in% M40502_joined_metadata$patid & bc_class == "HR" & strat2_recep == 1)

hr_pos_metadata <- M40502_joined_metadata |> 
  left_join(hr_pos_subset, by = "patid") |> 
  relocate(bestresp, .after = rna_decon_sampleid ) |> 
  select(patid:bestresp) |> 
  filter(!(is.na(bestresp)) & !(is.na(INVESTIGATOR_SAMPLENAME)))

# Making signature data tidy by removing extra columns and rows
tidy_signature_data <- signature_data |>  
  select(-c("GID","CLID","GWEIGHT")) |> 
  slice(-c(1,2)) |>
  column_to_rownames("NAME") 

# Removing leading X to match the metadata
colnames(tidy_signature_data) <- str_remove(colnames(tidy_signature_data),"^X")

# Removing columns not found in metadata
tidy_signature_data <- select(tidy_signature_data,any_of(hr_pos_metadata$INVESTIGATOR_SAMPLENAME))



transposed_signature_data <- as_tibble(t(tidy_signature_data),rownames = NA) |> 
  rownames_to_column("INVESTIGATOR_SAMPLENAME") |> 
  left_join(hr_pos_metadata |> select(INVESTIGATOR_SAMPLENAME,patid), by ="INVESTIGATOR_SAMPLENAME") |> 
  relocate(patid,.before=1) |> 
  select(-INVESTIGATOR_SAMPLENAME)

named_transposed_signature_data <- as_tibble(t(tidy_signature_data),rownames = NA)


# OR calculation ----------------------------------------------------------

OR_signature_data <-  named_transposed_signature_data |> 
   mutate(bestresp = hr_pos_metadata$bestresp[match(rownames(named_transposed_signature_data),hr_pos_metadata$INVESTIGATOR_SAMPLENAME)],.before = 1) |> 
   mutate(bestresp = ifelse(bestresp <= 2,1,0))

gene_signatures <- colnames(OR_signature_data |> select(-c(bestresp)))


OR_matrix <- as.matrix(OR_signature_data[, gene_signatures])

OR_matrix <- apply(OR_matrix, 2, as.numeric)

# Convert transposed_signature_data to big.matrix format
X <- as_FBM(OR_matrix)
y <- OR_signature_data$bestresp  # Response variable

# Perform univariate logistic regression
logistic_results <- big_univLogReg(X, y)

logistic_results$OR <- exp(logistic_results$estim)

logistic_results$CI_Lower <- exp(logistic_results$estim - 1.96 * logistic_results$std.err)
logistic_results$CI_Upper <- exp(logistic_results$estim + 1.96 * logistic_results$std.err)

logistic_results$or_P_value <- predict(logistic_results,log10 = FALSE)


logistic_results$Gene <- colnames(OR_matrix)
logistic_results$log2_OR <- log2(logistic_results$OR)

write_xlsx(logistic_results,"C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\Documents\\Research\\40502_rotation_project\\Results\\Outcomes\\LogReg_HR+_modules.xlsx")

# HR calculations ---------------------------------------------------------

hr_signature_data <- left_join(transposed_signature_data,hr_pos_subset |> select(patid,pfsmos,pfsstat),by = "patid") |> 
  relocate(c(pfsmos,pfsstat),.before = 1)
# Making character vectors numeric
hr_signature_data <- hr_signature_data %>%
  mutate(across(where(is.character), as.numeric))
gene_signatures <- colnames(hr_signature_data |> select(-c(pfsmos,pfsstat,patid)))


get_HR <- function(gene) {
  model <- coxph(Surv(pfsmos, pfsstat) ~ hr_signature_data[[gene]], data = hr_signature_data)
  HR <- exp(coef(model))
  CI <- exp(confint(model))
  hr_P_Value <- summary(model)$coefficients[, "Pr(>|z|)"]
  log2_HR <- log2(exp(coef(model)))
  return(data.frame(Gene = gene, HR = HR, CI_Lower = CI[1], CI_Upper = CI[2], hr_P_Value = hr_P_Value,log2_HR = log2_HR))
}

# Apply function to all gene signatures
cox_results <- do.call(rbind, pblapply(gene_signatures, get_HR))

write_xlsx(cox_results,"C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\Documents\\Research\\40502_rotation_project\\Results\\Outcomes\\Cox_HR+_modules.xlsx")


# Visualization -----------------------------------------------------------

#Combining data 
merged_results <- cox_results |> 
  full_join(logistic_results,by = "Gene") |> 
  column_to_rownames("Gene") |> 
  mutate(sig_stat = as.factor(case_when(
    hr_P_Value < 0.05 & or_P_value < 0.05 ~ 3,
    hr_P_Value < 0.05 & or_P_value >= 0.05 ~ 1,
    or_P_value < 0.05 & hr_P_Value >= 0.05 ~ 2,
    TRUE ~ 0
  ) )) |> 
  select(c(log2_HR,log2_OR,sig_stat,hr_P_Value,or_P_value))

#Selecting the signatures with most effect
top_genes_OR <- merged_results  |> 
  filter(sig_stat != 0) |> 
  arrange(abs(log2_OR))  |>   
  slice_tail(n = 10)  

top_genes_HR <- merged_results  |> 
  filter(sig_stat != 0) |> 
  arrange(abs(log2_HR))  |>   
  slice_tail(n = 10)  


combined_genes <- unique(bind_rows(top_genes_HR,top_genes_OR))

# Scatter plot
genome_scatter <- ggplot(merged_results,aes(log2_OR,log2_HR,color = sig_stat)) +
  geom_point() +
  geom_hline(aes(yintercept = 0),linetype = "dotted") +
  geom_vline(aes(xintercept = 0),linetype = "dotted") +
  scale_color_manual(
    values = c("0"="gray70","1"="#E41A1C","2" = "#377EB8","3"= "#4DAF4A"), 
    name = "Significant features (P < .05)",
    labels = c("0"="NS, 526 features","1"="PFS only, 384 features","2" = "Response only,101 features","3"= "PFS and Response, 10 features")) +
  labs(x = "Response(CR or PR) Log2 OR", y = "PFS Log2 HR") +
  geom_text_repel(
    data = combined_genes,
    aes(label = row.names(combined_genes)), 
    size = 2.5,  
    max.overlaps = 50,
    box.padding = .5,
    point.padding = 0.3,
    fontface = "bold",
    force = 10
  ) +
 theme_classic()
table(merged_results$sig_stat)
