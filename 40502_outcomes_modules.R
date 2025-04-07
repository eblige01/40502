library("tidyverse")
library("broom")
library("pbapply")
library("future.apply")
library("bigstatsr")
library("writexl")
library("survival")

# 40502 Outcome modules
M40502_joined_metadata <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv")
D40502_data <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv")
signature_data <- read_table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\cdt.txt")


# Removing rows from data that do not have patids in the metadata
D40502_data <- D40502_data |> 
  filter(patid %in% M40502_joined_metadata$patid)

M40502_joined_metadata <- M40502_joined_metadata |> 
  left_join(D40502_data, by = "patid") |> 
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
tidy_signature_data <- select(tidy_signature_data,any_of(M40502_joined_metadata$INVESTIGATOR_SAMPLENAME))


transposed_signature_data <- as_tibble(t(tidy_signature_data),rownames = NA) |> 
  rownames_to_column("INVESTIGATOR_SAMPLENAME") |> 
  left_join(M40502_joined_metadata |> select(INVESTIGATOR_SAMPLENAME,patid), by ="INVESTIGATOR_SAMPLENAME") |> 
  relocate(patid,.before=1) |> 
  select(-INVESTIGATOR_SAMPLENAME)

hr_signature_data <- left_join(transposed_signature_data,D40502_data |> select(patid,pfsmos,pfsstat),by = "patid") |> 
  relocate(c(pfsmos,pfsstat),.before = 1)
# For OR
# transposed_signature_data <-  transposed_signature_data |> 
#   mutate(bestresp = M40502_joined_metadata$bestresp[match(rownames(transposed_signature_data),M40502_joined_metadata$INVESTIGATOR_SAMPLENAME)],.before = 1) |> 
#   mutate(bestresp = ifelse(bestresp <= 2,1,0))

# OR calculation ----------------------------------------------------------

gene_signatures <- setdiff(colnames(transposed_signature_data), "bestresp")


new_matrix <- as.matrix(transposed_signature_data[, gene_signatures])

new_matrix <- apply(new_matrix, 2, as.numeric)

# Convert transposed_signature_data to big.matrix format
X <- as_FBM(new_matrix)
y <- transposed_signature_data$bestresp  # Response variable

# Perform univariate logistic regression
logistic_results <- big_univLogReg(X, y)

logistic_results$OR <- exp(logistic_results$estim)

logistic_results$CI_Lower <- exp(logistic_results$estim - 1.96 * logistic_results$std.err)
logistic_results$CI_Upper <- exp(logistic_results$estim + 1.96 * logistic_results$std.err)

logistic_results$P_value <- predict(logistic_results,log10 = FALSE)


logistic_results$gene_signatures <- colnames(new_matrix)
logistic_results$log2_OR <- log2(logistic_results$OR)

# HR calculations ---------------------------------------------------------

gene_signatures <- colnames(hr_signature_data |> select(-c(pfsmos,pfsstat,patid)))


get_HR <- function(gene) {
  model <- coxph(Surv(pfsmos, pfsstat) ~ hr_signature_data[[gene]], data = hr_signature_data)
  HR <- exp(coef(model))
  CI <- exp(confint(model))
  P_Value <- summary(model)$coefficients[, "Pr(>|z|)"]
  
  return(data.frame(Gene = gene, HR = HR, CI_Lower = CI[1], CI_Upper = CI[2], P_Value = P_Value))
}

# Apply function to all gene signatures
HR_results <- do.call(rbind, pblapply(gene_signatures, get_HR))

