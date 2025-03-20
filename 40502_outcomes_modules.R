library("dplyr")
library("broom")
library("ggplot2")
library("biglm") 
library("pbapply")

# 40502 Outcome script
M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv", dec=",")
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
mergeddata <- merge(D40502_data, M40502_joined_metadata[, c("patid","INVESTIGATOR_SAMPLENAME")], by = "patid")
mergeddata <- merge(resignature_data,mergeddata[, c("INVESTIGATOR_SAMPLENAME", "bestresp")],by = "INVESTIGATOR_SAMPLENAME")

# Removing values with no valesu for bestresp
mergeddata <-  mergeddata %>% filter(!is.na(bestresp))

#Adding binary outcome for pCR
mergeddata$pCR <- ifelse(mergeddata$bestresp %in% c(1,2), 1,0)

#Extracting signatures from data frame
features <- colnames(mergeddata)[!(colnames(mergeddata) %in% c ("INVESTIGATOR_SAMPLENAME","bestresp","pCR"))]

# logistic regression for all signatures

results <- pblapply(features, function(f) {
  formula <- as.formula(paste("pCR ~", f))
  model <- bigglm(formula, data = mergeddata, family = binomial())
  tidy_model <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  tidy_model$term <- f
  return(tidy_model[2,])
})

pcr_or_df <- do.call(rbind,results)

# Extract the signature you want to calculate OR for
feature <- "Up_regulated_upon_N_RAS_repression_Cell_Rep.2015_PMID.26166574"  # Replace this with the specific feature/column name you want

# Logistic regression for the chosen signature
formula <- as.formula(paste("pCR ~", feature))
model <- bigglm(formula, data = mergeddata, family = binomial())

# Tidy model to get OR, CI, and p-value
tidy_model <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)

# Extract and print OR (odds ratio) and its confidence interval for the chosen feature
tidy_model$term <- feature
print(tidy_model[2, ])  # The 2nd row corresponds to the OR for the feature
