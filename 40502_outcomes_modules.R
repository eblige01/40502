library("tidyverse")
library("broom")
library("biglm") 
library("pbapply")

# 40502 Outcome modules
M40502_joined_metadata <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv")
D40502_data <- read_csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv")
signature_data <- read_table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\cdt.txt")

# Removing rows with NA for Investigator_SAMPLENAME
M40502_joined_metadata <- M40502_joined_metadata |> 
  filter(!(is.na(INVESTIGATOR_SAMPLENAME)))

# Removing rows from data that do not have patids in the metadata
D40502_data <- D40502_data |> 
  filter(patid %in% M40502_joined_metadata$patid)

# Removing duplicate
# Making signature data tidy by removing extra columns and rows
tidy_signature_data <- signature_data |>  
  select(-c("GID","CLID","GWEIGHT")) |> 
  slice(-c(1,2)) |>
  column_to_rownames("NAME") 
  
# Removing leading X to match the metadata
colnames(tidy_signature_data) <- str_remove(colnames(tidy_signature_data),"^X")

# Removing columns not found in metadata
tidy_signature_data <- select(tidy_signature_data,any_of(M40502_joined_metadata$INVESTIGATOR_SAMPLENAME))

match(D40502_data$patid,M40502_joined_metadata$patid)

# Use mutate to create bestrep in M4 so we dont lose data due to duplicate paitent IDS



