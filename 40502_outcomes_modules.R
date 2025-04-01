library("tidyverse")


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






          

