library(ggplot2)
library("gplots")
library(dplyr)
# Loading in the data

M40502_joined_metadata <- read.csv("~/Desktop/Research/StoverLab_rotation/data/40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("~/Desktop/Research/StoverLab_rotation/data/D40502_data.csv", dec=",")
pam50 <-  read.table("~/Desktop/Research/StoverLab_rotation/data/PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE, sep = "\t") 
rna_seq_df <- read.csv("~/Desktop/Research/StoverLab_rotation/data/rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]


# Reformating and merging 
M40502_joined_metadata <- M40502_joined_metadata[!is.na(M40502_joined_metadata$Label.on.Curls),]
names(M40502_joined_metadata)[4] <- "slideid"
names(pam50)[1] <- "slideid"
mergeddata <- merge(D40502_data, M40502_joined_metadata[, c("patid","slideid", "sTILs","rna_decon_sampleid")], by = "patid")
mergeddata$rna_decon_sampleid <- gsub("Sample_", "sample",mergeddata$rna_decon_sampleid)
mergeddata <- merge(mergeddata,pam50,by="slideid")
mergeddata <- mergeddata[!is.na(mergeddata$sTILs),]

# Removing duplicate paitent IDs by selecting the highest TILS
mergeddata <- mergeddata[!is.na(mergeddata$sTILs), ]
sortedData <- mergeddata[order(mergeddata$patid),]
uniqueData <- sortedData %>%
  group_by(patid) %>%
  slice_max(order_by = sTILs, n = 1) %>%
  ungroup()
head(uniqueData)
uniqueData<- unique(uniqueData)
uniqueData <- uniqueData %>%
  group_by(rna_decon_sampleid) %>%
  filter(n() == 1) %>%
  ungroup()

sampleID <- uniqueData$rna_decon_sampleid
uniqueData$rna_decon_sampleid <- gsub("Sample_", "sample",sampleID)
uniqueData$sTILs <- ifelse(uniqueData$sTILs >= 5,"High","Low")
# Subsetting RNA SEQ to only containing high TILS samples

rna_seq_df <- rna_seq_df[,colnames(rna_seq_df) %in% sampleID]
ncol(rna_seq_df)

rownames(uniqueData) <- uniqueData$rna_decon_sampleid

trans_rna_seq_df <- t(rna_seq_df)
columnname2 <- c("B cell memory_CIBERSORT","T cell CD8+_CIBERSORT","NK cell activated_CIBERSORT","Macrophage M2_CIBERSORT","B cell_QUANTISEQ","Macrophage M1_QUANTISEQ","Macrophage M2_QUANTISEQ","NK cell_QUANTISEQ","T cell CD8+_QUANTISEQ","B cell_EPIC","T cell CD8+_EPIC","Macrophage_EPIC")
                 
cibersortdf <- as.data.frame(trans_rna_seq_df[,colnames(trans_rna_seq_df) %in% columnname2])
cibersortdf$rownames <- row.names(cibersortdf)
uniqueData$rownames <- row.names(uniqueData)
mergeCiber <- merge(cibersortdf,uniqueData, by="rownames")
mergeCiber <- mergeCiber[,-1]
columnmname1 <- c("B cell memory_CIBERSORT","T cell CD8+_CIBERSORT","NK cell activated_CIBERSORT","Macrophage M2_CIBERSORT","B cell_QUANTISEQ","Macrophage M1_QUANTISEQ","Macrophage M2_QUANTISEQ","NK cell_QUANTISEQ","T cell CD8+_QUANTISEQ","B cell_EPIC","T cell CD8+_EPIC","Macrophage_EPIC","Call","bc_class","sTILs")
mergeCiber <- mergeCiber[,colnames(mergeCiber) %in% columnmname1]

# # Spearman correlation 
# ## Convert columns to numeric 
# mergeCiber$sTILs <- as.numeric(as.character(mergeCiber$sTILs))
# corrList <- list()
# for (type in columnname2) {
#   # Convert the current column to numeric
#   mergeCiber[[type]] <- as.numeric(as.character(mergeCiber[[type]]))
#   
#   # Perform Spearman correlation
#   correlation <- cor(mergeCiber$sTILs, mergeCiber[[type]], method = "spearman")
#   
#   # Store the result in the list
#   corrList[[type]] <- correlation
# }
# 
# # Convert corrList to a dataframe
# corr_df <- data.frame(
#   Cell_Type = names(corrList),
#   Correlation = unlist(corrList)
# )
# 
# # Convert corrList to a dataframe
# corr_df <- data.frame(
#   Cell_Type = names(corrList),
#   Correlation = as.numeric(unlist(corrList))
# )
# 
# # Create a bar plot using ggplot2
# ggplot(corr_df, aes(x = Cell_Type, y = Correlation)) +
#   geom_bar(stat = "identity", fill = "skyblue") +  # Create bars
#   labs(y = "Spearman Correlation") +
#   theme_minimal() +
#   geom_segment(aes(x=0,xend=0,y=0,yend=0.4),color = "black",size = .5)+
#   geom_hline(yintercept = 0,color = "black",size=.5) +
#   theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 65, hjust = 1,size =7.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())  # Rotate x-axis labels for better readability

# Making box plot with colored dots
mergeCiber$sTILs <- as.factor(mergeCiber$sTILs)
mergeCiber$`T cell CD8+_CIBERSORT` <- as.numeric(mergeCiber$`T cell CD8+_CIBERSORT`)

ggplot(mergeCiber, aes(x = sTILs, y = `T cell CD8+_CIBERSORT`)) +
  geom_boxplot(fill = "black", color = "black", outlier.shape = NA) +  # Black boxes
  geom_jitter(aes(color = factor(Call)), width = 0.2, size = 2) +   # Colored dots for subtypes
  labs(title = "TILs Levels by Subtype",
       x = "TILs Group",
       y = "T cell CD8+_CIBERSORT") +
  scale_color_manual(values = c("red", "pink", "lightblue", "blue", "green")) +  # Assign colors to subtypes
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Customize x-axis label size
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) 

# making violin plot
mergeCiber$`Macrophage M2_CIBERSORT` <- as.numeric(mergeCiber$`Macrophage M2_CIBERSORT`)
ggplot(mergeCiber, aes(x = sTILs, y = `Macrophage M2_CIBERSORT`)) +
  geom_violin(fill = "black", color = "black") +  # Create black violin shapes
  geom_jitter(aes(color = factor(Call)), width = 0.2, size = 1) +   # Colored dots for subtypes
  labs(color = "PAM50 subtype",
       x = "TILs Group",
       y = "Macrophage M2_CIBERSORT") +
  scale_color_manual(values = c("red", "pink", "lightblue", "blue", "green")) +  # Assign colors to subtypes
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Customize x-axis label size
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

