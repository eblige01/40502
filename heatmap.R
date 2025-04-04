### NHS Heatmap Code ##
#######################
library(gplots)
library(dplyr)
library(tidyr)
library(writexl)
 


# Loading in the data 
M40502_joined_metadata <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\40502_joined_metadata_fixed.csv", dec=",")
D40502_data <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\NCTN-D3-recoded.csv", dec=",")
pam50 <-  read.table("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\PAM50scores_C40502_ZHAO4_AFM_09.16.21_pam50scores.txt", header = TRUE, sep = "\t") 
rna_seq_df <- read.csv("C:\\Users\\blig02\\OneDrive - The Ohio State University Wexner Medical Center\\40502_data\\Data\\rna_decon_matrix_40502.csv", dec=",")
rownames(rna_seq_df) <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]

# Subsetting based on Cell_Cat
#rna_seq_df <- subset(rna_seq_df, Cell_Cat == "Macrophage/Monocyte")

# Reordering RNA seq based on cell cat
rna_seq_df$Cell_Cat <- factor(rna_seq_df$Cell_Cat, levels = c("T cell","Macrophage/Monocyte","B cell","Dendritic cell","Granulocytes","NK cell","Mast cell","Score","Progenitor","Endothelial cell","Cancer associated fibroblast","Uncharacterized"))
rna_seq_df <- rna_seq_df %>% arrange(Cell_Cat)
# Making column colors for heatmap

cell_cat2 <- rna_seq_df$Cell_Cat
cell_cat2 <- gsub("B cell","yellow",cell_cat2)
cell_cat2 <- gsub("Dendritic cell","firebrick",cell_cat2)
cell_cat2 <- gsub("Granulocytes","goldenrod",cell_cat2)
cell_cat2 <- gsub("Mast cell","forestgreen",cell_cat2)
cell_cat2 <- gsub("Progenitor","darkorange",cell_cat2)
cell_cat2 <- gsub("T cell","deepskyblue",cell_cat2)
cell_cat2 <- gsub("Cancer associated fibroblast","mediumorchid",cell_cat2)
cell_cat2 <- gsub("Endothelial cell","chocolate",cell_cat2)
cell_cat2 <- gsub("Macrophage/Monocyte","yellowgreen",cell_cat2)
cell_cat2 <- gsub("NK cell","slateblue",cell_cat2)
cell_cat2 <- gsub("Score","tomato",cell_cat2)
cell_cat2 <- gsub("Uncharacterized","cyan4",cell_cat2)



# Removing Cell_Cat col
cell_cat <- rna_seq_df[,1]
rna_seq_df <- rna_seq_df[,-1]


table(rna_seq_df$Cell_Cat)

#Investigating size differences:


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
uniqueData$tils <- ifelse(uniqueData$sTILs >= 5,1,0)
# Subsetting RNA SEQ to only containing high TILS samples

rna_seq_df <- rna_seq_df[,colnames(rna_seq_df) %in% sampleID]
ncol(rna_seq_df)
rownames(uniqueData) <- uniqueData$rna_decon_sampleid

rna_seq_df <- rna_seq_df[, colSums(rna_seq_df != 0) > 0]
trans_rna <- t(rna_seq_df)


columnname2 <- colnames(trans_rna)
#columnname2 <- columnname2[grep("CIBERSORT",columnname2)]
cibersortdf <- as.data.frame(trans_rna[,colnames(trans_rna) %in% columnname2])
cibersortdf$rownames <- row.names(cibersortdf)
uniqueData$rownames <- row.names(uniqueData)
mergeCiber <- merge(cibersortdf,uniqueData, by="rownames")
row.names(mergeCiber) <- mergeCiber$rownames
mergeCiber <- mergeCiber[,-1]
columnmname1 <- c(columnname2,"Call","bc_class","tils")
mergeCiber <- mergeCiber[,colnames(mergeCiber) %in% columnmname1]

# Order by specific feature
mergeCiber <- mergeCiber %>% group_by(tils)  %>% arrange(tils)

# Formatting for heatmap color bars
PAM50vector=mergeCiber$Call
PAM50vector2=PAM50vector
PAM50vector2=gsub("Basal", "red", PAM50vector2)
PAM50vector2=gsub("LumA", "lightblue", PAM50vector2)
PAM50vector2=gsub("LumB", "blue", PAM50vector2)
PAM50vector2=gsub("Normal", "green", PAM50vector2)
PAM50vector2=gsub("Her2", "pink", PAM50vector2)

SubtypeVector= mergeCiber$bc_class
SubtypeVector2=SubtypeVector
SubtypeVector2=gsub("TNBC", "brown", SubtypeVector2)
SubtypeVector2=gsub("HR+", "orange", SubtypeVector2)
SubtypeVector2=gsub("HER2+", "yellow", SubtypeVector2)
SubtypeVector2=gsub('\\+', "", SubtypeVector2)

tilsVector <- mergeCiber$tils
tilsVector2 <- tilsVector
tilsVector2 <- gsub("1","grey",tilsVector2)
tilsVector2 <-  gsub("0","black",tilsVector2)
tils_cat <- tilsVector
tils_cat <-  gsub("1","High",tils_cat)
tils_cat <-  gsub("0","Low",tils_cat)
# Making tils categories for legend
mergeCiberdata <- apply(as.matrix(mergeCiber[,colnames(mergeCiber) %in% columnname2]),2,as.numeric)

# Making the heatmaps
heatmap.2(mergeCiberdata,
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          scale = "column",
          breaks = seq(-3, 3, length.out = 300),
          col = colorpanel(299, "blue", "white", "red"), 
          trace = "none", revC = F,
          key = TRUE,
          keysize = 1,  # Increased key size
          density.info = "none",
          key.title = NA,
          main = "Cell Type  Abundance Estimates",
          key.ylab = NULL,
          key.xlab = NULL,
          labCol  = FALSE,
          margins = c(5, 8),  # Increased right margin
          RowSideColors = tilsVector2,
          ColSideColors = cell_cat2 )  # tilsVector2 SubtypeVector2 PAM50vector2
          #cexRow = NULL, cexCol = 1, srtCol = -5, offsetCol = 0.5,adjCol = c(0,1))

# Note to self could not get legend to fit on heatmap so generate separately and screenshot
# legend("bottom",legend = unique(cell_cat),col = unique(cell_cat2), lty= 1,
#           lwd = 5,
#           cex=1.5)


dev.off()
heatmap_result <- heatmap.2(mergeCiberdata,
          Rowv=T,
          Colv=T,
          dendrogram = "none",
          scale="column",
          breaks = seq(-2, 2, length.out = 300),
          col=colorpanel(299, "blue", "white", "red"), 
          trace="none", revC=F,
          key=TRUE,
          keysize = 1,
          density.info="none",
          key.title = NA,
          key.xlab = NA,
          key.ylab = NA, margins=c(5,5),
          RowSideColors=tilsVector2, #PAM50vector2
          cexRow=NULL, 
          cexCol=.75,
          srtCol = 4,
          labCol = NULL)


dev.off()

row_hclust <- as.hclust(heatmap_result$rowDendrogram)
# k is how many clusters you want to assign
cluster_mem <- cutree(row_hclust, k = 20) 
colors_list <- rainbow(20)

# Creating a clustered heatmap with colored clusters
dend1 <- color_branches(heatmap_result$rowDendrogram,k=20, col = colors_list)
# Creating a cluster assignmnet tables used for subseting the original RNA decon estimates
cluster_data <- data.frame(Row=rownames(mergeCiber),Cluster = cluster_mem)
rownames(cluster_data) <- cluster_data$Row


col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

heatmap_result <- heatmap.2(mergeCiberdata,
                            Rowv=dend1,
                            Colv=F,
                            dendrogram = "row",
                            scale="column",
                            breaks = seq(-2, 2, length.out = 300),
                            col=colorpanel(299, "blue", "white", "red"), 
                            trace="none", revC=F,
                            key=TRUE,
                            keysize = 1,
                            density.info="none",
                            key.title = NA,
                            key.xlab = NA,
                            key.ylab = NA, margins=c(5,5),
                            RowSideColors=SubtypeVector2, #PAM50vector2
                            cexRow=NULL, cexCol=0.75,srtCol = 4,)

# Ensure rownames match between cluster_data and mergeCiberdata
if (!all(rownames(mergeCiber) == rownames(cluster_data))) {
  stop("Row names of mergeCiberdata and cluster_data do not match!")
}

# Subset rows in mergeCiberdata corresponding to cluster 1 in cluster_data
cluster_1_cibersort <- mergeCiber[cluster_data$Cluster== 1, ]
cluster_1_cibersort <- cluster_1_cibersort %>% select(-c(bc_class,Call,tils))
cluster_1_all <- trans_rna[cluster_data$Cluster== 1, ]

## Loops to interate and calculate the proportion of 0s for in each cluster
# cluster_list_cibersort <- c()
# ncluster_list <- c()
# for(i in 1:10){
#   cluster <- mergeCiber[cluster_data$Cluster== i, ]
#   cluster <- cluster %>% select(-c(bc_class,Call,tils))
#   ratio_cluster <- (sum(cluster==0) / (nrow(cluster)*ncol(cluster))) * 100
#   ncluster <- nrow(cluster)
#   cluster_list_cibersort[paste0("Cluster_",i)] <- ratio_cluster
#   ncluster_list[paste0("n_Cluster_",i)] <- ncluster
# }
# 
# cluster_list_all <- c()
# for(i in 1:10){
#   cluster <- trans_rna[cluster_data$Cluster== i, ]
#   ratio_cluster <- (sum(cluster==0) / (nrow(cluster)*ncol(cluster))) * 100
#   cluster_list_all[paste0("Cluster_",i)] <- ratio_cluster
# }


## Loops to interate and calculate the proportion of all 0 samples for in each cluster
cluster_list_cibersort <- c()
ncluster_list <- c()
for(i in 1:10){
  cluster <- mergeCiber[cluster_data$Cluster== i, ]
  cluster <- cluster %>% select(-c(bc_class,Call,tils))
  all_zero_rows <- sum(rowSums(cluster==0) == ncol(cluster))
  ncluster <- nrow(cluster)
  cluster_list_cibersort[paste0("Cluster_",i)] <- all_zero_rows
  ncluster_list[paste0("n_Cluster_",i)] <- ncluster
}

cluster_list_all <- c()
for(i in 1:10){
  cluster <- trans_rna[cluster_data$Cluster== i, ]
  all_zero_rows <- sum(rowSums(cluster==0) == ncol(cluster))
  cluster_list_all[paste0("Cluster_",i)] <- all_zero_rows
}

# Creating a summary table for visualization
summary_0 <- data.frame(cibersort = cluster_list_cibersort, all = cluster_list_all, n = ncluster_list)

