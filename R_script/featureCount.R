suppressPackageStartupMessages({
  library(Rsubread)
  library(tidyverse)
  library(ComplexHeatmap)
  library(cluster)
})

setwd("/home/toumi/Cedratvirus/bam/neff_mito_haplo2m2")

#### From BAM do quantification with featureCount ####
fileprefix <- c("A_1h", "A_3h", "A_6h", "A_MOCK", "B_1h", "B_3h", "B_6h", "B_MOCK", "C_1h",
                "C_3h", "C_6h", "C_MOCK", "C_1h30", "C_2h", "C_2h30", "C_4h", "C_8h")

paths = c()
for(i in fileprefix) {
 paths<- c(paths, sprintf("STAR_%s.Aligned.sortedByCoord.out.bam",i))
}

count<-featureCounts(files=paths,
              isPairedEnd = T,
              nthreads = 7,
              annot.ext = "/home/toumi/Cedratvirus/ref_kamtchatka/A_castellanii_HiC-Seq+mito+haplo2m2.gtf",
              isGTFAnnotationFile = TRUE,
              GTF.featureType = "exon")

#### Save the data & normalize count to TPM ####

stat <- data.frame(count$stat)
counts <- data.frame(count$count)

Counts_to_tpm <- function(counts, featureLength) {
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tpmatrix <- data.frame(Counts_to_tpm(count$counts, count$annotation$Length))
names(count2) <- c("A_1h", "A_3h", "A_6h", "A_MOCK", "B_1h", "B_3h", "B_6h", "B_MOCK", "C_1h",
                   "C_3h", "C_6h", "C_MOCK", "C_1h30", "C_2h", "C_2h30", "C_4h", "C_8h")
names(tpmatrix) <- c("A_1h", "A_3h", "A_6h", "A_MOCK", "B_1h", "B_3h", "B_6h", "B_MOCK", "C_1h",
                       "C_3h", "C_6h", "C_MOCK", "C_1h30", "C_2h", "C_2h30", "C_4h", "C_8h")
names(stat) <- c("Status", "A_1h", "A_3h", "A_6h", "A_MOCK", "B_1h", "B_3h", "B_6h", "B_MOCK", "C_1h",
                 "C_3h", "C_6h", "C_MOCK", "C_1h30", "C_2h", "C_2h30", "C_4h", "C_8h")

write.table(stat, "../../final_data/featureCounts_stat.csv",row.names = F, sep = "\t")
write.table(count2, "../../final_data/featureCounts_count.csv", sep = "\t")
write.table(tpmatrix, "../../final_data/featureCounts_tpm.csv", sep = "\t")

#### Read the data from before, and plot the TPM value after filtering ####
tpmatrix <- read.delim("/home/toumi/Cedratvirus/final_data/featureCounts_tpm.csv", sep = "\t", header = T)
counts <- read.delim("/home/toumi/Cedratvirus/final_data/featureCounts_count.csv", sep = "\t", header = T)
stat <- read.delim("/home/toumi/Cedratvirus/final_data/featureCounts_stat.csv", sep = "\t", header = T)


df_ck_fc <- tpmatrix %>% filter(str_detect(rownames(tpmatrix), "CK"))
##### Heatmap #### 
df_ck <- df_ck_fc[c(1,5,9,2,6,10,3,7,11)]
df_ck2 <- df_ck
df_ck2[df_ck2 < 1] <- 0
keep <- rowSums((df_ck2) >= 1) > 3
df_ck <- df_ck2[keep,]
df_ck <- log2(df_ck+1)
df_ck <- data.frame(t(scale(t(df_ck))))
df_ck_pam <- pam(df_ck, 3)
  
h<-Heatmap(as.matrix(df_ck),
           cluster_columns = F,
           clustering_distance_rows = "spearman",
           clustering_method_rows = "ward.D",
           row_split = df_ck_pam$clustering,
           show_row_names = F,
           name ="z-score",
           column_split = rep(1:3, each = 3),
           column_title = NULL, ) 
png("../../final_data/tpm_graph/heatmap_ck_hk.png", width =6.5, height=8, units="in", res=1200)
draw(h)
dev.off()


df_C_ck <- df_ck_fc[c(9,13:15,10,16,11,17)]
df_C_ck2 <- df_C_ck
df_C_ck2[df_C_ck2 < 1] <- 0
keep <- rowSums((df_C_ck2) >= 1) > 3
df_C_ck <- df_C_ck2[keep,]
colSums_tpm <- cbind(colSums(df_C_ck))
df_C_ck <- log2(df_C_ck+1)
df_C_ck <- data.frame(t(scale(t(df_C_ck))))
df_C_ck_pam <- pam(df_C_ck, 3)

h<-Heatmap(as.matrix(df_C_ck),
           cluster_columns = F,
           clustering_distance_rows = "spearman",
           clustering_method_rows = "ward.D",
           row_split = df_C_ck_pam$clustering,
           show_row_names = F,
           name ="z-score", 
           top_annotation = HeatmapAnnotation
           (TPM = anno_barplot(as.matrix(colSums_tpm))
           ))
png("/home/toumi/Cedratvirus/final_data/tpm_graph/heatmap_C_ck_pam.png", width =6.5, height=8.5, units="in", res=1200)
draw(h)
dev.off()


df_amibe_fc <- tpmatrix %>% filter(str_detect(rownames(tpmatrix), "BAE"))
#### dont run  ####
df_amibe <- df_amibe_fc[c(4,8,12,1,5,9,2,6,10,3,7,11)]
df_amibe2 <- df_amibe
df_amibe2[df_amibe2 < 1] <- 0
keep <- rowSums((df_amibe2) >= 1) > 4
df_amibe <- df_amibe2[keep,]
df_amibe <- log2(df_amibe+1)
df_amibe <- data.frame(t(scale(t(df_amibe))))

h<-Heatmap(as.matrix(df_amibe),
           cluster_columns = F,
           clustering_distance_rows = "pearson",
           clustering_method_rows = "ward.D",
           row_km = 5,
           show_row_names = F,
           row_title = NULL,
           name ="z-score",
           column_split = rep(1:4, each = 3),
           column_title = NULL,
           use_raster = T) 
png("../../final_data/tpm_graph/heatmap_amibe_hk.png", width =7.5, height=8, units="in", res=1200)
draw(h)
dev.off()

df_C_amibe <- df_amibe_fc[c(12,9,13:15,10,16,11,17)]
df_C_amibe <- df_C_amibe
df_C_amibe[df_C_amibe < 1] <- 0
keep <- rowSums((df_C_amibe) >= 1) > 4
df_C_amibe <- df_C_amibe[keep,]
colSums_tpm <- cbind(colSums(df_C_amibe))
df_C_amibe <- log2(df_C_amibe+1)
df_C_amibe <- data.frame(t(scale(t(df_C_amibe))))
df_C_amibe_pam <- pam(df_C_amibe, 5)

h<-Heatmap(as.matrix(df_C_amibe),
           cluster_columns = F,
           clustering_distance_rows = "spearman",
           clustering_method_rows = "ward.D",
           row_split = df_C_amibe_pam$clustering,
           show_row_names = F,
           name ="z-score",
           use_raster = T,
           top_annotation = HeatmapAnnotation
           (TPM = anno_barplot(as.matrix(colSums_tpm)))
           )

png("/home/toumi/Cedratvirus/final_data/tpm_graph/heatmap_C_amibe_pam.png", width =7.5, height=8, units="in", res=1200)
draw(h)
dev.off()


#### ACP ####
condition <- c("1h","3h","6h","1h","3h","6h","1h","3h","6h", "1h30", "2h", "2h30", "4h", "8h")
df_ck_acp <- df_ck_fc[c(1:3,5:7,9:11,13:17)]
df_pca <- as.data.frame(t(df_ck_acp))
df_pca <-data.frame(df_pca, condition)
res.pca <- PCA(df_pca[1:562])
acp1 <- fviz_pca_ind(res.pca,
                     col.ind = df_pca$condition,
                     addEllipses = F,
                     legend.title = "Time",
                     show.legend = F)+ theme_minimal()+
  scale_shape_manual(values=c(19,19,19,19,19,19,19,19))+ labs(title = "PCA - C. kamchatka")
png("/home/toumi/Cedratvirus/final_data/acp/acp_allck_new2.png", units = "in", width = 7, height = 6, res = 800)
print(acp1)
dev.off()


condition <- c("1h","3h","6h",'MOCK',"1h","3h","6h",'MOCK',"1h","3h","6h",'MOCK', "1h30", "2h", "2h30", "4h", "8h")
df_amibe_acp <- df_amibe_fc
df_pca <- as.data.frame(t(df_amibe_acp))
df_pca <-data.frame(df_pca, condition)
res.pca <- PCA(df_pca[1:15487])
acp <- fviz_pca_ind(res.pca,
                     col.ind = df_pca$condition,
                     addEllipses = F,
                     legend.title = "Time",
                     show.legend = F)+ theme_minimal()+
  scale_shape_manual(values=c(19,19,19,19,19,19,19,19,19))+ labs(title = "PCA - A. castellanii")
png("/home/toumi/Cedratvirus/final_data/acp/acp_allamibe_new.png", units = "in", width = 8, height = 6, res = 800)
print(acp)
dev.off()
