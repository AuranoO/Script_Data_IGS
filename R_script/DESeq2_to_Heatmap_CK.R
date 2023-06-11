suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggvenn)
  library(purrr)
  library(tidyverse)
  library(cluster)
})

#### Setup ####
setwd("/media/aurano/GPARTED-LIV/final_data/final_data/DEG/virus")
paths = c("DEG_ck_1hvs3h", "DEG_ck_3hvs6h")
names(paths) = c("1hvs3h", "3hvs6h")

for(i in 1:2) {
  p <- read.csv(paths[i], header =T, sep = "\t")
  p <- p[,c(1,3)]
  p <- p %>% as_tibble() %>% rename(!!names(paths[i]) := 2)
  assign(paste("deg",i, sep="_"),p)
}

degs <- list(deg_1,deg_2)
names(degs) <- c("1hvs3h", "3hvs6h")
k <- bind_rows(deg_1, deg_2)

unique_gene <- unique(k$gene)

library(ggvenn)
degs_to_venn <- list(deg_1$gene, deg_2$gene)
names(degs_to_venn) <- c("1hvs3h", "3hvs6h")
ggvenn(degs_to_venn)

ggsave("deg_not_all_ck.png", device = "png", dpi = "print", width = 5, height = 6)
paths = c("DES_ck_1hvs3h", "DES_ck_3hvs6h")
names(paths) = c("1hvs3h", "3hvs6h")

for(i in 1:2) {
  p <- read.csv(paths[i], header =T, sep = "\t")
  p <- p[,c(1,3)]
  p <- p %>% as_tibble() %>% rename(!!names(paths[i]) := 2)
  assign(paste("des",i, sep="_"),p)
}
list_des<-list(des_1, des_2)

df_des <- list_des %>% reduce(full_join, by="gene")
df_deg2 <- df_des[df_des$gene %in% unique_gene, ]
df_deg2 <- df_deg2 %>% column_to_rownames(var = "gene")

colSums_pos_neg <- cbind(colSums((df_deg2*(df_deg2<0))^!is.na(df_deg2))*-1, colSums((df_deg2*(df_deg2>0))^!is.na(df_deg2)))

h<-Heatmap(as.matrix(df_deg2),
           cluster_columns = F,
           cluster_rows = F,
           clustering_distance_rows = "spearman",
           clustering_method_rows = "ward.D2",
           show_column_names = T,
           show_row_names = F,
           # row_split = 2,
           row_names_gp = gpar(fontsize = 5),
           name ="log2FC",
           top_annotation = HeatmapAnnotation
           (log2FC = anno_barplot(as.matrix(colSums_pos_neg), 
                                  gp = gpar(fill = c("blue","red"), col = NA), beside = T))
           # col=col_fun,
           # heatmap_legend_param = list(at = c(-5,0,5,10,15)),
)

png("/heatmap_ck_test_deg.png",
    width = 6, height = 8, units = "in", res = 1200)
draw(h)
dev.off()

tpm <- read.delim("../../featureCounts_tpm.csv",
                  sep = "\t", header = T)
tpm_ck <- tpm %>% filter(str_detect(rownames(tpm), "CK"))
# tpm_ck <- tpm_ck[c(4,8,12,1,5,9,2,6,10,3,7,11)]
tpm_ck <- tpm_ck[c(1,5,9,2,6,10,3,7,11)]
tpm_ck <- tpm_ck[rownames(tpm_ck) %in% unique_gene,]

tpm_ck_pam <- pam(tpm_ck, 3)
tpm_ck <- log2(tpm_ck+1)
tpm_ck <- data.frame(t(scale(t(tpm_ck))))

h <- Heatmap(as.matrix(tpm_ck),
             cluster_columns = F,
             clustering_distance_rows = "spearman",
             clustering_method_rows = "ward.D2",
             row_split = 3,
             show_row_names = F,
             column_title = NULL,
             row_title = NULL,
             name ="z-score",
             use_raster = T,
             column_split = rep(1:3, each = 3),
             left_annotation = rowAnnotation(cluster_block = anno_block(gp = gpar(fill = c(2:4,6)),
                                                                        labels = c("1", "2", "3"),
                                                                        labels_gp = gpar(col = "white", fontsize = 14),)),
             # right_annotation = ha
)

png("/home/toumi/Cedratvirus/final_data/deseq_graph/heatmap_deg_ck_tpm_testvisual.png",
    width =6.8, height=8, units="in", res=1200)
h<-draw(h)
dev.off()

clusterlist <- row_order(h)

# for clustering classic
for (i in (1:length(clusterlist))) {
  clust <- df_deg[df_deg$gene %in% rownames(df_deg2)[clusterlist[[i]]],]
  assign(paste("cluster", i, sep="_"),  clust)
}
non_deg <- tpm_ck
tpm_ck$type <- "DEG"
deg <- tpm_ck
data_violin = bind_rows(deg, non_deg)
library(reshape2)
reshaped <- melt(data_violin, direction = "long")
reshaped <- reshaped %>% mutate(variable=factor(variable, levels = c("A_1h", "A_3h", "A_6h", "B_1h", "B_3h", "B_6h", "C_1h",
                                                                       "C_3h", "C_6h")))
p<- ggplot(reshaped, aes(x=variable, y=value, color = type))+ geom_violin()+
  xlab("TPM") + ylab("Sample")
p

ggsave(filename="violin_plot", dpi = "print", device = "png", width = 2000, height = 1000, units = "px") 

# cluster 1=78 2=39 3=301  