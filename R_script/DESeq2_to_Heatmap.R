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
setwd("/home/toumi/Cedratvirus/final_data/DEG/Amibe/")
paths = c("DEG_amibe_1hvsMOCK", "DEG_amibe_1hvs3h", "DEG_amibe_3hvs6h")
names(paths) = c("MOCKvs1h", "1hvs3h", "3hvs6h")


for(i in 1:3) {
  p <- read.csv(paths[i], header =T, sep = "\t")
  p <- p[,c(1,3)]
  p <- p %>% as_tibble() %>% rename(!!names(paths[i]) := 2)
  assign(paste("deg",i, sep="_"),p)
}

degs <- list(deg_1,deg_2, deg_3)
names(degs) = c("MOCKvs1h", "1hvs3h", "3hvs6h")
k <- bind_rows(deg_1, deg_2, deg_3)
# unique amibe = 3780
unique_gene <- unique(k$gene)

#### Venn Diagram ####

# library(ggvenn)
# degs_to_venn <- list(deg_1$gene, deg_2$gene, deg_3$gene)
# names(degs_to_venn) = c("MOCKvs1h", "1hvs3h", "3hvs6h")
# names(degs_to_venn) <- c("1hvs3h", "3hvd6h", "1hvs6h")
# ggvenn(degs_to_venn)

paths = c("DES_amibe_1hvsMOCK", "DES_amibe_1hvs3h", "DES_amibe_3hvs6h")
names(paths) = c("MOCKvs1h", "1hvs3h", "3hvs6h")

for(i in 1:3) {
  p <- read.csv(paths[i], header =T, sep = "\t")
  p <- p[,c(1,3)]
  p <- p %>% as_tibble() %>% rename(!!names(paths[i]) := 2)
  assign(paste("deg",i, sep="_"),p)
}

degs <- list(deg_1,deg_2, deg_3)
names(degs) <- c("MOCKvs1h", "1hvs3h", "3hvs6h")

df_deg <- degs %>% reduce(full_join, by="gene")
df_deg2 <- df_deg[df_deg$gene %in% unique_gene, ]
df_deg2 <- df_deg2 %>% column_to_rownames(var = "gene")

#### Heatmaps ####
colSums_pos_neg <- cbind(colSums((df_deg2*(df_deg2<0))^!is.na(df_deg2))*-1, colSums((df_deg2*(df_deg2>0))^!is.na(df_deg2)))
df_deg2_pam <- pam(df_deg2, 7)

h<-Heatmap(as.matrix(df_deg2),
           cluster_columns = F,
           cluster_rows = T,
           clustering_distance_rows = "pearson",
           clustering_method_rows = "ward.D2",
           show_column_names = T,
           show_row_names = F,
           row_split = df_deg2_pam$clustering,
           row_names_gp = gpar(fontsize = 5),
           name ="log2FC",
           top_annotation = HeatmapAnnotation
           (log2FC = anno_barplot(as.matrix(colSums_pos_neg), 
                            gp = gpar(fill = c("blue","red"), col = NA), beside = T))
           # col=col_fun,
           # heatmap_legend_param = list(at = c(-5,0,5,10,15)),
           )

png("/home/toumi/Cedratvirus/final_data/deseq_graph/heatmap_amine_test_deg.png",
    width = 6, height = 8, units = "in", res = 1200)
draw(h)
dev.off()


#### setup the data ####
tpm <- read.delim("/home/toumi/Cedratvirus/final_data/featureCounts_tpm.csv",
                     sep = "\t", header = T)
tpm_amibe <- tpm %>% filter(str_detect(rownames(tpm), "BAE"))
tpm_amibe <- tpm_amibe[c(4,8,12,1,5,9,2,6,10,3,7,11)]
tpm_amibe <- tpm_amibe[rownames(tpm_amibe) %in% unique_gene,]
tpm_amibe <- log2(tpm_amibe+1)
tpm_amibe <- data.frame(t(scale(t(tpm_amibe))))

ha = rowAnnotation(go = anno_empty(border = FALSE,
                   width = max_text_width(unlist(text_list)) + unit(4, "mm")))

h <- Heatmap(as.matrix(tpm_amibe),
             cluster_columns = F,
             clustering_distance_rows = "spearman",
             clustering_method_rows = "ward.D2",
             row_split = 4,
             show_row_names = F,
             column_title = NULL,
             row_title = NULL,
             name ="z-score",
             use_raster = T,
             column_split = rep(1:4, each = 3),
             left_annotation = rowAnnotation(cluster_block = anno_block(gp = gpar(fill = c(2:4,6)),
                labels = c("1", "2", "3", "4"),
                labels_gp = gpar(col = "white", fontsize = 14),)),
             right_annotation = ha
            )

png("/home/toumi/Cedratvirus/final_data/deseq_graph/heatmap_deg_amibe_tpm_testvisual.png",
    width =6.8, height=8, units="in", res=1200)
h<-draw(h)
dev.off()

for(i in 1:4) {
  decorate_annotation("go", slice = i, {
    grid.text(paste(text_list[[i]]), x = unit(3, "mm"), just = "left")
  })
}

####  topGO with gene cluster ####
paths = c("DES_amibe_1hvsMOCK", "DES_amibe_1hvs3h", "DES_amibe_3hvs6h")
names(paths) = c("MOCKvs1h", "1hvs3h", "3hvs6h")

for(i in 1:3) {
  p <- read.csv(paths[i], header =T, sep = "\t")
  p <- p[,c(1,6)]
  p <- p %>% as_tibble() %>% rename(!!names(paths[i]) := 2)
  assign(paste("des",i, sep="_"),p)
}

des <- list(des_1,des_2, des_3)
df_des <- des %>% reduce(full_join, by="gene")

clusterlist <- row_order(h)
# for clustering classic
for (i in (1:length(clusterlist))) {
  clust <- df_des[df_des$gene %in% rownames(df_deg2)[clusterlist[[i]]],]
  assign(paste("cluster", i, sep="_"),  clust)
}
# 1= 122 2= 1989 3=164 4=1505

# for pam clustering
# for (i in (1:length(unique(tpm_amibe_pam$clustering)))) {
#   clust <- df_deg[df_deg$gene %in% names(tpm_amibe_pam$clustering[tpm_amibe_pam$clustering == i]) ,]
#   assign(paste("cluster", i, sep="_"), clust)
# }

#### topgo graph ####
suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
})

path2map <- "/home/toumi/Cedratvirus/final_data/gene2GO_hic2.map"
geneID2GO <- readMappings(file = path2map)

list_cluster <- list(cluster_1, cluster_2, cluster_3, cluster_4)

allgenes <- df_des[]
list_gene <- apply(df_des[2:4],1, min)
names(list_gene) <- df_des$gene

for (i in (3)) {
  selection <- function(allScore){
    return(allScore < 1e-10)
  }
 list_gene2=list_gene[list_cluster[i][[1]]$gene]
 list_gene2=sort(list_gene2)

   GOdata <- new(
    "topGOdata",
    description="clust_enrich",
    ontology = "MF",
    allGenes = list_gene2,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO,
    nodeSize = 1,
    geneSel = selection,
    ) 
   
  resultsKS = runTest(GOdata, algorithm = "weight01", statistic = "ks")
  goEnrichment <- GenTable(GOdata, KS=resultsKS,
                           orderBy="KS", topNodes=10, numChar=1000)
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<=0.05,]
  # assign(paste("goEnrichmentMF", i, sep = ""), goEnrichment)
}


# text for heatmap annotation
text_list = list(
  text1="\n\nglycolytic process\nphosphatase activity\nphosphopyruvate hydratase activity\ntriose-phosphate isomerase activity",
  text2="acetyl-CoA biosynthetic from pyruvate process\nmalate metabolic process\nbarbed end actin filament capping\nammonium transmembrane transport\nphosphogluconate dehydrogenase\ntransaldolase activity",
  text3="cell redox homeostasis\nperoxiredoxin activity\n4 iron, 4 sulfur cluster binding, zinc ion binding\noxidoreductase activity, acting on NAD(P)H",
  text4="tricarboxylic acid cycle\nglycolytic process\nphosphoglycerate kinase activity\nglycine catabolic process\n3-hydroxyacyl-CoA dehydrogenase activity\n"
)

