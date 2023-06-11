suppressPackageStartupMessages({
  library(readr)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(tidyverse)
})

# setup the data ####
setwd("/media/aurano/GPARTED-LIV/final_data/final_data/DEG/virus")
counts <- read.delim("/home/toumi/Cedratvirus/final_data/featureCounts_count.csv",
                     sep = "\t", header = T)

setwd("/media/aurano/GPARTED-LIV/final_data/final_data/DEG/virus")
counts <- read.delim("../../featureCounts_count.csv",
                     sep = "\t", header = T)
# counts_ck <- counts %>% filter(str_detect(rownames(counts), ""))
counts_ck <- counts
counts_ck <- counts_ck[c(1,5,9,2,6,10,3,7,11)]
coldata <- data.frame(factor(rep(c("1h", "3h", "6h"), each=3)))
names(coldata) <- 'condition'
rownames(coldata) <- colnames(counts_ck)

# create the deseq2 object ####
dds <- DESeqDataSetFromMatrix(countData = counts_ck,
                              colData = coldata,
                              design = ~ condition)

# pre-filtering
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# dds$condition <- relevel(dds$condition, ref = "1h")

ddslist <- list(c("6h", "1h", "1hvs6h"),
                c("3h", "1h", "1hvs3h"),
                c("6h", "3h", "3hvs6h"))

pval = 1e-5
fc = 1

dds <- DESeq(dds)

for (data in ddslist) {
  res <- results(dds, contrast=c("condition", data[1], data[2]))
  resLFC <- lfcShrink(dds,res=res, type ="ashr")
  up <- resLFC[which(resLFC$log2FoldChange > fc & resLFC$padj < pval),]
  down <- resLFC[which(resLFC$log2FoldChange < -fc & resLFC$padj < pval),]
  up_dt <- data.frame(up)
  up_dt <- rownames_to_column(up_dt, var ="gene")
  down_dt <- data.frame(down)
  down_dt <- rownames_to_column(down_dt, var ="gene")
  assign(paste("up", data[3], sep="_"), up_dt)
  assign(paste("down", data[3], sep="_"), down_dt)
  resLFCorder <- resLFC[order(resLFC$padj),]
  assign(paste("resLFC", data[3], sep=""), resLFCorder)
  DES <- data.frame(resLFCorder)
  DES <- rownames_to_column(DES, var ="gene")
  DES <- DES %>% filter(str_detect(gene, "CK"))
  assign(paste("DES", data[3], sep=""),DES)
  DEG <- data.frame(rbind(as.data.frame(up), as.data.frame(down)))
  DEG <- rownames_to_column(DEG, var ="gene")
  DEG <- DEG %>% filter(str_detect(gene, "CK"))
  assign(paste("DEG", data[3], sep=""),DEG)
}


#### Volcano Plot ####
listDES<-list( DES1hvs3h, DES3hvs6h)
titleName<-c("1h vs 3h", "3h vs 6h")

volcanoplot <- function(results, titl, pval, fc){
  # Create columns from the column numbers specified
  dim(results)
  results <- results %>% mutate(fdr = .[[6]],
                                pvalue = .[[5]],
                                logfc = .[[3]],
                                labels = .[[1]])
  
  # Get names for legend
  down <- unlist(strsplit('Under,Not Sig,Over', split = ","))[1]
  notsig <- unlist(strsplit('Under,Not Sig,Over', split = ","))[2]
  up <- unlist(strsplit('Under,Not Sig,Over', split = ","))[3]
  
  # Set colours
  colours <- setNames(c("cornflowerblue", "grey", "firebrick"), c(down, notsig, up))
  
  # Create significant (sig) column
  results <- mutate(results, sig = case_when(
    fdr < pval & logfc > fc ~ up,
    fdr < pval & logfc < -fc ~ down,
    TRUE ~ notsig))
  print(table(results$sig))
  
  # Create plot -------------------------------------------------------------
  
  # Set up base plot
  p <- ggplot(data = results, aes(x = logfc, y = -log10(pvalue))) +
    geom_point(aes(colour = sig)) +
    scale_color_manual(values = colours) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key = element_blank(),
          legend.text = element_text(size=11),
          legend.title = element_blank()) +
    geom_vline(xintercept = c(-fc,fc),
               color="goldenrod2",
               linetype ="dashed") +
    labs(title = titl)
  
  #legend.key.size = unit(2,"cm"),
  #legend.key.width= unit(4, 'cm'),
  #legend.key.height = unit(3, 'cm'))
  
  p <- p + guides(color = guide_legend(override.aes = list(size = 6)))
  # # Set legend title
  # p <- p + theme()
  
  # Print plot
  return(p)
}

for (i in 1:2) {
  k <- data.frame(listDES[i])
  assign(paste("volcano", i, sep = ""), volcanoplot(k, titleName[i], pval,fc)) 
}

arrange<-ggarrange(volcano1,volcano2,
                   ncol=2,
                   nrow=1)

png("deg_ck_only_volcano.png",
    units="px", width=4300,
    height=1500, res=500)

arrange
dev.off()


listDEG <- list(DEG1hvs3h, DEG1hvs6h, DEG3hvs6h)
names(listDEG) <- c("1hvs3h", "1hvs6h", "3hvs6h")
listup <- list(up_1hvs3h, up_1hvs6h, up_3hvs6h)
names(listup) <- c("1hvs3h", "1hvs6h", "3hvs6h")
listdown <- list(down_1hvs3h,down_1hvs6h, down_3hvs6h)
names(listdown) <- c("1hvs3h", "1hvs6h", "3hvs6h")

for (i in 1:3) {
  k <- data.frame(listDEG[i])
  names(k)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  write_tsv(k, paste("/home/toumi/Cedratvirus/final_data/DEG/virus/DEG_ck",
                     names(listup[i]), sep="_"))
  
  top10 <- k %>% slice_max(k$log2FoldChange, n = 10)
  down10 <- k %>% slice_min(k$log2FoldChange, n = 10)
  extreme <- bind_rows(top10, down10)
  assign(paste("extreme_ck", names(listDEG[i]), sep="_"), extreme)
  
  k <- data.frame(listup[i])
  names(k)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  f <- data.frame(listdown[i])
  names(f)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  
  write_tsv(k, paste("/home/toumi/Cedratvirus/final_data/DEG/virus/DEG_up_ck",
                     names(listup[i]), sep="_"))
  write_tsv(f, paste("/home/toumi/Cedratvirus/final_data/DEG/virus/DEG_down_ck",
                     names(listdown[i]), sep="_"))
  
  k<-data.frame(listDES[i])
  write_tsv(k, paste("/home/toumi/Cedratvirus/final_data/DEG/virus/DES_ck",
                     names(listDEG[i]), sep="_"))
}
