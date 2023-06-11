suppressPackageStartupMessages({
  library(readr)
  library(tximeta)
  library(DESeq2)
  library(SummarizedExperiment)
  library(impute)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(tidyverse)
})

#### Setup & DESeq2 #####
setwd("/home/toumi/Cedratvirus/quant_hic_haplom1m2/")
setwd("/home/toumi/Cedratvirus/bam/neff_mito_haplo2m2")

a <- c("A_1h", "B_1h", "C_1h", "A_3h", "B_3h", "C_3h", "A_6h",
       "B_6h", "C_6h", "A_MOCK", "B_MOCK", "C_MOCK")

samples <- data.frame(a)
samples$condition <- factor(rep(c("1h", "3h", "6h", "MOCK"), each=3))
names(samples) <- c("run","condition")

path <- c()
for(i in a) {
  p <- sprintf("%s.genes.results", i)
  path <- c(path, p)
} 
names(path) <- samples$run
samples$files <- path
names(samples) <- c("names", "condition", "files")
 
se <- tximeta(samples, type="rsem", txIn=F, txOut=F, skipMeta=T)
# line with only 0 are remove
assays(se)$length[ assays(se)$length == 0] <- NA
idx <- rowSums(is.na(assays(se)$length)) >= 6
table(idx)
se <- se[!idx,]
length_imp <- impute.knn(assays(se)$length)
assays(se)$length <- length_imp$data

ddsSE <- DESeqDataSet(se, design = ~condition)
keep <- rowSums(counts(ddsSE) >= 10) >= 3
ddsSE <- ddsSE[keep,]
ddsSE$condition <- factor(ddsSE$condition, levels = c("3h","1h","6h",'MOCK'))
dds <- DESeq(ddsSE)
# p2<-rownames(ddsSE[grep("gene", rownames(ddsSE))])

#### Visualization with PCA ####
vsd <- vst(dds, blind=F)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))

# png("/home/toumi/Cedratvirus/graph/DESeq2/PCA.png",
#     units="px", width=3000,
#     height=3000, res=500)
ggplot(pcaData, aes(PC1, PC2, col=condition)) +
  geom_point(size = 4) +
  xlab(paste("PC1: ",percentVar[1], "% variance")) +
  ylab(paste("PC2: ",percentVar[2], "% variance"))
# dev.off()

#### Set each comparison  & create corresponding dataframe #####
ddslist <- list(c("1h", "MOCK", "1hvsMOCK"), c("3h", "MOCK", "3hvsMOCK"),
                c("6h", "MOCK", "6hvsMOCK"), c("3h", "1h", "1hvs3h"),
                c("6h", "1h", "1hvs6h"), c("6h", "3h", "3hvs6h"))

pval = 0.01
fc = 2

for (data in ddslist) {
  res <- results(dds, contrast=c("condition", data[1], data[2]), alpha = 0.01)
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
  DES$gene <- str_replace(DES$gene, "gene-", "ck")
  assign(paste("DES", data[3], sep=""),DES)
  DEG <- data.frame(rbind(as.data.frame(up), as.data.frame(down)))
  DEG <- rownames_to_column(DEG, var ="gene")
  DEG$gene <- str_replace(DEG$gene, "gene-", "ck")
  assign(paste("DEG", data[3], sep=""), DEG)
}

#### Volcano Plot ####
listDES<-list(DES1hvsMOCK, DES3hvsMOCK, DES6hvsMOCK,
              DES1hvs3h, DES1hvs6h, DES3hvs6h)
titleName<-c("1h vs MOCK", "3h vs MOCK", "6h vs MOCK",
             "1h vs 3h", "1h vs 6h", "3h vs 6h")

volcanoplot <- function(results, titl, pval, fc){
  # Create columns from the column numbers specified
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

for (i in 1:6) {
  k <- data.frame(listDES[i])
  assign(paste("volcano", i, sep = ""), volcanoplot(k, titleName[i], pval,fc)) 
  }

arrange<-ggarrange(volcano1,volcano2,volcano3,volcano4,volcano5,volcano6, 
                   ncol=3,
                   nrow=2)

png("/home/toumi/Cedratvirus/graph/DESeq2/multi_volcano.png",
    units="px", width=7000,
    height=5000, res=500)

arrange
dev.off()

#### Create a dataframe count of DEG/DES/up/down for each organisms ####
listDEG <- list(DEG1hvsMOCK, DEG3hvsMOCK, DEG6hvsMOCK,
                DEG1hvs3h, DEG1hvs6h, DEG3hvs6h)
names(listDEG) <- c("1hvsMOCK", "3hvsMOCK", "6hvsMOCK",
                    "1hvs3h", "1hvs6h", "3hvs6h")
listup <- list(up_1hvsMOCK, up_3hvsMOCK, up_6hvsMOCK,
               up_1hvs3h, up_1hvs6h, up_3hvs6h)
names(listup) <- c("1hvsMOCK", "3hvsMOCK", "6hvsMOCK",
                    "1hvs3h", "1hvs6h", "3hvs6h")
listdown <- list(down_1hvsMOCK, down_3hvsMOCK, down_6hvsMOCK,
                 down_1hvs3h, down_1hvs6h, down_3hvs6h)
names(listdown) <- c("1hvsMOCK", "3hvsMOCK", "6hvsMOCK",
                    "1hvs3h", "1hvs6h", "3hvs6h")

d <- tibble(Versus = character(), Amibe = numeric(),
            Virus = numeric(), Total = numeric())
up <- tibble(Versus = character(), Amibe = numeric(),
            Virus = numeric(), Total = numeric())
down <- tibble(Versus = character(), Amibe = numeric(),
              Virus = numeric(), Total = numeric())

for (i in 1:6) {
  k <- data.frame(listDEG[i])
  names(k)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  k_ck <- k %>% filter(str_detect(gene, "CKst1"))
  k_amibe <- k %>% filter(str_detect(gene, "BAESF"))
  d <- d %>% add_row(Versus = names(listDEG[i]), Amibe = nrow(k_amibe),
                     Virus = nrow(k_ck), Total = nrow(k))
  write_tsv(k_ck, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_ck",
                        names(listDEG[i]), sep="_"))
  write_tsv(k_amibe, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_amibe",
                           names(listDEG[i]), sep="_"))
  top10 <- k_amibe %>% slice_max(k_amibe$log2FoldChange, n = 10)
  down10 <- k_amibe %>% slice_min(k_amibe$log2FoldChange, n = 10)
  extreme <- bind_rows(top10, down10)
  assign(paste("extreme_amibe", names(listDEG[i]), sep="_"), extreme)
  top10 <- k_ck %>% slice_max(k_ck$log2FoldChange, n = 10)
  down10 <- k_ck %>% slice_min(k_ck$log2FoldChange, n = 10)
  extreme <- bind_rows(top10, down10)
  assign(paste("extreme_ck", names(listDEG[i]), sep="_"), extreme)
  k <- data.frame(listup[i])
  names(k)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  f <- data.frame(listdown[i])
  names(f)<-c("gene","baseMean","log2FoldChange","lfcSE","pvalue", "padj")
  up_ck <- k %>% filter(str_detect(gene, "CKst1"))
  up_amibe <- k %>% filter(str_detect(gene, "BAESF"))
  down_ck <- f %>% filter(str_detect(gene, "CKst1"))
  down_amibe <- f %>% filter(str_detect(gene, "BAESF"))
  up <- up %>% add_row(Versus = names(listup[i]), Amibe = nrow(up_amibe),
                       Virus = nrow(up_ck), Total = nrow(k))
  down <- down %>% add_row(Versus = names(listdown[i]), Amibe = nrow(down_amibe),
                           Virus = nrow(down_ck), Total = nrow(f))
  write_tsv(up_ck, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_up_ck",
                         names(listup[i]), sep="_"))
  write_tsv(up_amibe, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_up_amibe",
                         names(listup[i]), sep="_"))
  write_tsv(down_ck, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_down_ck",
                         names(listdown[i]), sep="_"))
  write_tsv(down_ck, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DEG_down_amibe",
                         names(listdown[i]), sep="_"))
  k<-data.frame(listDES[i])
  k_ck <- k %>% filter(str_detect(gene, "CKst1"))
  k_amibe <- k %>% filter(str_detect(gene, "BAESF"))
  write_tsv(k_ck, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DES_ck",
                        names(listDEG[i]), sep="_"))
  write_tsv(k_amibe, paste("/home/toumi/Cedratvirus/DEG_neff_HiC+haplo1m2/DES_amibe",
                           names(listDEG[i]), sep="_"))
}



#### get extreme gene in clipboard ####
# library(clipr)
# write_clip(extreme_amibe_3hvs6h$gene)


