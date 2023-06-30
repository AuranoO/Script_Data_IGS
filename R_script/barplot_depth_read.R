library(tidyverse)
library(ggplot2)
library(dplyr)


coverage <- read.delim("~/Cedratvirus/variant_calling/depth_coverage_consensus_4.samtools", header=FALSE)
repea <- read.delim("~/Cedratvirus/variant_calling/Red_repeat/rpt/Cedratvirus_consensus_r4.rpt", header=FALSE, sep = ":")
repea <- separate(data = repea, col = V2, into = c("start", "end"), sep = "-")

repea$start = as.numeric(repea$start)
repea$end = as.numeric(repea$end)
xmin = c(repea$start)
xmax = c(repea$end)

ggplot(data=coverage, aes(x=V2, y=V3)) +
  geom_bar(stat="identity") +
  annotate(geom = "rect", xmin=xmin, xmax=xmax,
            ymin = -Inf, ymax = Inf, alpha=0.1, fill="red") +
  ylim(0, 5000) +
  theme_minimal()

# lines = c(as.vector(repea$start), as.vector(repea$end))

# ggplot(data=coverage, aes(x=V2, y=V3))+
#   geom_bar(stat="identity")+
#   geom_vline(xintercept = lines, linetype = "dashed", colour = "red")+
#   theme_minimal()
# 
# ggsave(filename = "sequencing_depht_with_repeate_regions.png",
#        device = "png", dpi = "print", units = "px",
#        height =  1200, width = 5000)