library(ggplot2)

tools=rep(c("RSEM", "Salmon", "featureCounts"), each=3)
sample=rep(c("A_3h","B_3h","C_3h"), 3)
value=c(77.93, 79.31,76.45,
       67.4, 69.66, 69.01,
       87.05, 88.68, 87.98)
data = data.frame(tools, sample, value)

ggplot(data, aes(x=sample, y=value, fill=tools)) + 
    geom_bar(position="dodge", stat="identity")+
    ylab("Alignment rate")+
    xlab("Samples")+
    theme_minimal()

ggsave("Histogram.png",units = "px",
       width = 3000, height = 2500,  dpi = 800)

