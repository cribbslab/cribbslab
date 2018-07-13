setwd(".")

library("reshape2")

DMSO1 <- read.table("DMSO-H3K27ac-ChIP-1_S22_total.txt",header = TRUE, fill = T)
J41 <- read.table("CBP30-H3K27ac-ChIP-2_S5_total.txt",header = TRUE, fill=T)
input <- read.table("CBP30-input-1_S25_total.txt",header = TRUE, fill=T)
library(ggplot2)
DMSO1data <- DMSO1[1:length(DMSO1$Distance),c(1,2)]
colnames(DMSO1data) <- c("Distance","DMSO_Cov")

J41data <- J41[1:length(J41$Distance),c(1,2)]
colnames(J41data) <- c("Distance","CBP30_Cov")

inputdata <- input[1:length(input$Distance),c(1,2)]
colnames(inputdata) <- c("Distance","input")

donor1 <- merge(DMSO1data, J41data, by = "Distance")
donor1 <- merge(inputdata, donor1, by = "Distance")



donor1 <- melt(donor1, id="Distance")
setEPS()
postscript("Coverage_H3K27ac_1_ATAC_peaks.eps")

ggplot(donor1, aes(Distance, value, color=variable)) +
  scale_color_manual(values=c("#f44242","#BABABA", "#349ADE"))+geom_line(size=2) + theme_bw()+

  
  theme(axis.line = element_line(color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())


dev.off()




