setwd(".")

library("reshape2")
library(ggplot2)
library(optparse)

option_list <- list(
	    make_option(c("--control"), default="must_specify",
	    help="this specifiys the control output from annotatePeaksBed"),
	    make_option(c("--treatment"), default="must_specify",
	    help="this specifiys the treatment output from annotatePeaksBed"),
	    make_option(c("--input"), default="must_specify",
	    help="this specifiys the input output from annotatePeaksBed")
)

opt <- parse_args(OptionParser(option_list=option_list))


print("Running with following options:")
print(opt)


control <- read.table(opt$control,header = TRUE, fill = T)
treatment <- read.table(opt$treatment,header = TRUE, fill=T)
input <- read.table(opt$input,header = TRUE, fill=T)



control_data <- control[1:length(control$Distance),c(1,2)]
colnames(control_data) <- c("Distance","Control_Cov")

treatment_data <- treatment[1:length(treatment$Distance),c(1,2)]
colnames(treatment_data) <- c("Distance","Treatment_Cov")

inputdata <- input[1:length(input$Distance),c(1,2)]
colnames(inputdata) <- c("Distance","input")

donor1 <- merge(control_data, treatment_data, by = "Distance")
donor1 <- merge(inputdata, donor1, by = "Distance")
donor1 <- melt(donor1, id="Distance")

name = gsub(".txt", "", opt$control)
filename <- paste0(name, ".eps")

setEPS()
postscript(filename)

ggplot(donor1, aes(Distance, value, color=variable)) +
  scale_color_manual(values=c("#f44242","#BABABA", "#349ADE"))+geom_line(size=2) + theme_bw()+
  theme(axis.line = element_line(color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())


dev.off()




