
library(ggpubr)
#setwd("/Users/scottalexandermalec/Projects/enlil")
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots
confounderLists <- list.files(path = "sqlrawOriginal/")
confounderStats <- c()
for (cl in confounderLists) {
  try({confounderStat <- length(unlist(read.csv2(paste("sqlrawOriginal/", cl, sep = ""))))
  #print(paste(cl, confounderStat, sep = ""))
  #if (confounderStat < 5) {
    print(cl)
    print(length(confounderStats))
    confounderStats <- c(confounderStats, confounderStat)
  #}
  })
}

median(confounderStats)
mean(confounderStats)
min(confounderStats)
max(confounderStats)
postscript(file = "boxplots/confounderHistogram.eps", title = "confounders per drug / adverse event pair", width = 480, height = 480)
hist(confounderStats, breaks = 200)
dev.off()
confounderStats



