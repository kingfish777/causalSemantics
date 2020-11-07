library(sqldf)

txt = "/Users/scottalexandermalec/Projects/enlil/fakeConfounders.txt"

text = read.csv.sql(file = txt, header = FALSE, sep = ":")

dat <- read.csv("/Users/scottalexandermalec/Projects/enlil/fakeConfounders.txt", sep = ":")


