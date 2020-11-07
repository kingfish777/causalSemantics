


library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(compiler)

SQUELCH <- 999
# exp <- "naltrexone"; hoi <- "gastrointestinal_hemorrhage"; confounderSource <- "psi"

getCleanConfounders <- function(exp, hoi, confounderSource, dummy1, dummy2, dummy3) {
  try({
    print(paste(hoi, exp, sep = ":"))
  if (confounderSource == "sql") {
    confounders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM smdbprepreds 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'TREATS'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM smdbprepreds 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    confounders <- str_to_lower(subset(confounders, theta1*theta2 > 1)[,1])
  } else if (confounderSource == "psi") {
    print("GET COMORBIDITIES")
    cmd <- paste("./confounders.sh ", exp, " ", hoi, " 50 > ccs.txt", sep = "")
    system(cmd)
    cleanupResultsFile_compiled("ccs.txt")
    confounders <- as.data.frame(read.csv2(file = "ccs.txt", stringsAsFactors = FALSE), stringAsFactors = FALSE)
    print(confounders)
    confounders <- confounders[[1]]
  }
  print("RAW CONFOUNDERS")
  print(confounders)
  # 
  # colliders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
  #                       FROM smdbprepreds 
  #                       WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
  #                       AND predicate = 'CAUSES'
  #                       GROUP BY object_name), 
  #                       ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
  #                       FROM smdbprepreds 
  #                       WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
  #                       AND predicate IN ('CAUSES')
  #                       GROUP BY subject_name) 
  #                       SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
  # 
  # colliders <- str_to_lower(subset(colliders, theta1*theta2 > 1)[,1])
  # 
  # 
  # try({expSynonyms <-
  #   c(read.csv2(
  #     file = paste("variables/syns/synonyms/", exp, ".txt", sep = ""),
  #     stringsAsFactors = FALSE,
  #     header = FALSE
  #   )[1])})
  # try({hoiSynonyms <-
  #   c(read.csv2(
  #     file = paste("variables/syns/synonyms/", hoi, ".txt", sep = ""),
  #     stringsAsFactors = FALSE,
  #     header = FALSE
  #   )[1])})
  # try({stopwords <-
  #   #as.character(
  #   c(read.csv2(
  #     file = "variables/stoplist.txt",
  #     #   sep = "\n",
  #     stringsAsFactors = FALSE,
  #     header = FALSE
  #   )[1])})
  # #####
  # try({confounders <- setdiff(confounders, stopwords[[1]])
  # print("#### minus stopwords")
  # print(length(confounders))})
  # #
  # try({confounders <- setdiff(confounders, expSynonyms[[1]])
  # print("#### minus expSynonyms")
  # print(length(confounders))})
  # #
  # try({confounders <- setdiff(confounders, hoiSynonyms[[1]])
  # print("#### minus hoiSynonyms")
  # print(length(confounders))})
  # ##################
  # print(confounders)
  # 
  # confounders.v1 <- c()
  # for (confounder in confounders) {
  #   print(confounder)
  #   try({fName <- paste("utgz/", confounder, ".txt.gz", sep = "")})
  #   try({if (file.exists(fName)) { confounders.v1 <- c(confounders.v1, confounder); print(confounders.v1) }
  #     else if (!file.exists(fName)) { print(paste(confounder, "WAS NOT MEASURED OR RECORDED"))}})
  # }
  # confounders.v1
  confounders
  if (is.null(confounders)) { print("no confounders"); confounders <- "" }
  #system("mkdir sqlraw")
  cName = paste("sqlraw/", hoi, "==", exp, ".txt", sep = "")
  print(cName)
  writeVec(dat = confounders, fn = cName)
  })
}

gcc <- cmpfun(getCleanConfounders)

source("gcc.R")
 #####  exp <- "naltrexone"; hoi <- "gastrointestinal_hemorrhage"; confounderSource <- "sql"
#
######gcc(exp, hoi, confounderSource, '', '', '')
#  
#confounders <- gcc(exp, hoi, confounderSource, '', '', '')
#print(confounders)
# bmc('sql', 'zidovudine', 'liver_failure,_acute', 'liver_failure,_acute', '1', '10')


cons <- dbListConnections(pg)
for (con in cons) {
  print(con)
  dbDisconnect(con)
}
