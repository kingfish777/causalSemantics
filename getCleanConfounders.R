

library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)

#SQUELCH <- 10
 # exp <- "naltrexone"; hoi <- "gastrointestinal_hemorrhage"; confounderSource <- "psi"

getCleanConfounders <- function(exp, hoi, confounderSource) {

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
  
  colliders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM smdbprepreds 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'CAUSES'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM smdbprepreds 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
  
  colliders <- str_to_lower(subset(colliders, theta1*theta2 > 1)[,1])
  
  
  try({expSynonyms <-
   c(read.csv2(
      file = paste("variables/syns/synonyms/", exp, ".txt", sep = ""),
      stringsAsFactors = FALSE,
      header = FALSE
    )[1])})
  try({hoiSynonyms <-
    c(read.csv2(
      file = paste("variables/syns/synonyms/", hoi, ".txt", sep = ""),
      stringsAsFactors = FALSE,
      header = FALSE
    )[1])})
  try({stopwords <-
    #as.character(
   c(read.csv2(
      file = "variables/stoplist.txt",
      #   sep = "\n",
      stringsAsFactors = FALSE,
      header = FALSE
    )[1])})
   #####
   try({confounders <- setdiff(confounders, stopwords[[1]])
   print("#### minus stopwords")
   print(length(confounders))})
   #
   try({confounders <- setdiff(confounders, expSynonyms[[1]])
   print("#### minus expSynonyms")
   print(length(confounders))})
   #
   try({confounders <- setdiff(confounders, hoiSynonyms[[1]])
   print("#### minus hoiSynonyms")
   print(length(confounders))})
    ##################
   print(confounders)
  
   confounders.v1 <- c()
   for (confounder in confounders) {
     print(confounder)
     try({fName <- paste("utgz/", confounder, ".txt.gz", sep = "")})
     try({if (file.exists(fName)) { confounders.v1 <- c(confounders.v1, confounder); print(confounders.v1) }
     else if (!file.exists(fName)) { print(paste(confounder, "WAS NOT MEASURED OR RECORDED"))}})
   }
   confounders.v1
}

gcc <- cmpfun(getCleanConfounders)


getTrueConfounders <- function(confounders, adeData, drugData, SQUELCH) {
  trueConfounders.local <- c()
  for (cc in confounders) {   # cc = "confounder candidate"
    if (length(trueConfounders.local) >= SQUELCH)
    { break }
    print("DOOKIE")
    if (file.exists(paste("utgz/", cc, ".txt.gz", sep = ""))) {
      cc <- as.character(cc)
      print(cc)
      ccData <-
        read.csv2(file = paste("utgz/", cc, ".txt.gz", sep = ""),
                  stringsAsFactors = FALSE)
      cc <-
        str_replace_all(string = cc,
                        pattern = "-",
                        replacement = "HYPHEN")
      cc <- 
        str_replace_all(string = cc, 
                        pattern = ",",
                        replacement = "COMMA")
      dat <- data.frame(adeData, drugData, ccData)
      colnames(dat) <- c("ade", "medication", cc)
      colnames(dat) <-
        str_replace_all(string = colnames(dat),
                        pattern = "-",
                        replacement = "-")
      colnames(dat) <- 
        str_replace_all(string = colnames(dat),
                        pattern = ",",
                        replacement = ",")
      write.table(
        dat,
        file = "dat.txt",
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
      dat.names <- cc
      #print(dat.names)
      buildDAG <- function(dat.names) {
        bn.names.string <- ""
        bn.names.string.fin <- ""
        getNameString <- function(name.string) {
          n.str <- ""
          n.str <- paste("[", name.string, "]", sep = "")
        }
        bn.names.string <- sapply(dat.names, getNameString)
        bn.names.string.fin <- ""
        for (i in bn.names.string)
        {
          bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
        }
        print(bn.names.string.fin)
      }
      buildDAG(dat.names)
      buildDAGade <- function(dat.names) {
        bn.names.string <- ""
        bn.names.string.fin <- ""
        getNameString <- function(name.string) {
          n.str <- ""
          n.str <- paste(name.string, ":", sep = "")
        }
        bn.names.string <- sapply(dat.names, getNameString)
        bn.names.string.fin <- "ade|medication:"
        for (i in bn.names.string)
        {
          bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
        }
        print(bn.names.string.fin)
      }
      buildDAGdrug <- function(dat.names) {
        bn.names.string <- ""
        bn.names.string.fin <- ""
        getNameString <- function(name.string) {
          n.str <- ""
          n.str <- paste(name.string, ":", sep = "")
        }
        bn.names.string <- sapply(dat.names, getNameString)
        bn.names.string.fin <- "medication|"
        for (i in bn.names.string)
        {
          bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
        }
        print(bn.names.string.fin)
      }
      predSTR <- buildDAG(cc)
      adeSTR <-
        paste("[",
              gsubfn(buildDAGade(cc), pattern = ":$", replacement = ""),
              "]",
              sep = "")
      drugSTR <-
        paste("[",
              gsubfn(buildDAGdrug(cc), pattern = ":$", replacement = ""),
              "]",
              sep = "")
      print(predSTR)
      print(adeSTR)
      print(drugSTR)
      dagSTR <- paste(predSTR, drugSTR, adeSTR, sep = "")
      print("###################DAG #########################")
      #print(dagSTR)
      dag <- bnlearn::model2network(dagSTR)
      plot(dag)
      # convert data to factors
      cols.to.factor <- colnames(dat)
      dat[cols.to.factor] <- lapply(dat[cols.to.factor] , factor)
      conf.names <- colnames(dat)[3:ncol(dat)] #dat.names[3:(ncol(dat))]
      conf.names <-
        str_replace_all(string = conf.names,
                        pattern = "-",
                        replacement = "HYPHEN")
      conf.names <- 
        str_replace_all(string = conf.names,
                        pattern = ",",
                        replacement = "COMMA")
      buildWL <- function(cnames) {
        wl.pre <- c()
        for (i in cnames) {
          wl.pre <- c(wl.pre, i, "ade", i, "medication")
        }
        return(wl.pre)
      }
      buildBL <- function(cnames) {
        bl.pre <- c("ade", "medication")
        for (i in cnames) {
          bl.pre <- c(bl.pre, "ade", i, "medication", i)
        }
        return(bl.pre)
      }
      wl <- matrix(buildWL(conf.names), ncol = 2, byrow = TRUE)
      bl <- matrix(buildBL(conf.names), ncol = 2, byrow = TRUE)
      print("WHITELIST")
      print(wl)
      print("BLACKLIST")
      print(bl)
      bn <-
        hc(
          dat,
          score = "bic",
          start = dag,
          whitelist = wl,
          blacklist = bl
        )
      drug.ade.score = 0
      print("################## network topology ################")
      medSTR <- paste("[medication|", cc, "]", sep ="")
      adeSTR <- paste("[medication|", cc, "]", sep ="")
      bn.model <- modelstring(bn)
      grep(pattern = medSTR, x = bn.model)
      grep(pattern = adeSTR, x = bn.model)
      if (grep(pattern = medSTR, x = bn.model) &&
          grep(pattern = adeSTR, x = bn.model)) {
        trueConfounders.local <-
          c(trueConfounders.local, cc)
        print("GOT ONE!!! WOOHOO!!!")
      }
    }
  }
  return(trueConfounders.local)
}
getTrueConfounders_compiled <- cmpfun(getTrueConfounders)


# bmc('sql', 'zidovudine', 'liver_failure,_acute', 'liver_failure,_acute', '1', '10')
