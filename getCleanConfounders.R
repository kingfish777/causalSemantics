
library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(compiler)
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "smdb",
    host = "localhost",
    port = 5432,
    user = "smalec",
    password = "mandarin"
  )

# DEPRESSION
#######
#  Antidepressive_Agents                    | C0003289
#Antidepressive_Agents,_Second-Generation | C0242905
#Antidepressive_Agents,_Tricyclic         | C0003290
#Depressive_Symptoms                      | C0086132
#Depressive_disorder                      | C0011581
#Major_Depressive_Disorder                | C1269683

#
# Alzheimer's_Disease                      | C0002395
# Alzheimer_Disease,_Early_Onset           | C0750901
# Alzheimer_Disease,_Late_Onset            | C0494463
# Familial_Alzheimer's_disease             | C0276496
#PS2_protein_(alzheimer-associated)       | C0528480
#PS2_protein_(alzheimer-associated)|PSEN2 | C0528480|5664



# obesity C0028754

# Hypertension
#  Hypertensive_disease              | C002053
# Hypertensive                      | C0857121


#Post-traumatic_brain_syndrome | C0546983
#Traumatic_Brain_Injury        | C0876926

# select DISTINCT subject_name, subject_cui from smdbprepreds where LOWER(subject_name) LIKE LOWER('%alcohol%');
# 
# ALCOHOL_WITHDRAWAL_ACUTE  C0149821
# acute_alcohol_abuse   C0497318
# Acute_alcoholism      C0699757
# Alcohol-Related_Disorders C0236664
# Alcohol_Withdrawal_Delirium C0001957
# Alcohol_abuse         C00085762
# Alcoholic_Intoxication,_Chronic C0001973

# Conductive_hearing_loss      | C0018777
# HEARING_LOSS_CONGENITAL      | C0262435
# Hearing_Loss,_Bilateral      | C0018775
# Hearing_Loss,_Cochlear       | C1527369
# Hearing_Loss,_High-Frequency | C0018780
# Hearing_Loss,_Unilateral     | C0521785
# Neural_hearing_loss          | C0155550
# Sensorineural_Hearing_Loss   | C0018784



exp <- "%depressive%"; hoi <- "%alzheimer%"; confounderSource <- "rawConfounders"

getCleanConfounders <- function(exp, hoi, confounderSource) {
  
  if (confounderSource == "rawConfounders") {
    
    confounders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate IN ('TREATS', 'CAUSES', 'PREDISPOSES')
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    confounders <- str_to_lower(subset(confounders, theta1*theta2 > 1)[,1])
    
    covariates = confounders
    
  } else if (confounderSource == "cookedConfounders") {
    
    confounders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'TREATS'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    confounders <- str_to_lower(subset(confounders, theta1*theta2 > 1)[,1])
    
    colliders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'CAUSES'
                        GROUP BY object_name), 
                        ZY AS (SELECT object_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY object_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    colliders <- str_to_lower(subset(colliders, theta1*theta2 > 1)[,1])
    colliders
    
    mediators <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'CAUSES'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    mediators <- str_to_lower(subset(mediators, theta1*theta2 > 1)[,1])
    mediators
    
    try({mediators <- setdiff(mediators, confounders)
    mediators <- setdiff(mediators, colliders)
    print("#### minus stopwords")
    print(mediators)
    print(length(mediators))})             
    
    try({colliders <- setdiff(colliders, confounders)
    colliders <- setdiff(colliders, mediators)
    print("#### minus stopwords")
    print(colliders)
    print(length(colliders))})
    
    try({confounders <- setdiff(confounders, colliders)
    confounders <- setdiff(confounders, mediators)
    print("#### minus stopwords")
    print(confounders)
    print(length(confounders))})
    
    
    covariates = confounders
    
  } else if (confounderSource == "cookedConfoundersPlusPrecisionVars") {
    
    confounders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'TREATS'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    confounders <- str_to_lower(subset(confounders, theta1*theta2 > 1)[,1])
    
    
    colliders <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'CAUSES'
                        GROUP BY object_name), 
                        ZY AS (SELECT object_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY object_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    colliders <- str_to_lower(subset(colliders, theta1*theta2 > 1)[,1])
    colliders
    
    mediators <-dbGetQuery(con, paste("WITH ZX AS (SELECT object_name AS covar, COUNT(*) AS ocnt
                        FROM predication 
                        WHERE lower(subject_name) LIKE '%", tolower(exp),"%' 
                        AND predicate = 'CAUSES'
                        GROUP BY object_name), 
                        ZY AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.ocnt as theta1, zy.scnt AS theta2 FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta1 DESC, theta2 DESC", sep = ""))
    
    mediators <- str_to_lower(subset(mediators, theta1*theta2 > 1)[,1])
    mediators
    
    
    precisionVariables <-dbGetQuery(con, paste("WITH ZX AS (SELECT subject_name AS covar , COUNT(*) AS scnt
                        FROM predication 
                        WHERE lower(object_name) LIKE '%", tolower(hoi),"%' 
                        AND predicate IN ('CAUSES', 'PREDISPOSES')
                        GROUP BY subject_name) 
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1 FROM ZX AS zx ORDER BY theta1 DESC", sep = ""))
    
    precisionVariables <- str_to_lower(subset(precisionVariables, theta1 > 1)[,1])
    
    
    try({mediators <- setdiff(mediators, confounders)
    mediators <- setdiff(mediators, colliders)
    print("#### minus stopwords")
    print(mediators)
    print(length(mediators))})             
    
    try({colliders <- setdiff(colliders, confounders)
    colliders <- setdiff(colliders, mediators)
    print("#### minus stopwords")
    print(colliders)
    print(length(colliders))})
    
    try({confounders <- setdiff(confounders, colliders)
    confounders <- setdiff(confounders, mediators)
    print("#### minus stopwords")
    print(confounders)
    print(length(confounders))})
    
    try({precisionVariables <- setdiff(precisionVariables, confounders)
    precisionVariables <- setdiff(precisionVariables, colliders)
    precisionVariables <- setdiff(precisionVariables, mediators)
    print("#### minus stopwords")
    print(precisionVariables)
    print(length(precisionVariables))})
    
    covariates = c(confounders, precisionVariables)
  }
  
  # try({expSynonyms <- c(read.csv2(file = paste("variables/syns/synonyms/", exp, ".txt", sep = ""), stringsAsFactors = FALSE, header = FALSE)[1])})
  # try({hoiSynonyms <- c(read.csv2(file = paste("variables/syns/synonyms/", hoi, ".txt", sep = ""), stringsAsFactors = FALSE, header = FALSE)[1])})
  
  #try({stopwords <-
  #as.character(
  #  c(read.csv2(
  #    file = "variables/stoplist.txt",
  #    #   sep = "\n",
  #    stringsAsFactors = FALSE,
  #    header = FALSE
  #  )[1])})
  
  #####
  #try({confounders <- setdiff(confounders, stopwords[[1]])
  #print("#### minus stopwords")
  #print(length(confounders))})
  ###
  #try({confounders <- setdiff(confounders, expSynonyms[[1]])
  #print("#### minus expSynonyms")
  #print(length(confounders))})
  ###
  #try({confounders <- setdiff(confounders, hoiSynonyms[[1]])
  #print("#### minus hoiSynonyms")
  #print(length(confounders))})
  
  #confounders.v1 <- c()
  #for (confounder in confounders) {
  #  print(confounder)
  #  try({fName <- paste("utgz/", confounder, ".txt.gz", sep = "")})
  #  try({if (file.exists(fName)) { confounders.v1 <- c(confounders.v1, confounder); print(confounders.v1) }
  #    else if (!file.exists(fName)) { print(paste(confounder, "WAS NOT MEASURED OR RECORDED"))}})
  #}
  print("MINUS variables not in data")
  print(length(covariates))
  #print(length(confounders.v1))
  return(covariates)
}

gcc <- cmpfun(getCleanConfounders)

rawConfounders <- gcc(exp = exp, hoi = hoi, confounderSource = "rawConfounders")
length(rawConfounders)
rawConfounders

cookedConfounders <- gcc(exp = exp, hoi = hoi, confounderSource = "cookedConfounders")
length(cookedConfounders)
cookedConfounders

cookedConfoundersPlusPrecisionVars <- gcc(exp = exp, hoi = hoi, confounderSource = "cookedConfoundersPlusPrecisionVars")
length(cookedConfoundersPlusPrecisionVars)
cookedConfoundersPlusPrecisionVars



# 
# getTrueConfounders <- function(confounders, adeData, drugData, SQUELCH) {
#   trueConfounders.local <- c()
#   for (cc in confounders) {   # cc = "confounder candidate"
#     if (length(trueConfounders.local) >= SQUELCH)
#     { break }
#     print("DOOKIE")
#     if (file.exists(paste("utgz/", cc, ".txt.gz", sep = ""))) {
#       cc <- as.character(cc)
#       print(cc)
#       ccData <-
#         read.csv2(file = paste("utgz/", cc, ".txt.gz", sep = ""),
#                   stringsAsFactors = FALSE)
#       cc <-
#         str_replace_all(string = cc,
#                         pattern = "-",
#                         replacement = "HYPHEN")
#       cc <- 
#         str_replace_all(string = cc, 
#                         pattern = ",",
#                         replacement = "COMMA")
#       dat <- data.frame(adeData, drugData, ccData)
#       colnames(dat) <- c("ade", "medication", cc)
#       colnames(dat) <-
#         str_replace_all(string = colnames(dat),
#                         pattern = "-",
#                         replacement = "-")
#       colnames(dat) <- 
#         str_replace_all(string = colnames(dat),
#                         pattern = ",",
#                         replacement = ",")
#       write.table(
#         dat,
#         file = "dat.txt",
#         sep = "\t",
#         row.names = FALSE,
#         col.names = TRUE,
#         quote = FALSE
#       )
#       dat.names <- cc
#       #print(dat.names)
#       buildDAG <- function(dat.names) {
#         bn.names.string <- ""
#         bn.names.string.fin <- ""
#         getNameString <- function(name.string) {
#           n.str <- ""
#           n.str <- paste("[", name.string, "]", sep = "")
#         }
#         bn.names.string <- sapply(dat.names, getNameString)
#         bn.names.string.fin <- ""
#         for (i in bn.names.string)
#         {

#           bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
#         }
#         print(bn.names.string.fin)
#       }
#       buildDAG(dat.names)
#       buildDAGade <- function(dat.names) {
#         bn.names.string <- ""
#         bn.names.string.fin <- ""
#         getNameString <- function(name.string) {
#           n.str <- ""
#           n.str <- paste(name.string, ":", sep = "")
#         }
#         bn.names.string <- sapply(dat.names, getNameString)
#         bn.names.string.fin <- "ade|medication:"
#         for (i in bn.names.string)
#         {
#           bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
#         }
#         print(bn.names.string.fin)
#       }
#       buildDAGdrug <- function(dat.names) {
#         bn.names.string <- ""
#         bn.names.string.fin <- ""
#         getNameString <- function(name.string) {
#           n.str <- ""
#           n.str <- paste(name.string, ":", sep = "")
#         }
#         bn.names.string <- sapply(dat.names, getNameString)
#         bn.names.string.fin <- "medication|"
#         for (i in bn.names.string)
#         {
#           bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
#         }
#         print(bn.names.string.fin)
#       }
#       predSTR <- buildDAG(cc)
#       adeSTR <-
#         paste("[",
#               gsubfn(buildDAGade(cc), pattern = ":$", replacement = ""),
#               "]",
#               sep = "")
#       drugSTR <-
#         paste("[",
#               gsubfn(buildDAGdrug(cc), pattern = ":$", replacement = ""),
#               "]",
#               sep = "")
#       print(predSTR)
#       print(adeSTR)
#       print(drugSTR)
#       dagSTR <- paste(predSTR, drugSTR, adeSTR, sep = "")
#       print("###################DAG #########################")
#       #print(dagSTR)
#       dag <- bnlearn::model2network(dagSTR)
#       plot(dag)
#       # convert data to factors
#       cols.to.factor <- colnames(dat)
#       dat[cols.to.factor] <- lapply(dat[cols.to.factor] , factor)
#       conf.names <- colnames(dat)[3:ncol(dat)] #dat.names[3:(ncol(dat))]
#       conf.names <-
#         str_replace_all(string = conf.names,
#                         pattern = "-",
#                         replacement = "HYPHEN")
#       conf.names <- 
#         str_replace_all(string = conf.names,
#                         pattern = ",",
#                         replacement = "COMMA")
#       buildWL <- function(cnames) {
#         wl.pre <- c()
#         for (i in cnames) {
#           wl.pre <- c(wl.pre, i, "ade", i, "medication")
#         }
#         return(wl.pre)
#       }
#       buildBL <- function(cnames) {
#         bl.pre <- c("ade", "medication")
#         for (i in cnames) {
#           bl.pre <- c(bl.pre, "ade", i, "medication", i)
#         }
#         return(bl.pre)
#       }
#       wl <- matrix(buildWL(conf.names), ncol = 2, byrow = TRUE)
#       bl <- matrix(buildBL(conf.names), ncol = 2, byrow = TRUE)
#       print("WHITELIST")
#       print(wl)
#       print("BLACKLIST")
#       print(bl)
#       bn <-
#         hc(
#           dat,
#           score = "bic",
#           start = dag,
#           whitelist = wl,
#           blacklist = bl
#         )
#       drug.ade.score = 0
#       print("################## network topology ################")
#       medSTR <- paste("[medication|", cc, "]", sep ="")
#       adeSTR <- paste("[medication|", cc, "]", sep ="")
#       bn.model <- modelstring(bn)
#       grep(pattern = medSTR, x = bn.model)
#       grep(pattern = adeSTR, x = bn.model)
#       if (grep(pattern = medSTR, x = bn.model) &&
#           grep(pattern = adeSTR, x = bn.model)) {
#         trueConfounders.local <-
#           c(trueConfounders.local, cc)
#         print("GOT ONE!!! WOOHOO!!!")
#       }
#     }
#   }
#   return(trueConfounders.local)
# }
# getTrueConfounders_compiled <- cmpfun(getTrueConfounders)
# 
# 
# # stop <- read.csv("variables/stoplist.txt", stringsAsFactors = FALSE, col.names = "nork")
# # paste0(gsubfn(x = sort(data.frame(stop)[,1]), pattern = "_", replacement = " "), collapse = "; ")
# 
# 
# 
