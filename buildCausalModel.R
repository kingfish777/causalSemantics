library(tidyverse)
library(dplyr)
library(stringr)
library(stringi)
library(gsubfn)
library(Rgraphviz)
library(bnlearn)
library(igraph)
library(car)
library(Hmisc)
library(utils)
library(compiler)
library(tmle)
library(gRain)
setwd("/home/smalec/Projects/ea")
source("helperFunctions.R")
# setwd("/Users/smalec/Projects/confounders")

args <- commandArgs(trailingOnly = TRUE)
cl <- as.character(args[1])
drug <- as.character(args[2])
ade <- as.character(args[3])
adr <- as.character(args[4])
group <- as.character(args[5])
SQUELCH <- as.character(args[6])

buildModel <- function(cl, drug, ade, adr, group, SQUELCH) {
#  if (is.na(cl)) { 
  # cl <- 'sql'; drug <- 'zidovudine'; ade <- 'liver_failure,_acute'; adr <- 'liver_failure,_acute'; group <- '1'; SQUELCH <- '10';
  #  drug <- "ibuprofen"
  #  ade <- "liver_failure,_acute"
  #  adr <- "ali"
   # group <- "1" 
   # SQUELCH <- "10"
  #  cl <- "psi"
#  }
  
    print("#########################################")
    print("#########################################")
    print("### READ DRUG / ADE data #################")
    print("#########################################")
    print("#########################################")
  drugData <-
    read.csv2(file =  paste("utgz/", drug, ".txt.gz", sep = ""), stringsAsFactors = FALSE)
  adeData  <-
    read.csv2(
      file = paste("utgz/", ade, ".txt.gz", sep = ""),
      stringsAsFactors = FALSE
    )
    drug <- str_to_lower(drug)
  
  print("#########################################")
  print("#########################################")
  print("### GETTING CONFOUNDERS #################")
  print("#########################################")
  print("#########################################")
  source("getCleanConfoundersJJ.R")
  confounders <- gcc(drug, ade, cl)
  
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/rawConfounders.txt", sep = "")
  writeVec(confounders, resultsFolder)
  system(paste("cat ", resultsFolder, sep = ""))
  trueConfounders <- getTrueConfounders_compiled(confounders, adeData, drugData, SQUELCH)
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/trueConfounders.txt", sep = "")
  writeVec(confounders, resultsFolder)
  system(paste("cat ", resultsFolder, sep = ""))
  
  print("#########################################")
  print("#########################################")
  print("### LOADING CONFOUNDER DATA #############")
  print("#########################################")
  print("#########################################")
  cleanDat <- c()
  for (i in 1:length(trueConfounders)) {
    cleanDat <- c(cleanDat, trueConfounders[i])
    cleanDat <-
      str_replace_all(string = cleanDat,
                      pattern = "-",
                      replacement = "HYPHEN")
    cleanDat <-
      str_replace_all(string = cleanDat,
                      pattern = ",",
                      replacement = "COMMA")
  }
  cleanDat.length <- length(cleanDat)
  print("######### list of clean confounders")
  print(cleanDat)
  ccvar.concepts <- c()
  print("################### READING IN confounder data###################")
  for (i in 1:length(cleanDat)) {
    ccvar.concept <- cleanDat[i]
    cleanDat[i] <-
      str_replace_all(string = cleanDat[i],
                      pattern = "HYPHEN",
                      replacement = "-")
    cleanDat[i] <-
      str_replace_all(string = cleanDat[i],
                      pattern = "COMMA",
                      replacement = ",")
    ccvar.concepts[[i]] <-
      read.csv2(file = paste("utgz/", cleanDat[i], ".txt.gz", sep = ""),
                stringsAsFactors = FALSE)}
  dat <- data.frame(adeData, drugData, ccvar.concepts)
  colnames(dat) <- c("ade", "medication", cleanDat)
  print(colnames(dat))
  colnames(dat) <-
    str_replace_all(string = colnames(dat),
                    pattern = "-",
                    replacement = "HYPHEN")
  colnames(dat) <- 
    str_replace_all(string = colnames(dat),
                    pattern = ",",
                    replacement = "COMMA")
  dat = dat[,colSums(dat) > 0]
  cols.to.numeric <- colnames(dat)
  dat[cols.to.numeric] <- lapply(dat[cols.to.numeric] , as.numeric)
  
  print("#########################################")
  print("#########################################")
  print("### BUILDING GLM BASELINE ###############")
  print("#########################################")
  print("#########################################")
  
  dat.df <- data.frame(ftable(dat, col.vars = c("ade"), row.vars = names(dat)[2:length(names(dat))]))
  dat.df.raw <- dat.df[1:(nrow(dat.df)/2), 1:(ncol(dat.df) - 2)] 
  n.tot <- dat.df[1:(nrow(dat.df)/2), ncol(dat.df)]
  n.ade <- dat.df[(1 + nrow(dat.df)/2):nrow(dat.df), ncol(dat.df)] 
  ade.tbl <- cbind(n.ade, n.tot)
  dat.names <- colnames(dat[, 2:ncol(dat)])
  print(dat.names)
  ########################################
  buildGLM <- function(dnames) {
    glm.dnames <- ""
    #paste(dnames[1], "", sep = "")
    for (i in dnames[2:length(dnames)]) {
      glm.dnames <- paste(glm.dnames, "+", i, sep = "")
    }
    return(glm.dnames)
  }
  glm.preds <- buildGLM(dat.names[2:length(dat.names)])
  print(glm.preds)
  pv.glm <-
    paste(
      "pv.glm <- glm(",
      "ade.tbl",
      " ~ medication",
      glm.preds,
      ", data = dat.df.raw, family = 'binomial')",
      sep = ""
    )
  print(pv.glm)
  try(eval(parse(text = pv.glm)))
  med.coef <- 0
  med.coef.adj <- 0
  try(med.coef.adj <- coef(pv.glm)[2])
  
  try(pv.glm0 <-
        glm(ade.tbl ~ medication, data = dat.df.raw, family = 'binomial'))[2]
  nerd <- summary(pv.glm0)
  nerd.df <- data.frame(coef(nerd))
  
  intercept0 <- nerd.df[1,1]
  stdErrorIntercept0 <- nerd.df[1,2]
  zValueIntercept0 <- nerd.df[1,3]
  pvalueIntercept0 <- nerd.df[1,4]
  medCoef0 <- nerd.df[2,1]
  medStdError0 <- nerd.df[2,2]
  medZvalue0 <- nerd.df[2,3]
  medPvalue0 <- nerd.df[2,4]
  results0 <- c(intercept0, stdErrorIntercept0, zValueIntercept0, pvalueIntercept0, medCoef0, medStdError0, medZvalue0, medPvalue0)
  print(results0)
  
  if (length(med.coef) == 0) {
    med.coef <- 0
  } else {
    med.coef #OUTPUT
  }
  if (length(med.coef.adj) == 0) {
    med.coef.adj <- 0
  } else {
    med.coef.adj #OUTPUT
  }
  
  nerd <- summary(pv.glm)
  nerd.df <- data.frame(coef(nerd))
  
  intercept <- nerd.df[1,1]
  stdErrorIntercept <- nerd.df[1,2]
  zValueIntercept <- nerd.df[1,3]
  pvalueIntercept <- nerd.df[1,4]
  medCoef <- nerd.df[2,1]
  medStdError <- nerd.df[2,2]
  medZvalue <- nerd.df[2,3]
  medPvalue <- nerd.df[2,4]
  
  results <- c(intercept, stdErrorIntercept, zValueIntercept, pvalueIntercept, medCoef, medStdError, medZvalue, medPvalue)
  print(results)
  print(results0)
  
  resultsFolder0 <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/glm0.txt", sep = "")
  writeVec(results0, resultsFolder0)
  system(paste("cat ", resultsFolder0, sep = ""))
  
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/glm.txt", sep = "")
  writeVec(results, resultsFolder)
  system(paste("cat ", resultsFolder, sep = ""))
  
  a <- sum(adeData*drugData)
  b <- sum(drugData) - a
  c <- sum(adeData) - a
  d = nrow(drugData) - a - b - c
  ade.tbl <- matrix(c(a, b, c, d), byrow = TRUE, 2, 2)
  chisqResult <- chisq.test(ade.tbl)
  chisqStat <- chisqResult$statistic
  chisqPvalue <- chisqResult$p.value
  results <- c(chisqStat, chisqPvalue)
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/chisq.txt", sep = "")
  writeVec(results, resultsFolder)
  system(paste("cat ", resultsFolder, sep = ""))
  dat.names <- colnames(dat[, 3:ncol(dat)])
  print("################################################################")
  print("########### CONVERTING DATA BACK TO FACTOR FORMAT ##############")
  print("################################################################")
  # convert data to factors
  cols.to.factor <- colnames(dat)
  dat[cols.to.factor] <- lapply(dat[cols.to.factor] , factor)
  
  print("#########################################")
  print("#########################################")
  print("### BUILD GRAPHICAL CAUSAL MODELS #######")
  print("#########################################")
  print("#########################################")  
  
  
  print("learning model from raw data ...")
  bnRaw <- hc(dat, optimized = TRUE)
  mbAde <- mb(bnRaw, "ade")
  mbMedication <- mb(bnRaw, "medication")
  bestConfounders <- intersect(mbAde, mbMedication)
  
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/mbAde.txt", sep = "")
  writeVec(mbAde, resultsFolder)
  
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/mbMedication.txt", sep = "")
  writeVec(mbMedication, resultsFolder)
  
  resultsFolder <- paste(unlist(newResultsFolderRollout(cl, drug, ade, adr, group, SQUELCH)), "/bestConfounders.txt", sep = "")
  writeVec(bestConfounders, resultsFolder)
  system(paste("cat ", resultsFolder, sep = ""))
  
  print("#############################################################################")
  print("#############################################################################")
  print("now incorporating knowledge from the literature")
  print("#############################################################################")
  print("#############################################################################")
  dag <- constructDAG(dat.names)
  print(dag)
  #conf.names <- colnames(dat)[3:ncol(dat)] #dat.names[3:(ncol(dat))]
  
  buildWLwithEdge <- function(cnames) {
    wl.pre <- matrix(c("medication", "ade"), ncol = 2, byrow = TRUE)
    for (i in cnames) {
      wl.pre <- c(wl.pre, i, "ade", i, "medication")
    }
    return(wl.pre)
  }
  buildWLwithEdge(dat.names)
  
  buildWL <- function(cnames) {
    wl.pre <- c()
    for (i in cnames) {
      wl.pre <- c(wl.pre, i, "ade", i, "medication")
    }
    return(wl.pre)
  }
  
  buildBL <- function(cnames) {
    bl.pre <- c()
    for (i in cnames) {
      bl.pre <- c(bl.pre, "ade", i, "medication", i)
    }
    return(bl.pre)
  }
  wl <- matrix(buildWLwithEdge(dat.names), ncol = 2, byrow = TRUE)
  bl <- matrix(buildBL(dat.names), ncol = 2, byrow = TRUE)
  bl <- rbind(bl, tiers2blacklist(dat.names))
  print("WHITELIST")
  print(wl)
  print("BLACKLIST")
  print(bl)
  
  print("#########################################")
  print("#########################################")
  print("### BUILDING GRAPHICAL CAUSAL MODELS ####")
  print("######## NO EDGE ASSUMED ################")
  print("#########################################")  
  print("################# learning topology of bayesian network #######################")
  
  
  bn <-
      hc(
      dat,
      score = "bic",
      start = dag,
      blacklist = bl,
      optimized = TRUE#,
      #debug = TRUE
    )
  
  drug.ade.score = 0
  print("################## network topology ################")
  print(modelstring(bn))
  #bn <- bn.fit(data = dat, method = "mle", start = dag, whitelist = wl, blacklist = bl, iss = 20000)
  # result = tryCatch({
  power <- arc.strength(bn, data = dat, criterion = "bic")
  arc.score.coord <-
    intersect(which(power$from == "medication"),  which(power$to == "ade"))
  drug.ade.score0 <- power$strength[arc.score.coord]
  print("#############SCORE#################")
  if (isEmpty(drug.ade.score0)) { drug.ade.score0 <- 0 }
  print("drug ADE score")
  print(drug.ade.score0)
  print("#########################################################")
  
  
  print("########## GRAPHING NO EDGE ASSUMED ##################################")
  bnf <- bn.fit(bn, data = dat, method = "mle", keep.fitted = TRUE)
  print(bnf)
  # prop.table(table(dat[ , c("medication", "ade")]), margin = 2)
  bn.graphNEL <- bnlearn::as.graphNEL(bn)
  par(font.axis = 2)
  par(font.lab = 2)
  par(family = "mono")
  # save plot
  system(paste("mkdir ", cl, "eaBL", SQUELCH, sep = ""))
  png(paste(cl, "eaBL", SQUELCH, "/", adr, "_", group, "_", drug, ".png", sep = ""), width = 4, height = 4, units = 'in', res = 300)
  gg <-
    strength.plot(
      main = paste(drug, " ", ade, " ", group,  sep = ""),
      arc.strength(bn, data = dat, criterion = "bic"),
      x = bn,
      threshold = 0,
      highlight = list(
        nodes = c("ade", "medication"),
        fill = c("red"),
        textCol = c("black")
      ),
      shape = "ellipse"
    )
  dev.off(which = dev.cur())
  
  
  print("#########################################")
  print("#########################################")
  print("### BUILDING GRAPHICAL CAUSAL MODELS ####")
  print("########### EDGE ASSUMED ################")
  print("#########################################")  
  print("################# learning topology of bayesian network #######################")
  
  
  bn <-
    hc(
      dat,
      score = "bic",
      start = dag,
      whitelist = wl,
      blacklist = bl,
      optimized = TRUE#,
      #debug = TRUE
    )
  
  drug.ade.score = 0
  print("################## network topology ################")
  print(modelstring(bn))
  #bn <- bn.fit(data = dat, method = "mle", start = dag, whitelist = wl, blacklist = bl, iss = 20000)
  # result = tryCatch({
  power <- arc.strength(bn, data = dat, criterion = "bic")
  arc.score.coord <-
    intersect(which(power$from == "medication"),  which(power$to == "ade"))
  drug.ade.score <- power$strength[arc.score.coord]
  print("#############SCORE#################")
  if (isEmpty(drug.ade.score)) { drug.ade.score <- 0 }
  print("drug ADE score")
  print(drug.ade.score)
  print("#########################################################")
  

  bnf <- bn.fit(bn, data = dat, method = "mle", keep.fitted = TRUE)
  print(bnf)
  prop.table(table(dat[ , c("medication", "ade")]), margin = 2)
  bn.graphNEL <- bnlearn::as.graphNEL(bn)
  par(font.axis = 2)
  par(font.lab = 2)
  par(family = "mono")
  # save plot
  system(paste("mkdir ", cl, "eaWL", SQUELCH, sep = ""))
  png(paste(cl, "eaWL", SQUELCH, "/", adr, "_", group, "_", drug, ".png", sep = ""), width = 4, height = 4, units = 'in', res = 300)
  gg <-
    strength.plot(
      main = paste(drug, " ", ade, " ", group,  sep = ""),
      arc.strength(bn, data = dat, criterion = "bic"),
      x = bn,
      threshold = 0, #10^-10,
      highlight = list(
        nodes = c("ade", "medication"),
        fill = c("red"),
        textCol = c("black")
      ),
      shape = "ellipse"
    )
  dev.off(which = dev.cur())
  
  #print("############################")
  #print("############################")
  #print("### CALCULATING SCORE FROM MUTILATED NETWORK")
  #print("############################")
  #print("############################")
  try({
  bng <- as.graphNEL(bn)
  # check fit without edge
  bng0 <- removeEdge(from = "medication", to = "ade", graph = bng)
  bn0 <- as.bn(bng0)
  networkScore0 <- score(bn0, data = dat, type = "bic") # OUTPUT
  networkScore1 <- score(bn, data = dat, type = "bic") # OUTPUT
  print("networkScore0")
  print(networkScore0)
  print("networkScore1")
  print(networkScore1)
  networkScoreDiff <- networkScore1 - networkScore0
  print("networkScoreDiff")
  print(networkScoreDiff)
  })
  
  print("############################")
  print("############################")
  print("### get k-fold validated ATEs #######")
  print("### ---exact inference methods ----####")
  print("############################")
  print("############################")
  deltas <- c()
  dart <- bn.cv(data = dat, fit = "mle", bn = bn) #, method = "k-fold", k = 10, m = 10, runs = 10)
  for (i in 1:10) {
    #bnf.gg <- as.grain(bnf)
    bnf.gg <- as.grain(dart[[i]]$fitted)
    jtree <- compile(bnf.gg)
    y0 <- setEvidence(jtree, nodes = "medication", states = "0")
    y1 <- setEvidence(jtree, nodes = "medication", states = "1")
    yy0 <- querygrain(y0, nodes = "ade")$ade[2]
    yy1 <- querygrain(y1, nodes = "ade")$ade[2]
    delta <- yy1 - yy0
    deltas <- c(deltas, delta)
  }
  hist(deltas)
  exactATE <- mean(deltas)
  print("gRain y1 - y0")
  print(exactATE) #OUTPUT
  exactATEvar <- var(deltas) #OUTPUT
  
 # cols.to.integer <- colnames(dat)
 # dag[cols.to.integer] <- lapply(dat[cols.to.integer], as.integer)
  
 # tmleConfounders <- paste(dat.names, collapse = "+")
 # w <- dat[,3:length(dat.names)]
 # attach(w)
  #tmleString <- paste("result <- tmle(Y = ade, A = medication, W = dat, family = 'binomial', ",
  #                    "Qform = ade~ medication+", tmleConfounders, ", ", "gform = medication~ ", tmleConfounders, ")", sep = "")
 # print("######")
  #print(tmleString)
 # adeDat <- as.integer(dat[,1]) - 1
 # drugDat <- as.integer(dat[,2]) - 1
 # tmleStringprologue <- paste("result <- tmle(Y = adeDat, A = drugDat, W = w, family = 'binomial'", sep = "")
  
 # Qform <- paste("Qform = ade~ medication+", tmleConfounders, sep = "")
 # gform <- paste("gform = medication~ ", tmleConfounders, sep = "")
 # tmleString = paste(tmleStringprologue, Qform, gform, sep = ", ")
 # tmleString = paste(tmleString, ")", sep = "")
 # print("######")
 # print(tmleString)
 # try(eval(parse(text = tmleString)))
  
  
  
  
  #exposureName, hoiName, caseControl, confounderSource, pFinding0, pFinding1, approxATE, exactATE, medCoef0, medCoef1     
  #sqlStatement <- paste("INSERT INTO enkiModels (exposureName, hoiName, caseControl, confounderSource, squelch, pFinding0, pFinding1, approxATE, exactATE, medCoef0, medCoef1) VALUES('", drug, "', '", adr, "', '", group, "', '", cl, "', '", SQUELCH, "', '", pFinding0, "', '", pFinding1, "', '", approxATE, "', '", exactATE, "', '", med.coef, "', '", med.coef.adj, "');", sep = "")
  #print(sqlStatement)
  #dbGetQuery(con, sqlStatement)
  
  results <- data.frame(
    adr, #1
    cl, #2
    group,#3
    ade,#4
    drug,#5
    cl,#6
    SQUELCH,#7
    chisqPvalue,#8
    drug.ade.score0,#9
    drug.ade.score,#10
    exactATE,#11
    exactATEvar,#12
    medCoef0,#13
    medZvalue0,#14
    medPvalue0,#15
    medCoef,#16
    medZvalue,#17
    medPvalue,#18
    networkScore0,#19
    networkScore1,#20 
    networkScoreDiff #21,
    # att #,
  )
  print(paste("writing results: ", results))
  write.table(file = paste(cl, "ea", SQUELCH, ".txt", sep = ""), x = results, col.names = FALSE, sep = "\t", append = TRUE, row.names = FALSE)
}

bmc <- cmpfun(buildModel)
print("#########b ###M #######################")
print("########## u###o ######## E ###########")
print("############ i## d ######## N  ########")
print("############# d###e ######## K ########")
print("###################l ########  I ######")
print("#######################################")
print("bmc(cl, drug, ade, adr, group, SQUELCH)") 
print("############ compiling enki function ##")
print("#######################################")
print("#######################################")
print("#######################################")


