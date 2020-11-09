

library(MLmetrics)
library(pROC)
library(PRROC)
library(readr)
library(xtable)
library(ROCR)
library(openxlsx)
library(RPostgreSQL)
library(ggpubr)
#install.packages("extrafont")
library(extrafont)
font_import()

################
# LOAD DATABASE
################
# load database
# createdb research
# createdb confounders
# createuser research # create postgres user named 'research'
# # import confounders database
# psql confounders -U < confounders.dat
#

drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "confounders",
    host = "localhost",
    port = 5432,
    user = "research",
    password = "mandarin"
  )

setwd("~/Documents/causalSemantics")
#######################################################
# GET (rectangular) DATA(set) TO BE ANALYZED
#######################################################
dat <-
  dbGetQuery(
    con,
    "SELECT DISTINCT sql5.casecontrol, sql5.hoiname, sql5.exposurename, 
                          ref.a, ref.b, ref.c, ref.d, ref.ror, sql5.medcoef0,
                          sql5.medcoef1 AS sql5medcoef1,
                          psi5.medcoef1 AS psi5medcoef1, 
                          sql10.medcoef1 AS sql10medcoef1,
                          psi10.medcoef1 AS psi10medcoef1,
                          sql5.exactATE AS sql5ate, 
                          psi5.exactATE AS psi5ate, 
                          sql10.exactATE AS sql10ate, 
                          psi10.exactATE AS psi10ate
                FROM penkiModels sql5, penkiModels psi5, penkiModels sql10, 
                      penkiModels psi10, refsetbaselines ref
                 WHERE sql5.hoiname = psi5.hoiname AND sql5.hoiname = ref.hoiname AND
                        sql5.exposurename = psi5.exposurename AND sql5.exposurename = ref.exposurename AND
                        sql5.casecontrol = psi5.casecontrol AND CAST(sql5.casecontrol AS INT) = ref.casecontrol AND
                        sql10.hoiname = psi10.hoiname AND sql10.hoiname = ref.hoiname AND
                        sql10.exposurename = psi10.exposurename AND sql10.exposurename = ref.exposurename AND
                        sql10.casecontrol = psi10.casecontrol AND CAST(sql10.casecontrol AS INT) = ref.casecontrol AND
                        sql5.confoundersource = 'sql' AND psi5.confoundersource = 'psi' AND
                        sql10.confoundersource = 'sql' AND psi10.confoundersource = 'psi' AND
                        sql5.squelch = '5' AND psi5.squelch = '5' AND
                        sql10.squelch = '10' AND psi10.squelch = '10'"
  )

# close open psql connections
lapply(dbListConnections(drv = dbDriver("PostgreSQL")), function(x) {dbDisconnect(conn = x)})
on.exit(dbUnloadDriver(drv), add = TRUE)

#######################################################
# Calculate chisq - append to rectangular dataset
#######################################################
chisqList <- c()
for (i in 1:nrow(dat)) {
  chisq <-
    chisq.test(rbind(cbind(dat$a[i], dat$c[i]), cbind(dat$b[i], dat$d[i])))$p.value
  #print(chisq)
  chisqList <- c(chisqList, chisq)
}
dat$chisq <- chisqList
dat.orig <- dat
nrowTot <- nrow(dat)

#names(dat) <- c("adr", "cl", "group", "ade", "drug", "cl", "SQUELCH", "chisqPvalue", "drug.ade.score0", "drug.ade.score", "exactATE", "exactATEvar", "medCoef0", "medZvalue0", "medPvalue0", "medCoef", "medZvalue", "medPvalue", "networkScore0", "networkScore1", "networkScoreDiff")
#
names(dat) <- c("casecontrol", "hoiname", "exposurename", "a", "b", "c", "d", "ror", "medcoef0", "sql5medcoef1", "psi5medcoef1", "sql10medcoef1", "psi10medcoef1", "sql5ate", "psi5ate", "sql10ate", "psi10ate", "chisq")

#######################################################
# write data locally
#######################################################
write_delim("resultsData-dib.tsv",x = dat,delim = "\t")
system("head -n 3 resultsData-dib.tsv")

#######################################################
# Get number of +case and -control pairs for each adverse event type
#######################################################
getPosNegCounts <- function(dat2) {
  ades <-
    c(
      "kidney_failure,_acute",
      "liver_failure,_acute",
      "acute_myocardial_infarction",
      "gastrointestinal_hemorrhage"
    )
  summaries <- c()
  for (jk in 1:length(ades)) {
    dat <- subset(dat2, hoiname == ades[jk])
    nrowDat <- nrow(dat)
    pos <- sum(as.numeric(dat$casecontrol))
    neg <- nrowDat - pos
    summary <- c(ades[jk], nrowDat, pos, neg)
    #print(summary)
    summaries <- rbind(summaries, summary)
  }
  #cbind(summaries, summaries)
  summaries
  totals <-
    c("totals",
      sum(as.numeric(summaries[, 2])),
      sum(as.numeric(summaries[, 3])),
      sum(as.numeric(summaries[, 4])))
  rbind(summaries, totals)
}

counts <- getPosNegCounts(dat)
countsDF <- data.frame(counts, row.names = NULL)
names(countsDF) <- c("ADE_type", "Total", "Pos", "Neg")
print(countsDF)

# LaTeX output
xtable(countsDF)

#######################################################
# LOAD model performance statistics data from experiment using randomly selected covariates
#######################################################
randomAKI <- data.frame(read.csv("randomAKI.txt", header = FALSE, sep = "\t"))
randomALI <- data.frame(read.csv("randomALI.txt", header = FALSE, sep = "\t"))
randomAMI <- data.frame(read.csv("randomAMI.txt", header = FALSE, sep = "\t"))
randomGIB <- data.frame(read.csv("randomGIB.txt", header = FALSE, sep = "\t"))

columnNames <- c("hoiname", "cl", "casecontrol", "hoiname", "exposurename", "cl", "SQUELCH",  "chisqPvalue", "drug.ade.score0", "drug.ade.score", "exactATE", "exactATEvar",  "medCoef0", "medZvalue0", "medPvalue0", "medCoef", "medZvalue", "medPvalue",  "networkScore0", "networkScore1", "networkScoreDiff")
names(randomAKI) <- columnNames
names(randomALI) <- columnNames
names(randomAMI) <- columnNames
names(randomGIB) <- columnNames
randomAKI <- randomAKI[,unique(names(randomAKI))]
randomALI <- randomALI[,unique(names(randomALI))]
randomAMI <- randomAMI[,unique(names(randomAMI))]
randomGIB <- randomGIB[,unique(names(randomGIB))]

randomALL <- rbind(randomAKI, randomALI, randomAMI, randomGIB)

write_delim("resultsRandom-dib.tsv",x = randomALL,delim = "\t")
system("head -n 3 resultsRandom-dib.tsv")

# #######################################################
# # Calculate stats given:
# #######################################################
# # input: predictions + expected labels (1=causal or 0=non-causal)
# #######################################################
# getStats <- function(classifierPredictions, labels) {
#   #classifierPredictions <- as.numeric(dat$sql5medcoef1); labels <- as.factor(dat$casecontrol)
#   y = as.factor(labels)
#   pred <- prediction(classifierPredictions, y)
#   # Recall-Precision curve
#   RP.perf <- performance(pred, "prec", "rec")
#   plot(RP.perf)
#   # ROC curve
#   ROC.perf <- performance(pred, "tpr", "fpr")
#   plot (ROC.perf)
#   # ROC area under the curve
#   auc.tmp <- performance(pred, "auc")
#   auc <- as.numeric(auc.tmp@y.values)
#   # F1 score
#   print("F1 Score")
#   f1 <- performance(pred, "f")
#   plot(performance(pred, "f"))
#   # auc
#   print("AUC")
#   out <- roc(y, classifierPredictions, ci = TRUE)
#   round(as.numeric(out$auc), 4)
# }

#######################################################
# Calculate stats given:
#######################################################
# input: predictions + expected labels (1=causal or 0=non-causal)
#######################################################
getStats <- function(classifierPredictions, casecontrol) {
  y = as.factor(casecontrol)
  pred <- prediction(classifierPredictions, y)
  # Recall-Precision curve
  RP.perf <- performance(pred, "prec", "rec")
  print(RP.perf)
  plot(RP.perf)
  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr")
  plot (ROC.perf)
  # ROC area under the curve
  auc.tmp <- performance(pred, "auc")
  print("auc.tmp")
  print(auc.tmp)
  auc <- as.numeric(auc.tmp@y.values)
  # F1 score
  print("F1 Score")
  print(performance(pred, "f"))
  plot(performance(pred, "f"))
  # auc
  #print("AUC")
  out <- roc(y, classifierPredictions, ci = TRUE)
  round(as.numeric(out$auc), 4)
}
# 
# #######################################################
# # Calculate aggregated stats given:
# #######################################################
# # input: data + hoi name
# #######################################################
# aggAuc <- function(dat, hoi) {
#   dat <- subset(dat, hoiname == hoi)
#   aucRor <- getStats(as.numeric(dat$ror), dat$casecontrol)
#   aucChisq <- getStats(as.numeric(dat$chisq), dat$casecontrol)
#   #getStats(as.numeric(dat$medcoef0))
#   aucSql5MedCoef <-
#     getStats(as.numeric(dat$sql5medcoef1), dat$casecontrol)
#   aucPsi5MedCoef <-
#     getStats(as.numeric(dat$psi5medcoef1), dat$casecontrol)
#   aucSql10MedCoef <-
#     getStats(as.numeric(dat$sql10medcoef1), dat$casecontrol)
#   aucPsi10MedCoef <-
#     getStats(as.numeric(dat$psi10medcoef1), dat$casecontrol)
#   aucSql5ate <- getStats(as.numeric(dat$sql5ate), dat$casecontrol)
#   aucPsi5ate <- getStats(as.numeric(dat$psi5ate), dat$casecontrol)
#   aucSql10ate <- getStats(as.numeric(dat$sql10ate), dat$casecontrol)
#   aucPsi10ate <- getStats(as.numeric(dat$psi10ate)/as.numeric(dat$psi10ate), dat$casecontrol)
#   out <-
#     paste(
#       hoi,
#       aucChisq,
#       aucRor,
#       aucSql5MedCoef,
#       aucPsi5MedCoef,
#       aucSql10MedCoef,
#       aucPsi10MedCoef,
#       aucSql5ate,
#       aucPsi5ate,
#       aucSql10ate,
#       aucPsi10ate,
#       sep = " & "
#     )
#   out
# }
# aggAuc(dat, "kidney_failure,_acute")
# aggAuc(dat, "liver_failure,_acute")
# aggAuc(dat, "acute_myocardial_infarction")
# aggAuc(dat, "gastrointestinal_hemorrhage")

aggAuc2 <- function(hoi) {
  dat <- subset(dat, hoiname == hoi)
  aucRor <- getStats(as.numeric(dat$ror), dat$casecontrol)
  aucChisq <- getStats(as.numeric(dat$chisq), dat$casecontrol)
  #getStats(as.numeric(dat$medcoef0))
  aucSql5MedCoef <-
    getStats(as.numeric(dat$sql5medcoef1), dat$casecontrol)
  aucPsi5MedCoef <-
    getStats(as.numeric(dat$psi5medcoef1), dat$casecontrol)
  aucSql10MedCoef <-
    getStats(as.numeric(dat$sql10medcoef1), dat$casecontrol)
  aucPsi10MedCoef <-
    getStats(as.numeric(dat$psi10medcoef1), dat$casecontrol)
  aucSql5ate <- getStats(as.numeric(dat$sql5ate), dat$casecontrol)
  aucPsi5ate <- getStats(as.numeric(dat$psi5ate), dat$casecontrol)
  aucSql10ate <- getStats(as.numeric(dat$sql10ate), dat$casecontrol)
  aucPsi10ate <- getStats(as.numeric(dat$psi10ate), dat$casecontrol)
  nrowDat <- nrow(dat)
  out <-
    c(
      hoi,
      aucChisq,
      aucRor,
      aucSql5MedCoef,
      aucPsi5MedCoef,
      aucSql10MedCoef,
      aucPsi10MedCoef,
      aucSql5ate,
      aucPsi5ate,
      aucSql10ate,
      aucPsi10ate,
      nrowDat
    )
  out
}

#out <- lapply("kidney_failure,_acute", aggAuc2)

dat.orig <- dat


#dat2 <- read_delim("resultsData-dib.tsv", delim = "\t")
#system("head -n 3 resultsData-dib.tsv")
#dat.orig <- dat

weightedCombinedAuc <- function(dat) {
  ades <-
    c(
      "kidney_failure,_acute",
      #"liver_failure,_acute",
      "acute_myocardial_infarction",
      "gastrointestinal_hemorrhage"
    )
  dat <- subset(dat, hoiname %in% ades)
  out <- lapply(ades, aggAuc2)
  #print("out")
  #print(out)
  nrowTot <- nrow(dat)
  weightedCombinedScores <- c()
  for (k in 2:11) {
    weights <- c()
    scores <- c()
    for (a in 1:length(ades)) {
      weight <- as.numeric(unlist(out[a])[12]) / nrowTot #12
      #print(paste("WEIGHT: ", weight))
      score <- as.numeric(unlist(out[a])[k])
      weights <- c(weights, weight)
      scores <- c(scores, score)
    }
    #print(scores)
    #print(weights)
    weightedCombinedScore <- round(sum(weights * scores), 4)
    weightedCombinedScores <-
      c(weightedCombinedScores, weightedCombinedScore)
  }
  combinedAUC <- paste0(weightedCombinedScores, collapse = " & ")
  combinedAUC
}

nork <- weightedCombinedAuc(dat)

###############
##### PRAUC!!!!!!!!!
###############

aggPrauc <- function(dat, hoi) {
  dat <- subset(dat, hoiname == hoi)
  praucRor <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$ror)
    ), 4)
  praucChisq <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$chisq)
    ), 4)
  #MLmetrics::PRAUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$medcoef0))
  praucSql5MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql5medcoef1)
    ), 4)
  praucPsi5MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi5medcoef1)
    ), 4)
  praucSql10MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql10medcoef1)
    ), 4)
  praucPsi10MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi10medcoef1)
    ), 4)
  praucSql5ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql5ate)
    ), 4)
  praucPsi5ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi5ate)
    ), 4)
  praucSql10ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql10ate)
    ), 4)
  praucPsi10ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi10ate)
    ), 4)
  out <-
    data.frame(
      hoi,
      praucChisq,
      praucRor,
      praucSql5MedCoef,
      praucPsi5MedCoef,
      praucSql10MedCoef,
      praucPsi10MedCoef,
      praucSql5ate,
      praucPsi5ate,
      praucSql10ate,
      praucPsi10ate#,
      #sep = " & "
    )
  out
 # xtable(out)
}

aggPrauc(dat, "kidney_failure,_acute")
aggPrauc(dat, "liver_failure,_acute")
aggPrauc(dat, "acute_myocardial_infarction")
aggPrauc(dat, "gastrointestinal_hemorrhage")

aggPrauc2 <- function(hoi) {
  dat <- subset(dat, hoiname == hoi)
  praucRor <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$ror)
    ), 4)
  praucChisq <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$chisq)
    ), 4)
  #MLmetrics::PRAUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$medcoef0))
  praucSql5MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql5medcoef1)
    ), 4)
  praucPsi5MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi5medcoef1)
    ), 4)
  praucSql10MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql10medcoef1)
    ), 4)
  praucPsi10MedCoef <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi10medcoef1)
    ), 4)
  praucSql5ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql5ate)
    ), 4)
  praucPsi5ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi5ate)
    ), 4)
  praucSql10ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$sql10ate)
    ), 4)
  praucPsi10ate <-
    round(MLmetrics::PRAUC(
      y_true = as.factor(dat$casecontrol),
      y_pred = as.numeric(dat$psi10ate)
    ), 4)
  nrowDat <- nrow(dat)
  out <-
    c(
      hoi, #1
      praucChisq, #2
      praucRor, #3
      praucSql5MedCoef, #4
      praucPsi5MedCoef, #5
      praucSql10MedCoef, #8
      praucPsi10MedCoef, #9
      praucSql5ate, #6
      praucPsi5ate, #7
      praucSql10ate, #10
      praucPsi10ate, #11
      nrowDat #12
    )
  out
}

weightedCombinedPrauc <- function(dat) {
  ades <-
    c(
      "kidney_failure,_acute",
      "liver_failure,_acute",
      "acute_myocardial_infarction",
      "gastrointestinal_hemorrhage"
    )
  out <- lapply(ades, aggPrauc2)
  dat <- subset(dat, hoiname %in% ades)
  nrowTot <- nrow(dat)
  weightedCombinedScores <- c()
  for (k in 2:11) {
    weights <- c()
    scores <- c()
    for (a in 1:length(ades)) {
      weight <- as.numeric(unlist(out[a])[12]) / nrowTot #12
      print(paste("weight: ", weight))
      score <- as.numeric(unlist(out[a])[k])
      weights <- c(weights, weight)
      scores <- c(scores, score)
    }
    print(scores)
    print(weights)
    weightedCombinedScore <- round(sum(weights * scores), 4)
    weightedCombinedScores <-
      c(weightedCombinedScores, weightedCombinedScore)
  }
  paste0(weightedCombinedScores, collapse = " & ")
  #data.frame(weightedCombinedScores)
}

weightedCombinedPrauc(dat)

############
##### Mean Average Precision at K (MAP-K)
############

getMapk <- function(k, Yactual, Ypredicted) {
  Yactual <- as.numeric(Yactual)
  
  #if (sum(Yactual[1:k]) > 0) {
  if (k > length(Yactual)) {
    k = length(Yactual)
  }
  #  if (sum(as.numeric(dat$casecontrol[1:k])) > 0) {
  #print(Ypredicted[order(Ypredicted, decreasing = FALSE)][1:k])
  Yactual_sort <- na.omit(Yactual)
  # The actual Y values are sorted by predicted Y values in descending order
  Yactual_sort <- Yactual[order(Ypredicted, decreasing = TRUE)]
  Yactual_sort <- Yactual_sort[1:k]
  Yactual_sort <- na.omit(Yactual_sort)
  #print(Yactual_sort)
  mapk <-
    sum(cumsum(Yactual_sort) * Yactual_sort / seq_along(Yactual_sort)) / sum(Yactual_sort)
  out <- round(mapk, 4)
  #} else if (sum(Yactual[1:k]) == 0) { out <- 0 }
  out
}

#dat <- dat.orig
dat.orig <- dat

aggMapk <- function(dat, k, hoi) {
  dat <- data.frame(dat)
  # hoi = "gastrointestinal_hemorrhage"
  dat <- subset(dat, hoiname == hoi)
  #dat <- na.omit(dat)
  mapkRor <- getMapk(k, dat$casecontrol, dat$ror)
  mapkChisq <- getMapk(k, dat$casecontrol, dat$chisq)
  mapkSql5MedCoef <- getMapk(k, dat$casecontrol, dat$sql5medcoef1)
  mapkPsi5MedCoef <- getMapk(k, dat$casecontrol, dat$psi5medcoef1)
  mapkSql5ate <- getMapk(k, dat$casecontrol, dat$sql5ate)
  mapkPsi5ate <- getMapk(k, dat$casecontrol, dat$psi5ate)
  mapkSql10MedCoef <- getMapk(k, dat$casecontrol, dat$sql10medcoef1)
  mapkPsi10MedCoef <- getMapk(k, dat$casecontrol, dat$psi10medcoef1)
  mapkSql10ate <- getMapk(k, dat$casecontrol, dat$sql10ate)
  mapkPsi10ate <- getMapk(k, dat$casecontrol, dat$psi10ate)
  out <-
    paste(
      hoi,
      mapkChisq,
      mapkRor,
      mapkSql5MedCoef,
      mapkPsi5MedCoef,
      mapkSql10MedCoef,
      mapkPsi10MedCoef,
      mapkSql5ate,
      mapkPsi5ate,
      mapkSql10ate,
      mapkPsi10ate,
      sep = " & "
    )
  out
}

aggMapk(dat, 10, "kidney_failure,_acute")
aggMapk(dat, 10, "liver_failure,_acute")
aggMapk(dat, 10, "acute_myocardial_infarction")
aggMapk(dat, 10, "gastrointestinal_hemorrhage")

aggMapk(dat, 25, "kidney_failure,_acute")
aggMapk(dat, 25, "liver_failure,_acute")
aggMapk(dat, 25, "acute_myocardial_infarction")
aggMapk(dat, 25, "gastrointestinal_hemorrhage")




dat.orig <- dat

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
## png()
#par(mfrow = c(1, 2))
## dev.off()

dat.new <- dat
dat.new$ror <- remove_outliers(dat$ror)
dat.new$medcoef0 <- remove_outliers(dat$medcoef0)
dat.new$sql5medcoef1 <- remove_outliers(dat$sql5medcoef1)
dat.new$sql10medcoef1 <- remove_outliers(dat$sql10medcoef1)
dat.new$psi5medcoef1 <- remove_outliers(dat$psi5medcoef1)
dat.new$sql10psicoef1 <- remove_outliers(dat$psi10medcoef1)

dat.new$sql5ate <- remove_outliers(dat$sql5ate)
dat.new$sql10ate <- remove_outliers(dat$sql10ate)
dat.new$psi5ate <- remove_outliers(dat$psi5ate)
dat.new$psi10ate <- remove_outliers(dat$psi10ate)
dat.new.orig <- dat.new
dat.new <- na.omit(dat.new)

dat.new$sql5ate <- abs(remove_outliers(dat$sql5ate))
dat.new$sql10ate <- abs(remove_outliers(dat$sql10ate))
dat.new$psi5ate <- abs(remove_outliers(dat$psi5ate))
dat.new$psi10ate <- abs(remove_outliers(dat$psi10ate))
dat.new.orig1 <- dat.new
dat.new <- na.omit(dat.new)

dat.new$sql5ate <- abs(dat$sql5ate)
dat.new$sql10ate <- abs(dat$sql10ate)
dat.new$psi5ate <- abs(dat$psi5ate)
dat.new$psi10ate <- abs(dat$psi10ate)
dat.new.orig1 <- dat.new
dat.new <- na.omit(dat.new)

nrow(dat.new)
dat.new.good <- dat.new
dat.new <- dat
dat.new <- dat.new.good
dat.orig <- dat
dat <- dat.new
system("mkdir boxplots")

###################
# BOXPLOTS (literature-derived)
###################

# BASELINES (lit-derived)
getGGBoxplotBaseline <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_1_baselines.tiff", sep = "")  
  plotname <- paste("Baselines for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
            x = "casecontrol", y = c("ror", "medcoef0"), 
            color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
            ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
            combine = TRUE, add = c("mean_ci", "dotplot"), 
            fill = "light yellow", label = "exposurename",
            label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotBaseline(dat, "kidney_failure,_acute")
getGGBoxplotBaseline(dat, "liver_failure,_acute")
getGGBoxplotBaseline(dat, "acute_myocardial_infarction")
getGGBoxplotBaseline(dat, "gastrointestinal_hemorrhage")

# RAW REGRESSION (lit-derived)
getGGBoxplotRawcoefs <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_2_rawcoefs.tiff", sep = "")  
  plotname <- paste("Regression coefficients (β) using LBD for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"),
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Regression coefficients (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawcoefs(dat, "kidney_failure,_acute")
getGGBoxplotRawcoefs(dat, "liver_failure,_acute")
getGGBoxplotRawcoefs(dat, "acute_myocardial_infarction")
getGGBoxplotRawcoefs(dat, "gastrointestinal_hemorrhage")

# ATE (lit-derived)
getGGBoxplotRawATE<- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_3_ate.tiff", sep = "")  
  plotname <- paste("Average Treatment Effect (Δ) using LBD for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "ATE (Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawATE(dat, "kidney_failure,_acute")
getGGBoxplotRawATE(dat, "liver_failure,_acute")
getGGBoxplotRawATE(dat, "acute_myocardial_infarction")
getGGBoxplotRawATE(dat, "gastrointestinal_hemorrhage")

##################
##################
# BOXPLOTS (RANDOMLAND)
##################


# BASELINES  (RANDOMLAND)
getGGBoxplotBaseline <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_1_baselineRand.tiff", sep = "")  
  plotname <- paste("Baselines [random covar subset] for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("chisqPvalue", "medCoef0"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotBaseline(randomAKI, "kidney_failure,_acute")
getGGBoxplotBaseline(randomALI, "liver_failure,_acute")
getGGBoxplotBaseline(randomAMI, "acute_myocardial_infarction")
getGGBoxplotBaseline(randomGIB, "gastrointestinal_hemorrhage")

# RAW REGRESSION  (RANDOMLAND)
getGGBoxplotRawcoefs <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_1_rawcoefsRand.tiff", sep = "")  
  plotname <- paste("Regression coefficients (β) [random covar subset]  for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("medCoef0", "medCoef"), #, "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"),
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Exposure coefficients (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawcoefs(randomAKI, "kidney_failure,_acute")
getGGBoxplotRawcoefs(randomALI, "liver_failure,_acute")
getGGBoxplotRawcoefs(randomAMI, "acute_myocardial_infarction")
getGGBoxplotRawcoefs(randomGIB, "gastrointestinal_hemorrhage")

# ATE  (RANDOMLAND)
getGGBoxplotRawATE<- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_3_ateRand.tiff", sep = "")  
  plotname <- paste("Average Treatment Effect (Δ)  [random covar subset] for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("exactATE"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "ATE (Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawATE(randomAKI, "kidney_failure,_acute")
getGGBoxplotRawATE(randomALI, "liver_failure,_acute")
getGGBoxplotRawATE(randomAMI, "acute_myocardial_infarction")
getGGBoxplotRawATE(randomGIB, "gastrointestinal_hemorrhage")



###########################
###########################

##################
##################
# BOXPLOTS (COMBINED)
##################
combinedDat.orig <- merge(x = dat, y = randomALL, by.x = c("exposurename", "hoiname", "casecontrol"), by.y = c("exposurename", "hoiname", "casecontrol"))

combinedDat <- combinedDat.orig
combinedDat$ror <- remove_outliers(combinedDat$ror)
combinedDat$ror <- remove_outliers(combinedDat$medCoef0)
combinedDat <- na.omit(combinedDat)
# BASELINES  (combined)
getGGBoxplotBaseline <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_1_baselineCombined.tiff", sep = "")  
  plotname <- paste("Baselines [combined] for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("ror", "medCoef0"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotBaseline(combinedDat, "kidney_failure,_acute")
getGGBoxplotBaseline(combinedDat, "liver_failure,_acute")
getGGBoxplotBaseline(combinedDat, "acute_myocardial_infarction")
getGGBoxplotBaseline(combinedDat, "gastrointestinal_hemorrhage")


combinedDat <- combinedDat.orig
combinedDat$medCoef0 <- remove_outliers(combinedDat$medCoef0)
combinedDat$medCoef <- remove_outliers(combinedDat$medCoef)
combinedDat$sql5medcoef1 <- remove_outliers(combinedDat$sql5medcoef1)
combinedDat$sql10medcoef1 <- remove_outliers(combinedDat$sql10medcoef1)
combinedDat$psi5medcoef1 <- remove_outliers(combinedDat$psi5medcoef1)
combinedDat$psi10medcoef1 <- remove_outliers(combinedDat$psi10medcoef1)
combinedDat <- na.omit(combinedDat)

# RAW REGRESSION  (combined)
getGGBoxplotRawcoefs <- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_1_rawcoefsCombined.tiff", sep = "")  
  plotname <- paste("Regression coefficients (β) [combined]  for ", "GIB", sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y =  c("medCoef0", "medCoef", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), #, "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"),
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Exposure coefficients (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawcoefs(combinedDat, "kidney_failure,_acute")
getGGBoxplotRawcoefs(combinedDat, "liver_failure,_acute")
getGGBoxplotRawcoefs(combinedDat, "acute_myocardial_infarction")
getGGBoxplotRawcoefs(combinedDat, "gastrointestinal_hemorrhage")

# ATE  (combined)
combinedDat <- combinedDat.orig
#combinedDat$ <- remove_outliers(combinedDat$medCoef0)
combinedDat$exactATE <- remove_outliers(combinedDat$medCoef)
combinedDat$exactATE <- combinedDat$exactATE/500
combinedDat$sql5ate <- remove_outliers(combinedDat$sql5ate)
combinedDat$sql10ate <- remove_outliers(combinedDat$sql10ate)
combinedDat$psi5ate <- remove_outliers(combinedDat$psi5ate)
combinedDat$psi10ate <- remove_outliers(combinedDat$psi10ate)
combinedDat <- na.omit(combinedDat)

getGGBoxplotRawATE<- function(dat, hoiname) {
  fn <- paste("boxplots/", hoiname, "_3_ateCombined.tiff", sep = "")  
  plotname <- paste("Average Treatment Effect (Δ)  [combined] for ", hoiname, sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoiname), title = plotname, 
                 x = "casecontrol", y = c("exactATE", "sql5ate", "psi10ate"), #y = c("exactATE", "sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "ATE (Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), 
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawATE(combinedDat, "kidney_failure,_acute")
getGGBoxplotRawATE(combinedDat, "liver_failure,_acute")
getGGBoxplotRawATE(combinedDat, "acute_myocardial_infarction")
getGGBoxplotRawATE(combinedDat, "gastrointestinal_hemorrhage")

combinedDat <- combinedDat.orig
combinedDat$sql5medcoef1 <- remove_outliers(combinedDat$sql5medcoef1)
combinedDat$sql10medcoef1 <- remove_outliers(combinedDat$sql10medcoef1)
combinedDat$psi5medcoef1 <- remove_outliers(combinedDat$psi5medcoef1)
combinedDat$psi10medcoef1 <- remove_outliers(combinedDat$psi10medcoef1)
combinedDat <- na.omit(combinedDat)
combinedDat$coefDiffRandom <- combinedDat$medCoef0 - combinedDat$medCoef
combinedDat$coefDiffsql5 <- combinedDat$medCoef0 - combinedDat$sql5medcoef1
combinedDat$coefDiffsql10 <- combinedDat$medCoef0 - combinedDat$sql10medcoef1
combinedDat$coefDiffpsi5 <- combinedDat$medCoef0 - combinedDat$psi5medcoef1
combinedDat$coefDiffpsi10 <- combinedDat$medCoef0 - combinedDat$psi10medcoef1

getGGBoxplotRawATE<- function(dat, hoi) {
  fn <- paste("boxplots/", hoi, "_4_coefDiff.tiff", sep = "")  
  plotname <- paste("Difference in β [combinedCoefDiff] for ", "GIB", sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  p <- ggboxplot(subset(dat, dat$hoiname == hoi), title = plotname, 
                 x = "casecontrol", y = c("coefDiffRandom", "coefDiffsql5", "coefDiffpsi10"), # c("coefDiffRandom", "coefDiffsql5", "coefDiffsql10", "coefDiffpsi5", "coefDiffpsi10"), 
                 color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
                 ylab = "Difference in β", xlab = "casecontrol", bxp.errorbar = TRUE, 
                 combine = TRUE, add = c("mean_ci", "dotplot"), font.label = 5,
                 fill = "light yellow", label = "exposurename",
                 label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
  plot(p)
  dev.off()
}

getGGBoxplotRawATE(combinedDat, "kidney_failure,_acute")
getGGBoxplotRawATE(combinedDat, "liver_failure,_acute")
getGGBoxplotRawATE(combinedDat, "acute_myocardial_infarction")
getGGBoxplotRawATE(combinedDat, "gastrointestinal_hemorrhage")


#c("coefDiffRandom", "coefDiffsql5", "coefDiffsql10", "coefDiffpsi5", "coefDiffpsi10"), 

aggAuc2 <- function(hoi, dat) {
  dat <- subset(dat, hoiname == hoi)
  aucRor <- getStats(as.numeric(dat$ror), dat$casecontrol)
  aucChisq <- getStats(as.numeric(dat$chisq), dat$casecontrol)
  aucMedcoef0 <- getStats(as.numeric(dat$medcoef0), dat$casecontrol)
  #getStats(as.numeric(dat$medcoef0))
  aucCoefDiffRandom <- getStats(as.numeric(dat$coefDiffRandom), dat$casecontrol)
  aucCoefDiffSql5 <-
    getStats(as.numeric(dat$coefDiffsql5), dat$casecontrol)
  aucCoefDiffPsi5 <-
    getStats(as.numeric(dat$coefDiffpsi5), dat$casecontrol)
  aucCoefDiffSql10 <-
    getStats(as.numeric(dat$coefDiffsql10), dat$casecontrol)
  aucCoefDiffPsi10 <-
    getStats(as.numeric(dat$coefDiffpsi10), dat$casecontrol)
  nrowDat <- nrow(dat)
  out <-
    c(
      hoi,
      aucChisq,
      aucRor,
      aucMedcoef0,
      aucCoefDiffRandom,
      aucCoefDiffSql5,
      aucCoefDiffPsi5,
      aucCoefDiffSql10,
      aucCoefDiffPsi10,
      nrowDat
    )
  out
}

aggAuc2("kidney_failure,_acute", combinedDat)
aggAuc2("liver_failure,_acute", combinedDat)
aggAuc2("acute_myocardial_infarction", combinedDat)
aggAuc2("gastrointestinal_hemorrhage", combinedDat)




aggAuc2 <- function(dat) {
  #dat <- subset(dat, hoiname == hoi)
  aucRor <- getStats(as.numeric(dat$ror), dat$casecontrol)
  aucChisq <- getStats(as.numeric(dat$chisq), dat$casecontrol)
  aucMedcoef0 <- getStats(as.numeric(dat$medcoef0), dat$casecontrol)
  #getStats(as.numeric(dat$medcoef0))
  aucCoefDiffRandom <- getStats(as.numeric(dat$coefDiffRandom), dat$casecontrol)
  aucCoefDiffSql5 <-
    getStats(as.numeric(dat$coefDiffsql5), dat$casecontrol)
  aucCoefDiffPsi5 <-
    getStats(as.numeric(dat$coefDiffpsi5), dat$casecontrol)
  aucCoefDiffSql10 <-
    getStats(as.numeric(dat$coefDiffsql10), dat$casecontrol)
  aucCoefDiffPsi10 <-
    getStats(as.numeric(dat$coefDiffpsi10), dat$casecontrol)
  nrowDat <- nrow(dat)
  out <-
    c(
      #hoi,
      aucChisq,
      aucRor,
      aucMedcoef0,
      aucCoefDiffRandom,
      aucCoefDiffSql5,
      aucCoefDiffSql10,
      aucCoefDiffPsi5,
      aucCoefDiffPsi10,
      nrowDat
    )
  out
}

aggAuc2(combinedDat)
aggAuc2("liver_failure,_acute", combinedDat)
aggAuc2("acute_myocardial_infarction", combinedDat)
aggAuc2("gastrointestinal_hemorrhage", combinedDat)

#c("coefDiffRandom", "coefDiffsql5", "coefDiffsql10", "coefDiffpsi5", "coefDiffpsi10"),
#baselineColor <- "#663500"; glmColor <- "#63a11b"; gcmColor <- "#ff845d" # https://medialab.github.io/iwanthue/

drawCombinedDiffAUC <- function(combinedDat, hoi) {
  combinedDat <- subset(combinedDat, hoiname == hoi)
  baselineColor <- "#495100" # dark green
  randomColor <-  "#590035" # dark purple
  sql5Color <- "#dbd976" # light yellowish green
  sql10Color <-  "#b453c4" # light purple
  psi5Color <-  "#e5c57d" # puke yellow
  psi10Color <-  "#ff736e" # orange/pink
  fn <- paste(hoi, "_aucCombinedDIFF.tiff", sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  plot(roc(combinedDat$casecontrol, combinedDat$ror, smooth = FALSE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffRandom, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
                         col = randomColor, print.auc.y = .4, add = TRUE, colorize = FALSE, pch=18, type = "b")
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffsql5, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
     col = sql5Color, print.auc.y = .3, add = TRUE, colorize = TRUE)

  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffsql10, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
     col = sql10Color, print.auc.y = .2, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffpsi5, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
     col = psi5Color, print.auc.y = .1, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffpsi10, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
                        col = psi10Color, print.auc.y = .0, add = TRUE, colorize = TRUE, pch=10, type = "b")
#  lines(combinedDat$coefDiffpsi10, y2, pch=18, col="blue", type="b", lty=2)
  legendColors <- c("#495100", "#590035", "#dbd976", "#b453c4", "#e5c57d", "#ff736e")
  legendLabels <- c("baseline", "randomCovariates", "sql5", "sql10", "psi5", "psi10")
  legend(1, 95, legend=legendLabels,
         col=legendColors, lty=1:2, cex=0.8)
  
  dev.off()
}


drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "kidney_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "liver_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "acute_myocardial_infarction")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "gastrointestinal_hemorrhage")




drawCombinedDiffAUC <- function(combinedDat, hoi) {
  combinedDat <- subset(combinedDat, hoiname == hoi)
  baselineColor <- "#495100" # dark green
  randomColor <-  "#590035" # dark purple
  psi10Color <-  "#ff736e" # orange/pink
  fn <- paste("randomDiffs/", hoi, "_aucCombinedDIFF.tiff", sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  plot(roc(combinedDat$casecontrol, combinedDat$ror, smooth = FALSE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffRandom, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = randomColor, print.auc.y = .4, add = TRUE, colorize = FALSE, pch=18, type = "b")
  #plot(roc(combinedDat$casecontrol, combinedDat$coefDiffsql5, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
  #     col = sql5Color, print.auc.y = .3, add = TRUE, colorize = TRUE)
  #plot(roc(combinedDat$casecontrol, combinedDat$coefDiffsql10, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
  #     col = sql10Color, print.auc.y = .2, add = TRUE, colorize = TRUE)
  #plot(roc(combinedDat$casecontrol, combinedDat$coefDiffpsi5, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
  #     col = psi5Color, print.auc.y = .1, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$coefDiffpsi10, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = psi10Color, print.auc.y = .3, add = TRUE, colorize = TRUE, pch=10, type = "b")
  #  lines(combinedDat$coefDiffpsi10, y2, pch=18, col="blue", type="b", lty=2)
  legendColors <- c("#495100", "#590035", "#ff736e")
  legendLabels <- c("baseline", "randomCovariates", "psi10")
  legend(1, 95, legend=legendLabels,
         col=legendColors, lty=1:2, cex=0.8)
  
  dev.off()
}


drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "kidney_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "liver_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "acute_myocardial_infarction")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "gastrointestinal_hemorrhage")



drawCombinedDiffAUC <- function(combinedDat) {
  #combinedDat <- subset(combinedDat, hoiname == hoi)
  baselineColor <- "#495100" # dark green
  randomColor <-  "#590035" # dark purple
  psi10Color <-  "#ff736e" # orange/pink
  fn <- paste("all", "_aucCombinedDIFF.tiff", sep = "")
  tiff(file = fn, height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
  plot(roc(combinedDat$casecontrol, combinedDat$ror, smooth = FALSE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$exactATE, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = randomColor, print.auc.y = .4, add = TRUE, colorize = FALSE, pch=18, type = "b")
  plot(roc(combinedDat$casecontrol, combinedDat$sql5medcoef1, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = sql5Color, print.auc.y = .3, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$sql10medcoef1, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = sql10Color, print.auc.y = .2, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$psi5medcoef1, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = psi5Color, print.auc.y = .1, add = TRUE, colorize = TRUE)
  plot(roc(combinedDat$casecontrol, combinedDat$psi10medcoef1, ci = TRUE, algorithm = 4, smooth = FALSE), print.auc = TRUE,
       col = psi10Color, print.auc.y = .0, add = TRUE, colorize = TRUE, pch=10, type = "b")
  #  lines(combinedDat$coefDiffpsi10, y2, pch=18, col="blue", type="b", lty=2)
  #legendColors <- c("#495100", "#590035", "#ff736e")
  #legendLabels <- c("baseline", "randomCovariates", "psi10")
  #legend(1, 95, legend=legendLabels,
  #       col=legendColors, lty=1:2, cex=0.8)
  
  dev.off()
}


drawCombinedDiffAUC(combinedDat = combinedDat)#, hoi = "kidney_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "liver_failure,_acute")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "acute_myocardial_infarction")
drawCombinedDiffAUC(combinedDat = combinedDat, hoi = "gastrointestinal_hemorrhage")






ggboxplot(data.frame(randomAKI), x = "casecontrol", y = "chisqPvalue")


#randomAKI$chisqPvalue <- remove_outliers(randomAKI$chisqPvalue)
#randomAKI$exactATE <- remove_outliers(randomAKI$exactATE)
randomAKI$medCoef0 <- remove_outliers(randomAKI$medcoef0)
randomAKI$medCoef <- remove_outliers(randomAKI$medCoef)
randomAKI <- na.omit(randomAKI)
ggboxplot(data.frame(randomAKI), 
          title = "random (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("medCoef0", "medCoef"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 


randomALI$medCoef0 <- remove_outliers(randomALI$medcoef0)
randomALI$medCoef <- remove_outliers(randomALI$medCoef)
randomALI <- na.omit(randomALI)
ggboxplot(data.frame(randomALI), 
          title = "random (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("medCoef0", "medCoef"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomAMI$medCoef0 <- remove_outliers(randomAMI$medcoef0)
randomAMI$medCoef <- remove_outliers(randomAMI$medCoef)
randomAMI <- na.omit(randomAMI)
ggboxplot(data.frame(randomAMI), 
          title = "random (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("medCoef0", "medCoef"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 


randomGIB$medCoef0 <- remove_outliers(randomGIB$medcoef0)
randomGIB$medCoef <- remove_outliers(randomGIB$medCoef)
randomGIB <- na.omit(randomGIB)
ggboxplot(data.frame(randomGIB), 
          title = "random (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("medCoef0", "medCoef"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 


### exact ATE with random

randomAKI$exactATE <- remove_outliers(randomAKI$exactATE)
randomAKI <- na.omit(randomAKI)
ggboxplot(data.frame(randomAKI), x = "casecontrol", y = "medCoef0")
ggboxplot(data.frame(randomAKI), 
          title = "ATE using randomly selected covariates for AKI",
          x = "casecontrol", y = c("exactATE"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "ATE", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomALI$exactATE <- remove_outliers(randomALI$exactATE)
randomALI <- na.omit(randomALI)
ggboxplot(data.frame(randomALI), x = "casecontrol", y = "medCoef0")
ggboxplot(data.frame(randomALI), 
          title = "ATE using randomly selected covariates for ALI",
          x = "casecontrol", y = c("exactATE"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "ATE", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomAMI$exactATE <- remove_outliers(randomAMI$exactATE)
randomAMI <- na.omit(randomAMI)
ggboxplot(data.frame(randomAMI), x = "casecontrol", y = "medCoef0")
ggboxplot(data.frame(randomAMI), 
          title = "ATE using randomly selected covariates for AMI",
          x = "casecontrol", y = c("exactATE"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "ATE", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomGIB$exactATE <- remove_outliers(randomGIB$exactATE)
randomGIB <- na.omit(randomGIB)
ggboxplot(data.frame(randomGIB), x = "casecontrol", y = "medCoef0")
ggboxplot(data.frame(randomGIB), 
          title = "ATE using randomly selected covariates for GIB",
          x = "casecontrol", y = c("medCoef0"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "ATE", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE)
ggboxplot(data.frame(randomGIB), 
          title = "ATE using randomly selected covariates for GIB",
          x = "casecontrol", y = c("exactATE"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "ATE", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomGIB$medCoef <- remove_outliers(randomGIB$medCoef)
randomGIB <- na.omit(randomGIB)
ggboxplot(data.frame(randomAMI), 
          title = "random (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("exactATE"), #, "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 


ggboxplot(data.frame(randomAKI), 
          x = "casecontrol", y = "exactATE", 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 1, top.down = 1)) 

ggboxplot(data.frame(randomALI), 
          x = "casecontrol", y = "exactATE", 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 1, top.down = 1)) 

ggboxplot(data.frame(randomAMI), 
          x = "casecontrol", y = "exactATE", 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 1, top.down = 1)) 

ggboxplot(data.frame(randomGIB), 
          x = "casecontrol", y = "exactATE", 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 1, top.down = 1)) 

randomAE <- rbind(randomAKI, randomALI, randomAMI, randomGIB)



ggboxplot(data.frame(randomAKI), 
          x = "casecontrol", y = c("medCoef0", "medCoef", "exactATE"), #, "networkScoreDiff"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "score", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(data.frame(randomAKI), 
          x = "casecontrol", y = c("networkScoreDiff"), #, "networkScoreDiff"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "score", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(data.frame(randomALI), 
          x = "casecontrol", y = c("medCoef0", "medCoef", "exactATE"), #, "networkScoreDiff"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "score", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(data.frame(randomAMI), 
          x = "casecontrol", y = c("medCoef0", "medCoef", "exactATE"), #, "networkScoreDiff"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "score", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(data.frame(randomGIB), 
          x = "casecontrol", y = c("medCoef0", "medCoef", "exactATE"), #, "networkScoreDiff"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "score", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

randomAE$medCoef0


ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 


ggboxplot(dat.new, x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(dat.new, x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(dat.new, x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 




# 'kidney_failure,_acute'
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = "median_iqr", fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = "median_iqr", fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = "median_iqr", fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
# 'liver_failure,_acute'
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
# 'acute_myocardial_infarction'
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
# 'gastrointestinal_hemorrhage'
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 





ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 1, top.down = 1)) 

label.select = list(top.up = 10, top.down = 4)

# 'kidney_failure,_acute'
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = "median_iqr")
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
ggboxplot(subset(dat, dat$hoiname == 'kidney_failure,_acute'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
# 'liver_failure,_acute'
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE)
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
ggboxplot(subset(dat, dat$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
# 'acute_myocardial_infarction'
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE)
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
ggboxplot(subset(dat, dat$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
# 'gastrointestinal_hemorrhage'
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE)
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
ggboxplot(subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, fill = "yellow", facet.by = "hoiname") 


ggboxplot(dat, x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", 
          palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", 
          bxp.errorbar = TRUE, combine = TRUE, facet.by = "hoiname")
ggboxplot(dat, x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE) 
ggboxplot(dat, x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", 
          bxp.errorbar = TRUE, combine = TRUE, fill = "yellow", facet.by = "hoiname") 

ggboxplot(dat, x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", 
          bxp.errorbar = TRUE, combine = TRUE, fill = "yellow", facet.by = "hoiname", 
          select = subset(dat, dat$sql5ate > .1)$exposurename)


ggboxplot(dat, 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", 
          bxp.errorbar = TRUE, combine = TRUE, fill = "yellow")#, facet.by = "hoiname") 




############################


round(MLmetrics::PRAUC(
  y_true = as.factor(randomAKI$casecontrol),
  y_pred = as.numeric(randomAKI$medCoef0)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomALI$casecontrol),
  y_pred = as.numeric(randomALI$medCoef)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomAKI$casecontrol),
  y_pred = as.numeric(randomAKI$exactATE)
), 4)

round(MLmetrics::PRAUC(
  y_true = as.factor(randomALI$casecontrol),
  y_pred = as.numeric(randomALI$medCoef0)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomALI$casecontrol),
  y_pred = as.numeric(randomALI$medCoef)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomALI$casecontrol),
  y_pred = as.numeric(randomALI$exactATE)
), 4)

round(MLmetrics::PRAUC(
  y_true = as.factor(randomAMI$casecontrol),
  y_pred = as.numeric(randomAMI$medCoef0)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomAMI$casecontrol),
  y_pred = as.numeric(randomAMI$medCoef)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomAMI$casecontrol),
  y_pred = as.numeric(randomAMI$exactATE)
), 4)



round(MLmetrics::PRAUC(
  y_true = as.factor(randomGIB$casecontrol),
  y_pred = as.numeric(randomGIB$medCoef0)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomGIB$casecontrol),
  y_pred = as.numeric(randomGIB$medCoef)
), 4)
round(MLmetrics::PRAUC(
  y_true = as.factor(randomGIB$casecontrol),
  y_pred = as.numeric(randomGIB$exactATE)
), 4)




