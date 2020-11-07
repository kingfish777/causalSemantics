 #setwd("/Users/scottalexandermalec/Projects/enlil")
library(MLmetrics)
library(pROC)
library(PRROC)
library(readr)
library(xtable)
library(ROCR)
library(openxlsx)
library(RPostgreSQL)
library(ggpubr)
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

#######################################################
# GET (rectangular) DATA(set) TO BE ANALYZED
#######################################################
dat <-
  dbGetQuery(
    con,
    "SELECT DISTINCT sql5.casecontrol, sql5.hoiname, sql5.exposurename, 
                          ref.a, ref.b, ref.c, ref.d, ref.ror, sql5.medcoef0,
                          sql5.exactATE AS sql5ate, sql5.medcoef1 AS sql5medcoef1,
                          psi5.exactATE AS psi5ate, psi5.medcoef1 AS psi5medcoef1,
                          sql10.exactATE AS sql10ate, sql10.medcoef1 AS sql10medcoef1,
                          psi10.exactATE AS psi10ate, psi10.medcoef1 AS psi10medcoef1
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
names(dat) <- c("casecontrol", "hoiname", "exposurename", "a", "b", "c", "d", "ror", "medcoef0", "sql5ate", "sql5medcoef1", "psi5ate", "psi5medcoef1", "sql10ate", "sql10medcoef1", "psi10ate", "psi10medcoef1", "chisq")

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


randomAKI <- data.frame(read.csv("randomAKI.txt", header = FALSE, sep = "\t"))
randomALI <- data.frame(read.csv("randomALI.txt", header = FALSE, sep = "\t"))
randomAMI <- data.frame(read.csv("randomAMI.txt", header = FALSE, sep = "\t"))
randomGIB <- data.frame(read.csv("randomGIB.txt", header = FALSE, sep = "\t"))

columnNames <- c("hoiname", "cl", "casecontrol", "hoiname", "exposurename", "cl", "SQUELCH",  "chisqPvalue", "drug.ade.score0", "drug.ade.score", "exactATE", "exactATEvar",  "medCoef0", "medZvalue0", "medPvalue0", "medCoef", "medZvalue", "medPvalue",  "networkScore0", "networkScore1", "networkScoreDiff")
names(randomAKI) <- columnNames
names(randomALI) <- columnNames
names(randomAMI) <- columnNames
names(randomGIB) <- columnNames


#######################################################
# Calculate stats given:
#######################################################
# input: predictions + expected labels (1=causal or 0=non-causal)
#######################################################
getStats <- function(classifierPredictions, labels) {
  #classifierPredictions <- as.numeric(dat$sql5medcoef1); labels <- as.factor(dat$casecontrol)
  y = as.factor(labels)
  pred <- prediction(classifierPredictions, y)
  # Recall-Precision curve
  RP.perf <- performance(pred, "prec", "rec")
  plot(RP.perf)
  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr")
  plot (ROC.perf)
  # ROC area under the curve
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  # F1 score
  print("F1 Score")
  f1 <- performance(pred, "f")
  plot(performance(pred, "f"))
  # auc
  print("AUC")
  out <- roc(y, classifierPredictions, ci = TRUE)
  round(as.numeric(out$auc), 4)
}

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

#######################################################
# Calculate aggregated stats given:
#######################################################
# input: data + hoi name
#######################################################
aggAuc <- function(dat, hoi) {
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
  out <-
    paste(
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
      sep = " & "
    )
  out
}
aggAuc(dat, "kidney_failure,_acute")
aggAuc(dat, "liver_failure,_acute")
aggAuc(dat, "acute_myocardial_infarction")
aggAuc(dat, "gastrointestinal_hemorrhage")

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
      aucSql5ate,
      aucPsi5ate,
      aucSql10MedCoef,
      aucPsi10MedCoef,
      aucSql10ate,
      aucPsi10ate,
      nrowDat
    )
  out
}

#out <- lapply("kidney_failure,_acute", aggAuc2)

dat.orig <- dat


dat2 <- read_delim("resultsData-dib.tsv", delim = "\t")
#system("head -n 3 resultsData-dib.tsv")


dat.orig <- dat


dat.aki <- subset(dat, dat$hoiname == 'kidney_failure,_acute')
dat.ali <- subset(dat, dat$hoiname == 'liver_failure,_acute')
dat.ami <- subset(dat, dat$hoiname == 'acute_myocardial_infarction')
dat.gib <- subset(dat, dat$hoiname == 'gastrointestinal_hemorrhage')

# weighted tables
dat.aki <- subset(dat, dat$hoiname == 'kidney_failure,_acute')
dat.aki.orig <- dat.aki
dat.aki$ror <- dat.aki$ror/mean(dat.aki$ror)
dat.aki$medcoef0 <- dat.aki$medcoef0/mean(dat.aki$medcoef0)
dat.aki$sql5medcoef1 <- dat.aki$sql5medcoef1/mean(dat.aki$sql5medcoef1)
dat.aki$sql10medcoef1 <- dat.aki$sql10medcoef1/mean(dat.aki$sql10medcoef1)
dat.aki$psi5medcoef1 <- dat.aki$psi5medcoef1/mean(dat.aki$psi5medcoef1)
dat.aki$psi10medcoef1 <- dat.aki$psi10medcoef1/mean(dat.aki$psi10medcoef1)

dat.aki$sql5ate <- dat.aki$sql5ate/mean(dat.aki$sql5ate)
dat.aki$sql10ate <- dat.aki$sql10ate/mean(dat.aki$sql10ate)
dat.aki$psi5ate <- dat.aki$psi5ate/mean(dat.aki$psi5ate)
dat.aki$psi10ate <- dat.aki$psi10ate/mean(dat.aki$psi10ate)

## ali
dat.ali.orig <- dat.aki
dat.ali$ror <- dat.ali$ror/mean(dat.ali$ror)
dat.ali$medcoef0 <- dat.ali$medcoef0/mean(dat.ali$medcoef0)
dat.ali$sql5medcoef1 <- dat.ali$sql5medcoef1/mean(dat.ali$sql5medcoef1)
dat.ali$sql10medcoef1 <- dat.ali$sql10medcoef1/mean(dat.ali$sql10medcoef1)
dat.ali$psi5medcoef1 <- dat.ali$psi5medcoef1/mean(dat.ali$psi5medcoef1)
dat.ali$psi10medcoef1 <- dat.ali$psi10medcoef1/mean(dat.ali$psi10medcoef1)

dat.ali$sql5ate <- dat.ali$sql5ate/mean(dat.ali$sql5ate)
dat.ali$sql10ate <- dat.ali$sql10ate/mean(dat.ali$sql10ate)
dat.ali$psi5ate <- dat.ali$psi5ate/mean(dat.ali$psi5ate)
dat.ali$psi10ate <- dat.ali$psi10ate/mean(dat.ali$psi10ate)

## ali
dat.ami.orig <- dat.ami
dat.ami$ror <- dat.ami$ror/mean(dat.ami$ror)
dat.ami$medcoef0 <- dat.ami$medcoef0/mean(dat.ami$medcoef0)
dat.ami$sql5medcoef1 <- dat.ami$sql5medcoef1/mean(dat.ami$sql5medcoef1)
dat.ami$sql10medcoef1 <- dat.ami$sql10medcoef1/mean(dat.ami$sql10medcoef1)
dat.ami$psi5medcoef1 <- dat.ami$psi5medcoef1/mean(dat.ami$psi5medcoef1)
dat.ami$psi10medcoef1 <- dat.ami$psi10medcoef1/mean(dat.ami$psi10medcoef1)

dat.ami$sql5ate <- dat.ami$sql5ate/mean(dat.ami$sql5ate)
dat.ami$sql10ate <- dat.ami$sql10ate/mean(dat.ami$sql10ate)
dat.ami$psi5ate <- dat.ami$psi5ate/mean(dat.ami$psi5ate)
dat.ami$psi10ate <- dat.ami$psi10ate/mean(dat.ami$psi10ate)

## ali
dat.gib.orig <- dat.gib
dat.gib$ror <- dat.gib$ror/mean(dat.gib$ror)
dat.gib$medcoef0 <- dat.gib$medcoef0/mean(dat.gib$medcoef0)
dat.gib$sql5medcoef1 <- dat.gib$sql5medcoef1/mean(dat.gib$sql5medcoef1)
dat.gib$sql10medcoef1 <- dat.gib$sql10medcoef1/mean(dat.gib$sql10medcoef1)
dat.gib$psi5medcoef1 <- dat.gib$psi5medcoef1/mean(dat.gib$psi5medcoef1)
dat.gib$psi10medcoef1 <- dat.gib$psi10medcoef1/mean(dat.gib$psi10medcoef1)

dat.gib$sql5ate <- dat.gib$sql5ate/mean(dat.gib$sql5ate)
dat.gib$sql10ate <- dat.gib$sql10ate/mean(dat.gib$sql10ate)
dat.gib$psi5ate <- dat.gib$psi5ate/mean(dat.gib$psi5ate)
dat.gib$psi10ate <- dat.gib$psi10ate/mean(dat.gib$psi10ate)




### AKI 
dat.aki.summary <- c()
dat.aki.summary$ae <- "kidney_failure,_acute"
dat.aki.summary$ror <- auc(dat.aki$casecontrol, dat.aki$ror)
dat.aki.summary$chisq <- auc(dat.aki$casecontrol, dat.aki$chisq)
dat.aki.summary$medcoef0 <- auc(dat.aki$casecontrol, dat.aki$medcoef0)
dat.aki.summary$sql5medcoef1 <- auc(dat.aki$casecontrol, dat.aki$sql5medcoef1)
dat.aki.summary$sql10medcoef1 <- auc(dat.aki$casecontrol, dat.aki$sql10medcoef1)
dat.aki.summary$psi5medcoef1 <- auc(dat.aki$casecontrol, dat.aki$psi5medcoef1)
dat.aki.summary$psi10medcoef1 <- auc(dat.aki$casecontrol, dat.aki$psi10medcoef1)

dat.aki.summary$sql5ate <- auc(dat.aki$casecontrol, dat.aki$sql5ate)
dat.aki.summary$sql10ate <- auc(dat.aki$casecontrol, dat.aki$sql10ate)
dat.aki.summary$psi5ate <- auc(dat.aki$casecontrol, dat.aki$psi5ate)
dat.aki.summary$psi10ate <- auc(dat.aki$casecontrol, dat.aki$psi10ate)
dat.aki.summary <- as.data.frame(dat.aki.summary)
print(dat.aki.summary)

### ALI
dat.ali.summary <- c()
dat.ali.summary$ae <- "liver_failure,_acute"
dat.ali.summary$ror <- auc(dat.ali$casecontrol, dat.ali$ror)
dat.ali.summary$chisq <- auc(dat.ali$casecontrol, dat.ali$chisq)
dat.ali.summary$medcoef0 <- auc(dat.ali$casecontrol, dat.ali$medcoef0)
dat.ali.summary$sql5medcoef1 <- auc(dat.ali$casecontrol, dat.ali$sql5medcoef1)
dat.ali.summary$sql10medcoef1 <- auc(dat.ali$casecontrol, dat.ali$sql10medcoef1)
dat.ali.summary$psi5medcoef1 <- auc(dat.ali$casecontrol, dat.ali$psi5medcoef1)
dat.ali.summary$psi10medcoef1 <- auc(dat.ali$casecontrol, dat.ali$psi10medcoef1)

dat.ali.summary$sql5ate <- auc(dat.ali$casecontrol, dat.ali$sql5ate)
dat.ali.summary$sql10ate <- auc(dat.ali$casecontrol, dat.ali$sql10ate)
dat.ali.summary$psi5ate <- auc(dat.ali$casecontrol, dat.ali$psi5ate)
dat.ali.summary$psi10ate <- auc(dat.ali$casecontrol, dat.ali$psi10ate)
dat.ali.summary <- as.data.frame(dat.ali.summary)
print(dat.ali.summary)

### AMI
dat.ami.summary <- c()
dat.ami.summary$ae <- "acute_myocardial_infarction"
dat.ami.summary$ror <- auc(dat.ami$casecontrol, dat.ami$ror)
dat.ami.summary$chisq <- auc(dat.ami$casecontrol, dat.ami$chisq)
dat.ami.summary$medcoef0 <- auc(dat.ami$casecontrol, dat.ami$medcoef0)
dat.ami.summary$sql5medcoef1 <- auc(dat.ami$casecontrol, dat.ami$sql5medcoef1)
dat.ami.summary$sql10medcoef1 <- auc(dat.ami$casecontrol, dat.ami$sql10medcoef1)
dat.ami.summary$psi5medcoef1 <- auc(dat.ami$casecontrol, dat.ami$psi5medcoef1)
dat.ami.summary$psi10medcoef1 <- auc(dat.ami$casecontrol, dat.ami$psi10medcoef1)

dat.ami.summary$sql5ate <- auc(dat.ami$casecontrol, dat.ami$sql5ate)
dat.ami.summary$sql10ate <- auc(dat.ami$casecontrol, dat.ami$sql10ate)
dat.ami.summary$psi5ate <- auc(dat.ami$casecontrol, dat.ami$psi5ate)
dat.ami.summary$psi10ate <- auc(dat.ami$casecontrol, dat.ami$psi10ate)
dat.ami.summary <- as.data.frame(dat.ami.summary)
print(dat.ali.summary)

### GIB
dat.gib.summary <- c()
dat.gib.summary$ae <- "gastrointestinal_hemorrhage"
dat.gib.summary$ror <- auc(dat.gib$casecontrol, dat.gib$ror)
dat.gib.summary$chisq <- auc(dat.gib$casecontrol, dat.gib$chisq)
dat.gib.summary$medcoef0 <- auc(dat.gib$casecontrol, dat.gib$medcoef0)
dat.gib.summary$sql5medcoef1 <- auc(dat.gib$casecontrol, dat.gib$sql5medcoef1)
dat.gib.summary$sql10medcoef1 <- auc(dat.gib$casecontrol, dat.gib$sql10medcoef1)
dat.gib.summary$psi5medcoef1 <- auc(dat.gib$casecontrol, dat.gib$psi5medcoef1)
dat.gib.summary$psi10medcoef1 <- auc(dat.gib$casecontrol, dat.gib$psi10medcoef1)

dat.gib.summary$sql5ate <- auc(dat.gib$casecontrol, dat.gib$sql5ate)
dat.gib.summary$sql10ate <- auc(dat.gib$casecontrol, dat.gib$sql10ate)
dat.gib.summary$psi5ate <- auc(dat.gib$casecontrol, dat.gib$psi5ate)
dat.gib.summary$psi10ate <- auc(dat.gib$casecontrol, dat.gib$psi10ate)
dat.gib.summary <- as.data.frame(dat.gib.summary)
print(dat.gib.summary)


###############
### AGGREGATE
dat.agg <- rbind(dat.aki, dat.ali, dat.ami, dat.gib)
dat.agg.summary <- c()
dat.agg.summary$ae <- "total"
dat.agg.summary$ror <- auc(dat.agg$casecontrol, dat.agg$ror)
dat.agg.summary$chisq <- auc(dat.agg$casecontrol, dat.agg$chisq)
dat.agg.summary$medcoef0 <- auc(dat.agg$casecontrol, dat.agg$medcoef0)
dat.agg.summary$sql5medcoef1 <- auc(dat.agg$casecontrol, dat.agg$sql5medcoef1)
dat.agg.summary$sql10medcoef1 <- auc(dat.agg$casecontrol, dat.agg$sql10medcoef1)
dat.agg.summary$psi5medcoef1 <- auc(dat.agg$casecontrol, dat.agg$psi5medcoef1)
dat.agg.summary$psi10medcoef1 <- auc(dat.agg$casecontrol, dat.agg$psi10medcoef1)

dat.agg.summary$sql5ate <- auc(dat.agg$casecontrol, dat.agg$sql5ate)
dat.agg.summary$sql10ate <- auc(dat.agg$casecontrol, dat.agg$sql10ate)
dat.agg.summary$psi5ate <- auc(dat.agg$casecontrol, dat.agg$psi5ate)
dat.agg.summary$psi10ate <- auc(dat.agg$casecontrol, dat.agg$psi10ate)
dat.agg.summary <- as.data.frame(dat.agg.summary)
print(dat.agg.summary)
dat.agg <- rbind(dat.aki.summary, dat.ali.summary, dat.ami.summary, dat.gib.summary, dat.agg.summary)
library(xtable)
xtable(dat.agg, digits = 4)


weightedCombinedAuc <- function(dat) {
  ades <-
    c(
      "kidney_failure,_acute",
      "liver_failure,_acute",
      "acute_myocardial_infarction",
      "gastrointestinal_hemorrhage"
    )
  out <- lapply(ades, aggAuc2)
  print("out")
  print(out)
  nrowTot <- nrow(dat)
  weightedCombinedScores <- c()
  for (k in 2:11) {
    weights <- c()
    scores <- c()
    for (a in 1:length(ades)) {
      weight <- as.numeric(unlist(out[a])[12]) / nrowTot #12
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
    paste(
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
      praucPsi10ate,
      sep = " & "
    )
  out
}

aggPrauc(dat, "kidney_failure,_acute")
aggPrauc(dat, "liver_failure,_acute")
aggPrauc(dat, "acute_myocardial_infarction")
aggPrauc(dat, "gastrointestinal_hemorrhage")


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

aggMapk2 <- function(hoi) {
  k = 15
  dat <- data.frame(dat)
  # hoi = "gastrointestinal_hemorrhage"
  dat <- subset(dat, hoiname == hoi)
  dat <- na.omit(dat)
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
  nrowDat <- nrow(dat)
  out <-
    c(
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
      nrowDat
    )
  out
}

weightedCombinedMapk <- function(dat) {
  ades <-
    c(
      "kidney_failure,_acute",
      "liver_failure,_acute",
      "acute_myocardial_infarction",
      "gastrointestinal_hemorrhage"
    )
  out <- lapply(ades, aggMapk2)
  #print(out)
  nrowTot <- nrow(dat)
  weightedCombinedScores <- c()
  for (j in 2:11) {
    weights <- c()
    scores <- c()
    for (a in 1:length(ades)) {
      weight <- as.numeric(unlist(out[a])[j]) / nrowTot #j = 12
      score <- round(as.numeric(unlist(out[a])[j]), 4)
      weights <- c(weights, weight)
      scores <- c(scores, score)
    }
    #print(scores)
    #print(weights)
    weightedCombinedScore <- round(sum(weights * scores), 4)
    weightedCombinedScores <-
      c(weightedCombinedScores, weightedCombinedScore)
  }
  paste0(weightedCombinedScores, collapse = " & ")
}

weightedCombinedMapk(dat)

##############################################################
##############################################################
##### CONFOUNDERUPDOWN ##### (setup) #########################
### idea: what goes down, improves, reducing false positive rate and the like
### both arms of the journey #################################
##############################################################

dat.test <-
  dat #, dat$medcoef0, dat$medcoef1, dat$exposurename, dat$hoiname, dat$exactate)
dat.test <- dat.test[order(dat.test$chisq, decreasing = FALSE),]
dat.test$orderchisq <- order(dat.test$chisq, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$chisq, decreasing = TRUE),]
#head(dat.test, n = 10)
dat.test <- dat.test[order(dat.test$ror, decreasing = FALSE),]
dat.test$orderror <- order(dat.test$ror, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$ror, decreasing = TRUE),]
#head(dat.test, n = 10)
dat.test <- dat.test[order(dat.test$medcoef0, decreasing = FALSE),]
dat.test$ordermedcoef0 <-
  order(dat.test$medcoef0, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$medcoef0, decreasing = TRUE),]

###################
dat.test <-
  dat.test[order(dat.test$sql5medcoef1, decreasing = FALSE),]
dat.test$orderpsi5medcoef1 <-
  order(dat.test$sql5medcoef1, decreasing = TRUE)
dat.test <-
  dat.test[order(dat.test$sql5medcoef1, decreasing = TRUE),]
#---
dat.test <-
  dat.test[order(dat.test$psi5medcoef1, decreasing = FALSE),]
dat.test$orderpsi5medcoef1 <-
  order(dat.test$psi5medcoef1, decreasing = TRUE)
dat.test <-
  dat.test[order(dat.test$psi5medcoef1, decreasing = TRUE),]
#---
dat.test <-
  dat.test[order(dat.test$sql10medcoef1, decreasing = FALSE),]
dat.test$orderpsi10medcoef1 <-
  order(dat.test$sql10medcoef1, decreasing = TRUE)
dat.test <-
  dat.test[order(dat.test$sql10medcoef1, decreasing = TRUE),]
#---
dat.test <-
  dat.test[order(dat.test$psi10medcoef1, decreasing = FALSE),]
dat.test$orderpsi10medcoef1 <-
  order(dat.test$psi10medcoef1, decreasing = TRUE)
dat.test <-
  dat.test[order(dat.test$psi10medcoef1, decreasing = TRUE),]

#####################
dat.test <- dat.test[order(dat.test$sql5ate, decreasing = FALSE),]
dat.test$orderpsi5ate <- order(dat.test$sql5ate, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$sql5ate, decreasing = TRUE),]
#---
dat.test <- dat.test[order(dat.test$psi5ate, decreasing = FALSE),]
dat.test$orderpsi5ate <- order(dat.test$psi5ate, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$psi5ate, decreasing = TRUE),]
#---
dat.test <- dat.test[order(dat.test$sql10ate, decreasing = FALSE),]
dat.test$orderpsi10ate <-
  order(dat.test$sql10ate, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$sql10ate, decreasing = TRUE),]
#---
dat.test <- dat.test[order(dat.test$psi10ate, decreasing = FALSE),]
dat.test$orderpsi10ate <-
  order(dat.test$psi10ate, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$psi10ate, decreasing = TRUE),]

dat.test

#### how much movement in ordered rank?
dat.derived <-
  data.frame(
    dat.test$hoiname,
    dat.test$exposurename,
    as.integer(dat.test$casecontrol),-1 * as.numeric(dat.test$orderror - dat.test$orderpsi10ate)
  )
names(dat.derived) <-
  c("hoiname", "exposurename", "casecontrol", "diff")

dat.derived <-
  dat.derived[order(dat.derived$diff, decreasing = FALSE), ]
dat.derived$orderDD <- order(dat.derived$diff, decreasing = TRUE)
dat.derived <- dat.derived[order(dat.derived$orderDD),]
dat.derived

getScores <- function(dat, hoi, i) {
  hoi <- "acute_myocardial_infarction"
  dat <- subset(dat, hoiname == hoi)
  labs <- as.logical(as.integer(dat$casecontrol))
  score.ror <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$ror)[i], 4)
  score.chisq <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$chisq)[i],
          4)
  score.sql5medcoef <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$sql5medcoef1)[i],
          4)
  score.psi5medcoef <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$psi5medcoef1)[i],
          4)
  score.sql10medcoef <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$sql10medcoef1)[i],
          4)
  score.psi10medcoef <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$psi10medcoef1)[i],
          4)
  score.sql5ate <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$sql5ate)[i],
          4)
  score.psi5ate <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$psi5ate)[i],
          4)
  score.sql10ate <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$sql10ate)[i],
          4)
  score.psi10ate <-
    round(PerfMeas::F.measure.single(labels = labs, pred = dat$psi10ate)[i],
          4)
  scores <-
    c(
      score.ror,
      score.chisq,
      score.sql5medcoef,
      score.psi5medcoef,
      score.sql10medcoef,
      score.psi10medcoef,
      score.sql5ate,
      score.psi5ate,
      score.sql5ate,
      score.psi5ate
    )
  scores
}

getDoobie <- function(jkMin, perfIndex, hoi) {
  dat <- subset(dat, hoiname == hoi)
  #dat <- subset(dat, )
  dat <- subset(dat, a >= jkMin)
  getScores(dat, "kidney_failure,_acute", perfIndex)
  getScores(dat, "liver_failure,_acute", perfIndex)
  getScores(dat, "gastrointestinal_hemorrhage", perfIndex)
  getScores(dat, "acute_myocardial_infarction", perfIndex)
}

getDoobie(50, 4, "gastrointestinal_hemorrhage")


getStats <- function(classifierPredictions, Y_actual) {
  y = as.factor(Y_actual)
  pred <- prediction(classifierPredictions, y)
  # Recall-Precision curve
  RP.perf <- performance(pred, "prec", "rec")
  #plot(RP.perf)
  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr")
  #plot (ROC.perf)
  # ROC area under the curve
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  # F1 score
  print("F1 Score")
  performance(pred, "f")
  plot(performance(pred, "f"))
  # auc
  print("AUC")
  roc(y, classifierPredictions, ci = TRUE)
}

dat.derived$diff2 <-
  log(1 + min(dat.derived$diff) * -1 + dat.derived$diff)
getStats(dat.derived$diff2, dat.derived$casecontrol)

library(ROCit)
logistic.model <-
  glm(
    as.factor(dat.derived$casecontrol) ~ dat.derived$diff,
    data = dat.derived,
    family = "binomial"
  )
class <- logistic.model$y
score <- logistic.model$fitted.values
measure <- measureit(
  score = score,
  class = class,
  measure = c("ACC", "SENS", "FSCR")
)
names(measure)
#> [1] "Cutoff" "Depth"  "TP"     "FP"     "TN"     "FN"     "ACC"    "SENS"
#> [9] "FSCR"
plot(measure$ACC ~ measure$Cutoff, type = "l")



###################
############ JUNK

##########
getMapk(10, as.numeric(dat.derived$casecontrol), 1 + dat.derived$diff)
getMapk(25, dat.derived$casecontrol, dat$ror)
getMapk(25, dat.derived$casecontrol, dat$sql5medcoef1)
getMapk(25, dat.derived$casecontrol, dat$psi5medcoef1)
getMapk(25, dat.derived$casecontrol, dat$sql5ate)
getMapk(25, dat.derived$casecontrol, dat$psi5ate)
getMapk(25, dat.derived$casecontrol, dat$sql10medcoef1)
getMapk(25, dat.derived$casecontrol, dat$psi10medcoef1)
getMapk(25, dat.derived$casecontrol, dat$sql10ate)
getMapk(25, dat.derived$casecontrol, dat$psi10ate)


dat.test$orderRor <- order(dat.test$ror, decreasing = FALSE)
dat.test$orderRor <- order(dat.test$ror)


dat.test <- dat.test[order(dat.test$chisq),]
dat.test$orderRor <- order(dat.test$ror)
dat.test <- dat.test[order(dat.test$ror),]

dat.test$orderMedcoef0 <-
  order(dat.test$medcoef0, decreasing = TRUE)
dat.test <- dat.test[order(dat.test$medcoef0),]
dat.test

dat.test$orderRor <- order(dat.test$dat.ror, decreasing = TRUE)
dat.test$orderChisq <- order(dat.test$dat.chisq, decreasing = TRUE)
dat.test$orderMedcoef0 <-
  order(dat.test$dat.medcoef0, decreasing = TRUE)
dat.test$orderMedcoef1 <-
  order(dat.test$dat.medcoef1, decreasing = TRUE)
dat.test$orderExactate <-
  order(dat.test$dat.exactate, decreasing = TRUE)
dat.test[order(dat.test$orderRor),]
dat.test[order(dat.test$orderChisq),]

dat.test[order(dat.test$orderChisq, decreasing = TRUE),]

narg <- dat.test[order(dat.test$dat.ror, decreasing = TRUE),]

dat.test <- data.frame(dat$casecontrol, dat$medcoef1)
dat.test$order <- order(dat.test$dat.medcoef1, decreasing = FALSE)

#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$chisq))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$ror))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$medcoef0))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$sql5medcoef1))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$psi5medcoef1))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$sql10medcoef1))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$psi10medcoef1))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$sql5ate))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$psi5ate))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$sql10ate))
#MLmetrics::AUC(y_true = as.factor(dat$casecontrol), y_pred = as.numeric(dat$psi10ate))
#############

predictions = as.numeric(dat$psi10ate)
y = as.factor(dat$casecontrol)
pred <- prediction(predictions, y)

# Recall-Precision curve
RP.perf <- performance(pred, "prec", "rec")

plot(RP.perf)

# ROC curve
ROC.perf <- performance(pred, "tpr", "fpr")

plot (ROC.perf)

# ROC area under the curve
auc.tmp <- performance(pred, "auc")

auc <- as.numeric(auc.tmp@y.values)
# F1 score
print("F1 Score")
performance(pred, "f")
plot(performance(pred, "f"))
# auc
print("AUC")
roc(y, predictions, ci = TRUE)

logreg <-
  glm(
    formula = as.numeric(casecontrol) ~ as.numeric(sql5ate),
    family = binomial(link = "logit"),
    data = dat
  )
pred <- ifelse(logreg$fitted.values < 0.4, 0, 1)
F1_Score(
  y_pred = pred,
  y_true = as.numeric(dat$casecontrol),
  positive = "1"
)

optimalF1 <- 0
highestScoreHIGH <- 0
for (jk in 1:10000) {
  spin <- .0001 * jk
  logreg <-
    glm(
      formula = as.numeric(casecontrol) ~ as.numeric(psi10ate),
      family = binomial(link = "logit"),
      data = subset(dat, hoiname = "gastrointestinal_hemorrhage")
    )
  pred <- ifelse(logreg$fitted.values < 0.4, 0, 1)
  F1score <-
    F1_Score(
      y_pred = pred,
      y_true = as.numeric(dat$casecontrol),
      positive = "1"
    )
  if (F1score > highestScoreHIGH) {
    optimalF1 <- spin
    highestScoreHIGH <- F1score
    
  }
}
c(optimalF1, highestScoreHIGH)



pred <-
  prediction(predictions = as.numeric(dat$sql5ate),
             labels = as.numeric(dat$casecontrol))
fmeasure <- performance(pred, "f")
fmeasure = fmeasure@y.values[[1]]

f1_score <- function(predicted, expected, positive.class = "1") {
  predicted <-
    factor(as.character(predicted), levels = unique(as.character(expected)))
  expected  <- as.factor(expected)
  cm = as.matrix(table(expected, predicted))
  print(cm)
  precision <- diag(cm) / colSums(cm)
  recall <- diag(cm) / rowSums(cm)
  f1 <-
    ifelse(precision + recall == 0,
           0,
           2 * precision * recall / (precision + recall))
  
  #Assuming that F1 is zero when it's not possible compute it
  f1[is.na(f1)] <- 0
  
  #Binary F1 or Multi-class macro-averaged F1
  ifelse(nlevels(expected) == 2, f1[positive.class], mean(f1))
}
dat2 <- subset(dat, hoiname == "kidney_failure,_acute")
dat2 <- subset(dat, hoiname == "liver_failure,_acute")
dat2 <- subset(dat, hoiname == "gastrointestinal_hemorrhage")
dat2 <- subset(dat, hoiname == "acute_myocardial_infarction")
f1_score(as.numeric(dat2$sql5ate),
         as.factor(dat2$casecontrol),
         positive.class = 1)

measurePrecisionRecall <- function(hoi) {
  dat <- subset(dat, hoiname == hoi)
  predict <- as.numeric(dat$sql5ate)
  actual_labels <- as.integer(dat$casecontrol)
  precision <- sum(predict & actual_labels) / sum(predict)
  recall <- sum(predict & actual_labels) / sum(actual_labels)
  fmeasure <- 2 * precision * recall / (precision + recall)
  
  cat('precision:  ')
  cat(precision)
  cat('%')
  cat('\n')
  
  cat('recall:     ')
  cat(recall)
  cat('%')
  cat('\n')
  
  cat('f-measure:  ')
  cat(fmeasure)
  cat('%')
  cat('\n')
}

measurePrecisionRecall("kidney_failure,_acute") 
measurePrecisionRecall("liver_failure,_acute") 
measurePrecisionRecall("acute_myocardial_infarction") 
measurePrecisionRecall("gastrointestinal_hemorrhage") 

dat.test[order(dat.test$order),]


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

system("mkdir boxplots")

# BASELINES

ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), title = "Baselines for kidney_failure,_acute", 
          x = "casecontrol", y = c("ror", "medcoef0"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), title = "Baselines for liver_failure,_acute", 
          x = "casecontrol", y = c("ror", "medcoef0"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), title = "Baselines for acute_myocardial_infarction", 
          x = "casecontrol", y = c("ror", "medcoef0"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), title = "Baselines for gastrointestinal_hemorrhage",
          x = "casecontrol", y = c("ror", "medcoef0"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 




# medcoef
ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), title = "Exposure Coefficients (β) for kidney_failure,_acute",
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'),
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'),
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'),
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Exposure Coefficient (β)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

# ATE
ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), title = "Average Treatment Effect (Δ) for kidney_failure,_acute",
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Average Treatment Effect (ATE Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), title = "Average Treatment Effect (Δ) for liver_failure,_acute",
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Average Treatment Effect (ATE Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), title = "Average Treatment Effect (Δ) for acute_myocardial_infarction",
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Average Treatment Effect (ATE Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), title = "Average Treatment Effect (Δ) for gastrointestinal_hemorrhage",
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Average Treatment Effect (ATE Δ)", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 



# 'kidney_failure,_acute'
ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), x = "casecontrol", 
          y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow", 
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), 
          font.label = list(size = 10, style = "bold", color = "black"), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), 
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "coefficients", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow",
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'kidney_failure,_acute'), 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#'liver_failure,_acute'
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), x = "casecontrol", 
          y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow", 
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), 
          font.label = list(size = 10, style = "bold", color = "black"), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), 
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "coefficients", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow",
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 



ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'liver_failure,_acute'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
# 'acute_myocardial_infarction'


ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", 
          y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow", 
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), 
          font.label = list(size = 10, style = "bold", color = "black"), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), 
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "coefficients", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow",
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 



ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("ror", "medcoef0"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "medcoef1", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction'), x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, combine = TRUE, add = c("dotplot", "median_iqr"), fill = "light yellow", label = "exposurename", label.select = list(top.up = 1, top.down = 1)) 
# 'gastrointestinal_hemorrhage'



ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), x = "casecontrol", 
          y = c("ror"), color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow", 
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), 
          font.label = list(size = 10, style = "bold", color = "black"), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), 
          x = "casecontrol", y = c("medcoef0", "sql5medcoef1", "sql10medcoef1", "psi5medcoef1", "psi10medcoef1"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "coefficients", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), fill = "light yellow",
          label = "exposurename", label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
ggboxplot(subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage'), 
          x = "casecontrol", y = c("sql5ate", "sql10ate", "psi5ate", "psi10ate"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "sql5ate", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", 
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 








gib <- subset(dat.new, dat.new$hoiname == 'gastrointestinal_hemorrhage')
gib0 <- subset(gib, gib$casecontrol == '0')
gib1 <- subset(gib, gib$casecontrol == '1')
gib0median <- median(gib0$psi5ate)
subset(gib1, gib1$psi10ate > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$psi5ate > gib0median))

gib0median <- median(gib0$ror)
subset(gib1, gib1$ror > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$ror > gib0median))

gib0median <- median(gib0$ror)
subset(gib0, gib0$ror > gib0median)
nrow(gib0)
nrow(subset(gib0, gib0$ror > gib0median))

gib0median <- median(gib0$sql5ate)
nrow(gib0)
nrow(subset(gib0, gib0$sql5ate > gib0median))



gib <- subset(dat.new, dat.new$hoiname == 'liver_failure,_acute')
gib0 <- subset(gib, gib$casecontrol == '0')
gib1 <- subset(gib, gib$casecontrol == '1')
gib0median <- median(gib0$psi5ate)
subset(gib1, gib1$psi10ate > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$psi5ate > gib0median))

gib0median <- median(gib0$ror)
subset(gib1, gib1$ror > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$ror > gib0median))

gib0median <- median(gib0$ror)
subset(gib0, gib0$ror > gib0median)
nrow(gib0)
nrow(subset(gib0, gib0$ror > gib0median))

gib0median <- median(gib0$sql5ate)
nrow(gib0)
nrow(subset(gib0, gib0$sql5ate > gib0median))

gib <- subset(dat.new, dat.new$hoiname == 'acute_myocardial_infarction')
gib0 <- subset(gib, gib$casecontrol == '0')
gib1 <- subset(gib, gib$casecontrol == '1')
gib0median <- median(gib0$psi5ate)
subset(gib1, gib1$psi10ate > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$psi5ate > gib0median))

gib0median <- median(gib0$ror)
subset(gib1, gib1$ror > gib0median)
nrow(gib1)
nrow(subset(gib1, gib1$ror > gib0median))

gib0median <- median(gib0$ror)
subset(gib0, gib0$ror > gib0median)
nrow(gib0)
nrow(subset(gib0, gib0$ror > gib0median))


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


