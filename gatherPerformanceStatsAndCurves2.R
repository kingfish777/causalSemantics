
library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
 #setwd("/Users/scottalexandermalec/Projects/causalSemantics/")
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
#baselineColor <- "#ffb64e"; glmColor <- "#18c3ff"; gcmColor <- "#37b937" # https://medialab.github.io/iwanthue/
baselineColor <- "#663500"; glmColor <- "#63a11b"; gcmColor <- "#ff845d" # https://medialab.github.io/iwanthue/

randomAKI <- data.frame(read.csv("randomAKI.txt", header = FALSE, sep = "\t"))
randomALI <- data.frame(read.csv("randomALI.txt", header = FALSE, sep = "\t"))
randomAMI <- data.frame(read.csv("randomAMI.txt", header = FALSE, sep = "\t"))
randomGIB <- data.frame(read.csv("randomGIB.txt", header = FALSE, sep = "\t"))


drawAUPRC <- function(squelch, confounderSource, hoiname) {
  dat <-
    dbGetQuery(
      con,
      paste(
        "SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a,
                                j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue,
                                j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND
                  j.hoiname LIKE '",
        hoiname,
        "' AND j.casecontrol = r.casecontrol AND a >= 10 AND (a+b) >= 100 AND
                  j.squelch = ",
        squelch,
        " AND confoundersource = '",
        confounderSource,
        "';",
        sep = ""
      )
    )
  
  pos <- print(sum(dat$casecontrol))
  tot <- nrow(dat)
  neg <- tot - pos
  dat <- na.omit(dat)
  pred0 <- prediction(dat$ror, dat$casecontrol)
  perf0 <- performance (pred0, "prec", "rec")
  pred1 <- prediction(dat$medcoef, dat$casecontrol)
  perf1 <- performance (pred1, "prec", "rec")
  pred2 <- prediction(dat$exactate, dat$casecontrol)
  perf2 <- performance (pred2, "prec", "rec")
  nerding <-
    plot(
      main = paste(
        confounderSource,
        squelch,
        "(pos: ",
        pos,
        ", neg: ",
        neg,
        "total: ",
        tot,
        ")"
      ),
      perf0,
      col = baselineColor,
      print.auc.y = .3
    )
  nerding <- plot(perf1,
                  add = TRUE,
                  col = glmColor,
                  print.auc.y = .3)
  nerding <-
    plot(perf2,
         add = TRUE,
         col = gcmColor,
         print.auc.y = .3)
  #print(F1_Score(y_pred = ifelse(as.numeric(dat$exactate) > 0.1, 0, 1), y_true = as.numeric(dat$casecontrol), positive = "1"))
}

drawAUPRC(5, "sql", "kidney_failure,_acute")
drawAUPRC(10, "sql", "kidney_failure,_acute")
drawAUPRC(5, "psi", "kidney_failure,_acute")
drawAUPRC(10, "psi", "kidney_failure,_acute")

drawAUPRC(5, "sql", "liver_failure,_acute")
drawAUPRC(10, "sql", "liver_failure,_acute")
drawAUPRC(5, "psi", "liver_failure,_acute")
drawAUPRC(10, "psi", "liver_failure,_acute")

drawAUPRC(5, "sql", "acute_myocardial_infarction")
drawAUPRC(10, "sql", "acute_myocardial_infarction")
drawAUPRC(5, "psi", "acute_myocardial_infarction")
drawAUPRC(10, "psi", "acute_myocardial_infarction")

drawAUPRC(5, "sql", "gastrointestinal_hemorrhage")
drawAUPRC(10, "sql", "gastrointestinal_hemorrhage")
drawAUPRC(5, "psi", "gastrointestinal_hemorrhage")
drawAUPRC(10, "psi", "gastrointestinal_hemorrhage")

drawAUPRC(5, "sql", "%")
drawAUPRC(10, "sql", "%")
drawAUPRC(5, "psi", "%")
drawAUPRC(10, "psi", "%")

##############
# BY ADE
##############

drawAUROC <- function(squelch, confounderSource, hoiname) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue, j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND  
                  j.hoiname LIKE '", hoiname, "' AND j.casecontrol = r.casecontrol AND
                  j.squelch = ", squelch, " AND confoundersource = '", confounderSource, "';", sep = ""))
  
  pos <- print(sum(dat$casecontrol))
  tot <- nrow(dat)
  neg <- tot - pos
  # rocSQL5Chisq <- plot(main = paste(hoiname,confounderSource, squelch), smooth(roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark grey")
  # rocSQL5Medcoef <- plot(smooth(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE)
  # rocSQL5exactate <- plot(smooth(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
  rocSQL5Chisq <- plot(main = paste(confounderSource, squelch, 
                                    "(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)
  
}

drawAUROC(5, "sql", "kidney_failure,_acute") # g > d > r
drawAUROC(10, "sql", "kidney_failure,_acute") # g > d r
drawAUROC(5, "psi", "kidney_failure,_acute") # d > g > r
drawAUROC(10, "psi", "kidney_failure,_acute") # d g r

drawAUROC(5, "sql", "liver_failure,_acute") # r g > d
drawAUROC(10, "sql", "liver_failure,_acute") # r g > d
drawAUROC(5, "psi", "liver_failure,_acute") # d g > r
drawAUROC(10, "psi", "liver_failure,_acute") # d g > r

drawAUROC(5, "sql", "acute_myocardial_infarction") # r > g d
drawAUROC(10, "sql", "acute_myocardial_infarction") # r  g > d
drawAUROC(5, "psi", "acute_myocardial_infarction") # r > g d
drawAUROC(10, "psi", "acute_myocardial_infarction") # r > g d

drawAUROC(5, "sql", "gastrointestinal_hemorrhage") # r g d
drawAUROC(10, "sql", "gastrointestinal_hemorrhage") # r > d g
drawAUROC(5, "psi", "gastrointestinal_hemorrhage") # d > r > g
drawAUROC(10, "psi", "gastrointestinal_hemorrhage") # d > r > g

pos <- sum(randomAKI$casecontrol)
neg <- nrow(randomAKI) - pos
tot = nrow(randomAKI)
#plot(main = paste("randomly selected covariates for AKI ",
#                  "(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), 


####### calculate chisqvalues

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


#### AKI AKI AKI AKI AKI random vs non random
overlaps <- intersect(randomAKI$exposurename, dat$exposurename)

nonrandom <- subset(dat, hoiname == 'kidney_failure,_acute' & exposurename %in% overlaps)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(randomAKI$casecontrol, randomAKI$medCoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(randomAKI$casecontrol, randomAKI$exactATE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

akiCombined <- data.frame(merge(x = randomAKI, y = nonrandom, by = "exposurename"))
ggboxplot(akiCombined, title = "Baselines for kidney_failure,_acute", 
          x = "casecontrol.y", y = c("medCoef0", "medCoef", "psi5medcoef1", "psi10medcoef1"),  #c("exactATE", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(akiCombined, title = "Baselines for kidney_failure,_acute", 
          x = "casecontrol.y", y = c("exactATE", "psi5ate", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

akiCombined$medCoefDiffRandom <- akiCombined$medCoef0 - akiCombined$medCoef
akiCombined$medCoefDiffLit <- akiCombined$medCoef0 - akiCombined$psi10medcoef1

ggboxplot(akiCombined, title = "Random vs Lit-derived confounders for kidney_failure,_acute", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#### GOOD
###### GOOD
akiCombined <- data.frame(merge(x = randomAKI, y = nonrandom, by = "exposurename"))
akiCombined$medCoefDiff_Random <- akiCombined$medCoef0 - akiCombined$medCoef
akiCombined$medCoefDiffLit_sql5 <- akiCombined$medCoef0 - akiCombined$sql5medcoef1
akiCombined$medCoefDiffLit_sql10 <- akiCombined$medCoef0 - akiCombined$sql10medcoef1
akiCombined$medCoefDiffLit_psi5 <- akiCombined$medCoef0 - akiCombined$psi10medcoef1
akiCombined$medCoefDiffLit_psi10 <- akiCombined$medCoef0 - akiCombined$psi10medcoef1
names(akiCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiff_Random', 'medCoefDiffLit_sql5', 'medCoefDiffLit_sql10', 'medCoefDiffLit_psi5', 'medCoefDiffLit_psi10')

ggboxplot(akiCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiff_Random", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

auc(akiCombined$casecontrol, akiCombined$medCoefDiffLit_sql5)
auc(akiCombined$casecontrol, akiCombined$medCoefDiffLit_sql10)
auc(akiCombined$casecontrol, akiCombined$medCoefDiffLit_psi5)
auc(akiCombined$casecontrol, akiCombined$medCoefDiffLit_psi10)

#### ALI ALI ALI ALI ALI random vs non random
overlaps <- intersect(randomALI$exposurename, dat$exposurename)

nonrandom <- subset(dat, hoiname == 'liver_failure,_acute' & exposurename %in% overlaps)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(randomALI$casecontrol, randomALI$medCoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(randomALI$casecontrol, randomALI$exactATE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

aliCombined <- data.frame(merge(x = randomALI, y = nonrandom, by = "exposurename"))
ggboxplot(aliCombined, title = "Baselines for liver_failure,_acute", 
          x = "casecontrol.y", y = c("medCoef0", "medCoef", "psi5medcoef1", "psi10medcoef1"), #c("exactATE", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(aliCombined, title = "Baselines for liver_failure,_acute", 
          x = "casecontrol.y", y = c("exactATE", "psi5ate", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

aliCombined$medCoefDiff_Random <- aliCombined$medCoef0 - aliCombined$medCoef
aliCombined$medCoefDiffLit_psi10 <- aliCombined$medCoef0 - aliCombined$psi10medcoef1
aliCombined$medCoefDiffLit_sql5 <- aliCombined$medCoef0 - aliCombined$sql5medcoef1


ggboxplot(aliCombined, title = "Random vs Lit-derived confounders for liver_failure,_acute", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#### GOOD
###### GOOD
aliCombined <- data.frame(merge(x = randomAKI, y = nonrandom, by = "exposurename"))
aliCombined$medCoefDiff_Random <- aliCombined$medCoef0 - aliCombined$medCoef
aliCombined$medCoefDiffLit_sql5 <- aliCombined$medCoef0 - aliCombined$sql5medcoef1
aliCombined$medCoefDiffLit_sql10 <- aliCombined$medCoef0 - aliCombined$sql10medcoef1
aliCombined$medCoefDiffLit_psi5 <- aliCombined$medCoef0 - aliCombined$psi10medcoef1
aliCombined$medCoefDiffLit_psi10 <- aliCombined$medCoef0 - aliCombined$psi10medcoef1
names(aliCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiff_Random', 'medCoefDiffLit_sql5', 'medCoefDiffLit_sql10', 'medCoefDiffLit_psi5', 'medCoefDiffLit_psi10')

ggboxplot(aliCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiff_Random", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

auc(aliCombined$casecontrol, aliCombined$medCoefDiffLit_sql5)
auc(aliCombined$casecontrol, aliCombined$medCoefDiffLit_sql10)
auc(aliCombined$casecontrol, aliCombined$medCoefDiffLit_psi5)
auc(aliCombined$casecontrol, aliCombined$medCoefDiffLit_psi10)


#### AMI AMI AMI AMI AMI random vs non random

overlaps <- intersect(randomAMI$exposurename, dat$exposurename)

nonrandom <- subset(dat, hoiname == 'acute_myocardial_infarction' & exposurename %in% overlaps)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(randomAMI$casecontrol, randomAMI$medCoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(randomAMI$casecontrol, randomAMI$exactATE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

amiCombined <- data.frame(merge(x = randomAMI, y = nonrandom, by = "exposurename", no.dups = FALSE))
unique(amiCombined)
ggboxplot(amiCombined, title = "Regression coefficients for acute_myocardial_infarction", 
          x = "casecontrol.y", y = c("medCoef0", "medCoef", "psi5medcoef1", "psi10medcoef1"),  #c("exactATE", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(amiCombined, title = "Effect Estimates for acute_myocardial_infarction", 
          x = "casecontrol.y", y = c("exactATE", "psi5ate", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

amiCombined$medCoefDiffRandom <- amiCombined$medCoef0 - amiCombined$medCoef
amiCombined$medCoefDiffLit <- amiCombined$medCoef0 - amiCombined$psi10medcoef1
names(amiCombined$casecontrol.x) <- "casecontrol"
ggboxplot(amiCombined, title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#### GOOD
###### GOOD
amiCombined <- data.frame(merge(x = randomAMI, y = nonrandom, by = "exposurename"))
amiCombined$medCoefDiff_Random <- amiCombined$medCoef0 - amiCombined$medCoef
amiCombined$medCoefDiffLit_sql5 <- amiCombined$medCoef0 - amiCombined$sql5medcoef1
amiCombined$medCoefDiffLit_sql10 <- amiCombined$medCoef0 - amiCombined$sql10medcoef1
amiCombined$medCoefDiffLit_psi5 <- amiCombined$medCoef0 - amiCombined$psi10medcoef1
amiCombined$medCoefDiffLit_psi10 <- amiCombined$medCoef0 - amiCombined$psi10medcoef1
names(amiCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiff_Random', 'medCoefDiffLit_sql5', 'medCoefDiffLit_sql10', 'medCoefDiffLit_psi5', 'medCoefDiffLit_psi10')

ggboxplot(amiCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiff_Random", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

auc(amiCombined$casecontrol, amiCombined$medCoefDiffLit_sql5)
auc(amiCombined$casecontrol, amiCombined$medCoefDiffLit_sql10)
auc(amiCombined$casecontrol, amiCombined$medCoefDiffLit_psi5)
auc(amiCombined$casecontrol, amiCombined$medCoefDiffLit_psi10)



#ggboxplot(amiCombined, x = "casecontrol.y", y = c("exactATE", "psi5ate"), combine = TRUE) #, label = "exposurename", label.select = list(top.up = 3, top.down = 3))

amiCombined <- data.frame(merge(x = randomAMI, y = nonrandom, by = "exposurename"))

amiCombined$medCoefDiffRandom <- amiCombined$medCoef0 - amiCombined$medCoef
amiCombined$medCoefDiffGood <- amiCombined$medCoef0 - amiCombined$psi10medcoef1

ggboxplot(amiCombined, title = "Baselines for acute_myocardial_infarction", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffGood"), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#### GIB GIB GIB GIB GIB random vs non random

overlaps <- intersect(randomGIB$exposurename, dat$exposurename)

nonrandom <- subset(dat, hoiname == 'gastrointestinal_hemorrhage' & exposurename %in% overlaps)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(randomGIB$casecontrol, randomGIB$medCoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(randomGIB$casecontrol, randomGIB$exactATE, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$sql10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi5ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

plot(roc(nonrandom$casecontrol, nonrandom$ror, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, print.auc.y = .1, colorize = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10medcoef1, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .2, add = TRUE)
plot(roc(nonrandom$casecontrol, nonrandom$psi10ate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)

gibCombined <- data.frame(merge(x = randomGIB, y = nonrandom, by = "exposurename"))
ggboxplot(gibCombined, title = "Baselines for gastrointestinal_hemorrhage", 
          x = "casecontrol.y", y = c("medCoef0", "medCoef", "psi5medcoef1", "psi10medcoef1"),  #c("exactATE", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

ggboxplot(gibCombined, title = "Baselines for gastrointestinal_hemorrhage", 
          x = "casecontrol.y", y = c("exactATE", "psi5ate", "psi10ate" ), 
          color = "casecontrol.y", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

gibCombined$medCoefDiffRandom <- gibCombined$medCoef0 - gibCombined$medCoef
gibCombined$medCoefDiffLit <- gibCombined$medCoef0 - gibCombined$psi10medcoef1

ggboxplot(gibCombined, title = "Random vs Lit-derived confounders for gastrointestinal_hemorrhage", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

#### GOOD
###### GOOD
gibCombined <- data.frame(merge(x = randomGIB, y = nonrandom, by = "exposurename"))
gibCombined$medCoefDiff_Random <- gibCombined$medCoef0 - gibCombined$medCoef
gibCombined$medCoefDiffLit_sql5 <- gibCombined$medCoef0 - gibCombined$sql5medcoef1
gibCombined$medCoefDiffLit_sql10 <- gibCombined$medCoef0 - gibCombined$sql10medcoef1
gibCombined$medCoefDiffLit_psi5 <- gibCombined$medCoef0 - gibCombined$psi10medcoef1
gibCombined$medCoefDiffLit_psi10 <- gibCombined$medCoef0 - gibCombined$psi10medcoef1
names(gibCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiff_Random', 'medCoefDiffLit_sql5', 'medCoefDiffLit_sql10', 'medCoefDiffLit_psi5', 'medCoefDiffLit_psi10')

ggboxplot(gibCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiff_Random", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

auc(gibCombined$casecontrol, gibCombined$medCoefDiffLit_sql5)
auc(gibCombined$casecontrol, gibCombined$medCoefDiffLit_sql10)
auc(gibCombined$casecontrol, gibCombined$medCoefDiffLit_psi5)
auc(gibCombined$casecontrol, gibCombined$medCoefDiffLit_psi10)



###############


#### GOOD
###### GOOD
datCombined <- dat
datCombinedall <- dat
#datCombined$medCoefDiff_Random <- datCombined$medcoef0 - datCombined$medCoef
datCombined$medCoefDiffLit_sql5 <- datCombined$medcoef0 - datCombined$sql5medcoef1
datCombined$medCoefDiffLit_sql10 <- datCombined$medcoef0 - datCombined$sql10medcoef1
datCombined$medCoefDiffLit_psi5 <- datCombined$medcoef0 - datCombined$psi10medcoef1
datCombined$medCoefDiffLit_psi10 <- datCombined$medcoef0 - datCombined$psi10medcoef1
#names(datCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiffLit_sql5', 'medCoefDiffLit_sql10', 'medCoefDiffLit_psi5', 'medCoefDiffLit_psi10')


ggboxplot(datCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

datCombined <- subset(datCombined, hoiname == 'acute_myocardial_infarction')
datCombined <- subset(datCombined, hoiname == 'kidney_failure,_acute')
datCombined <- subset(datCombined, hoiname == 'liver_failure,_acute')
datCombined <- subset(datCombined, hoiname == 'gastrointestinal_hemorrhage')

auc(datCombined$casecontrol, datCombined$medCoefDiffLit_sql5)
auc(datCombined$casecontrol, datCombined$medCoefDiffLit_sql10)
auc(datCombined$casecontrol, datCombined$medCoefDiffLit_psi5)
auc(datCombined$casecontrol, datCombined$medCoefDiffLit_psi10)

######## all combined AUCS
# sql5 = 0.4732 AUROC
# sql10 = 0.4961 AUROC
# psi5 = 0.6508
# psi10 = 0.6508

plot(roc(datCombined$casecontrol, datCombined$medCoefDiffLit_psi10))


allCombined <- rbind(akiCombined, aliCombined, akiCombined, gibCombined)
allCombined$medCoefDiffRandom <- allCombined$medCoef0 - allCombined$medCoef
allCombined$medCoefDiffLit <- allCombined$medCoef0 - allCombined$psi10medcoef1

ggboxplot(allCombined, title = "Random vs Lit-derived confounders for aLL", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 
#### GOOD
###### GOOD
allCombined$medCoefDiff_Random <- allCombined$medCoef0 - allCombined$medCoef
allCombined$medCoefDiffLit_psi10 <- allCombined$medCoef0 - allCombined$psi10medcoef1
allCombined$medCoefDiffLit_sql5 <- allCombined$medCoef0 - allCombined$sql5medcoef1
names(allCombined) <- c('exposurename', 'hoiname', 'cl', 'casecontrol', 'hoiname.x.1', 'cl.1', 'SQUELCH', 'chisqPvalue', 'drug.ade.score0', 'drug.ade.score', 'exactATE', 'exactATEvar', 'medCoef0', 'medZvalue0', 'medPvalue0', 'medCoef', 'medZvalue', 'medPvalue', 'networkScore0', 'networkScore1', 'networkScoreDiff', 'casecontrol.y', 'hoiname.y', 'a', 'b', 'c', 'd', 'ror', 'medcoef0', 'sql5ate', 'sql5medcoef1', 'psi5ate', 'psi5medcoef1', 'sql10ate', 'sql10medcoef1', 'psi10ate', 'psi10medcoef1', 'chisq', 'medCoefDiffRandom', 'medCoefDiffLit', 'medCoefDiff_Random', 'medCoefDiffLit_psi10', 'medCoefDiffLit_sql5')

ggboxplot(allCombined, #title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol", y = c("medCoefDiffRandom", "medCoefDiffLit_sql5", "medCoefDiffLit_psi10"), 
          color = "casecontrol", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "case control labels", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename", font.label = 15,
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 

auc(allCombined$casecontrol, allCombined$medCoefDiffLit_sql5)

######################


datCombined <- dat
datCombinedall <- dat
#datCombined$medCoefDiff_Random <- datCombined$medcoef0 - datCombined$medCoef
datCombined$medCoefDiffLit_sql5 <- datCombined$medcoef0 - datCombined$sql5medcoef1
datCombined$medCoefDiffLit_sql10 <- datCombined$medcoef0 - datCombined$sql10medcoef1
datCombined$medCoefDiffLit_psi5 <- datCombined$medcoef0 - datCombined$psi10medcoef1
datCombined$medCoefDiffLit_psi10 <- datCombined$medcoef0 - datCombined$psi10medcoef1

datCombined.orig <- datCombined

datCombined <- datCombined.orig
datCombined <- subset(datCombined, hoiname == 'kidney_failure,_acute')
datCombined <- subset(datCombined, hoiname == 'liver_failure,_acute')
datCombined <- subset(datCombined, hoiname == 'acute_myocardial_infarction')
datCombined <- subset(datCombined, hoiname == 'gastrointestinal_hemorrhage')
datCombined <- subset(datCombined, hoiname != 'liver_failure,_acute')

getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medcoef0)
getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql5)
getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql10)
getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi5)
getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi10)


getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medcoef0)
getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql5)
getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql10)
getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi5)
getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi10)

#[1] 0.916
#> getMapk(k = 10, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi10)
#[1] 0.916
#> getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi10)
#[1] 0.822
#> getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi10)
#[1] 0.822
#> getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_psi5)
#[1] 0.822
#> getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql5)
#[1] 0.6112
#> getMapk(k = 25, Yactual = datCombined$casecontrol, Ypredicted = datCombined$medCoefDiffLit_sql10)
#[1] 0.5713
#> nrows(allCombined)
#Error in nrows(allCombined) : could not find function "nrows"
#> nrow(allCombined)
#[1] 48
#> getMapk(k = 25, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_sql10)
#[1] 0.7804
#> getMapk(k = 25, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_psi10)
#[1] 0.9342
#> getMapk(k = 25, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiff_Random)
#[1] 0.4859
#> getMapk(k = 25, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_sql5)
#[1] 0.8872
#> getMapk(k = 25, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_sql10)
#[1] 0.7804
#> getMapk(k = 10, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_psi10)
#[1] 1
#> getMapk(k = 10, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_psi5)
#[1] 1
#> getMapk(k = 10, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_sql5)
#[1] 1
#> getMapk(k = 10, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiffLit_sql10)
#[1] 0.8386
#> getMapk(k = 10, Yactual = allCombined$casecontrol, Ypredicted = allCombined$medCoefDiff_Random)
#[1] 0.3048
#
