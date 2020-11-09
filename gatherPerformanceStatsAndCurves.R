library(ggpubr)
library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
setwd("/Users/scottalexandermalec/Projects/causalSemantics/")
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "confounders",
    host = "localhost",
    port = 5432,
    user = "smalec",
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

aliCombined$medCoefDiffRandom <- aliCombined$medCoef0 - aliCombined$medCoef
aliCombined$medCoefDiffLit <- aliCombined$medCoef0 - aliCombined$psi10medcoef1

ggboxplot(aliCombined, title = "Random vs Lit-derived confounders for liver_failure,_acute", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 


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

amiCombined <- data.frame(merge(x = randomAMI, y = nonrandom, by = "exposurename"))
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

ggboxplot(amiCombined, title = "Random vs Lit-derived confounders for acute_myocardial_infarction", 
          x = "casecontrol.x", y = c("medCoefDiffRandom", "medCoefDiffLit"), 
          color = "casecontrol.x", palette = c("#00AFBB", "#E7B800"), 
          ylab = "Baselines", xlab = "casecontrol.x.x", bxp.errorbar = TRUE, 
          combine = TRUE, add = c("mean_ci", "dotplot"), 
          fill = "light yellow", label = "exposurename",
          label.select = list(top.up = 3, top.down = 3), repel = TRUE) 



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

