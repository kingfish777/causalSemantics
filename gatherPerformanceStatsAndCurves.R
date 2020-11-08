
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




