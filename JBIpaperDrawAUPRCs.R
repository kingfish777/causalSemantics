

library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
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
baselineColor <- "#ffb64e"; glmColor <- "#18c3ff"; gcmColor <- "#37b937" # https://medialab.github.io/iwanthue/
baselineColor <- "#663500"; glmColor <- "#63a11b"; gcmColor <- "#ff845d" # https://medialab.github.io/iwanthue/

##############
# OVERALL
##############

drawAUROC <- function(squelch, confounderSource) {
  squelch <- 5
  confounderSource <- "sql"
  dat <-
    dbGetQuery(
      con,
      paste(
        "SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue,
                                       j.medcoef0, j.medpvalue0, j.medcoef, j.medpvalue,
                                       j.exactate, j.exactatevar
                                         FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND a >= 100 AND -- j.medpvalue0 < .05 AND
                    j.casecontrol = r.casecontrol AND j.squelch = ",
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
  
  #rocSQL5Chisq <- plot(main = paste(confounderSource, squelch), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4),
  #print.auc = TRUE, col = "dark grey", colorize = TRUE)
  #rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE,
  #                       col = "dark red", print.auc.y = .4, add = TRUE, colorize = TRUE)
  #rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE,
  #                        col = "red", print.auc.y = .3, add = TRUE)
}
par(mar=c(1,1,1,1))
drawAUROC(5, "sql")
drawAUROC(10, "sql")
drawAUROC(5, "psi")
drawAUROC(10, "psi")

##############
# BY ADE
##############

drawAUROC <- function(squelch, confounderSource, hoiname) {
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

drawAUROC(5, "sql", "kidney_failure,_acute")
drawAUROC(10, "sql", "kidney_failure,_acute")
drawAUROC(5, "psi", "kidney_failure,_acute")
drawAUROC(10, "psi", "kidney_failure,_acute")

drawAUROC(5, "sql", "liver_failure,_acute")
drawAUROC(10, "sql", "liver_failure,_acute")
drawAUROC(5, "psi", "liver_failure,_acute")
drawAUROC(10, "psi", "liver_failure,_acute")

drawAUROC(5, "sql", "acute_myocardial_infarction")
drawAUROC(10, "sql", "acute_myocardial_infarction")
drawAUROC(5, "psi", "acute_myocardial_infarction")
drawAUROC(10, "psi", "acute_myocardial_infarction")

drawAUROC(5, "sql", "gastrointestinal_hemorrhage")
drawAUROC(10, "sql", "gastrointestinal_hemorrhage")
drawAUROC(5, "psi", "gastrointestinal_hemorrhage")
drawAUROC(10, "psi", "gastrointestinal_hemorrhage")

drawAUROC(5, "sql", "%")
drawAUROC(10, "sql", "%")
drawAUROC(5, "psi", "%")
drawAUROC(10, "psi", "%")

