
library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
setwd("/Users/scottalexandermalec/Projects/causalSemantics/")
drv <- dbDriver('PostgreSQL')
drv = RPostgreSQL::PostgreSQL(max.con=20)
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
#baselineColor <- "#663500"; glmColor <- "#63a11b"; gcmColor <- "#ff845d" # https://medialab.github.io/iwanthue/

rorColor <- "#e75e39"
chi2Color <- "#832000"
coefSql5Color <- "#d2870a"
coefSql10Color <- "#8c6800"
coefPsi5Color <- "#91ac19"
coefPsi10Color <- "#24af4b"
ateSql5Color <- "#019adb"
ateSql10Color <- "#a26de5"
atePsi5Color <- "#540050"
atePsi10Color <- "#a9006f"

drawAUPRC <- function(hoiname) {
  hoiname = "kidney_failure,_acute"
  dat <- dbGetQuery(
    con,
    paste("SELECT DISTINCT sql5.casecontrol, sql5.hoiname, sql5.exposurename, 
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
                        sql5.squelch = '5' AND psi5.squelch = '5' AND psi10.hoiname = '", hoiname, "' AND 
                        sql10.squelch = '10' AND psi10.squelch = '10'", sep = "")
  )
  print(dat)
  pos <- print(sum(as.numeric(dat$casecontrol)))
  tot <- nrow(dat)
  neg <- tot - pos
  dat <- na.omit(dat)
  dat$ror <- dat$a*dat$c/(dat$c*dat$d)
  predRor <- prediction(dat$ror, dat$casecontrol)
  perfRor <- performance (predRor, "prec", "rec")
# perfRor perfSql5medcoef perfSql10medcoef perfPsi5medcoef perfPsi10medcoef perfSql5ate perfSql10ate perfPsi5ate perfPsi10ate
  #  predChisq <- prediction(dat$chisq, dat$casecontrol)
#  perfChisq <- performance (predChisq, "prec", "rec")
  predSql5medcoef <- prediction(dat$sql5medcoef1, dat$casecontrol)
  perfSql5medcoef <- performance (predSql5medcoef, "prec", "rec")
  predSql10medcoef <- prediction(dat$sql10medcoef1, dat$casecontrol)
  perfSql10medcoef <- performance (predSql10medcoef, "prec", "rec")
  predPsi5medcoef <- prediction(dat$sql5medcoef1, dat$casecontrol)
  perfPsi5medcoef <- performance (predPsi5medcoef, "prec", "rec")
  predPsi10medcoef <- prediction(dat$psi10medcoef1, dat$casecontrol)
  perfPsi10medcoef <- performance(predPsi10medcoef, "prec", "rec")
  predSql5ate <- prediction(dat$sql5ate, dat$casecontrol)
  perfSql5ate <- performance (predSql5ate, "prec", "rec")
  predSql10ate <- prediction(dat$sql10ate, dat$casecontrol)
  perfSql10ate <- performance (predSql10ate, "prec", "rec")
  predPsi5ate <- prediction(dat$sql5ate, dat$casecontrol)
  perfPsi5ate <- performance (predPsi5ate, "prec", "rec")
  predPsi10ate <- prediction(dat$psi10ate, dat$casecontrol)
  perfPsi10ate <- performance(predPsi10ate, "prec", "rec")
  #### 
  statName = "AUPRC"
  fn <- paste("PerformanceStatsAndCurves/", hoiname, "-", statName, "_", squelch,  ".png", sep = "")
  #confounderSource,squelch, "confounders ","(pos: ",pos,", neg: ",neg,"total: ",tot,")"),
  # perfRor # perfSql5medcoef perfSql10medcoef perfPsi5medcoef perfPsi10medcoef 
  # perfSql5ate perfSql10ate perfPsi5ate perfPsi10ate  
  # rorColor chi2Color coefSql5Color coefSql10Color coefPsi5Color coefPsi10Color
  # ateSql5Color ateSql10Color atePsi5Color atePsi10Color
  plot(main = paste(
        perfRor,col = rorColor,print.auc.y = .3,legend = TRUE))
  plot(perfSql5medcoef, add = TRUE, col = coefSql5Color, print.auc.y = .3)
  plot(perfSql10medcoef, add = TRUE, col = coefSql10Color, print.auc.y = .3)
  plot(perfPsi5medcoef, add = TRUE, col = coefPsi5Color, print.auc.y = .3)
  plot(perfPsi10medcoef, add = TRUE, col = coefPsi10Color, print.auc.y = .3)
  
  
  plot(perfSql5ate, add = TRUE, col = ateSql5Color, print.auc.y = .3)
  plot(perfSql10ate, add = TRUE, col = ateSql10Color, print.auc.y = .3)
  plot(perfPsi5ate, add = TRUE, col = atePsi5Color, print.auc.y = .3)
  plot(perfPsi10ate, add = TRUE, col = atePsi10Color, print.auc.y = .3)
  
  print(xtable(dat[, c(1,3,8,9,10,11,14,15,12,13,16,17)]))
  xtable(data.frame(dat$casecontrol, dat$exposurename, round(dat[,"sql5ate"],4), round(dat[,"sql10ate"],4)), digits = 4)
}
drawAUPRC('kidney_failure,_acute')


dbGetQuery(con, "SELECT * FROM penkiModels LIMIT 1")

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
  statName = "AUPRC"
  fn <- paste("PerformanceStatsAndCurves/", hoiname, "-", statName, "_", squelch,  ".png", sep = "")
  png(fn)
  nerding <-
    plot(
      main = paste(
        #confounderSource,
        squelch, "confounders ",
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
  dev.off()
  
}
drawAUPRC('kidney_failure,_acute')


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
        "' AND j.casecontrol = r.casecontrol AND -- a >= 10 AND (a+b) >= 100 AND
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
  statName = "AUPRC"
  fn <- paste("PerformanceStatsAndCurves/", hoiname, "-", statName, "_", squelch,  ".png", sep = "")
  ##png(fn)
  nerding <-
    plot(
      main = paste(
        #confounderSource,
        squelch, "confounders ",
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
  ##dev.off()
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

##############
# BY ADE
##############

drawAUROC <- function(squelch, confounderSource, hoiname) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue, j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND  
                  j.hoiname LIKE '", hoiname, "' AND j.casecontrol = r.casecontrol AND
                  j.squelch = ", squelch, " AND confoundersource = '", confounderSource, "';", sep = ""))
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
  
  pos <- print(sum(dat$casecontrol))
  tot <- nrow(dat)
  neg <- tot - pos
  # rocSQL5Chisq <- plot(main = paste(hoiname,confounderSource, squelch), smooth(roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark grey")
  # rocSQL5Medcoef <- plot(smooth(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE)
  # rocSQL5exactate <- plot(smooth(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
  statName = "AUROC"
  fn <- paste("PerformanceStatsAndCurves/", hoiname, "-", statName, "_", squelch,  ".png", sep = "")
  png(fn)
  rocSQL5Chisq <- plot(main = paste(
    #confounderSource,
    squelch, "confounders ",
    "(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = baselineColor, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = glmColor, print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = gcmColor, print.auc.y = .3, add = TRUE)
  dev.off()
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




