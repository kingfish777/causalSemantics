library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
drv <- dbDriver('PostgreSQL')
con <- dbConnect(drv, dbname = "confounders", host = "localhost", port = 5432, user = "smalec", password = "mandarin")



##############
# OVERALL
##############


drawAUROC <- function(squelch, confounderSource) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue, j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND -- a >= 100 AND -- j.medpvalue0 < .05 AND 
                    j.casecontrol = r.casecontrol AND j.squelch = ", squelch, " AND confoundersource = '", confounderSource, "';", sep = ""))
  
  
  rocSQL5Chisq <- plot(main = paste(confounderSource, squelch), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark grey", colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                         col = "dark red", print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "red", print.auc.y = .3, add = TRUE)
}
drawAUROC(5, "sql")
drawAUROC(10, "sql")
drawAUROC(5, "psi")
drawAUROC(10, "psi")

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
  rocSQL5Chisq <- plot(main = paste(confounderSource, squelch, "(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark grey", colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
  
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
drawAUROC(10,"psi", "%")

##############
# BY ADE [PRIMO]
##############

drawAUROC <- function(squelch, confounderSource, hoiname) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue, j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND (r.a+r.b) > 100 AND a >= 50 AND -- j.medpvalue0 < .05 AND 
                  j.hoiname LIKE '", hoiname, "' AND j.casecontrol = r.casecontrol AND 
                  j.squelch = ", squelch, " AND confoundersource = '", confounderSource, "';", sep = ""))
  
  pos <- print(sum(dat$casecontrol))
  tot <- nrow(dat)
  neg <- tot - pos
  # rocSQL5Chisq <- plot(main = paste(hoiname,confounderSource, squelch), smooth(roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark grey")
  # rocSQL5Medcoef <- plot(smooth(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE)
  # rocSQL5exactate <- plot(smooth(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
  rocSQL5Chisq <- plot(main = paste(confounderSource, squelch, "(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark grey", colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
  
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

drawAUROC(10, "psi", "%")












ggroc(list("sql 5 chisq" = rocSQL5Chisq, "sql 5 coef" = rocSQL5Medcoef, "sql 5 ATE" = rocSQL5exactate))


##############
# SQL 10 OVERALL
##############
dat <- dbGetQuery(con, "SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, j.chisqpvalue, j.medcoef0, j.medpvalue, j.medcoef, j.medpvalue, j.exactate, j.exactate, j.exactatevar FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND -- a >= 100 AND -- j.medpvalue0 < .05 AND 
                    j.casecontrol = r.casecontrol AND j.squelch = 10 AND confoundersource = 'sql';")


rocSQL5 <- plot(roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), print.auc = TRUE, col = "dark grey")
rocSQL5 <- plot(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                col = "dark red", print.auc.y = .4, add = TRUE)
rocSQL5 <- plot(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                col = "red", print.auc.y = .3, add = TRUE)
rocSQL5 <- plot(roc(dat$casecontrol, (dat$exactate - 100000*dat$exactatevar, algorithm = 4), ci = TRUE, algorithm = 4), print.auc = TRUE, 
                col = "orange", print.auc.y = .2, add = TRUE)








roc_rose <- plot(roc(hacide.test$cls, pred_rose[,2]), print.auc = TRUE, col = "blue")

And for the second ROC curve you need to change the y position of the AUC and use add the plot the two curves on the same plot:
  
  roc_rose <- plot(roc(hacide.test$cls, pred_both[,2]), print.auc = TRUE, 
                   col = "green", print.auc.y = .4, add = TRUE)




preds <- cbind(ror = dat$ror, chisq = dat$chisqpvlue, exactate = dat$exactate)
pred.mat <- prediction(preds, labels = matrix(dat$casecontrol, nrow = length(dat$casecontrol), ncol = 2))


data(ROCR.simple)
preds <- cbind(p1 = ROCR.simple$predictions, 
               p2 = abs(ROCR.simple$predictions + 
                          rnorm(length(ROCR.simple$predictions), 0, 0.1)))

pred.mat <- prediction(preds, labels = matrix(ROCR.simple$labels, 
                                              nrow = length(ROCR.simple$labels), ncol = 2) )
perf.mat <- performance(pred.mat, "tpr", "fpr")
plot(perf.mat)
df <- data.frame(Curve=as.factor(rep(c(1,2), each=length(perf.mat@x.values[[1]]))), 
                 FalsePositive=c(perf.mat@x.values[[1]],perf.mat@x.values[[2]]),
                 TruePositive=c(perf.mat@y.values[[1]],perf.mat@y.values[[2]]))
plt <- ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Curve)) + geom_line()
print(plt)

  
dat <- dbGetQuery(con, "SELECT * FROM jbiscores j, refsetbaselines r
                  WHERE j.hoiname = r.hoiname AND j.exposurename = r.exposurename AND 
                    j.casecontrol = r.casecontrol AND j.squelch = 5 AND confoundersource = 'sql';")

aucOR <- auc(dat$casecontrol, dat$chisqpvalue)
auc(dat$casecontrol, dat$ror)
auc(dat$casecontrol, dat$medcoef0)
auc(dat$casecontrol, dat$exactate)
library(MLmetrics)
subset(dat, subset = hoiname == "acute_myocardial_infarction")







dat <- dbGetQuery(con, "SELECT DISTINCT e1.casecontrol, e1.hoiname, e1.exposurename, e1.exactATE, e1.medcoef0, e1.medcoef1, 
                              ref.a, ref.b, ref.c, ref.ror
                FROM  jbiscores enkiModels e2, refsetbaselines ref
                 WHERE e1.hoiname = e2.hoiname AND e1.hoiname = ref.hoiname AND ref.a >= 100 AND ref.ror >= 3 AND
                        e1.exposurename = e2.exposurename AND e1.exposurename = ref.exposurename AND
                        e1.casecontrol = e2.casecontrol AND CAST(e1.casecontrol AS INT) = ref.casecontrol AND
                        e1.confoundersource = 'psi' AND e2.confoundersource = 'sql';")



dat <- dbGetQuery(con, "SELECT DISTINCT e1.casecontrol, e1.hoiname, e1.exposurename, e1.exactATE, e1.medcoef0, e1.medcoef1, ref.a, ref.b, ref.c, ref.ror
                FROM enkiModels e1, enkiModels e2, refsetbaselines ref
                 WHERE e1.hoiname = e2.hoiname AND e1.hoiname = ref.hoiname AND ref.a >= 100 AND ref.ror >= 3 AND
                        e1.exposurename = e2.exposurename AND e1.exposurename = ref.exposurename AND
                        e1.casecontrol = e2.casecontrol AND CAST(e1.casecontrol AS INT) = ref.casecontrol AND
                        e1.confoundersource = 'psi' AND e2.confoundersource = 'sql';")

datpsi <- dbGetQuery(con, "SELECT DISTINCT e2.casecontrol, e2.hoiname, e2.exposurename, e2.exactATE, e2.medcoef0, e2.medcoef1, ref.a, ref.b, ref.c, ref.ror
                FROM enkiModels e1, enkiModels e2, refsetbaselines ref
                  WHERE e1.hoiname = e2.hoiname AND e1.hoiname = ref.hoiname AND ref.a >= 100 AND ref.ror >= 3 AND 
                  e1.exposurename = e2.exposurename AND e1.exposurename = ref.exposurename AND
                  e1.casecontrol = e2.casecontrol AND CAST(e1.casecontrol AS INT) = ref.casecontrol AND
                  e1.confoundersource = 'psi' AND e2.confoundersource = 'sql' AND ref.refsetid = 'omop';")

nDat <- nrow(dat); 
nCase <- sum(as.numeric(dat$casecontrol)); 
