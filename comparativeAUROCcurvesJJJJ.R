
library(pROC)
library(PRROC)
library(ROCR)
library(ggplot2)
library(RPostgreSQL)
drv <- dbDriver('PostgreSQL')
con <- dbConnect(drv, dbname = "confounders", host = "localhost", port = 5432, user = "research", password = "mandarin")

dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND (a+b) >= 100 AND a >= 10 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))

rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                     print.auc = TRUE, col = "dark grey", colorize = TRUE)

rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql5medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "gold", print.auc.y = .4, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql10medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "dark orange", print.auc.y = .3, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi5medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "red", print.auc.y = .2, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi10medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "dark red", print.auc.y = .1, add = TRUE, colorize = TRUE)


rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                     print.auc = TRUE, col = "dark grey", colorize = TRUE)

rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "light green", print.auc.y = .4, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "green", print.auc.y = .3, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "blue", print.auc.y = 0.2, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "dark blue", print.auc.y = .1, add = TRUE)
##############
# BY ADE
##############

getCPERF <- function(hoiname) {
dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND
                                          r.hoiname = '", hoiname, "' AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND --(a+b) >= 100 AND a >= 10 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))

rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC) ", hoiname), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                     print.auc = TRUE, col = "dark grey", colorize = TRUE)

rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql5medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "gold", print.auc.y = .4, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql10medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "dark orange", print.auc.y = .3, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi5medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "red", print.auc.y = .2, add = TRUE, colorize = TRUE)
rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi10medcoef, ci = TRUE), print.auc = TRUE, 
                       col = "dark red", print.auc.y = .1, add = TRUE, colorize = TRUE)


#rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
#                     print.auc = TRUE, col = "dark grey", colorize = TRUE)

rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "light green", print.auc.y = .45, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "green", print.auc.y = .35, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "blue", print.auc.y = 0.25, add = TRUE)
rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                        col = "dark blue", print.auc.y = .15, add = TRUE)
}

getCPERF('kidney_failure,_acute')
getCPERF('liver_failure,_acute')
getCPERF('acute_myocardial_infarction')
getCPERF('gastrointestinal_hemorrhage')



getCPERF <- function(hoiname, venn) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND
                                          r.hoiname = '", hoiname, "' AND r.a >= ", venn," AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND --(a+b) >= 100 AND a >= 10 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))
  
  rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC) ", hoiname), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                       print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql5medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "gold", print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql10medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "dark orange", print.auc.y = .3, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi5medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "red", print.auc.y = .2, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi10medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "dark red", print.auc.y = .1, add = TRUE, colorize = TRUE)
  
  
  #rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
  #                     print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "light green", print.auc.y = .45, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "green", print.auc.y = .35, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "blue", print.auc.y = 0.25, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "dark blue", print.auc.y = .15, add = TRUE)
}

getCPERF('kidney_failure,_acute', 25)
getCPERF('liver_failure,_acute', 25)
getCPERF('acute_myocardial_infarction', 25)
getCPERF('gastrointestinal_hemorrhage', 25)




getCPERF <- function(hoiname, venn) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND
                                          r.hoiname = '", hoiname, "' AND r.a >= ", venn," AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND --(a+b) >= 100 AND a >= 10 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))
  
  rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC) ", hoiname), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                       print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql5medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "gold", print.auc.y = .4, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql10medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "dark orange", print.auc.y = .3, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi5medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "red", print.auc.y = .2, add = TRUE, colorize = TRUE)
  rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi10medcoef, ci = TRUE), print.auc = TRUE, 
                         col = "dark red", print.auc.y = .1, add = TRUE, colorize = TRUE)
  
  
  #rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
  #                     print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  #rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
  #                        col = "light green", print.auc.y = .45, add = TRUE)
  #rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
  #                        col = "green", print.auc.y = .35, add = TRUE)
  #rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
  #                        col = "blue", print.auc.y = 0.25, add = TRUE)
  #rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
  #                        col = "dark blue", print.auc.y = .15, add = TRUE)
}

getCPERF('kidney_failure,_acute', 25)
getCPERF('liver_failure,_acute', 25)
getCPERF('acute_myocardial_infarction', 25)
getCPERF('gastrointestinal_hemorrhage', 25)


getCPERF <- function(hoiname, venn) {
  dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND
                                          r.hoiname = '", hoiname, "' AND r.a >= ", venn," AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND --(a+b) >= 100 AND a >= 10 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))
  
  rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC) ", hoiname), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
                       print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  #rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql5medcoef, ci = TRUE), print.auc = TRUE, 
  #                       col = "gold", print.auc.y = .4, add = TRUE, colorize = TRUE)
  #rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$sql10medcoef, ci = TRUE), print.auc = TRUE, 
  #                       col = "dark orange", print.auc.y = .3, add = TRUE, colorize = TRUE)
  #rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi5medcoef, ci = TRUE), print.auc = TRUE, 
  #                       col = "red", print.auc.y = .2, add = TRUE, colorize = TRUE)
  #rocSQL5Medcoef <- plot(roc(dat$casecontrol, dat$psi10medcoef, ci = TRUE), print.auc = TRUE, 
  #                       col = "dark red", print.auc.y = .1, add = TRUE, colorize = TRUE)
  
  
  #rocSQL5Chisq <- plot(main = paste("Comparative Performance (AUROC)"), roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4), 
  #                     print.auc = TRUE, col = "dark grey", colorize = TRUE)
  
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "light green", print.auc.y = .45, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$sql10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "green", print.auc.y = .35, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi5exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "blue", print.auc.y = 0.25, add = TRUE)
  rocSQL5exactate <- plot(roc(dat$casecontrol, dat$psi10exactate, ci = TRUE, algorithm = 4), print.auc = TRUE, 
                          col = "dark blue", print.auc.y = .15, add = TRUE)
}

getCPERF('kidney_failure,_acute', 25)
getCPERF('liver_failure,_acute', 25)
getCPERF('acute_myocardial_infarction', 25)
getCPERF('gastrointestinal_hemorrhage', 25)

