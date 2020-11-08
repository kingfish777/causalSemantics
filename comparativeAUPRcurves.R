


dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND r.ror >= 1 AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND (a+b) >= 100 AND a >= 25 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))


pos <- print(sum(dat$casecontrol))
tot <- nrow(dat)
neg <- tot - pos
pred0 <- prediction(dat$chisqpvalue, dat$casecontrol)
perf0 <- performance (pred0, "prec", "rec")
pred1 <- prediction(dat$sql5medcoef, dat$casecontrol)
perf1 <- performance (pred1, "prec", "rec")
pred2 <- prediction(dat$psi10exactatevar, dat$casecontrol)
perf2 <- performance (pred2, "prec", "rec")
nerding <- plot(main = paste("(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), perf0, col = "dark grey", print.auc.y = .3)
nerding <- plot(perf1, add = TRUE, col = "red", print.auc.y = .3)
nerding <- plot(perf2, add = TRUE, col = "dark red", print.auc.y = .3)


dd1 <- subset(dat, select = casecontrol == 1)
dd0 <- subset(dat, select = casecontrol == 0)



predBaseline <- prediction(dd$chisqpvalue, dd$casecontrol)
perfBaseline <- performance(predBaseline, "prec", "rec")
predSQL5coef <- prediction(dd$sql5medcoef, dd$casecontrol)
perfSQL5coef <- performance(predSQL5coef, "prec", "rec")
predSQL10coef <- prediction(dd$sql10medcoef, dd$casecontrol)
perfSQL10coef <- performance(predSQL10coef, "prec", "rec")
predPSI5coef <- prediction(dd$psi5medcoef, dd$casecontrol)
perfPSI5coef <- performance(predPSI5coef, "prec", "rec")
predPSI10coef <- prediction(dd$psi10medcoef, dd$casecontrol)
perfPSI10coef <- performance(predPSI10coef, "prec", "rec")


predBaseline <- prediction(dd$chisqpvalue, dd$casecontrol)
perfBaseline <- performance(predBaseline, "prec", "rec")
predSQL5ate <- prediction(dd$sql5exactate, dd$casecontrol)
perfSQL5ate <- performance(predSQL5ate, "prec", "rec")
predSQL10ate <- prediction(dd$sql10exactate, dd$casecontrol)
perfSQL10ate <- performance(predSQL10ate, "prec", "rec")
predPSI5ate <- prediction(dd$psi5exactate, dd$casecontrol)
perfPSI5ate <- performance(predPSI5coef, "prec", "rec")
predPSI10ate <- prediction(dd$psi10exactate, dd$casecontrol)
perfPSI10ate <- performance(predPSI10ate, "prec", "rec")


plot(perfBaseline, col = 'dark grey')
plot(perfSQL5coef, col = 'yellow', add = TRUE)
plot(perfSQL10coef, col = 'orange', add = TRUE)
plot(perfPSI5coef, col = 'red', add = TRUE)
plot(perfPSI10coef, col = 'dark red', add = TRUE)

plot(perfBaseline, col = 'dark grey')
plot(perfSQL5ate, col = 'yellow', add = TRUE)
plot(perfSQL10ate, col = 'orange', add = TRUE)
plot(perfPSI5ate, col = 'red', add = TRUE)
plot(perfPSI10ate, col = 'dark red', add = TRUE)


###########################


getAUPRC <- function(hoiname, squelch) {

dat <- dbGetQuery(con, paste("SELECT DISTINCT r.exposurename, r.hoiname, r.casecontrol, r.ror, r.a, 
                                      sql5.chisqpvalue, sql5.medcoef0 AS sql5medcoef0, sql5.medpvalue0 AS sql5medpvalue0, sql5.medcoef AS sql5medcoef, sql5.medpvalue AS sql5medpvalue, sql5.exactate AS sql5exactate, sql5.exactatevar AS sql5exactatevar,
                                      sql10.medcoef0 AS sql10medcoef0, sql10.medpvalue0 AS sql10medpvalue0, sql10.medcoef AS sql10medcoef, sql10.medpvalue AS sql10medpvalue, sql10.exactate AS sql10exactate, sql10.exactatevar AS sql10exactatevar,
                                      psi5.medcoef0 AS psi5medcoef0, psi5.medpvalue0 AS psi5medpvalue0, psi5.medcoef AS psi5medcoef, psi5.medpvalue AS psi5medpvalue, psi5.exactate AS psi5exactate, psi5.exactatevar AS psi5exactatevar,
                                      psi10.medcoef0 AS psi10medcoef0, psi10.medpvalue AS psi10medpvalue, psi10.medcoef AS psi10medcoef, psi10.medpvalue AS psi10medpvalue, psi10.exactate AS psi10exactate, psi10.exactatevar AS psi10exactatevar
                                    FROM refsetbaselines r, jbiscores sql5, jbiscores sql10, jbiscores psi5, jbiscores psi10
                                    WHERE r.hoiname = sql5.hoiname AND r.exposurename = sql5.exposurename AND r.casecontrol = sql5.casecontrol 
                                          AND sql5.squelch = 5 AND sql5.confoundersource = 'sql' AND 
                                          r.hoiname LIKE '", hoiname, "' AND 
                                          r.hoiname = sql10.hoiname AND r.exposurename = sql10.exposurename AND r.casecontrol = sql10.casecontrol 
                                          AND sql10.squelch = 10 AND sql10.confoundersource = 'sql' AND r.ror >= 1 AND
                                          r.hoiname = psi5.hoiname AND r.exposurename = psi5.exposurename AND r.casecontrol = psi5.casecontrol 
                                          AND psi5.squelch = 5 AND psi5.confoundersource = 'psi' AND -- (a+b) >= 100 AND a >= 25 AND
                                          r.hoiname = psi10.hoiname AND r.exposurename = psi10.exposurename AND r.casecontrol = psi10.casecontrol 
                                          AND psi10.squelch = 10 AND psi10.confoundersource = 'psi';", sep = ""))


pos <- print(sum(dat$casecontrol))
tot <- nrow(dat)
neg <- tot - pos
pred0 <- prediction(dat$chisqpvalue, dat$casecontrol)
perf0 <- performance (pred0, "prec", "rec")
pred1 <- prediction(dat$sql5medcoef, dat$casecontrol)
perf1 <- performance (pred1, "prec", "rec")
pred2 <- prediction(dat$psi10exactatevar, dat$casecontrol)
perf2 <- performance (pred2, "prec", "rec")
nerding <- plot(main = paste("(pos: ", pos, ", neg: ", neg, "total: ", tot, ")"), perf0, col = "dark grey", print.auc.y = .3)
nerding <- plot(perf1, add = TRUE, col = "red", print.auc.y = .3)
nerding <- plot(perf2, add = TRUE, col = "dark red", print.auc.y = .2)

dd <- dat # subset(dat, select = psi5exactate > squelch)
predBaseline <- prediction(dd$chisqpvalue, dd$casecontrol)
perfBaseline <- performance(predBaseline, "prec", "rec")
predSQL5coef <- prediction(dd$sql5medcoef, dd$casecontrol)
perfSQL5coef <- performance(predSQL5coef, "prec", "rec")
predSQL10coef <- prediction(dd$sql10medcoef, dd$casecontrol)
perfSQL10coef <- performance(predSQL10coef, "prec", "rec")
predPSI5coef <- prediction(dd$psi5medcoef, dd$casecontrol)
perfPSI5coef <- performance(predPSI5coef, "prec", "rec")
predPSI10coef <- prediction(dd$psi10medcoef, dd$casecontrol)
perfPSI10coef <- performance(predPSI10coef, "prec", "rec")


pos <- print(sum(dat$casecontrol))
tot <- nrow(dat)
neg <- tot - pos
# rocSQL5Chisq <- plot(main = paste(hoiname,confounderSource, squelch), smooth(roc(dat$casecontrol, dat$chisqpvalue, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark grey")
# rocSQL5Medcoef <- plot(smooth(roc(dat$casecontrol, dat$medcoef, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "dark red", print.auc.y = .4, add = TRUE)
# rocSQL5exactate <- plot(smooth(roc(dat$casecontrol, dat$exactate, ci = TRUE, algorithm = 4)), print.auc = TRUE, col = "red", print.auc.y = .3, add = TRUE)
mainTitle = paste("(pos: ", pos, ", neg: ", neg, " total: ", tot, ")", sep = "")


predBaseline <- prediction(dd$chisqpvalue, dd$casecontrol)
perfBaseline <- performance(predBaseline, "prec", "rec")
predSQL5ate <- prediction(dd$sql5exactate, dd$casecontrol)
perfSQL5ate <- performance(predSQL5ate, "prec", "rec")
predSQL10ate <- prediction(dd$sql10exactate, dd$casecontrol)
perfSQL10ate <- performance(predSQL10ate, "prec", "rec")
predPSI5ate <- prediction(dd$psi5exactate, dd$casecontrol)
perfPSI5ate <- performance(predPSI5coef, "prec", "rec")
predPSI10ate <- prediction(dd$psi10exactate, dd$casecontrol)
perfPSI10ate <- performance(predPSI10ate, "prec", "rec")


plot(perfBaseline, col = 'grey', main = mainTitle)
#plot(perfSQL5coef, col = 'dark grey', add = TRUE)
#plot(perfSQL10coef, col = 'brown', add = TRUE)
#plot(perfPSI5coef, col = 'orange', add = TRUE)
#plot(perfPSI10coef, col = 'dark orange', add = TRUE)

#plot(perfBaseline, col = 'grey', add = TRUE)
plot(perfSQL5ate, col = 'dark grey', add = TRUE)
plot(perfSQL10ate, col = 'brown', add = TRUE)
plot(perfPSI5ate, col = 'dark red', add = TRUE)
plot(perfPSI10ate, col = 'red', add = TRUE)

}

getAUPRC('kidney_failure,_acute', .02)
getAUPRC('liver_failure,_acute', .02) 
getAUPRC('acute_myocardial_infarction', .02) 
getAUPRC('gastrointestinal_hemorrhage', .02) 

getAUPRC('liver_failure,_acute', .3)
getAUPRC('liver_failure,_acute', .2) 
getAUPRC('liver_failure,_acute', .1) 
getAUPRC('liver_failure,_acute', .05) 

getAUPRC('acute_myocardial_infarction', .3)
getAUPRC('acute_myocardial_infarction', .2) 
getAUPRC('acute_myocardial_infarction', .1) 
getAUPRC('acute_myocardial_infarction', .05)
