library(wordcloud)
library(tm)
library(stringr)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

source("helperFunctions.R")

drv <- dbDriver('PostgreSQL')
con <- dbConnect(drv, dbname = "confounders", host = "localhost", port = 5432, user = "smalec", password = "mandarin")


dat <- dbGetQuery(con, "SELECT exposurename, hoiname, casecontrol 
                          FROM refsetbaselines 
                          WHERE a >= 5 AND refsetid = 'omop' 
                            AND hoiname = 'gastrointestinal_hemorrhage'
                  ORDER BY casecontrol DESC, hoiname DESC, exposurename DESC;")
confounderList <- c()
for (i in 1:nrow(dat)) {
  drow <- dat[i,]
  cl <- "sql"
  drug <- as.character(drow[1])
  ade <- as.character(drow[2])
  adr <- as.character(drow[2])
  group <- as.character(drow[3])
  SQUELCH <- "10" #as.character(10)
  try({ fn <- paste(outputResultsFolder(cl, drug, ade, adr, group, SQUELCH), "/rawConfounders.txt", sep = "")
  print(fn)
  confounders <- unlist(read.csv2(fn, header = FALSE, stringsAsFactors = FALSE))
  confounders <- paste0(confounders, collapse = ", ")
  print(confounders)
  confounderList <- c(confounders, confounderList)})
  # try({bmc(cl, drug, ade, adr, group, SQUELCH)})
  #} else if (doneAlready != 0) { print("DID THAT ALREADY") }
}

text <- confounderList # readLines(filePath)
docs <- Corpus(VectorSource(text))
inspect(docs)
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
# docs <- tm_map(docs, stemDocument)
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Most frequent words",
        ylab = "Word frequencies")





###################
##
###################

getWordcloud <- function(hoiname) {
  
  dat <- dbGetQuery(con, paste("SELECT exposurename, hoiname, casecontrol 
                          FROM refsetbaselines 
                          WHERE a >= 5 AND refsetid = 'omop' 
                            AND hoiname = '", hoiname, "' 
                  ORDER BY casecontrol DESC, hoiname DESC, exposurename DESC;", sep = ""))
  confounderList <- c()
  for (i in 1:nrow(dat)) {
    drow <- dat[i,]
    cl <- "sql"
    drug <- as.character(drow[1])
    ade <- as.character(drow[2])
    adr <- as.character(drow[2])
    group <- as.character(drow[3])
    SQUELCH <- "10" #as.character(10)
    try({ fn <- paste(outputResultsFolder(cl, drug, ade, adr, group, SQUELCH), "/rawConfounders.txt", sep = "")
    print(fn)
    confounders <- unlist(read.csv2(fn, header = FALSE, stringsAsFactors = FALSE))
    confounders <- paste0(confounders, collapse = ", ")
    print(confounders)
    confounderList <- c(confounders, confounderList)})
    # try({bmc(cl, drug, ade, adr, group, SQUELCH)})
    #} else if (doneAlready != 0) { print("DID THAT ALREADY") }
  }
  
  text <- confounderList # readLines(filePath)
  docs <- Corpus(VectorSource(text))
  inspect(docs)
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  head(d, 10)
  set.seed(1232)
  wordcloud(words = d$word, freq = d$freq, min.freq = 5,
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))
  
 # barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
 #         col ="lightblue", main ="Most frequent words",
 #         ylab = "Word frequencies")
}

getWordcloud("kidney_failure,_acute")
getWordcloud("liver_failure,_acute")
getWordcloud("acute_myocardial_infarction")
getWordcloud("gastrointestinal_hemorrhage")

getWordcloud <- function(hoiname, casecontrol) {
  
  dat <- dbGetQuery(con, paste("SELECT exposurename, hoiname, casecontrol 
                          FROM refsetbaselines 
                          WHERE a >= 5 AND refsetid = 'omop' AND casecontrol = '", casecontrol, "' 
                            AND hoiname = '", hoiname, "' 
                  ORDER BY casecontrol DESC, hoiname DESC, exposurename DESC;", sep = ""))
  confounderList <- c()
  for (i in 1:nrow(dat)) {
    drow <- dat[i,]
    cl <- "sql"
    drug <- as.character(drow[1])
    ade <- as.character(drow[2])
    adr <- as.character(drow[2])
    group <- as.character(drow[3])
    SQUELCH <- "10" #as.character(10)
    try({ fn <- paste(outputResultsFolder(cl, drug, ade, adr, group, SQUELCH), "/rawConfounders.txt", sep = "")
    print(fn)
    confounders <- unlist(read.csv2(fn, header = FALSE, stringsAsFactors = FALSE))
    confounders <- paste0(confounders, collapse = ", ")
    confounders <- str_replace_all(string = confounders, pattern = "adverse event", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "effects", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "adverse", replacement = "")
    confounderList <- c(confounders, confounderList)})
    # try({bmc(cl, drug, ade, adr, group, SQUELCH)})
    #} else if (doneAlready != 0) { print("DID THAT ALREADY") }
  }
  
  text <- confounderList # readLines(filePath)
  docs <- Corpus(VectorSource(text))
  inspect(docs)
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  head(d, 10)
  set.seed(1232)
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))
  
  # barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
  #         col ="lightblue", main ="Most frequent words",
  #         ylab = "Word frequencies")
}

getWordcloud("kidney_failure,_acute", 0)
getWordcloud("kidney_failure,_acute", 1)
getWordcloud("liver_failure,_acute", 0)
getWordcloud("liver_failure,_acute", 1)
getWordcloud("acute_myocardial_infarction", 0)
getWordcloud("acute_myocardial_infarction", 1)
getWordcloud("gastrointestinal_hemorrhage", 0)
getWordcloud("gastrointestinal_hemorrhage", 1)


getWordcloud <- function(hoiname, casecontrol) {
  
  dat <- dbGetQuery(con, paste("SELECT exposurename, hoiname, casecontrol 
                          FROM refsetbaselines 
                          WHERE a >= 5 AND refsetid = 'omop' AND casecontrol = '", casecontrol, "' 
                            AND hoiname = '", hoiname, "' 
                  ORDER BY casecontrol DESC, hoiname DESC, exposurename DESC;", sep = ""))
  confounderList <- c()
  for (i in 1:nrow(dat)) {
    drow <- dat[i,]
    cl <- "sql"
    drug <- as.character(drow[1])
    ade <- as.character(drow[2])
    adr <- as.character(drow[2])
    group <- as.character(drow[3])
    SQUELCH <- "10" #as.character(10)
    try({ fn <- paste(outputResultsFolder(cl, drug, ade, adr, group, SQUELCH), "/trueConfounders.txt", sep = "")
    print(fn)
    confounders <- unlist(read.csv2(fn, header = FALSE, stringsAsFactors = FALSE))
    confounders <- paste0(confounders, collapse = ", ")
    confounders <- str_replace_all(string = confounders, pattern = "adverse event", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "effects", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "adverse", replacement = "")
    confounderList <- c(confounders, confounderList)})
    # try({bmc(cl, drug, ade, adr, group, SQUELCH)})
    #} else if (doneAlready != 0) { print("DID THAT ALREADY") }
  }
  
  text <- confounderList # readLines(filePath)
  docs <- Corpus(VectorSource(text))
  inspect(docs)
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  head(d, 10)
  set.seed(1232)
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))
  
  # barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
  #         col ="lightblue", main ="Most frequent words",
  #         ylab = "Word frequencies")
}

getWordcloud("kidney_failure,_acute", 0)
getWordcloud("kidney_failure,_acute", 1)
getWordcloud("liver_failure,_acute", 0)
getWordcloud("liver_failure,_acute", 1)
getWordcloud("acute_myocardial_infarction", 0)
getWordcloud("acute_myocardial_infarction", 1)
getWordcloud("gastrointestinal_hemorrhage", 0)
getWordcloud("gastrointestinal_hemorrhage", 1)


getWordcloud <- function(hoiname, casecontrol) {
  
  dat <- dbGetQuery(con, paste("SELECT exposurename, hoiname, casecontrol 
                          FROM refsetbaselines 
                          WHERE a >= 5 AND refsetid = 'omop' AND casecontrol = '", casecontrol, "' 
                            AND hoiname = '", hoiname, "' 
                  ORDER BY casecontrol DESC, hoiname DESC, exposurename DESC;", sep = ""))
  confounderList <- c()
  for (i in 1:nrow(dat)) {
    drow <- dat[i,]
    cl <- "sql"
    drug <- as.character(drow[1])
    ade <- as.character(drow[2])
    adr <- as.character(drow[2])
    group <- as.character(drow[3])
    SQUELCH <- "10" #as.character(10)
    try({ fn <- paste(outputResultsFolder(cl, drug, ade, adr, group, SQUELCH), "/bestConfounders.txt", sep = "")
    print(fn)
    confounders <- unlist(read.csv2(fn, header = FALSE, stringsAsFactors = FALSE))
    confounders <- paste0(confounders, collapse = ", ")
    confounders <- str_replace_all(string = confounders, pattern = "adverse event", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "effects", replacement = "")
    confounders <- str_replace_all(string = confounders, pattern = "adverse", replacement = "")
    confounderList <- c(confounders, confounderList)})
    # try({bmc(cl, drug, ade, adr, group, SQUELCH)})
    #} else if (doneAlready != 0) { print("DID THAT ALREADY") }
  }
  
  text <- confounderList # readLines(filePath)
  docs <- Corpus(VectorSource(text))
  inspect(docs)
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Text stemming
  # docs <- tm_map(docs, stemDocument)
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  head(d, 10)
  #set.seed(1232)
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))
  
  # barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
  #         col ="lightblue", main ="Most frequent words",
  #         ylab = "Word frequencies")
}

getWordcloud("kidney_failure,_acute", 0)
getWordcloud("kidney_failure,_acute", 1)
getWordcloud("liver_failure,_acute", 0)
getWordcloud("liver_failure,_acute", 1)
getWordcloud("acute_myocardial_infarction", 0)
getWordcloud("acute_myocardial_infarction", 1)
getWordcloud("gastrointestinal_hemorrhage", 0)
getWordcloud("gastrointestinal_hemorrhage", 1)
