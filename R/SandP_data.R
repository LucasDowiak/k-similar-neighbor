setwd("~/Git/k-similar-neighbor/")
source("R/alpha_adv_api_key.R")
source("R/auto_marginal.R")
library(rvest)
library(RCurl)
library(jsonlite)

# Pull List of S&P 500 companies and their stock tickers from wikipedia
u <- "https://en.wikipedia.org/wiki/List_of_S&P_500_companies"
doc <- read_html(u)
tables <- html_nodes(doc, "table")
dtfCos <- html_table(tables[1])[[1]]
names(dtfCos) <- c("ticker", "security", "GICS_sector", "GICS_sub_ind",
                   "hq", "date_added", "cik", "founded")


# Use stock ticker to pull data from alpha advantage api
hit_api <- function(tickers, outputsize=c("compact", "full"), saveJSON=FALSE)
{
  
  outputsize <- match.arg(outputsize)
  
  save_file <- function(json, tick)
  {
    dir_ <- "/Users/lucasdowiak/Git/k-similar-neighbor/data/raw_json/"
    file_ <- paste0(tick, ".json")
    write(minify(json), file=paste0(dir_, file_))
  }
  
  baseuri <- "https://www.alphavantage.co/query"
  parstring <- "function=TIME_SERIES_DAILY_ADJUSTED&symbol=%s&outputsize=%s&apikey=%s"
  uri <- paste(baseuri, 
               sprintf(parstring, tickers, outputsize, .APIKEY),
               sep="?")
  
  if (length(uri) > 1) {
    jsonText <- getURIAsynchronous(uri)
  } else {
    jsonText <- getForm(uri)
  }
  
  if (saveJSON) {
    mapply(save_file, jsonText, tickers)
  }
  
  jsonText
}


# Use keywords to search for available data
ticker_search <- function(keywords)
{
  baseuri <- "https://www.alphavantage.co/query"
  parstring <- "function=SYMBOL_SEARCH&keywords=%s&apikey=%s"
  uri <- paste(baseuri, 
               sprintf(parstring, keywords, .APIKEY),
               sep="?")
  
  if (length(uri) > 1) {
    jsonText <- getURIAsynchronous(uri)
  } else {
    jsonText <- getForm(uri)
  }
  return(jsonText)
}


list_all_available <- function()
{
  baseuri <- "https://www.alphavantage.co/query"
  parstring <- "function=LISTING_STATUS&apikey=%s"
  uri <- paste(baseuri, 
               sprintf(parstring, .APIKEY),
               sep="?")
  
  jsonText <- getForm(uri)
  return(read.csv(textConnection(rawToChar(aa))))
}


# Run through list of stock ticks and obtain data
# ------------------------------------------------------------------------------

for (tick in dtfCos$ticker) {
  cat('\r', sprintf("Tick: %s - %d / %d", tick, which(tick == dtfCos$ticker), nrow(dtfCos)))
  hit_api(tick, outputsize="full", saveJSON=TRUE)
  Sys.sleep(2.1)
}


# Create single data.frame with stock tick values as columns
lst_dtfs <- lapply(dtfCos$ticker, parse_stock_data, 5) # 5 represents adjusted closing price
dtf <- rbindlist(lapply(lst_dtfs, function(DT) {DT[, "tick"] <- attr(DT, "ticker"); return(DT)}))

# Order stocks by length of history
dtfMinDate <- dtf[, min(Date), by=tick][order(V1)]
dtf <- dcast(dtf, Date ~ tick, value.var="adjusted close")
dtf <- dtf[, .SD, .SDcols=c(dtfMinDate$tick, "Date")]
fwrite(dtf, file="data/SandP_tick_history.csv")

# Missing map
Amelia::missmap(dtf, ylab="Feb 2021 - Nov 1999", xlab="Tick")


files <- dir("data/raw_json/", pattern="*.json")

# Check for failed data pulls and re-run data pull from AlphaAdvantage API
for (f_ in files) {
  dir_ <- paste0("data/raw_json/", f_)
  fs <- file.size(dir_)
  if (fs < 300) {
    print(sprintf("Re-running %s", dir_))
    tick <- unlist(strsplit(f_, "\\."))[1]
    hit_api(tick, outputsize="full", saveJSON=TRUE)
    fs <- file.size(dir_)
    if (fs > 300) {
      print(sprintf("Cured %s. New file size %s", dir_, fs))
    } else {
      print(sprintf("Not Cured %s", dir_))
    }
    Sys.sleep(2.01)
  }
}


