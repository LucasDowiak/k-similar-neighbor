library(jsonlite)
library(data.table)
library(rugarch)
setwd("~/Git/k-similar-neighbor")
source("R/auto_marginal.R")



# ------------------------------------------------------------------------
# code snippet to run through all S&P 500 stock tickers

SPEC_PATH <- "data/marginal_specifications_20190201_20200214.json"
st_date <- "2019-02-01"
ed_date <- "2020-02-14"

# Pull the tickers from the json files
tickers <- sapply(strsplit(dir("data/raw_json", pattern="json$"), "\\."), `[[`, 1)

# Auto_fit each ticker, saving after each tick
SPECS <- vector("list", length(tickers))
names(SPECS) <- tickers
for (tick in tickers) {
  print(sprintf("TICK %s at %s", tick, Sys.time()))
  dtfU <- try(data.table(parse_json(tick)))
  if (nrow(dtfU) == 0) {
    next
  } else {
    u <- dtfU[Date >= st_date & Date <= ed_date, diff(log(close))]
    out <- try(auto_fit(u))
    if (inherits(out, "try-error"))
      out <- as.character(out)
    SPECS[[tick]] <- out
    write_json(SPECS, path=SPEC_PATH)
  }
}


# Create single data.frame with stock tick values as columns
lst_dtfs <- lapply(tickers, parse_json)
dtf <- rbindlist(lapply(lst_dtfs, function(DT) {DT[, "tick"] <- attr(DT, "ticker"); return(DT)}))
# Order stocks by length of history
dtfMinDate <- dtf[, min(Date), by=tick][order(V1)]
dtf <- dcast(dtf, Date ~ tick, value.var="close")
dtf <- dtf[, .SD, .SDcols=c(dtfMinDate$tick, "Date")]
fwrite(dtf, file="data/SandP_tick_history.csv")
# Missing map
# Amelia::missmap(dtf)


tst <- specs_to_resid_matrix(SPEC_PATH)

# For the ticks that failed to produce a proper fit, check if forcing the ar component to 7 does the trick
# ------------------------------------------------------------------------
SPEC_PATH <- "data/marginal_specifications_20200201_20210215.json"
SPECIFICATIONS <- read_json(SPEC_PATH)
SPECS <- SPECIFICATIONS[-length(SPECIFICATIONS)]

ticks <- good_vs_bad_symbols(SPEC_PATH)

out_spec <- run_seasonal_modifications(badticks, SPECS)

write_json(out_spec, path="data/marginal_specifications_seasonal.json")

