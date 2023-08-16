setwd("~/Git/k-similar-neighbor")
library(jsonlite)
library(data.table)
library(rugarch)
source("R/auto_marginal.R")



# ------------------------------------------------------------------------
# code snippet to run through all S&P 500 stock tickers for a given year

run_arima_garch_by_years <- function(years, rerun=TRUE)
{
  if (length(years) > 1) {
    year_range <- paste0(years, collapse="_")
  } else {
    year_range <- as.character(years)
  }
  tickers <- sapply(strsplit(dir("data/raw_json", pattern="json$"), "\\."), `[[`, 1)
  spec_path <- sprintf("data/marginal_specifications_%s.json", year_range)
  
  if (rerun) {
    spec_list <- interpret_spec_list(spec_path)
    tmp <- good_vs_bad_symbols(spec_path)
    tickers <- c(tmp$fail_ticks, tmp$not_run_ticks)
    
  } else {
    spec_list <- vector("list", length(tickers))
    names(spec_list) <- tickers
  }
  
  for (tick in tickers) {
    print(sprintf("TICK %s at %s", tick, Sys.time()))
    dtfU <- try(data.table(parse_json(tick)))
    if (any(nrow(dtfU) == 0, length(dtfU) == 0, is(dtfU, "NULL"))) {
      next
    } else {
      u <- dtfU[year(Date) %in% years, diff(log(close))]
      out <- try(auto_fit(u))
      if (inherits(out, "try-error"))
        out <- as.character(out)
      spec_list[[tick]] <- out
      write_json(spec_list, path=spec_path)
    }
  }
}

years <- c(2007, 2008)

run_arima_garch_by_years(c(2007, 2008), rerun=FALSE)

for (yr in years) {
  run_arima_garch_by_year(yr, rerun = FALSE)
}

lst_files <- list.files('data', pattern = "marginal_specifications_[0-9]{4}.json",
                        full.names = TRUE)
gvb <- lapply(lst_files, good_vs_bad_symbols)
o <- lapply(gvb, `[[`, 'fail_ticks')
names(o) <- lst_files

all_fails <- Reduce(union, o)
univ_intersect <- Reduce(intersect, lapply(gvb, `[[`, 'pass_ticks'))
gvb_sums <- t(sapply(gvb, function(x) c(pass_ticks=length(x$pass_ticks), fail_ticks=length(x$fail_ticks), not_run=length(x$not_run_ticks))))


SPEC_PATH <- "data/marginal_specifications_20200323_20210219.json"
st_date <- "2020-03-23"
ed_date <- "2021-02-19"

tmp <- good_vs_bad_symbols(SPECS)
tickers <- c(tmp$fail_ticks, tmp$not_run_ticks)
# Pull the tickers from the json files
tickers <- sapply(strsplit(dir("data/raw_json", pattern="json$"), "\\."), `[[`, 1)


# Create single data.frame with stock tick values as columns
lst_dtfs <- lapply(tickers, parse_json)
dtf <- rbindlist(lapply(lst_dtfs, function(DT) {DT[, "tick"] <- attr(DT, "ticker"); return(DT)}))
# Order stocks by length of history
dtfMinDate <- dtf[, min(Date), by=tick][order(V1)]
dtf <- dcast(dtf, Date ~ tick, value.var="close")
dtf <- dtf[, .SD, .SDcols=c(dtfMinDate$tick, "Date")]
fwrite(dtf, file="data/SandP_tick_history.csv")
# Missing map
Amelia::missmap(dtf, ylab="Feb 2021 - Nov 1999",
                xlab="Tick")


tst <- specs_to_resid_matrix(SPEC_PATH)

# For the ticks that failed to produce a proper fit, check if forcing the ar component to 7 does the trick
# ------------------------------------------------------------------------
SPEC_PATH <- "data/marginal_specifications_20190201_20200214.json"
SPECIFICATIONS <- read_json(SPEC_PATH)
st_date <- "2019-02-01"
ed_date <- "2020-02-14"

ticks <- good_vs_bad_symbols(SPEC_PATH)
tickers <- ticks$bad_ticks
# Auto_fit each ticker, saving after each tick
SPECS <- vector("list", length(tickers))
names(SPECS) <- tickers
for (tick in tickers) {
  print(sprintf("TICK %s at %s", tick, Sys.time()))
  dtfU <- try(data.table(parse_json(tick)))
  if (any(nrow(dtfU) == 0, length(dtfU) == 0, is(dtfU, "NULL"))) {
    next
  } else {
    u <- dtfU[Date >= st_date & Date <= ed_date, diff(log(close))]
    out <- try(auto_fit(u))
    if (inherits(out, "try-error"))
      out <- as.character(out)
    SPECS[[tick]] <- out
    #write_json(SPECS, path=SPEC_PATH)
  }
}

xx <- sapply(SPECS, length)
for (tick in names(xx[xx==1])) {
  dtfU <- data.table(parse_json(tick))
  print(dtfU[, range(Date)])
}


produce_rGARCHfit_from_spec("ABC", SPECIFICATIONS)

out_spec <- run_seasonal_modifications(badticks, SPECS)

write_json(out_spec, path="data/marginal_specifications_seasonal.json")

tmp1 <- fread("data/std_resids_20190201_20200214.csv")
tmp2 <- fread("data/std_resids_20200323_20210219.csv")

