setwd("~/Git/k-similar-neighbor")
library(jsonlite)
library(data.table)
library(rugarch)
source("R/auto_marginal.R")




# ------------------------------------------------------------------------
# code snippet to run through all S&P 500 stock tickers for a given year

# Function to run ARIMA-GARCH models by year
#
#   years - list of years to include as time series domain
#
#   rerun - boolean indicating if the process should start from the ticker where
#           the previous run left off
#
debugonce(arima_garch_optimization)
aa <- arima_garch_optimization(start_date="2020-01-01",
                               end_date="2020-12-31",
                               label=year_label,
                               tickers=tickers,
                               write_output=TRUE,
                               save_models=TRUE,
                               rerun=FALSE,
                               ignore_nyblom=TRUE)

# Run model optimization in batches ---------------------------------------------
dtfSP <- fread("data/SandP_tick_history.csv")

year_label <- "2000"
tickers <- names(Filter(function(x) sum(is.na(x)) == 0, dtfSP[year(Date) == as.integer(year_label)]))
tickers <- setdiff(tickers, "Date")
aa <- arima_garch_optimization(start_date=sprintf("%s-01-01", year_label),
                               end_date=sprintf("%s-12-31", year_label),
                               label=year_label,
                               tickers=tickers,
                               write_output=TRUE,
                               save_models=TRUE,
                               rerun=FALSE,
                               ignore_nyblom=TRUE)


files_ <- list.files(sprintf("~/Git/k-similar-neighbor/data/model_objects/%s/", year_label),
                     full.names = T)
ftick <- list.files(sprintf("~/Git/k-similar-neighbor/data/model_objects/%s/", year_label))
ftick <- sapply(strsplit(ftick, "\\.|_"), `[[`, 2)
models <- lapply(files_, readRDS); names(models) <- ftick
results <- lapply(models, function(x) verify_marginal_test(marginal_tests(x$model), "0.05", ignore_nyblom=TRUE))
for (ii in names(results)) {
  results[[ii]][, tick := ii]
}
results <- rbindlist(results)
failed_ticks <- results[, all(pass_test), by=tick][!as.logical(V1), tick]
failed_ticks


diagnostic <- vector("list")
for (yr in as.character(2000:2022)) {
  model_diag <- run_label_analysis(yr,
                                   start_date=sprintf("%s-01-01", yr),
                                   end_date=sprintf("%s-12-31", yr),
                                   write_output = TRUE)
  diagnostic[[yr]] <- model_diag
}


# QA on the model selection process --------------------------------------------
# DEPRECATED - Uses old system where the model specification, and not a ugarch obj,
#              is serialized in json
lst_files <- list.files('data', pattern = "marginal_specifications_[0-9]{4}.json",
                        full.names = TRUE)
gvb <- lapply(lst_files, good_vs_bad_symbols)
o <- lapply(gvb, `[[`, 'fail_ticks')
names(o) <- lst_files

all_fails <- Reduce(union, o)
univ_intersect <- Reduce(intersect, lapply(gvb, `[[`, 'pass_ticks'))
gvb_sums <- t(sapply(gvb,
                     function(x) c(pass_ticks=length(x$pass_ticks),
                                   fail_ticks=length(x$fail_ticks),
                                   not_run=length(x$not_run_ticks))))


SPEC_PATH <- "data/marginal_specifications_20200323_20210219.json"
st_date <- "2020-03-23"
ed_date <- "2021-02-19"

tmp <- good_vs_bad_symbols(SPECS)
tickers <- c(tmp$fail_ticks, tmp$not_run_ticks)
# Pull the tickers from the json files
tickers <- sapply(strsplit(dir("data/raw_json", pattern="json$"), "\\."), `[[`, 1)
