setwd("~/Git/dtw-in-finance")
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


files_ <- list.files(sprintf("~/Git/dtw-in-finance/data/model_objects/%s/", year_label),
                     full.names = T)
ftick <- list.files(sprintf("~/Git/dtw-in-finance/data/model_objects/%s/", year_label))
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
