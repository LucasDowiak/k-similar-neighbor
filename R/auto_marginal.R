library(jsonlite)
library(data.table)
library(rugarch)
source("R/normality_tests.R")


LABEL_PAT_4_TICK <- "(%s_)([A-Z]+\\.?[A-Z]*)(\\.rds)"

# Check to see which stocks had a successful auto_fit run
check_spec <- function(lst)
{
  norun <- length(lst) <= 1
  if (norun)
    return(NA)
  else
    return(isTRUE(unlist(lst$pass)))
}


# Given a path containing lists of model specifications, read in and return list
# containing 1) good and 2) bad vector of stock symbol
good_vs_bad_symbols <- function(spec.path)
{
  specs <- interpret_spec_list(spec.path)
  nms <- names(specs)
  goodspecs <- sapply(specs, check_spec)
  return(list(pass_ticks=nms[goodspecs & !is.na(goodspecs)],
              fail_ticks=nms[!goodspecs & !is.na(goodspecs)],
              not_run_ticks=nms[is.na(goodspecs)]))
}


# Given a ticker value, this function reads in the corresponding
# json file and pulls out the values you want (e.g. "open", "close", "high")
#
#   Input:
#       tick   - stock tick symbol
#       values - integer values between 1 and 8 are accepted
#                   1. open
#                   2. high
#                   3. low
#                   4. close
#                   5. adjusted close
#                   6. dividend amount
#                   7. volume
#                   8. split coefficient
#
#   Output:
#       data.table
#
parse_stock_data <- function(tick, values=1:8)
{
  values <- unique(values)
  stopifnot(all(values %in% 1:8))
  
  file_ <- paste0(tick, ".json")
  full_path <- paste0("data/raw_json/", file_)
  if (file.exists(full_path)) {
    fs <- file.size(full_path)
    if (fs < 300) {
      print(sprintf("Tick %s has no data. Returning NULL data.frame", tick))
      return(data.frame())
    }
    lst <- read_json(full_path)
  } else {
    lst <- fromJSON(tick)
  }
  
  prices <- do.call(rbind, lst[[2]])
  dates <- row.names(prices)
  cols <- vapply(strsplit(colnames(prices), "\\. "), `[[`, 2, FUN.VALUE=character(1))
  dtf <- as.data.table(matrix(as.numeric(prices),
                              ncol=length(cols),
                              dimnames=list(NULL, cols)))
  dtf <- dtf[, values, with=FALSE]
  dtf[, Date := as.Date(dates)]
  attr(dtf, "ticker") <- lst[[1]]$`2. Symbol`
  attr(dtf, "last_refreshed") <- lst[[1]]$`3. Last Refreshed`
  attr(dtf, "date_range") <- range(dtf$Date)
  return(dtf[order(Date)])
}


# Given a series x, automatically try various conditional mean
# and condition variance models. check whether any pass all the 
# verification tests (pass). Choose the best model and store the
# model specifications
#
# auto_fit(u)
#
auto_fit <- function(x, max_arch=2, max_garch=2, ignore_nyblom=FALSE)
{
  "ar ma sar sma period i si"
  aa <- forecast::auto.arima(x, max.p=10, allowmean=TRUE, allowdrift=FALSE)
  spec <- aa$arma
  arma <- spec[1:2]
  if (spec[3] > 0)
    arma[1] <- 7 * arma[3]
  if (spec[4] > 0)
    arma[2] <- 7 * arma[4]
  
  arch <- 1:max_arch
  garch <- 1:max_garch
  garchmodels <- c("sGARCH", "gjrGARCH", "csGARCH") # "apARCH")
  distributions <- c("std", "sstd")
  outerpaste <- function(a, b) as.vector(outer(a, b, paste, sep="_"))
  rnames <- Reduce(outerpaste, list(garchmodels, arch, garch), init=distributions)
  lik <- matrix(NA_real_,
                nrow=prod(length(arch), length(garch),
                          length(garchmodels), length(distributions)),
                ncol=3,
                dimnames=list(rnames, c("bic", "likelihood", "pass_marginal")))
  all_pass_test <- FALSE
  lowest_bic_model <- NULL
  llh_bic <- Inf
  lowest_pass_model <- NULL
  llh_pass <- Inf
  for (mod in garchmodels) {
    for (d in distributions) {
      for (g in garch) {
        for (a in arch) {
          nn <- paste(d, mod, g, a, sep="_")
          print(nn)
          spec_mod <- ugarchspec(
            variance.model = list(model=mod, garchOrder=c(a,g)),
            mean.model = list(armaOrder=arma, include.mean=TRUE),
            distribution.model = d
          )
          
          fit_ <- try(ugarchfit(spec=spec_mod, data=x))

          # If the fitting process failed
          if (inherits(fit_, "try-error") || fit_@fit$convergence != 0) {
            lik[nn,] <- c(NA_real_, NA_real_, FALSE)

          # If the fitting process succeeded
          } else {
            marg_tests <- try(marginal_tests(fit_, PRINT=FALSE, PLOT=FALSE))
            marg_check <- try(verify_marginal_test(marg_tests, alpha="0.05", ignore_nyblom=ignore_nyblom))
            l <- try(infocriteria(fit_)['Bayes', 1])
            LLH <- try(fit_@fit$LLH)
            
            if (inherits(l, "try-error") || inherits(marg_check, "try-error") || inherits(marg_tests, "try-error")) {
              lik[nn,] <- c(NA_real_, NA_real_, FALSE)
              
            } else {
              all_pass_test <- marg_check[, all(pass_test, na.rm=TRUE)]
              lik[nn,] <- c(l, LLH, all_pass_test)
              
              # Check if we have a new lower bic model
              if (l < llh_bic) {
                lowest_bic_model <- fit_
                llh_bic <- l
              }
              # Check if we have a new lower bic model that also passes marginal checks
              if (l < llh_pass & all_pass_test) {
                lowest_pass_model <- fit_
                llh_pass <- l
              }
            }
          }
          # Stop looking if we have a model specification that passes all tests
          if (all_pass_test) {
            break
          }
        }
        if (all_pass_test) {
          break
        }
      }
      if (all_pass_test) {
        break
      }
    }
    if (all_pass_test) {
      break
    }
  }
  # Default choice that minimizes BIC across specifications
  choice <- unlist(strsplit(names(which.min(lik[, "bic"])), "_"))

  # Prefer to consider only specifications that pass
  pass <- any(lik[, "pass_marginal"] == 1)
  if (pass) {
    num_pass_specs <- sum(lik[, "pass_marginal"], na.rm=TRUE)
    model <- lowest_pass_model
    
    if (num_pass_specs == 1) {
      choice <- unlist(strsplit(names(which(lik[, "pass_marginal"] == 1)), "_"))
      
    } else if (num_pass_specs > 1) {
      lik2 <- lik[which(lik[, "pass_marginal"] > 0), ]
      choice <- unlist(strsplit(names(which.min(lik2[, "bic"])), "_"))
      
    }
  } else {
    model <- lowest_bic_model
  }
  return(list(ar=arma[1], ma=arma[2], arch=choice[4], garch=choice[3],
              garchmod=choice[2], distr=choice[1], pass=pass, lik=lik,
              model=model))
}


# For a tick that has already been run through the auto_fit function but
# failed to find a fit that passes all the marginal tests,
# try setting the the AR component to 7

# batch process produce_rGarchfit_s over all failed marginal models
# >>> Not a actually working
run_seasonal_modification <- function(badticks, specification, st_date, ed_date)
{
  bool <- c()
  names(bool) <- badticks
  for (tt in seq_along(badticks)) {
    mod <- try(produce_rGARCHfit_s(badticks[tt], st_date, ed_date))
    if (inherits(mod, "try-error"))
      verify <- FALSE
    else
      verify <-  try(verify_marginal_test(marginal_tests(mod)))
    if (inherits(verify, "try-error"))
      verify <- FALSE
    bool <- c(bool, verify)
  }
  
  for (tick in names(bool[bool])) {
    specification[[tick]][["ar"]] <- list(7L)
    specification[[tick]][["pass"]] <- list(TRUE)
  }
  return(specification)
}


# Takes a stock ticker `tick`
# Grabs the raw stock values:
#   a) Transforms into log returns
#   b) Returns subset of dates from `start_date` to `end_date`
# Fits a model from it's specification in the list object `specs` 
# Returns:
#   a) Fitted ugarch model
#   b) Data.table with log returns
produce_rGARCHfit_from_spec <- function(tick, specs, start_date, end_date)
{
  
  file_ <- paste0(toupper(tick), ".json")
  
  # Grab data from raw json files
  # Subset dates
  # calculate log returns
  if (file.exists(paste0("data/raw_json/", file_))) {
    pj <- parse_stock_data(tick)
    dtfU <- data.table(pj)
    # No data for this tick
    if (nrow(dtfU) == 0) {
      print(sprintf("No data for ticker [%s] data/raw_json/", tick))
      return(list(fit=NULL, dates=NULL))
    }
    dtfU <- dtfU[, x := c(NA_real_, diff(log(close)))][Date >= start_date & Date <= end_date]
    x <- dtfU[, x]
    if (nrow(dtfU) == 0) {
      print(sprintf("No data for ticker [%s] from %s to %s", tick, start_date, end_date))
      return(list(fit=NULL, dates=NULL))
    }
  } else {
    print(sprintf("JSON file doesn't exist for ticker [%s] in directory data/raw_json/", tick))
    return(list(fit=NULL, dates=NULL))
  }
  
  # Extract model specification from specs
  spec <- specs[[toupper(tick)]]
  if (is(spec, "NULL")) {
    stop("Tick does not exist in `specs`")
  }
  
  if (!unlist(spec[["pass"]])) {
    warning(sprintf("Tick %s did not pass all the marginal tests", tick))
  }
  # Fit garch model and return along with dates of analysis
  spec_mod <- ugarchspec(
    variance.model = list(model = unlist(spec$garchmod),
                          garchOrder = c(as.integer(unlist(spec$arch)), as.integer(unlist(spec$garch)))),
    mean.model = list(armaOrder = c(unlist(spec$ar), unlist(spec$ma)), include.mean = TRUE),
    distribution.model = unlist(spec$distr)
  )
  fit_ <- try(ugarchfit(spec=spec_mod, data=x))
  return(list(fit=fit_, dates=dtfU[, Date]))
}


# Function that performs univariate ARIMA-GARCH modeling and dynamic time warping for a efine function that takes a spec_path and creates a matrix of each
#  Formally 'specs_to_resid_matrix'
#
#   spec_path - file path to specifications
arima_garch_optimization <- function(start_date,
                                     end_date,
                                     tickers,
                                     label="default",
                                     write_output=FALSE,
                                     save_models=FALSE,
                                     rerun=FALSE,
                                     ...)
{
  model_dir <- sprintf("~/Git/dtw-in-finance/data/model_objects/%s/", label)
  if (rerun) {
    if (!dir.exists(model_dir))
      stop(sprintf("The directory does not exists: %s", model_dir))
    pat <- sprintf(LABEL_PAT_4_TICK, label)
    existing_ticks <- gsub(pat, "\\2", list.files(model_dir, pattern=label))
    tickers <- setdiff(tickers, existing_ticks)
  }
  resid_list_U <- resid_list_Z <- resid_list_I <- model_struct <- model_params <- vector("list", length(tickers))
  names(resid_list_U) <- tickers
  names(resid_list_Z) <- tickers
  names(resid_list_I) <- tickers
  names(model_params) <- tickers
  names(model_struct) <- tickers
  
  for (tick in tickers) {
    ii <- which(tick == tickers)
    print(sprintf("Tick: %s - %d out of %d", tick, ii, length(tickers)))
    dtfU <- parse_stock_data(tick)
    setnames(dtfU, "adjusted close", "adjclose")
    dtfU <- dtfU[, logR := c(NA_real_, diff(log(adjclose)))][Date >= start_date & Date <= end_date]
    obj <- auto_fit(dtfU$logR, ...)

    # If successfully fit from the specification
    if (inherits(obj$model, "uGARCHfit")) {

      nms <- c("Date", tick)
      
      # Standardized residuals
      tmpZ <- data.table(Date=dtfU$Date, residuals(obj$model, standardize=TRUE))
      setnames(tmpZ, names(tmpZ), nms)
      
      # Uniform marginals by way of the probability integral transform
      tmpU <- data.table(Date=dtfU$Date, rugarch::pit(obj$model))
      setnames(tmpU, names(tmpU), nms)
      
      # Standardize time series - Divide all prices by the price on the first day of the series
      tmpI <- data.table(Date=dtfU$Date, dtfU$adjclose / dtfU$adjclose[1])
      setnames(tmpI, names(tmpI), nms)
      
      # Model Summary
      tmpP <- as.data.table(c(as.list(obj$model@model$modelinc), obj$model@model$modeldesc))
      tmpP[, tick := tick]
      
      # ARIMA-GARCH Coefficients
      tmpM <- as.data.table(t(coef(obj$model)))
      tmpM[, tick := tick]
      
    } else {
      # No ARIMA-GARCH was fitted
      
      tmpU <- data.table(Date=dtfU$Date, NA_real_)
      setnames(tmpU, "V1", tick)
      
      tmpZ <- data.table(Date=dtfU$Date, NA_real_)
      setnames(tmpZ, "V1", tick)
      
      tmpI <- data.table(Date=dtfU$Date, NA_real_)
      setnames(tmpI, "V1", tick)
      
      tmpP <- null_model_spec()[['modelsum']]
      
      tmpM <- null_model_spec()[['modelcoef']]

    }
    resid_list_U[[tick]] <- tmpU
    resid_list_Z[[tick]] <- tmpZ
    resid_list_I[[tick]] <- tmpI
    model_struct[[tick]] <- tmpP
    model_params[[tick]] <- tmpM
    
    if (save_models) {
      if (!dir.exists(model_dir)) {
        dir.create(model_dir)
      }
      model_file <- paste(c(model_dir, sprintf("%s_%s.rds", label, tick)), collapse = "/")
      saveRDS(obj, file=model_file)
    }
  }
  outZ <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_Z)
  outU <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_U)
  outI <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_I)
  outP <- rbindlist(model_struct, use.names=TRUE, fill=TRUE)
  outM <- rbindlist(model_params, use.names=TRUE, fill=TRUE)
  
  if (write_output) {
    base_path <- "~/Git/dtw-in-finance/data/label_analysis/"
    create_name <- function(x) sprintf("%s%s_%s.csv", base_path, label, x)
    fwrite(outU, file=create_name("uni_marg"))
    fwrite(outZ, file=create_name("std_resids"))
    fwrite(outI, file=create_name("std_price"))
    fwrite(outP, file=create_name("model_summary"))
    fwrite(outM, file=create_name("model_coef"))
  }
  return(list(uni_marg=outU, std_resids=outZ, std_price=outI, model_summary=outP, model_coef=outM))
}


# Re-run label analysis. Use existing model objects for the label to re-calculate
# standard residuals, standard prices, and model summaries
# 
# tst <- run_label_analysis(label="DA_2015_2016", start_date="2015-01-01", end_date="2016-12-31", write=TRUE)
#
run_label_analysis <- function(label, start_date, end_date, write_output=TRUE)
{
  
  base_dir <- "~/Git/dtw-in-finance/data/"
  dtfSP <- fread(paste0(base_dir, "SandP_tick_history.csv"))
  # Buffer the start date by one trading day so log returns can be taken for
  # the first day
  si <- dtfSP[Date >= start_date, which(dtfSP$Date == min(Date))]
  ei <- dtfSP[Date <= end_date, which(dtfSP$Date == max(Date))]
  dtfU <- dtfSP[(si-1):ei]
  
  mobjs <- list.files(paste0(base_dir, sprintf("model_objects/%s/", label)))
  pat <- sprintf(LABEL_PAT_4_TICK, label)
  
  resid_list_U <- resid_list_Z <- resid_list_I <- model_struct <- model_params <- vector("list", length(mobjs))

  for (i in seq_along(mobjs)) {
    tick <- gsub(pat, "\\2", mobjs[i])
    print(sprintf("Start ticker [%d / %d]: %s", i, length(mobjs), tick))
    obj <- readRDS(paste0(base_dir, sprintf("model_objects/%s/%s", label, mobjs[i])))
    
    if (inherits(obj$model, "uGARCHfit")) {
      nms <- c("Date", tick)
      
      # Standardized residuals
      tmpZ <- data.table(Date=dtfU[-1]$Date, residuals(obj$model, standardize=TRUE))
      setnames(tmpZ, names(tmpZ), nms)
      
      # Uniform marginals by way of the probability integral transform
      tmpU <- data.table(Date=dtfU[-1]$Date, rugarch::pit(obj$model))
      setnames(tmpU, names(tmpU), nms)
      
      # Standardize time series with cumulative return transformation
      logR <- diff(log(dtfU[[tick]]))
      compretn <- cumprod(1 + logR) - 1
      tmpI <- data.table(Date=dtfU[-1]$Date, cumprod(1 + logR) - 1)
      setnames(tmpI, names(tmpI), nms)
      
      # Model Summary
      tmpP <- as.data.table(c(as.list(obj$model@model$modelinc), obj$model@model$modeldesc))
      tmpP[, tick := tick]
      
      # ARIMA-GARCH Coefficients
      tmpM <- as.data.table(t(coef(obj$model)))
      tmpM[, tick := tick]
    } else {
      warning("File %s is not a uGARCHfit object", mobjs[i])
      tmpZ <- tmpU <- tmpI <- tmpP <- tmpM <- data.table()
    }
    resid_list_U[[i]] <- tmpU
    resid_list_Z[[i]] <- tmpZ
    resid_list_I[[i]] <- tmpI
    model_struct[[i]] <- tmpP
    model_params[[i]] <- tmpM
  }
  outZ <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_Z)
  outU <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_U)
  outI <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_I)
  outP <- rbindlist(model_struct, use.names=TRUE, fill=TRUE)
  outM <- rbindlist(model_params, use.names=TRUE, fill=TRUE)
  for (dtbl in list(outZ, outU, outI, outP, outM)) {
    dtbl[, label := label]
  }
  if (write_output) {
    base_path <- "~/Git/dtw-in-finance/data/label_analysis/"
    create_name <- function(x) sprintf("%s%s_%s.csv", base_path, label, x)

    fwrite(outU, file=create_name("uni_marg"))
    fwrite(outZ, file=create_name("std_resids"))
    fwrite(outI, file=create_name("std_price"))
    fwrite(outP, file=create_name("model_summary"))
    fwrite(outM, file=create_name("model_coef"))
  }
  
  return(list(outZ=outZ, outU=outU, outI=outI, outP=outP, outM=outM))
}


# Define NULL model spec summary and NULL model coefficients
null_model_spec <- function()
{
  modelsum <-  data.table(
     mu=          NA_integer_,
     ar=          NA_integer_,
     arfima=      NA_integer_,
     archm=       NA_integer_,
     mxreg=       NA_integer_,
     omega=       NA_integer_,
     alpha=       NA_integer_,
     beta=        NA_integer_,
     gamma=       NA_integer_,
     eta1=        NA_integer_,
     eta2=        NA_integer_,
     delta=       NA_integer_,
     lambda=      NA_integer_,
     vxreg=       NA_integer_,
     skew=        NA_integer_,
     shape=       NA_integer_,
     ghlambda=    NA_integer_,
     xi=          NA_integer_,
     aux=         NA_integer_,
     aux=         NA_integer_,
     aux=         NA_integer_,
     distribution=NA_character_,
     distno=      NA_integer_,
     vmodel=      NA_character_
  )
  modelcoef <- data.table(mu=NA_real_, omega=NA_real_, shape=NA_real_)
  return(list(modelsummary=modelsum, modelcoef=modelcoef))
}


# Takes a specification file and produces a table with ARIMA-GARCH specifications
#
#   spec_list - either a list of ARIMA-GARCH specifications or a file path to the json object
#
#   dtfSpec <- summarize_specs("data/marginal_specifications_20200323_20210219.json")
#
summarize_model_spec <- function(spec_list)
{
  spec_list <- interpret_spec_list(spec_list)

  # Extract model specifications
  grab_specifications <- function(x) {
    if (length(x) < 2) {
      out <- data.table(ar=NA, ma=NA, arch=NA, garch=NA, garchmod=NA, distr=NA, pass=NA)
    } else {
      out <- data.table(ar=x$ar, ma=x$ma, arch=x$arch, garch=x$garch,
                        garchmod=x$garchmod, distr=x$distr, pass=x$pass)
      out <- out[, lapply(.SD, unlist)]
    }
    return(out)
  }
  dtfSpec <- rbindlist(lapply(spec_list, grab_specifications))
  dtfSpec[, ticker := names(spec_list)]
  return(dtfSpec)
}

