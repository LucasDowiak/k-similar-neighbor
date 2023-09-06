library(jsonlite)
library(data.table)
library(rugarch)
source("R/normality_tests.R")

# source("R/auto_marginal.R")
# -------------------------------------------------
# dtfTmp <- parse_stock_data("A")
# tfTmp[, x := c(NA_real_, diff(log(close)))]
# x <- dtfTmp[Date >= "2019-01-01" & Date <= "2019-12-31", x]
# aa <- forecast::auto.arima(x, max.p=10, allowmean=TRUE, allowdrift=FALSE)
# afit <- auto_fit(x)
# spec_mod <- ugarchspec(
#   variance.model = list(model="gjrGARCH", garchOrder=c(2,2)),
#   mean.model = list(armaOrder=aa$arma[1:2], include.mean=TRUE),
#   distribution.model = "std"
# )

# fit_ <- try(ugarchfit(spec=spec_mod, data=x))



# Model specification lists can be given either as:
# i) a file path to the specification list stored as a json object
# ii) an R list object
# This function can handle both inputs and adds some sanity checks
#
#   spec_list - either a path to a specification file or the R list object itself
#
interpret_spec_list <- function(spec_list)
{
  if (is(spec_list, "character")) {
    if (file.exists(spec_list)) {
      spec_list <- read_json(spec_list)
    } else {
      stop(sprintf("JSON file doesn't exist for specification file located at [%s] ", spec_list))
    }
  }
  if (!is(spec_list, "list")) {
    stop("spec_list input is not a list.")
  }
  return(spec_list)
}


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
# dtfU <- data.table(parse_json("GOOG"))
# u <- dtfU[Date >= st_date & Date <= ed_date, diff(log(close))]
#
parse_stock_data <- function(tick, values="close")
{
  values <- unique(values)
  value_opts <- c("open", "high", "low", "close", "volume")
  stopifnot(all(values %in% value_opts))
  idx <- which(value_opts %in% values)
  
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
  prices <- t(vapply(lst[[2]], function(x) unlist(x[idx]),
                     FUN.VALUE=character(length(values))))
  
  if (length(values) == 1) {
    dates <- colnames(prices)
    cols <- values
  } else {
    cols <- vapply(strsplit(colnames(prices), " "), `[[`, 2, FUN.VALUE=character(1))
    dates <- row.names(prices)
  }
  dtf <- as.data.frame(matrix(as.numeric(prices), ncol=length(cols),
                              dimnames=list(NULL, cols)))
  dtf['Date'] <- as.Date(dates)
  attr(dtf, "ticker") <- lst[[1]]$`2. Symbol`
  attr(dtf, "last_refreshed") <- lst[[1]]$`3. Last Refreshed`
  attr(dtf, "date_range") <- range(dtf$Date)
  return(data.table(dtf[order(dtf$Date), ]))
}


# Given a series x, automatically try various conditional mean
# and condition variance models. check whether any pass all the 
# verification tests (pass). Choose the best model and store the
# model specifications
#
# auto_fit(u)
#
auto_fit <- function(x, max_arch=2, max_garch=2)
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
  garchmodels <- c("sGARCH", "gjrGARCH", "csGARCH")
  distributions <- c("norm", "std", "sstd")
  outerpaste <- function(a, b) as.vector(outer(a, b, paste, sep="_"))
  rnames <- Reduce(outerpaste, list(garchmodels, arch, garch), init=distributions)
  lik <- matrix(NA_real_,
                nrow=prod(length(arch), length(garch),
                          length(garchmodels), length(distributions)),
                ncol=3,
                dimnames=list(rnames, c("bic", "likelihood", "pass_marginal")))
  lowest_bic_model <- NULL
  llh_bic <- Inf
  lowest_pass_model <- NULL
  llh_pass <- Inf
  for (a in arch) {
    for (g in garch) {
      for (mod in garchmodels) {
        for (d in distributions) {
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
            marg_tests <- marginal_tests(fit_, print=FALSE, plot=FALSE)
            marg_check <- try(verify_marginal_test(marg_tests, 0.05))
            l <- try(infocriteria(fit_)['Bayes', 1])
            LLH <- try(fit_@fit$LLH)
            
            if (inherits(l, "try-error") || inherits(marg_check, "try-error")) {
              lik[nn,] <- c(NA_real_, NA_real_, FALSE)
            } else {
              lik[nn,] <- c(l, LLH, marg_check)
              # Check if we have a new lower bic model
              if (l < llh_bic) {
                lowest_bic_model <- fit_
                llh_bic <- l
              }
              # Check if we have a new lower bic model that also passes marginal checks
              if (l < llh_pass & marg_check) {
                lowest_pass_model <- fit_
                llh_pass <- l
              }
            }
          }
        }
      }
    }
  }
  # Default choice that minimizes BIC across specifications
  choice <- unlist(strsplit(names(which.min(lik[, "bic"])), "_"))

  # Prefer to consider only specifications that pass
  pass <- any(lik[, "pass_marginal"] == 1)
  if (pass) {
    num_pass_specs <- sum(lik[, "pass_marginal"])
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
                                     ...)
{
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
    dtfU <- dtfU[, logR := c(NA_real_, diff(log(close)))][Date >= start_date & Date <= end_date]
    obj <- auto_fit(dtfU$logR, ...)

    # If successfully fit from the specification
    if (inherits(obj$model, "uGARCHfit")) {

      nms <- c("Date", tick)
      
      # Standardized residuals
      tmpZ <- data.table(Date=dtfU$Date, residuals(obj$model, standardize=TRUE))
      setnames(tmpZ, names(tmpZ), nms)
      
      # Uniform marginals by way of the probability integral transform
      tmpU <- data.table(Date=dtfU$Date, PIT(obj$model))
      setnames(tmpU, names(tmpU), nms)
      
      # Standardize time series - Divide all prices by the price on the first day of the series
      tmpI <- data.table(Date=dtfU$Date, dtfU$close / dtfU$close[1])
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
      model_file <- sprintf("~/Git/k-similar-neighbor/data/model_objects/%s_%s.rds", label, tick)
      saveRDS(obj, file=model_file)
    }
  }
  outZ <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_Z)
  outU <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_U)
  outI <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list_I)
  outP <- rbindlist(model_struct, use.names=TRUE, fill=TRUE)
  outM <- rbindlist(model_params, use.names=TRUE, fill=TRUE)
  
  if (write_output) {
    base_path <- "~/Git/k-similar-neighbor/data/label_analysis/"
    create_name <- function(x) sprintf("%s%s_%s.csv", base_path, label, x)
    
    fwrite(outU, file=create_name("uni_marg"))
    
    fwrite(outZ, file=create_name("std_resids"))
    
    fwrite(outI, file=create_name("std_price"))
    
    fwrite(outP, file=create_name("model_summary"))
    
    fwrite(outM, file=create_name("model_coef"))
  }
  return(list(uni_marg=outU, std_resids=outZ, std_price=outI, model_summary=outP, model_coef=outM))
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

