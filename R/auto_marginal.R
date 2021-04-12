library(jsonlite)
library(data.table)
library(rugarch)
source("R/normality_tests.R")
source("R/likelihood_models.R")


# Check to see which stocks had a successful auto_fit run
check_spec <- function(lst)
{
  norun <- length(lst) == 1
  if (norun)
    return(FALSE)
  else
    return(isTRUE(unlist(lst$pass)))
}


# Given a path containing lists of model specifications, read in and return list containing 1) good and 2) bad vector of stock symbol
good_vs_bad_symbols <- function(spec.path)
{
  specs <- read_json(spec.path)
  nms <- names(specs)
  goodspecs <- sapply(specs, check_spec)
  return(list(good_ticks=nms[goodspecs], bad_ticks=nms[!goodspecs]))
}


# Given a ticker value, this function reads in the corresponding
# json file and pulls out the values you want (e.g. "open", "close", "high")
#
# dtfU <- data.table(parse_json("GOOG"))
# u <- dtfU[Date >= st_date & Date <= ed_date, diff(log(close))]
#
parse_json <- function(tick, values="close")
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
  dtf[order(dtf$Date), ]
}


# Given a series x, automatically try various conditional mean
# and condition variance models. check whether any pass all the 
# verification tests (pass). Choose the best model and store the
# model specifications
#
# auto_fit(u)
#
auto_fit <- function(x)
{
  "ar ma sar sma period i si"
  aa <- forecast::auto.arima(x, max.p=7, allowmean=TRUE, allowdrift=FALSE)
  spec <- aa$arma
  arma <- spec[1:2]
  if (spec[3] > 0)
    arma[1] <- 7 * arma[3]
  if (spec[4] > 0)
    arma[2] <- 7 * arma[4]
  
  arch <- 1:2
  garch <- 1:2
  garchmodels <- c("sGARCH", "gjrGARCH", "csGARCH")
  distributions <- c("norm", "std", "sstd")
  outerpaste <- function(a, b) as.vector(outer(a, b, paste, sep="_"))
  rnames <- Reduce(outerpaste, list(garchmodels, arch, garch), init=distributions)
  lik <- matrix(NA_real_,
                nrow=prod(length(arch), length(garch),
                          length(garchmodels), length(distributions)),
                ncol=3,
                dimnames=list(rnames, c("bic", "likelihood", "pass_marginal")))
  
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
          
          if (inherits(fit_, "try-error")) {
            lik[nn,] <- c(NA_real_, NA_real_, FALSE)
            
          } else {
            marg_check <- try(verify_marginal_test(
              marginal_tests(fit_, print=FALSE, plot=FALSE), 0.05))
            l <- try(infocriteria(fit_)['Bayes', 1])
            LLH <- try(fit_@fit$LLH)
            
            if (inherits(l, "try-error") || inherits(marg_check, "try-error"))
              lik[nn,] <- c(NA_real_, NA_real_, FALSE)
            else
              lik[nn,] <- c(l, LLH, marg_check)
          }
        }
      }
    }
  }
  pass <- any(lik[, "pass_marginal"] == 1)
  if (pass) {
    lik2 <- lik[which(lik[, "pass_marginal"] > 0), ]
    choice <- unlist(strsplit(names(which.min(lik2[, "bic"])), "_"))
  } else {
    choice <- c("Error", "Error", "Error", "Error")
  }
  return(list(ar=arma[1], ma=arma[2], arch=choice[4], garch=choice[3],
              garchmod=choice[2], distr=choice[1], pass=pass, lik=lik))
}


# For a tick that has already been run through the auto_fit function but
# failed to find a fit that passes all the marginal tests,
# try setting the the AR component to 7
produce_rGARCHfit_s <- function(tick, start_date, end_date)
{
  
  file_ <- paste0(toupper(tick), ".json")
  
  if (file.exists(paste0("data/raw_json/", file_))) {
    pj <- parse_json(tick)
    dtfU <- data.table(pj)
    x <- dtfU[Date >= start_date & Date <= end_date, diff(log(close))]
  } else {
    stop(sprintf("JSON file doesn't exist for ticker [%s] in directory data/", tick))
  }
  
  if (!exists("SPECIFICATIONS"))
    SPECIFICATIONS <<- read_json(SPEC_PATH)
  
  spec <- SPECIFICATIONS[[toupper(tick)]]
  if (is(spec, "NULL")) {
    stop("Tick does not exist in SPECIFICATIONS.")
  }
  
  spec_mod <- ugarchspec(
    variance.model = list(model = unlist(spec$garchmod),
                          garchOrder = c(as.integer(unlist(spec$arch)), as.integer(unlist(spec$garch)))),
    mean.model = list(armaOrder = c(7, unlist(spec$ma)), include.mean = TRUE),
    distribution.model = unlist(spec$distr)
  )
  fit_norm <- ugarchfit(spec=spec_mod, data=x)
  return(fit_norm)
}


# batch process produce_rGarchfit_s over all failed marginal models
run_seasonal_modification(ticks, specification)
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
    pj <- parse_json(tick)
    dtfU <- data.table(pj)
    # No data for this tick
    if (nrow(dtfU) == 0) {
      print(sprintf("No data for ticker [%s] data/raw_json/", tick))
      return(list(fit=NULL, dates=NULL))
    }
    dtfU <- dtfU[Date >= start_date & Date <= end_date]
    x <- dtfU[, diff(log(close))]
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
  
  if (unlist(spec[["pass"]])) {
    # Fit garch model and return along with dates of analysis
    spec_mod <- ugarchspec(
      variance.model = list(model = unlist(spec$garchmod),
                            garchOrder = c(as.integer(unlist(spec$arch)), as.integer(unlist(spec$garch)))),
      mean.model = list(armaOrder = c(unlist(spec$ar), unlist(spec$ma)), include.mean = TRUE),
      distribution.model = unlist(spec$distr)
    )
    fit_ <- try(ugarchfit(spec=spec_mod, data=x))
    return(list(fit=fit_, dates=dtfU[2:.N, Date]))
  } else {
    print(sprintf("No ARIMA-GARCH specification passed for %s", tick))
    return(list(fit=NULL, dates=dtfU[2:.N, Date]))
  }
  
}


# define function that takes a spec_path and creates a matrix of each
# stocks standardized residuals
specs_to_resid_matrix <- function(path, write_to_file=TRUE)
{
  
  # grab dates from the path
  m <- gregexpr("[0-9]{8}", path)
  dts <- unlist(regmatches(path, m))
  dts2 <- sapply(dts, function(x) as.character(as.IDate(x, format="%Y%m%d")))
  
  spec_list <- read_json(path)
  good_vs_bad <- good_vs_bad_symbols(path)
  warning(sprintf("Removed the following tickers with bad specifications: [%s]",
                  paste0(good_vs_bad[['bad_ticks']], collapse=", ")))
  for (bad_tick in good_vs_bad[['bad_ticks']])
    spec_list[[bad_tick]] <- NULL
  
  tickers <- names(spec_list)
  resid_list <- vector("list", length(spec_list))
  names(resid_list) <- tickers
  
  for (tick in tickers) {
    ii <- which(tick == tickers)
    print(sprintf("Tick: %s - %d out of %d", tick, ii, length(tickers)))
    obj <- produce_rGARCHfit_from_spec(tick, spec_list, dts2[1], dts2[2])
    
    # If successfully fit from the specification
    if (inherits(obj$fit, "uGARCHfit")) {
      
      z <- PIT(obj$fit)
      tmp <- data.table(Date=obj[['dates']], Z=z)
      setnames(tmp, "Z", tick)
      
    } else {
      
      if (is.null(obj[['dates']])) {
        # No data for this tick
        print(sprintf("Removing %s from output list.", tick))
        resid_list[[tick]] <- NULL
        
      } else {
        # No ARIMA-GARCH specification passed all the marginal tests
        tmp <- data.table(Date=obj[['dates']], Z=NA_real_)
        setnames(tmp, "Z", tick) 
      }
    }
    resid_list[[tick]] <- tmp
  }
  out <- Reduce(function(x, y) merge(x, y, on="Date", all=TRUE), resid_list)
  if (write_to_file) {
    file_name <- paste0("std_resids_", dts[1], "_", dts[2], ".csv")
    file_path <- paste0("~/Git/k-similar-neighbor/data/", file_name)
    fwrite(out, file=file_path)
  }
  return(out)
}
