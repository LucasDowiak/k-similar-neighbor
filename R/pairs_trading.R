# Pairs-Trading ----------------------------------------------------------------
# Find the N most closely associated pairs
find_pairs <- function(M, k=10L, si=1, dist.method=c("cor", "dtw", "l2"), replace=TRUE,
                       omit_ticks=NULL)
{
  # M - association matrix
  dist.method <- match.arg(dist.method)
  cnames <- colnames(M)
  if (!is.null(omit_ticks)) {
    cnames <- cnames[!cnames %in% omit_ticks]
    M <- M[cnames, cnames]
  }
  M[upper.tri(M, diag=TRUE)] <- NA_real_
  M <- as.data.table(reshape2::melt(M))[!is.na(value)]
  M <- M[, `:=`(Var1=as.character(Var1), Var2=as.character(Var2))]
  #> vals <- sort(M[lower.tri(M)])
  #> svals <- summary(vals)
  if (dist.method == "cor") {
    M <- M[order(-value)]
  } else {
    M <- M[order(value)]
  }
  if (replace) {
    # take top k values
    topn <- M[si:(si + k - 1L), value]
    out <- lapply(split(M[si:(si + k - 1L), .SD], 1:k), function(x) c(x$Var1, x$Var2))
    out <- unname(out)
    #out <- M[si:(si + k - 1L), .(Var1, Var2)][, apply(.SD, 1, c)]
    #out <- as.list(out)
    svals <- M[, summary(value)]
    #> topn <- vals[si:(si + k - 1L)]
    #> out <- lapply(topn, function(x) row.names(which(M == x, arr.ind=TRUE)))
  } else {
    # while loop that rejects a pair if one of it's stocks is already part of the portfolio
    MM <- vector("list", floor(length(cnames) / 2))
    used_ticks <- c()
    for (ii in seq_along(MM)) { #length(avail_ticks) > 0) {
      used_ticks <- c(used_ticks, M[1, c(Var1, Var2)])
      MM[[ii]] <- M[1, ]
      M <- M[!(Var1 %in% used_ticks | Var2 %in% used_ticks)]
    }
    MM <- rbindlist(MM)
    topn <- MM[si:(si + k - 1L), value]
    out <- lapply(split(MM[si:(si + k - 1L), .SD], 1:k), function(x) c(x$Var1, x$Var2))
    out <- unname(out)
    svals <- MM[, summary(value)]
    #> out <- vector("list", k)
    #> topn <- c()
    #> ii <- 1
    #> while (any(sapply(out, is.null))) {
    #>   p <- vals[1]
    #>   ticks <- row.names(which(M == p, arr.ind=TRUE))
    #>   if (!any(ticks %in% unlist(out))) {
    #>     out[[ii]] <- ticks
    #>     topn <- c(topn, p)
    #>     ii <- ii + 1
    #> }
    #> vals <- vals[-1]
  }
  attr(out, 'summary') <- svals
  attr(out, 'values') <- topn
  attr(out, "dist.method") <- dist.method
  return(out)
}


.pull_ticks_and_melt <- function(pair, DT)
{
  out <- copy(DT[, .SD, .SDcols=c("Date", pair)])
  out <- melt(out, id.vars="Date")
  out[, pair := paste(pair, collapse="-")]
  setnames(out, "variable", "tick")
  return(out)
}


plot_list_of_pairs <- function(pairs, DT)
{
  require(ggplot2)
  dtfplot <- rbindlist(lapply(pairs, .pull_ticks_and_melt, DT=DT))
  fpairs <- sapply(pairs, function(x) paste(x, collapse="-"))
  p <- ggplot() +
    geom_line(data=dtfplot, aes(x=Date, y=value, group=tick)) +
    facet_wrap(~factor(pair, levels=fpairs), scales = "free")
  plot(p)
}


open_position_signal <- function(d, Z, threshold=2)
{
  z <- (abs(d) - attr(Z, "scaled:center")) / attr(Z, "scaled:scale")
  return(z > threshold)
}


close_position_signal <- function(d)
{
  return(c(NA, diff(sign(d)) != 0))
}


.calculate_return <- function(DT, short_premium=0.025)
{
  t1 <- DT[, unique(tick1)]; t2 <- DT[, unique(tick2)]
  strike_prices <- DT[action=="open", unlist(.SD), .SDcols=c(t1, t2)]
  spot_prices <- DT[action=="close", unlist(.SD), .SDcols=c(t1, t2)]
  long_pos <- names(which.min(strike_prices))
  short_pos <- names(which.max(strike_prices))
  
  long_ret <- spot_prices[long_pos] - strike_prices[long_pos]
  short_ret <- strike_prices[short_pos] - spot_prices[short_pos]
  return(data.table(pair=DT[1, pair],
                    open_date=DT[action=="open", Date],
                    close_date=DT[action=="close", Date],
                    long_pos=long_pos,
                    long_strike=strike_prices[long_pos],
                    long_spot=spot_prices[long_pos],
                    long_ret=long_ret,
                    short_pos=short_pos,
                    short_strike=strike_prices[short_pos],
                    short_spot=spot_prices[short_pos],
                    short_ret=short_ret,
                    pair_return=long_ret + short_ret))
}


back_test_pair <- function(pairs, dt_form, dt_trade, threshold=2)
{
  t1 <- pairs[1]; t2 <- pairs[2]
  rf1 <- dt_form[[t1]]; rf2 <- dt_form[[t2]]
  rt1 <- dt_trade[[t1]]; rt2 <- dt_trade[[t2]]
  open_index <- which(open_position_signal(rt1 - rt2, scale(abs(rf1 - rf2)),
                                           threshold=threshold))
  close_index <- which(close_position_signal(rt1 - rt2))
  pname <- paste(pairs, collapse="-")
  
  .close_position <- function(op, cls)
  {
    if (any(cls[cls > op$timet])) {
      ci <- cls[cls > op$timet][1]
    } else {
      ci <- Tf
    }
    cp <- list(action="close", timet=ci, trade_id=op[["trade_id"]])
    return(rbindlist(list(op, cp)))
  }
  
  operations <- list()
  Tf <- dt_trade[, .N]
  tt <- 0
  while (tt < Tf) {
    # If there is an open position signal
    if (length(open_index) > 0) {
      # Open the position
      oi <- open_index[1]
      open_position <- list(action="open", timet=oi, trade_id=as.character(round(runif(1), 15) * 1e15))
      
      # Close the position
      trade <- .close_position(open_position, close_index)
      # aa <- excess_return(trade, dtfTradeCumRet)
      operations[[length(operations) + 1]] <- trade
      
      # Increment loop logic
      tt <- trade[action=="close", timet] + 1
      open_index <- open_index[open_index >= tt]
    } else {
      break
    }
  }
  if (length(operations) == 0L) {
    return(data.table())
  } else {
    pair_trades <- rbindlist(operations)
    pair_trades[, `:=`(pair=pname, tick1=t1, tick2=t2)]
    tmp <- dt_trade[pair_trades$timet, .SD, .SDcols=c('Date', t1, t2)]
    pair_trades <- cbind(pair_trades, tmp)
    revenues <- pair_trades[, .calculate_return(.SD), by=trade_id]
  }
  
  
  return(revenues)
}


back_test_strategy <- function(label, k=10L, start_index=1, threshold=2, replace=TRUE,
                               buy_signal=c("unadjusted_cor", "model_cor", "dtw", "l2"))
{
  buy_signal <- match.arg(buy_signal)
  trade_year <- as.character(as.integer(label) + 1)
  
  ## Given buy-signal read in association matrix
  assoc_file <- sprintf("data/association_results/%s_%s.tsv", label, buy_signal)
  M <- read_association_table(assoc_file)
  dtfSP <- fread("data/SandP_companies.csv")
  dtfFormCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", label))
  dtfTradeCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", trade_year))
  Tf <- dtfTradeCumRet[, .N]
  bad_ticks <- colnames(M)[!colnames(M) %in% names(dtfTradeCumRet)]
  
  ## Find the top K pairs
  dm <- if (grepl("cor", buy_signal)) "cor" else buy_signal
  top_pairs <- find_pairs(M, k=k, si=start_index, dist.method=dm, replace=replace,
                          omit_ticks=bad_ticks)
  pnames <- unlist(lapply(top_pairs, paste, collapse="-"))
  
  # Execute strategy on following year's returns
  ## Read in standard price in 2013
  revn <- vector("list", length=k); names(revn) <- pnames
  for (ii in seq_len(k)) {
    pn <- paste(top_pairs[[ii]], collapse="-")
    revn[[pn]] <- back_test_pair(top_pairs[[ii]], dtfFormCumRet, dtfTradeCumRet,
                                 threshold=threshold)
  }
  
  # What annual return do we get by going long in the portfolio for the year
  bl_long <- data.table(pair=pnames,
                        baseline_return=sapply(top_pairs, base_line_long_return,
                                               dt_trade=dtfTradeCumRet))
  
  # data.table to record pair metadata
  dtfPairEval <- data.table(pair=pnames)
  dtfPairEval[, trade_year := trade_year]
  dtfPairEval[, threshold := threshold]
  dtfPairEval[, start_index := start_index]
  dtfPairEval[, n_pairs := k]
  dtfPairEval[, buy_signal := buy_signal]
  
  # Return metrics by pair
  if (all(sapply(revn, length) == 0)) {
    revenues <- data.table()
    dtfPairEval[, pair_return := NA_real_]
    commited_return <- invested_return <- 1
  } else {
    revenues <- rbindlist(revn)[order(pair, close_date)]
    rev_by_pair <- revenues[, .(pair_return=prod(1 + pair_return)), by=pair]
    dtfPairEval <- merge(dtfPairEval, rev_by_pair, all.x=TRUE)
    
    # Provide invested and committed return values
    invested_return <- dtfPairEval[, mean(pair_return, na.rm=TRUE)]
    commited_return <- dtfPairEval[, (sum(pair_return, na.rm=TRUE) + sum(is.na(pair_return))) / k]
  }
  dtfPairEval[, c("pair_1", "pair_2") := tstrsplit(pair, "-", fixed=TRUE)]
  dtfPairEval <- merge(dtfPairEval, dtfSP[, .(ticker, sector)], by.x="pair_1", by.y="ticker")
  dtfPairEval <- merge(dtfPairEval, dtfSP[, .(ticker, sector)], by.x="pair_2", by.y="ticker")
  dtfPairEval <- merge(dtfPairEval, bl_long, by="pair")
  setnames(dtfPairEval, c("sector.x", "sector.y"), c("sector_1", "sector_2"))
  
  dtfBenchmarks = data.table(
             trade_year=trade_year,
             threshold=threshold,
             start_index=start_index,
             n_pairs=k,
             buy_signal=buy_signal,
             baseline_long=bl_long[, mean(baseline_return)],
             committed=commited_return,
             invested=invested_return)

  revenues[, trade_year := trade_year]
  revenues[, threshold := threshold]
  revenues[, start_index := start_index]
  revenues[, n_pairs := k]
  revenues[, buy_signal := buy_signal]
  if ("pair" %in% names(revenues)) {
    revenues[, c("pair_1", "pair_2") := tstrsplit(pair, "-", fixed=TRUE)]
    revenues <- merge(revenues, dtfSP[, .(ticker, sector)], by.x="pair_1", by.y="ticker")
    revenues <- merge(revenues, dtfSP[, .(ticker, sector)], by.x="pair_2", by.y="ticker")
    setnames(revenues, c("sector.x", "sector.y"), c("sector_1", "sector_2"))
  }

  
  return(list(revenues=revenues, annual_returns=dtfBenchmarks, pairs=dtfPairEval))
}


base_line_long_return <- function(portfolio_pairs, dt_trade)
{
  pos <- unlist(lapply(portfolio_pairs, strsplit, "-"))
  wts <- table(pos)
  rtn <- dt_trade[.N, unlist(.SD), .SDcols=unique(pos)]
  ret <- (rtn[names(wts)] %*% wts)[1, 1] / sum(wts)
  return(1 + ret)
}
