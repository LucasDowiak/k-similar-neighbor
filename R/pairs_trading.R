source("R/similarity_functions.R")

# Pairs-Trading ----------------------------------------------------------------
# Find the N most closely associated pairs
find_pairs <- function(M, k=10L, dist.method=c("cor", "dtw"))
{
  # M - association matrix
  dist.method <- match.arg(dist.method)
  vals <- sort(M[lower.tri(M)])
  if (dist.method == "cor") {
    topn <- rev(tail(vals, n=k))
  } else {
    topn <- head(vals, n=k) 
  }
  out <- lapply(topn, function(x) row.names(which(M == x, arr.ind=TRUE)))
  attr(out, 'summary') <- summary(vals)
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

dtfFormCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", "2012"))
dtfTradeCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", "2013"))

aa_dtw <- read_association_table("data/association_results/2012_dtw.tsv")
aa_udj <- read_association_table("data/association_results/2012_unadjusted_cor.tsv")
aa_adj <- read_association_table("data/association_results/2012_model_cor.tsv")

pairs_udj <- find_pairs(aa_udj, n=10L)
pairs_adj <- find_pairs(aa_adj, n=10L)
pairs_dtw <- find_pairs(aa_dtw, n=10L, dist.method = "dtw")



back_test_strategy <- function(label, k=10L, buy_signal=c("unadjusted_cor", "model_cor", "dtw"))
{
  buy_signal <- match.arg(buy_signal)
  trade_year <- as.character(as.integer(label) + 1)
  
  # Stat and train on 2012 data
  ## Given buy-signal read in association matrix
  assoc_file <- sprintf("data/association_results/%s_%s.tsv", label, buy_signal)
  M <- read_association_table(assoc_file)
  dtfFormCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", label))
  dtfTradeCumRet <- fread(sprintf("data/label_analysis/%s_std_price.csv", trade_year))
  Tf <- dtfFormCumRet[, .N]
  
  ## Find the top K pairs
  dm <- if (buy_signal == "dtw") "dtw" else "cor"
  top_pairs <- find_pairs(M, k=k, dist.method=dm)
  
  # Execute strategy on following year's returns
  ## Read in standard price in 2013
  out <- vector("list", length=k)
  names(out) <- sapply(top_pairs, paste, collapse="-")
  
  .open_position()
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
  
  for (ii in seq_len(k)) {
    t1 <- top_pairs[[ii]][1]
    t2 <- top_pairs[[ii]][2]
    pname <- paste(top_pairs[[ii]], collapse="-")
    
    rf1 <- dtfFormCumRet[year(Date) == as.integer(label)][[t1]]
    rf2 <- dtfFormCumRet[year(Date) == as.integer(label)][[t2]]
    rt1 <- dtfTradeCumRet[year(Date) == as.integer(trade_year)][[t1]]
    rt2 <- dtfTradeCumRet[year(Date) == as.integer(trade_year)][[t2]]
    
    open_index <- which(open_position_signal(rt1 - rt2, scale(abs(rf1 - rf2))))
    close_index <- which(close_position_signal(rt1 - rt2))
    

    operations <- list()
    tt <- 0
    while (tt < Tf) {
      # If there is an open position signal
      if (length(open_index) > 0) {
        # Open the position
        oi <- open_index[1]
        open_position <- list(action="open", timet=oi, trade_id=as.character(round(runif(1), 15) * 1e15))
        
        # Close the position
        trade <- .close_position(open_position, close_index)
        operations[[length(operations) + 1]] <- trade
        
        # Increment loop logic
        tt <- trade[action=="close", timet] + 1
        open_index <- open_index[open_index >= tt]
      } else {
        break
      }
    }
    pair_trades <- rbindlist(operations)
    pair_trades[, `:=`(pair = pname, tick1=t1, tick2=t2)]
    
    tmp <- dtfFormCumRet[pair_trades$timet, .SD, .SDcols=c('Date', t1, t2)]
    pair_trades <- cbind(pair_trades, tmp)
    out[[pname]] <- pair_trades
  }
  return(out)
}

.position <- function(type=c("open", "close"), prices)

calculate_returns <- function(bts, returns)
{
  
  position <- NULL
  nms <- names(bts)
  for (k in seq_along(bts)) {
    t1 <- unlist(strsplit(nms[[k]], "-"))[1]
    t2 <- unlist(strsplit(nms[[k]], "-"))[2]
    
    for (a in seq_along(bts[[k]])) {
      p <- dtfTradeCumRet[bts[[k]][[a]]$timet, .SD, .SDcols=c(t1, t2)]
      position <- .position(bts[[k]][[a]]$action, prices, position)
      wmax <- which.max(prices)
    }


  }
}



plot_list_of_pairs(pairs_dtw, dtfFormCumRet)
plot_list_of_pairs(pairs_dtw, dtfTradeCumRet)

tst <- back_test_strategy("2012", k=10, "dtw")
#  - given buy-signal find open and close intervals for each pair through out the
#    trading period
