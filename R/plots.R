library(ggplot2)
library(data.table)
setwd("~/Git/k-similar-neighbor/")


# Read in the raw S&P 500 data file
#
#   ticks - optional subset of ticks to return
#
read_SandP_data <- function(ticks=NULL)
{
  out <- fread("data/SandP_companies.csv")
  out <- out[order(sector, subindustry, ticker)]
  if (!is(ticks, "NULL"))
    out <- out[ticker %in% ticks]
  return(out)
}


# Plot the similarity matrices
#
#   assoc_file - file path for stored similarity matrix
#
#   type - type of similarity matrix, which effects certain plot parameters
#
#   tits- optional title for the plot
#
plot_sim_matrix <- function(assoc_file, tickers=NULL, tits="")
{
  require(ggplot2)
  if (is.matrix(assoc_file)) {
    X <- assoc_file
  } else if (is.character(assoc_file)) {
    X <- read.table(assoc_file, header=TRUE, row.names=1)
  } else {
    stop("assoc_file must be a matrix or a file path to an association file")
  }
  x <- c(X[upper.tri(X)])
  midpoint <- median(x)
  limit <- c(min(x), quantile(x, 0.75) + (IQR(x) * 1.5))

  dtfStock <- read_SandP_data(names(X))
  
  if (is.null(tickers)) {
    X <- X[dtfStock$ticker, dtfStock$ticker]
  } else {
    X <- X[tickers, tickers]
  }
  dtfStock[, idx := 1:.N]
  X$ticks <- names(X)
  X$ticks <- factor(X$ticks, levels=X$ticks)
  
  X <- as.data.table(X)
  idxsec <- which(dtfStock[1:(.N - 1), sector] != dtfStock[2:.N, sector])
  idxlab <- dtfStock[, unique(sector)]
  idxloc <- dtfStock[, round(mean(idx)), by=sector][, V1]
  
  X <- melt(X, id.vars="ticks", measure.vars=setdiff(names(X), "ticks"))
  
  p <- ggplot(data=X, aes(x=ticks, y=variable, fill=value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = midpoint, limit = limit, space = "Lab", 
                         name="") +
    theme_minimal() +
    #theme(axis.text.x=element_text(angle = -15)) +
    theme(axis.text=element_blank()) +
    annotate(geom="text", x=idxloc, y=idxloc, label=idxlab) +
    scale_x_discrete(breaks=dtfStock$ticker[idxloc], labels=idxlab) + 
    scale_y_discrete(breaks=dtfStock$ticker[idxloc], labels=idxlab) + 
    geom_vline(xintercept=idxsec, size=0.5, lty=3) + 
    geom_hline(yintercept=idxsec, size=0.5, lty=3) +
    labs(x="", y="", title=tits)
  print(p)
}



# Plot S&P time series by group
#
#   D - distance matrix
# 
#   DT - data.table with price data
#
#   k - number of clusters
#

WW <- W_SCB
hc <- agnes(as.dist(WW), method="ward")
cluster_and_plot_series <- function(D, DT, k=3, return_melted=FALSE, ...)
{
  require(ggplot2)
  hc <- agnes(as.dist(D), ...)
  dtfsp <- read_SandP_data()[, .(ticker, security, sector, subindustry)]
  setnames(dtfsp, "ticker", "tick")
  
  dtfSubGrp <- merge(dtfsp,
                     data.table(group=cutree(tree=hc, k=k), tick=colnames(D)))
  M <- melt(DT, id.vars="Date", variable.name="tick", value.name="std.price")
  M <- merge(dtfSubGrp, M)
  
  p <- ggplot() +
    geom_line(data=M, aes(x = Date, y=std.price, group=tick, color=sector), alpha=0.75) +
    facet_wrap(~group, scales = "free")
  plot(p)
  if (return_melted)
    return(M)
}
cluster_and_plot_series(W_SCB, dtfStdPrc1920, k=12, method="ward")
cbind(colnames(W_0), hc$order.lab)