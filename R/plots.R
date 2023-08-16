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
plot_sim_matrix <- function(assoc_file, type=c("cor", "dtw"), tits="")
{
  require(ggplot2)
  type <- match.arg(type)
  if (type == "cor") {
    midpoint <- 0
    limit <- c(-1, 1)
  } else if (type == "dtw") {
    midpoint <- 100
    limit <- c(0, 400)
  }
  X <- read.table(assoc_file, header=TRUE, row.names=1)
  dtfStock <- read_SandP_data(names(X))
  dtfStock[, idx := 1:.N]
  
  X <- X[dtfStock$ticker, dtfStock$ticker]
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



# Plot an example of the DTW algorithm at work
#
#   x, y (numeric) - time series to plot
#
#   standardize (boolean) - should the series be standardized by dividing values by setting x[1] = 1
#
#   ... - additional parameters passed to the DTW algorithm
#
plot_dtw_ts <- function(x, y, standardize=TRUE, ...)
{
  if (standardize) {
    x <- x / x[1]; y <- y / y[1]
  }
  pdtw <- dtw(x, y, keep.internals=TRUE, ...)
  
  yrange <- range(x, y)
  padding <- diff(yrange) * 0.05
  plot(x, type="l", ylim=c(yrange[1] - padding, yrange[2] + padding),
       main="DTW Index Alignment")
  lines(y, col="red")
  grid()
  segments(pdtw$index1, pdtw$query[pdtw$index1],
           pdtw$index2, pdtw$reference[pdtw$index2],
           col="gray80")
}
tmp <- dtfSP[Date >= "2019-01-01" & Date <= "2019-06-30", .(Date, A, GE)]
tmp[, A := A / A[1]]
tmp[, GE := GE / GE[1]]
pdtw <- dtw(tmp$A, tmp$GE, keep.internals = TRUE,
            window.type = "sakoechiba", window.size=15)

plot(tmp$GE, type="l", ylim=c(0.9, 1.5))
lines(tmp$A, col="red")
segments(pdtw$index1, pdtw$query[pdtw$index1],
         pdtw$index2, pdtw$reference[pdtw$index2],
         col="gray80")

