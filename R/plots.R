library(ggplot2)
library(data.table)
setwd("~/Git/dtw-in-finance/")


# Missing map
# Amelia::missmap(dtfSP, ylab="Feb 2021 - Nov 1999", xlab="Tick")

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
cluster_and_plot_series <- function(D, Price, k=3, return_melted=FALSE, ...)
{
  require(ggplot2)
  hc <- agnes(as.dist(D), ...)
  dtfsp <- read_SandP_data()[, .(ticker, security, sector, subindustry)]
  setnames(dtfsp, "ticker", "tick")
  
  dtfSubGrp <- merge(dtfsp,
                     data.table(group=cutree(tree=hc, k=k), tick=colnames(D)))
  M <- melt(Price, id.vars="Date", variable.name="tick", value.name="std.price")
  M <- merge(dtfSubGrp, M)
  
  p <- ggplot() +
    geom_line(data=M, aes(x = Date, y=std.price, group=tick, color=sector), alpha=0.75) +
    facet_wrap(~group, scales = "free")
  plot(p)
  if (return_melted)
    return(M)
}



# Helper-function for plotting line segments in the dynamic time warping graphs
draw_square <- function(x, y)
{
  o <- 1/2
  segments(x0 = c(x - o, x - o, x + o, x + o),
           y0 = c(y - o, y - o, y + o, y + o),
           x1 = c(x - o, x + o, x + o, x - o),
           y1 = c(y + o, y - o, y - o, y + o),
           col="white")
}


# Modified plot from the dtw package
dtwPlotTwoWay2 <- function (d, xts = NULL, yts = NULL, offset = 0, ts.type = "l", 
                            pch = 21, match.indices = NULL, match.col = "gray70", match.lty = 3, 
                            xlab = "Index", ylab = "Query value", ...) 
{
  if (is.null(xts) || is.null(yts)) {
    xts <- d$query
    yts <- d$reference
  }
  if (is.null(xts) || is.null(yts)) 
    stop("Original timeseries are required")
  ytso <- yts + offset
  maxlen <- max(length(xts), length(ytso))
  length(xts) <- maxlen
  length(ytso) <- maxlen
  # def.par <- par(no.readonly = TRUE)
  if (offset != 0) {
    par(mar = c(5, 4, 4, 4) + 0.1)
  }
  matplot(cbind(xts, ytso), type = ts.type, pch = pch, xlab = xlab, 
          ylab = ylab, axes = FALSE, ...)
  box()
  axis(1)
  axis(2, at = pretty(xts))
  if (offset != 0) {
    rightTicks <- pretty(yts)
    axis(4, at = rightTicks + offset, labels = rightTicks)
  }
  if (is.null(match.indices)) {
    ml <- length(d$index1)
    idx <- 1:ml
  }
  else if (length(match.indices) == 1) {
    idx <- seq(from = 1, to = length(d$index1), length.out = match.indices)
  }
  else {
    idx <- match.indices
  }
  # For the diss only show every-other line to make it more visually appealing
  idx <- idx[idx %% 2 == 0]
  segments(d$index1[idx], xts[d$index1[idx]], d$index2[idx], 
           ytso[d$index2[idx]], col = match.col, lty = match.lty)
  # par(def.par)
}

# Same as dtwPlotTwoWay2 but with line segments showing normal synchroneous time alignment
plotTwoWayEuclid <- function (d, offset = 0, ts.type = "l", 
                               pch = 21, match.indices = NULL, match.col = "gray70", match.lty = 3, 
                               xlab = "Index", ylab = "Query value", ...) 
{
  xts <- d$query
  yts <- d$reference
  ytso <- yts + offset
  maxlen <- max(length(xts), length(ytso))
  length(xts) <- maxlen
  length(ytso) <- maxlen
  # def.par <- par(no.readonly = TRUE)
  if (offset != 0) {
    par(mar = c(5, 4, 4, 4) + 0.1)
  }
  matplot(cbind(xts, ytso), type = ts.type, pch = pch, xlab = xlab, 
          ylab = ylab, axes = FALSE, ...)
  box()
  axis(1)
  axis(2, at = pretty(xts))
  if (offset != 0) {
    rightTicks <- pretty(yts)
    axis(4, at = rightTicks + offset, labels = rightTicks)
  }
  if (is.null(match.indices)) {
    ml <- length(d$index1)
    idx <- 1:ml
  }
  else if (length(match.indices) == 1) {
    idx <- seq(from = 1, to = length(d$index1), length.out = match.indices)
  }
  else {
    idx <- match.indices
  }
  # For the diss only show every-other line to make it more visually appealing
  idx <- idx[idx %% 2 == 0]
  segments(1:length(xts), xts[1:length(xts)], 1:length(ytso), ytso[1:length(ytso)], 
           col = match.col, lty = match.lty)
  # par(def.par)
}



# Missing data map used in the Data chapter
dtf <- fread("data/SandP_log_return_history.csv")
dtf[, first_day_of_year := Date == min(Date), by=year(Date)]
numNAs <- dtf[, sapply(.SD, function(x) sum(is.na(x))), .SDcols=setdiff(names(dtf), "Date")]
numNAs <- sort(numNAs, decreasing=TRUE)
dtf[, Amelia::missmap(.SD, y.labels=2000:2023, y.at=which(first_day_of_year)[-1],
                      ylab="", xlab="Stock Symbol"),
    .SDcols=names(numNAs)]







# Vanilla DTW vs Weighted DTW
t1 <- "TER"; t2 <- "LRCX"
yr_ <- 2021
dtfReturns <- fread("data/SandP_log_return_history.csv")
m1 <- readRDS(sprintf("data/model_objects/%d/%d_%s.rds", yr_, yr_, t1))
m2 <- readRDS(sprintf("data/model_objects/%d/%d_%s.rds", yr_, yr_, t2))

P <- dtfReturns[year(Date) == yr_, t1, with=F]
Q <- dtfReturns[year(Date) == yr_, t2, with=F]
p <- cumprod(1 + P)[[1]]
q <- cumprod(1 + Q)[[1]]
N <- length(p)

unc_dtw <- dtw(p, q, keep.internals=TRUE, window.type="none")
scb_dtw <- dtw(p, q, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=15)


# 3x2 plot to compare Euclidean Distance, Constrained DTW, and Unconstrained DTW
par(mfrow=c(3, 2))

plotTwoWayEuclid(scb_dtw, xlab=as.character(yr_), match.lty=2, match.col="gray70",
                 ylab="Standardized Price",
                 main="A. Standard Squared Distance")
legend("topleft", legend=c(t1, t2), col=c("black", "red"), lty=c(1,2), bty="n")
legend("bottomright", legend=sprintf("D=%.*f", 2, sum((p - q)**2)), bty="n")

plot(seq_along(p), seq_along(q), type="l",
     main="A. Optimal Warping Path",
     xlab=sprintf("%s Index", t1),
     ylab=sprintf("%s Index", t2))
abline(15, 1, lty=3)
abline(-15, 1, lty=3)


dtwPlotTwoWay2(scb_dtw, xlab=as.character(yr_), match.lty=2, match.col="gray70",
               ylab="Standardized Price",
               main="B. Constrained DTW")
legend("topleft", legend=c(t1, t2), col=c("black", "red"), lty=c(1,2), bty="n")
legend("bottomright", legend=sprintf("D=%.*f", 2, scb_dtw$distance), bty="n")

plot(scb_dtw, type="alignment",
     main="B. Constrained Optimal Warping Path",
     xlab=sprintf("%s Index", t1), ylab=sprintf("%s Index", t2))
abline(15, 1, lty=3)
abline(-15, 1, lty=3)


dtwPlotTwoWay2(unc_dtw, xlab=as.character(yr_), match.lty=2, match.col="gray90",
               ylab="Standardized Price",
               main="C. Unconstrained DTW")
legend("topleft", legend=c(t1, t2), col=c("black", "red"), lty=c(1,2), bty="n")
legend("bottomright", legend=sprintf("D=%.*f", 2, unc_dtw$distance), bty="n")

plot(unc_dtw, type="alignment",
     main="C. Unconstrained Optimal Warping Path",
     xlab=sprintf("%s Index", t1), ylab=sprintf("%s Index", t2))
abline(15, 1, lty=3)
abline(-15, 1, lty=3)


# 2x2 plot used in the Distance Measure chapter
X <- scale(cbind(P, Q))

par(mfrow=c(2,2))

# Top-left
levs <- c(0.0001, 0.001, 0.01, 0.05, 0.25, 0.5)
f1 <- fitdist(distribution="std", X[, 1])
u1 <- pdist("std", X[,1], mu=f1$pars['mu'], sigma=f1$pars['sigma'], shape=f1$pars['shape'])
f2 <- fitdist(distribution="std", X[, 2])
u2 <- pdist("std", X[,2], mu=f2$pars['mu'], sigma=f2$pars['sigma'], shape=f2$pars['shape'])
estTCop_nomodel <- BiCopEst(u1, u2, family=2, se=TRUE)
plot(estTCop_nomodel, "contour", levels=levs, drawlabels=FALSE,
     xlim=c(-4,4), ylim=c(-4,4), xlab=t1, ylab=t2,
     main="A. Unadjusted Cor")
points(X, cex=0.2)
legend("topleft", expression(rho==0.86, nu==4.79), bty="n")

# Top-right
estTCop_garch <- BiCopEst(pit(m1$model), pit(m2$model), family=2, se=TRUE)
plot(estTCop_garch, "contour", levels=levs, drawlabels=FALSE,
     xlim=c(-4,4), ylim=c(-4,4), xlab=t1, ylab=t2,
     main="B. Modeled Cor")
points(X, cex=0.2)
legend("topleft", expression(rho==0.84, nu==9.01), bty="n")


# Bottom-Left
SIG <- cbind(m1$model@fit$sigma**2, m2$model@fit$sigma**2)
plot(SIG[,1], type="l", ylim=c(0, 1.05 * max(SIG)), xlab=as.character(yr_),
     ylab="Variance",
     main="C. Fitted Variance")
lines(SIG[,2], col="red", lty=2)
legend("bottomleft", legend=c(t1, t2), col=c("black", "red"), lty=c(1,2), bty="n")


# Bottom-Right
MU <- cbind(m1$model@fit$fitted.values, m2$model@fit$fitted.values)
plot(MU[,2], type="l", ylim=c(0.95 * min(MU), 1.05 * max(MU)), xlab=as.character(yr_),
     ylab="Expected Return", col="red", lty=2,
     main="D. Fitted Mean")
lines(MU[,1], col="black", lty=1, lwd=2)
legend("topright", legend=c(t1, t2), col=c("black", "red"), lty=c(1,2), bty="n")



# Cost Matrix and Accumulative Cost Matrix used in the Distance Measure chapter
par(mfrow=c(1,2))
tmpC <- round(scb_dtw$localCostMatrix[1:10, 1:10], 3)
tmpA <- round(scb_dtw$costMatrix[1:10, 1:10], 3)
par(mfrow=c(1, 2))
image(x=1:10, y=1:10, z=tmpC, col=gray.colors(50),
      main="Cost Matrix",
      xlab="Teradyne Time Index", ylab="Lam Research Time Index",
      bty="n")
for (i in 1:10) {
  for (j in 1:10) {
    text(i, j, labels=as.character(tmpC[i,j]), cex=0.75,
         col=ifelse(tmpC[i,j] >= 0.02, "black", "gray"))
  }
}
for (i in 1:13) {
  draw_square(scb_dtw$index1[i], scb_dtw$index2[i])
}

image(x=1:10, y=1:10, z=tmpA, col=gray.colors(50),
      main="Accumulative Cost Matrix",
      xlab="Teradyne Time Index", ylab="Lam Research Time Index",
      bty="n")
for (i in 1:10) {
  for (j in 1:10) {
    text(i, j, labels=as.character(tmpA[i,j]), cex=0.75,
         col=ifelse(tmpA[i,j] >= 0.1, "black", "gray"))
  }
}
for (i in 1:13) {
  draw_square(scb_dtw$index1[i], scb_dtw$index2[i])
}

