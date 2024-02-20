# pairs-trading script
setwd("~/Git/dtw-in-finance/")
source("R/similarity_functions.R")
source("R/pairs_trading.R")
library(data.table)

yrs.int <- 2000:2021
yr_labels <- as.character(yrs.int)
signals <- c("unadjusted", "model", "dtw", "l2")
K <- c(25, 50, 100)
sidx <- c(1, 1, 1)
thresh <- 2

# 1
outerpaste <- function(a, b) as.vector(outer(a, b, paste, sep="_"))
rnames <- Reduce(outerpaste, list(signals, K, sidx, thresh), init=yr_labels)

lst_ann_trade_ret <- vector("list", length(rnames))
names(lst_ann_trade_ret) <- rnames

for (yr_lab in yr_labels) {
  print(sprintf("Start year %s at %s", yr_lab, Sys.time()))
  for (sig in signals) {
    for (k in seq_along(K)) {
      for (thsh in thresh) {
        tyr <- paste(c(yr_lab, sig, K[k], sidx[k], threshold=thsh), collapse="_")
        tmp <- back_test_strategy(yr_lab, k=K[k], threshold=thsh, start_index=sidx[k], buy_signal=sig, replace=TRUE)
        lst_ann_trade_ret[[tyr]] <- tmp
      }
    }
  }
}


lst_ <- lst_ann_trade_ret


revenue_streams <- rbindlist(lapply(lst_, `[[`, 1), fill=TRUE)
revenue_streams[, group := paste(start_index, n_pairs, sep=" : ")]

benchmarks <- rbindlist(lapply(lst_, `[[`, 2), fill=TRUE)
benchmarks[, group := paste(start_index, n_pairs, sep=" : ")]
benchmarks[, threshold := as.character(threshold)]
benchmarks[, buy_signal:=factor(buy_signal, levels=c("unadjusted_cor", "model_cor", "dtw", "l2"))]
benchmarks[, group:=factor(group, levels=c("1 : 25", "1 : 50", "1 : 100"))]
benchmarks[, .(baseline_long=mean(baseline_long - 1),
               committed=mean(committed - 1),
               invested=mean(invested - 1)), by=c("buy_signal", "group")]
benchmarks[, trade_year := as.integer(trade_year )]

pairs <- rbindlist(lapply(lst_, `[[`, 3), fill=TRUE)
pairs[, group := paste(start_index, n_pairs, sep=" : ")]
pairs[, threshold := as.character(threshold)]

levs <- pairs[, sort(unique(c(sector_1, sector_2)))]
pairs[, sector_1 := factor(sector_1, levels=levs)]
pairs[, sector_2 := factor(sector_2, levels=levs)]

# Percent of formed pairs that are traded
tmpP <- pairs[, as.list(table(is.na(pair_return))), by=c("buy_signal", "group")]
tmpP[, percent :=  `FALSE` / rowSums(.SD), .SDcols=c("FALSE", "TRUE")]

# Standard Errors for the average can be found via regression
lm1 <- lm((baseline_long - 1) ~ -1 + buy_signal : group, data=benchmarks[threshold=="2"])
W <- sandwich::NeweyWest(lm1, lag=3, prewhite=FALSE, order.by=benchmarks[threshold=="2"]$trade_year)
round(lmtest::coeftest(lm1, df=Inf, vcov. = W), 4)


# Calculate annual returns of the portfolio
# -------------------------------------------------------------

distribution_stats <- function(x)
{
  data.table(mean=mean(x),
             sd=sd(x),
             median=median(x),
             skew=moments::skewness(x),
             kurtosis=moments::kurtosis(x),
             mininum=min(x),
             maximum=max(x))
}

dtfRbs <- benchmarks[threshold=="2", distribution_stats(baseline_long - 1), by=c("buy_signal", "group")]
dtfRbs[, metric := "baseline_long"]
dtfRcom <- benchmarks[threshold=="2", distribution_stats(committed - 1), by=c("buy_signal", "group")]
dtfRcom[, metric := "committed"]
dtfRinv <- benchmarks[threshold=="2", distribution_stats(invested - 1), by=c("buy_signal", "group")]
dtfRinv[, metric := "invested"]
dtfR <- rbindlist(list(dtfRbs, dtfRcom, dtfRinv))
dtfR[order(group)][metric == "invested"]
dtfR[order(group)][metric == "invested", round(.SD, 3), .SDcols=4:9]

# -------------------------------------------------------------


# Plot annual returns by group and buy signal
# -------------------------------------------------------------
library(ggplot2); library(colorBlindness); library(viridisLite)
benchmarks[, portfolio_returns := committed - 1]
benchmarks[, trade_year_int := as.integer(trade_year)]
y_lab <- "Return on Committed Capital"
# Annual returns on invested capital by buy-signal
p1 <- ggplot(data=benchmarks[group=="1 : 25"], aes(x=trade_year, y=portfolio_returns, fill=buy_signal)) +
  geom_col(width=0.5, position=position_dodge()) +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.title=element_text(size=8)) + 
  scale_fill_viridis_d() +
  ylim(c(-0.1, .12)) +
  ggtitle("Top 25 Pairs") +
  xlab(label=NULL) + ylab(label=y_lab)


# Annual returns on committed capital by buy-signal
p2 <- ggplot(data=benchmarks[group=="1 : 50"]) +
  geom_col(aes(x=trade_year, y=portfolio_returns, fill=buy_signal), width=0.5, position=position_dodge()) +
  theme(legend.position="none", legend.title=element_blank(),
        axis.title=element_text(size=8)) + 
  scale_fill_viridis_d() +
  ylim(c(-0.1, .12)) +
  ggtitle("Top 50 Pairs") +
  xlab(label=NULL) + ylab(label=y_lab)


# Annual returns on invested capital by buy-signal
p3 <- ggplot(data=benchmarks[group=="1 : 100"]) +
  geom_col(aes(x=trade_year, y=portfolio_returns, fill=buy_signal), width=0.5, position=position_dodge()) +
  theme(legend.position="none", legend.title=element_blank(),
        axis.title = element_text(size=8)) +
  scale_fill_viridis_d() +
  ylim(c(-0.1, .12)) +
  ggtitle("Top 100 Pairs") +
  xlab(label=NULL) + ylab(label=y_lab)


ggpubr::ggarrange(p1, p2, p3,
                  nrow=3, ncol=1,
                  legend="bottom", common.legend=TRUE)

# -------------------------------------------------------------


# Summarize pair formation between buy-signals
# -------------------------------------------------------------

# Pair Matrix
dcast(pairs[group=="1 : 25" & buy_signal=="unadjusted_cor" & threshold=="2"],
      sector_1 ~ sector_2,
      value.var="pair",
      fun.aggregate=length)
M <- dcast(pairs[group=="1 : 25" & buy_signal=="unadjusted_cor" & threshold=="2"],
      sector_1 ~ sector_2,
      value.var="pair",
      fun.aggregate=length)[, -1]
M <- as.matrix(M)
cbind(colSums(M), round(colSums(M) / sum(M), 3))
M <- (M + t(M)) - diag(diag(M))
M[upper.tri(M)] <- 0
unname(M)

# Marginal Sector Counts
get_counts <- function(x)
{
  cnts <- table(x)
  cbind(cnts, 100 * round(cnts / sum(cnts), 3))
}
pairs[group=="1 : 25" & buy_signal=="l2" & threshold=="2", get_counts(c(sector_1, sector_2))]


UNA <- M
ADJ <- M
DTW <- M
# EUC <- M
# -------------------------------------------------------------
