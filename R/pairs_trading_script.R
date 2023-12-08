# pairs-trading script
setwd("~/Git/dtw-in-finance/")
source("R/similarity_functions.R")
source("R/pairs_trading.R")
library(data.table)

yrs.int <- 2000:2021
yr_labels <- as.character(yrs.int)
signals <- c("unadjusted", "model", "dtw")
K <- c(25, 50, 100, 100)
sidx <- c(1, 1, 101, 1001)
thresh <- seq(1.5, 3., 0.1)

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
        tmp <- back_test_strategy(yr_lab, k=K[k], threshold=thsh, start_index=sidx[k], buy_signal=sig)
        lst_ann_trade_ret[[tyr]] <- tmp
      }
    }
  }
}

# 2
outerpaste <- function(a, b) as.vector(outer(a, b, paste, sep="_"))
rnames <- Reduce(outerpaste, list(thresh, K, sidx), init=yr_labels)

lst_ann_ret_dtw <- vector("list", length(rnames))
names(lst_ann_ret_dtw) <- rnames

for (yr_lab in yr_labels) {
  print(sprintf("Start year %s at %s", yr_lab, Sys.time()))
  for (thsh in thresh) {
    for (k in seq_along(K)) {
      tyr <- paste(c(yr_lab, thsh, K[k], sidx[k]), collapse="_")
      tmp <- back_test_strategy(yr_lab, k=K[k], threshold=thsh, start_index=sidx[k], buy_signal="dtw")
      lst_ann_ret_dtw[[tyr]] <- tmp
    }
  }
}


lst_ <- lst_ann_trade_ret


revenue_streams <- rbindlist(lapply(lst_, `[[`, 1), fill=TRUE)
revenue_streams[, group := paste(start_index, n_pairs, sep=" : ")]

benchmarks <- rbindlist(lapply(lst_, `[[`, 2), fill=TRUE)
benchmarks[, group := paste(start_index, n_pairs, sep=" : ")]
benchmarks[, threshold := as.character(threshold)]
benchmarks[, buy_signal:=factor(buy_signal, levels=c("unadjusted_cor", "model_cor", "dtw"))]
benchmarks[, group:=factor(group, levels=c("1 : 25", "1 : 50", "101 : 100", "1001 : 100"))]

pairs <- rbindlist(lapply(lst_, `[[`, 3), fill=TRUE)
pairs[, group := paste(start_index, n_pairs, sep=" : ")]
pairs[, threshold := as.character(threshold)]

# Standard Errors for the average can be found via regression
lm1 <- lm((invested - 1) ~ -1 + buy_signal : group, data=benchmarks[threshold=="2"])
W <- sandwich::NeweyWest(lm1, lag=10, prewhite=FALSE)
lmtest::coeftest(lm1, df=Inf, vcov. = W)


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

dtfRbs <- benchmarks[threshold=="2", distribution_stats(base_line_long - 1), by=c("buy_signal", "group")]
dtfRbs[, metric := "baseline_long"]
dtfRcom <- benchmarks[threshold=="2", distribution_stats(committed - 1), by=c("buy_signal", "group")]
dtfRcom[, metric := "committed"]
dtfRinv <- benchmarks[threshold=="2", distribution_stats(invested - 1), by=c("buy_signal", "group")]
dtfRinv[, metric := "invested"]
dtfR <- rbindlist(list(dtfRbs, dtfRcom, dtfRinv))
dtfR[order(group)][metric == "baseline_long"]
dtfR[order(group)][metric == "baseline_long", round(.SD, 4), .SDcols=3:9]

# -------------------------------------------------------------


# Plot annual returns by group and buy signal
# -------------------------------------------------------------

# Annual returns on invested capital by buy-signal
p1 <- ggplot(data=benchmarks[group=="Top 5"],
            aes(x=trade_year, y=committed, fill=buy_signal)) +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) + 
  scale_fill_grey() +
  ggtitle("Top 5 Pairs") +
  ylab("Return on Committed Capital") +
  xlab("")


# Annual returns on committed capital by buy-signal
p2 <- ggplot(data=benchmarks[group=="Top 20"],
            aes(x=trade_year, y=committed, fill=buy_signal)) +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) + 
  scale_fill_grey() +
  ggtitle("Top 20 Pairs") +
  ylab("Return on Committed Capital") +
  xlab("")


# Annual returns on invested capital by buy-signal
p3 <- ggplot(data=benchmarks[group=="101-120"],
             aes(x=trade_year, y=committed, fill=buy_signal)) +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  scale_fill_grey() +
  ggtitle("Pairs 101-120") +
  ylab("Return on Committed Capital") +
  xlab("Trade Year")


ggpubr::ggarrange(p1, p2, p3, nrow=3, ncol=1)

# -------------------------------------------------------------


# Summarize pair formation between buy-signals
# -------------------------------------------------------------
dtfSP <- fread("data/SandP_companies.csv")
revenue_streams <- merge(revenue_streams, dtfSP[, .(ticker, sector)], by.x="long_pos", by.y="ticker")
setnames(revenue_streams, "sector", "long_sector")
revenue_streams <- merge(revenue_streams, dtfSP[, .(ticker, sector)], by.x="short_pos", by.y="ticker")
setnames(revenue_streams, "sector", "short_sector")
revenue_streams[, same_industry := long_sector == short_sector]
dcast(revenue_streams[group=="Top 5"], long_sector ~ short_sector, value.var="pair", fun.aggregate=length)

dtfPairs <- rbindlist(lapply(lst_ann_trade_ret, `[[`, 3))
dtfPairs[, c("pair_1", "pair_2") := tstrsplit(pair, "-", fixed=TRUE)]
dtfPairs[, c("sector_1", "sector_2") := tstrsplit(sector_pair, "-", fixed=TRUE)]
dtfPairs[start_index==1 & n_pairs==5, group := "Top 5"]
dtfPairs[start_index==1 & n_pairs==20, group := "Top 20"]
dtfPairs[start_index==101 & n_pairs==20, group := "101-120"]

dtfPairs <- merge(dtfPairs, dtfSP[, .(ticker, sector)], by.x="pair1", by.y="ticker")
setnames(dtfPairs, "sector", "sector_1")
dtfPairs <- merge(dtfPairs, dtfSP[, .(ticker, sector)], by.x="pair2", by.y="ticker")
setnames(dtfPairs, "sector", "sector_2")


dcast(dtfPairs[group=="Top 20" & buy_signal == "dtw"],
      sector_1 ~ sector_2,
      value.var="pair",
      fun.aggregate=length)[, 2:12]
M <- dcast(dtfPairs[group=="Top 20" & buy_signal == "dtw"],
      sector_1 ~ sector_2,
      value.var="pair",
      fun.aggregate=length)[, 2:12]
M <- as.matrix(M)
M <- (M + t(M)) - diag(diag(M))
M[upper.tri(M)] <- 0
dtfPairs[, .N, by=c("sector_pair", "buy_signal", "group")]

# -------------------------------------------------------------
