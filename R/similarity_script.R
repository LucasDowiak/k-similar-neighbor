setwd("~/Git/k-similar-neighbor")
library(data.table)
source("R/similarity_functions.R")
source("R/auto_marginal.R")
source("R/plots.R")

# Standard scale the flows
dtf <- fread("data/SandP_tick_history.csv")
dtfTmp <- dtf[Date >= "2019-02-01" & Date <= "2021-02-19"]
pt1 <- dtfTmp[, .(Date)]
pt2 <- Filter(is.numeric, dtfTmp)
pt2 <- pt2[, lapply(.SD, function(x) x / x[1])]
dtfTmp2 <- cbind(pt1, pt2)
fwrite(dtfTmp2, "data/std_scale_20190201_20210219.csv")



dticks <- c('BA', 'GD', 'HWM', 'HII', 'LHX', 'LMT', 'NOC', 'RTX', 'TDY', 'TXT', 'TDG')


# Marginal specifications to standardized residuals
tmp <- specs_to_resid_matrix("data/marginal_specifications_20200323_20210219.json",
                             write_to_file=FALSE)


# Pearson Correlation ---------------------------------------------------------- 

# Data Spec
dtf <- fread("data/std_resids_20200323_20210219.csv")
dtf1920 <- Filter(function(x) sum(is.na(x)) == 0, dtf)
tickers <- setdiff(names(dtf1920), "Date")

file_nm <- "data/association_results/corr_20200323_20210219"
res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, correlation, ref_file=file_nm, overwrite=TRUE, ncores=4L)

# Correlation Matrix
plot_sim_matrix("data/association_results/corr_20190201_20210219", type="cor", tits="Correlation: 2019-02-01 to 2021-02-19")
plot_sim_matrix("data/association_results/corr_20190201_20200214", type="cor", tits="Correlation: 2019-02-01 to 2020-02-14")
plot_sim_matrix("data/association_results/corr_20200323_20210219", type="cor", tits="Correlation: 2020-03-23 to 2021-02-19")


# Kullback-Leibler Divergence --------------------------------------------------

# Data Spec
dtfModelCoef <- tmp$model_coef # "data/model_coef_20200323_20210219.csv"
dtfModelSums <- tmp$model_summary # "data/model_coef_20200323_20210219.csv"
if (!"lambda" %in% dtfModelCoef)
  dtfModelCoef[, lambda := NA_real_]

dtfModelStats <- merge(dtfModelCoef[, .SD, .SDcols=c("tick", "shape", "skew", "lambda")],
                       dtfModelSums[, .SD, .SDcols=c("tick", "distribution")],
                       on="tick", how="outer")


# Distribution factory
create_ddist <- function(ticker, modelstats)
{
  na_to_null <- function(x)
  {
    if (is.na(x) || length(x) == 0)
      return(NULL)
    else
      return(x)
  }
  shape <- na_to_null(modelstats[tick==ticker, shape])
  skew <- na_to_null(modelstats[tick==ticker, skew])
  lambda <- na_to_null(modelstats[tick==ticker, lambda])
  distr <- na_to_null(modelstats[tick==ticker, distribution])
  f <- function(x)
  {
    ddist(distribution=distr, x, lambda=lambda, skew=skew, shape=shape)
  }
  return(f)
}

# Create matrix of theoretical density values for each stocks probability distribution
support <- seq(-5, 5, by=0.01)
Qs <- lapply(dtfModelStats$tick, create_ddist, dtfModelStats)
Qs <- as.data.table(lapply(Qs, function(Q) Q(support)))
setnames(Qs, dtfModelStats$tick)

debugonce(kullback_leibler)

tst <- kullback_leibler("A",
                        c("EXR", "MS", "AAL", "CMG", "UAL", "FRC", "UA", "KO", "HIG", "MHK"),
                        Qs)

file_nm <- "data/association_results/kld_20200323_20210219"
res <- matrix(NA_real_, nrow=dtfModelStats[, .N], ncol=dtfModelStats[, .N],
              dimnames=list(dtfModelStats$tick, dtfModelStats$tick))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(Qs, kullback_leibler, ref_file=file_nm,
                            overwrite=FALSE, ncores=1L, epsilon=1e-4)


# Dynamic Time Warping ---------------------------------------------------------

# Data Spec
dtf <- fread("data/std_scale_20200323_20210219.csv")

dtf1920 <- dtf[, tickers, with=FALSE]

file_nm <- "data/association_results/dtw_20200323_20210219"

res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, dynamic_time_warp, ref_file=file_nm, overwrite=TRUE, ncores=4L)


# DTW Matrix
plot_sim_matrix("data/association_results/dtw_20190201_20210219", type="dtw", tits="DTW: 2019-02-01 to 2021-02-19")
plot_sim_matrix("data/association_results/dtw_20190201_20200214", type="dtw", tits="DTW: 2019-02-01 to 2020-02-14")
plot_sim_matrix("data/association_results/dtw_20200323_20210219", type="dtw", tits="DTW: 2020-03-23 to 2021-02-19")


# Plot DTW distance vs Correlation ---------------------------------------------

read_and_melt <- function(file, column_name=NULL)
{
  X <- fread(file)
  ticks <- X$V1
  stopifnot(identical(ticks, names(X)[-1]))
  XX <- as.matrix(X[, -c("V1")])
  XX[lower.tri(XX, diag=TRUE)] <- NA_real_
  out <- data.table(XX)
  out[, TICK1 := ticks]
  out <- melt(out, id.vars="TICK1", variable.factor=FALSE)
  setnames(out, "variable", "TICK2")
  if (!is.null(column_name)) {
    setnames(out, "value", column_name)
  }
  return(out)
}

RHO_all <- read_and_melt("data/association_results/corr_20190201_20210219", "RHO_all")
RHO_pre <- read_and_melt("data/association_results/corr_20190201_20200214", "RHO_pre")
RHO_post <- read_and_melt("data/association_results/corr_20200323_20210219", "RHO_post")
RHO <- merge.data.table(RHO_post, RHO_pre, by=c("TICK1", "TICK2"), all=TRUE)
RHO <- merge.data.table(RHO, RHO_all, by=c("TICK1", "TICK2"), all=TRUE)

DTW_all <- read_and_melt("data/association_results/dtw_20190201_20210219", "DTW_all")
DTW_pre <- read_and_melt("data/association_results/dtw_20190201_20200214", "DTW_pre")
DTW_post <- read_and_melt("data/association_results/dtw_20200323_20210219", "DTW_post")
DTW <- merge.data.table(DTW_post, DTW_pre, by=c("TICK1", "TICK2"), all=TRUE)
DTW <- merge.data.table(DTW, DTW_all, by=c("TICK1", "TICK2"), all=TRUE)

DT <- merge(DTW, RHO, by=c("TICK1", "TICK2"), all=TRUE)
DT[, DTW_diff := DTW_post - DTW_pre]
DT[, RHO_diff := RHO_post - RHO_pre]
DT2 <- copy(DT[complete.cases(DT)])

ggplot(DT2, aes(x=DTW_all, y=RHO_all)) +
  geom_point(alpha=0.1)

DT2[, plot(DTW_pre, RHO_pre, cex=0.25, main="Pre-Covid", xlim=c(-10, 800), ylim=c(-1,1))]
DT2[, plot(DTW_post, RHO_post, cex=0.25, main="Post-Covid", xlim=c(-10, 800), ylim=c(-1,1))]

# Similar shape but negatively correlated
tmp <- DT2[DTW_pre < 50 & RHO_pre < -0.25][order(TICK1, TICK2)]

# No correlation but various degrees 
tmp <- DT2[RHO_pre < 0.1 & RHO_pre > -0.1]



# Cluster similarity matrices --------------------------------------------------

# Given a cluster of stock ticks, run a batch of summary statistics
#
# grps - named vector of groupings (possible from call to cutree)
#
# DT - data.table of standardized stock prices
#
calculate_group_stats <- function(grps, DT)
{
  DT <- copy(DT)
  dtf_grp <- data.table(grps, ticker=names(grps))
  DT <- merge(melt(DT, id.vars="Date", variable.name="ticker"), dtf_grp, on="ticker")
  DT_avg <- DT[, mean(value), by=c("Date", "grps")]
  DT_avg[, ticker := "zzz_avg"]
  setnames(DT_avg, "V1", "value")
  DT <- rbindlist(list(DT, DT_avg), use.names=TRUE)
  DT[, `:=`(color="gray", alpha=0.5)]
  DT[ticker=="zzz_avg", `:=`(color="red", alpha=1)]
  
  dtfSandP <- read_SandP_data()
  DT <- merge(DT, dtfSandP[, .SD,
                           .SDcols=c("ticker", "sector", "subindustry", "security")],
              on="ticker")
  
  print(table(grps))
  DT[!duplicated(DT$ticker), print(table(grps, sector))]
  
  p <- ggplot(DT, aes(x=Date, y=value, group=ticker)) +
    geom_line(aes(color=DT$color, alpha=DT$alpha)) +
    scale_color_manual(values=c("gray"="gray", "red"="red")) +
    scale_alpha_identity() + 
    facet_wrap(~ grps)
  print(p)
  
  
}

dtf <- fread("data/std_scale_20190201_20210219.csv")
dtf <- Filter(function(x) all(!is.na(x)), dtf)
X <- read.table("data/association_results/corr_20190201_20210219", header=T, row.names=1)
ii <- which(names(X) == "ENPH")
X <- X[-ii, -ii]

tree <- cluster_similarity(X, type="cor", visualize = FALSE)

groups <- cutree(tree, k=9)
table(groups)
calculate_group_stats(groups, dtf)
