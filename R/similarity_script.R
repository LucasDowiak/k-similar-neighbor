setwd("~/Git/k-similar-neighbor")
library(data.table)
source("R/similarity_functions.R")
source("R/auto_marginal.R")

# Standard scale the flows
dtf <- fread("data/SandP_tick_history.csv")
dtfTmp <- dtf[Date >= "2020-03-23" & Date <= "2021-02-19"]
pt1 <- dtfTmp[, .(Date)]
pt2 <- Filter(function(x) is(x, "numeric"), dtfTmp)
pt2 <- pt2[, lapply(.SD, function(x) x / x[1])]
dtfTmp2 <- cbind(pt1, pt2)
fwrite(dtfTmp2, "data/std_scale_20200323_20210219.csv")



dticks <- c('BA', 'GD', 'HWM', 'HII', 'LHX', 'LMT', 'NOC', 'RTX', 'TDY', 'TXT', 'TDG')


# Marginal specifications to standardized residuals
tmp <- specs_to_resid_matrix("data/marginal_specifications_20200323_20210219.json")


# Pearson Correlation ---------------------------------------------------------- 

# Data Spec
dtf <- fread("data/std_resids_20200323_20210219.csv")
dtf1920 <- Filter(function(x) sum(is.na(x)) == 0, dtf)
tickers <- setdiff(names(dtf1920), "Date")

file_nm <- "data/association_results/corr_20200323_20210219_485"
res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, correlation, ref_file=file_nm, overwrite=TRUE, ncores=4L)

# Correlation Matrix
plot_sim_matrix("data/association_results/corr_20190201_20200214_483", type="cor", tits="Correlation: 2019-02-01 to 2020-02-14")
plot_sim_matrix("data/association_results/corr_20200323_20210219_485", type="cor", tits="Correlation: 2020-03-23 to 2021-02-19")


# Dynamic Time Warping ---------------------------------------------------------

# Data Spec
dtf <- fread("data/std_scale_20200323_20210219.csv")
dtf1920 <- dtf[, tickers, with=FALSE]

file_nm <- "data/association_results/dtw_20200323_20210219_486"

res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, dynamic_time_warp, ref_file=file_nm, overwrite=TRUE, ncores=4L)


# DTW Matrix
plot_sim_matrix("data/association_results/dtw_20190201_20200214_483", type="dtw", tits="DTW: 2019-02-01 to 2020-02-14")
plot_sim_matrix("data/association_results/dtw_20200323_20210219_486", type="dtw", tits="DTW: 2020-03-23 to 2021-02-19")


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

RHO_pre <- read_and_melt("data/association_results/corr_20190201_20200214_483", "RHO_pre")
RHO_post <- read_and_melt("data/association_results/corr_20200323_20210219_485", "RHO_post")
RHO <- merge.data.table(RHO_pre, RHO_post, by=c("TICK1", "TICK2"))

DTW_pre <- read_and_melt("data/association_results/dtw_20190201_20200214_483", "DTW_pre")
DTW_post <- read_and_melt("data/association_results/dtw_20200323_20210219_486", "DTW_post")
DTW <- merge.data.table(DTW_pre, DTW_post, by=c("TICK1", "TICK2"))

DT <- merge(DTW, RHO, by=c("TICK1", "TICK2"))
DT[, DTW_diff := DTW_post - DTW_pre]
DT[, RHO_diff := RHO_post - RHO_pre]
DT2 <- copy(DT[complete.cases(DT)])

ggplot(DT2, aes(x=DTW_pre, y=RHO_pre)) +
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

dtf <- fread("data/std_scale_20200323_20210219.csv")
dtf <- Filter(function(x) all(!is.na(x)), dtf)
X <- read.table("data/association_results/dtw_20200323_20210219_486", header=T, row.names=1)
tree <- cluster_similarity(X, type="dtw", visualize = TRUE)

groups <- cutree(tree, k=8)
table(groups)
calculate_group_stats(groups, dtf)
