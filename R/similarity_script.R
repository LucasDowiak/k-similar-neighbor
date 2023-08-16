setwd("~/Git/k-similar-neighbor")
library(data.table)
source("R/similarity_functions.R")
source("R/auto_marginal.R")
source("R/plots.R")

# Standard scale the flows
dtf <- fread("data/SandP_tick_history.csv")
dtfTmp <- dtf[Date >= "2019-01-01" & Date <= "2020-12-31"]
pt1 <- dtfTmp[, .(Date)]
base_date <- pt1[, which(Date == "2020-01-02")]
pt2 <- Filter(is.numeric, dtfTmp)
pt2 <- pt2[, lapply(.SD, function(x) x / x[base_date])]
dtfTmp2 <- cbind(pt1, pt2)
fwrite(dtfTmp2, "data/std_scale_2019_2020.csv")



dticks <- c('BA', 'GD', 'HWM', 'HII', 'LHX', 'LMT', 'NOC', 'RTX', 'TDY', 'TXT', 'TDG')


# Marginal specifications to standardized residuals
tmp <- specs_to_resid_matrix("data/marginal_specifications_20200323_20210219.json",
                             write_to_file=FALSE)


# Pearson Correlation ---------------------------------------------------------- 

# Data Spec
dtf <- fread("data/std_resids_19_20.csv")
dtf1920 <- Filter(function(x) sum(is.na(x)) == 0, dtf)
tickers <- setdiff(names(dtf1920), "Date")

file_nm <- "data/association_results/corr_2019_2020"
res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, correlation, ref_file=file_nm, overwrite=TRUE, ncores=3L)

# Correlation Matrix
plot_sim_matrix("data/association_results/corr_2007_2008", type="cor", tits="Correlation: 2007 to 2008")
plot_sim_matrix("data/association_results/corr_2015_2016", type="cor", tits="Correlation: 2015 to 2016")
plot_sim_matrix("data/association_results/corr_2019_2020", type="cor", tits="Correlation: 2019 to 2020")


# Kullback-Leibler Divergence --------------------------------------------------

# Data Spec
yr <- 2016
dtfModelCoef <- fread(sprintf("data/model_coef_%d.csv", yr))
dtfModelSums <- fread(sprintf("data/model_summary_%d.csv", yr))
file_nm <- sprintf("data/association_results/kld_%d", yr)
if (!"lambda" %in% dtfModelCoef)
  dtfModelCoef[, lambda := NA_real_]

dtfModelStats <- merge(dtfModelCoef[, .SD, .SDcols=c("tick", "shape", "skew", "lambda")],
                       dtfModelSums[, .SD, .SDcols=c("tick", "distribution")],
                       on="tick", how="outer")

Qs <- lapply(dtfModelStats$tick, create_ddist, DT=dtfModelStats)
names(Qs) <- dtfModelStats$tick

M_spaceholder = t(dtfModelStats)
tickers <- M_spaceholder[1,]
M_spaceholder  <- data.table(M_spaceholder[2:5, ])
names(M_spaceholder) <- tickers

res <- matrix(NA_real_, nrow=dtfModelStats[, .N], ncol=dtfModelStats[, .N],
              dimnames=list(dtfModelStats$tick, dtfModelStats$tick))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(M_spaceholder, kullback_leibler, ref_file=file_nm,
                            overwrite=T, ncores=3L, df_model_distr=dtfModelStats)


# Dynamic Time Warping ---------------------------------------------------------

# Data Spec
dtf <- fread("data/std_scale_2019_2020.csv")
dtf1920 <- Filter(function(x) sum(is.na(x)) == 0, dtf)
tickers <- setdiff(names(dtf1920), "Date")
dtf1920 <- dtf[, tickers, with=FALSE]

file_nm <- "data/association_results/dtw_2019_2020"

res <- matrix(NA_real_, nrow=length(tickers), ncol=length(tickers),
              dimnames=list(tickers, tickers))
write.table(res, file=file_nm, sep="\t", col.names=NA)

tst <- calculate_similarity(dtf1920, dynamic_time_warp, ref_file=file_nm, overwrite=TRUE, ncores=3L)


# DTW Matrix
plot_sim_matrix("data/association_results/dtw_20190201_20210219", type="dtw", tits="DTW: 2019-02-01 to 2021-02-19")
plot_sim_matrix("data/association_results/dtw_20190201_20200214", type="dtw", tits="DTW: 2019-02-01 to 2020-02-14")
plot_sim_matrix("data/association_results/dtw_20200323_20210219", type="dtw", tits="DTW: 2020-03-23 to 2021-02-19")


# Plot Evolution of Metric over Periods ---------------------------------------------
rm(list=ls())
dtfSNP <- fread("data/SandP_companies.csv")
bpalette <- c('#c62828','#f44336','#9c27b0','#673ab7','#3f51b5','#2196f3','#29b6f6','#006064','#009688','#4caf50','#8bc34a','#BEC7C7')
names(bpalette) <- dtfSNP[, c(unique(sector), "ADefault")]

read_and_melt <- function(file, column_name=NULL, unique_pairs=TRUE)
{
  X <- fread(file)
  ticks <- X$V1
  stopifnot(identical(ticks, names(X)[-1]))
  XX <- as.matrix(X[, -c("V1")])
  if (unique_pairs) {
    XX[lower.tri(XX, diag=TRUE)] <- NA_real_
  }
  out <- data.table(XX)
  out[, TICK1 := ticks]
  out <- melt(out, id.vars="TICK1", variable.factor=FALSE)
  setnames(out, "variable", "TICK2")
  if (!is.null(column_name)) {
    setnames(out, "value", column_name)
  }
  return(out)
}

patt <- "dtw_[0-9]{4}_[0-9]{4}|corr_[0-9]{4}_[0-9]{4}|kld_[0-9]{4}_[0-9]{4}"
files_ <- list.files("data/association_results/", pattern=patt)
files_ <- files_[!files_ %in% c("kld_2007", "kld_2015", "kld_2019")]

lst_dtf <- lapply(files_, function(x) read_and_melt(sprintf("data/association_results/%s", x)))
names(lst_dtf) <- files_
lapply(files_, function(x) lst_dtf[[x]][, vtype := x])
DT <- rbindlist(lst_dtf)
DT <- DT[!is.na(value)]
DT <- merge(DT, dtfSNP[, .SD, .SDcols=c("ticker", "sector", "subindustry")],
                       by.x="TICK1", by.y="ticker", all.x=TRUE)
DT <- merge(DT, dtfSNP[, .SD, .SDcols=c("ticker", "sector", "subindustry")],
                       by.x="TICK2", by.y="ticker", all.x=TRUE, suffixes = c("_t1", "_t2"))

mtype = "kld" # "corr", "kld", "dtw"
for (industry in dtfSNP[, unique(sector)]) {
  set.seed(522193)
  dtDTWYoy <- copy(DT[complete.cases(DT)])[sample.int(.N, .N * .3)]
  dtDTWYoy <- dtDTWYoy[grepl(mtype, vtype)]
  dtDTWYoy[, sect_colr := ifelse(sector_t1==sector_t2 & sector_t1==industry, industry, "ADefault")]
  dtDTWYoy <- dtDTWYoy[!is.na(sect_colr)][order(sect_colr)]
  dtDTWYoy[vtype %in% c("kld_2008", "kld_2016"),
           p1_pairs := paste(TICK1, TICK2, sep="-")]
  dtDTWYoy[vtype %in% c("kld_2016", "kld_2020"),
           p2_pairs := paste(TICK1, TICK2, sep="-")]
  
  
  p <- ggplot(dtDTWYoy, aes(x=factor(vtype), y=value)) +
    geom_line(data=dtDTWYoy[sect_colr != industry], aes(group=factor(p1_pairs)), color=bpalette['ADefault'], alpha=0.1) +
    geom_line(data=dtDTWYoy[sect_colr != industry], aes(group=factor(p2_pairs)), color=bpalette['ADefault'], alpha=0.1) +
    geom_line(data=dtDTWYoy[sect_colr == industry], aes(group=factor(p1_pairs)), color=bpalette[industry], alpha=0.5) +
    geom_line(data=dtDTWYoy[sect_colr == industry], aes(group=factor(p2_pairs)), color=bpalette[industry], alpha=0.5) +
    geom_point(data=dtDTWYoy[sect_colr != industry], alpha=0.1, color=bpalette['ADefault']) +
    geom_point(data=dtDTWYoy[sect_colr == industry], alpha=0.5, color=bpalette[industry]) +
    labs(y="Metric Value", x=NULL, title=industry)
  print(p)
  ggsave(sprintf("article/images/%s_%s.png", mtype, industry))
}



# Plot DTW distance vs Correlation ---------------------------------------------
files_ <- list.files("data/association_results/", pattern="dtw|corr|kld")
files_ <- files_[!files_ %in% c("kld_2007", "kld_2015", "kld_2019")]

lst_dtf <- lapply(files_, function(x) read_and_melt(sprintf("data/association_results/%s", x), unique_pairs = FALSE))
lst_dtf
names(lst_dtf) <- files_
lapply(files_, function(x) lst_dtf[[x]][, vtype := x])
DT <- rbindlist(lst_dtf)

DT <- merge(DT, dtfSNP[, .SD, .SDcols=c("ticker", "sector", "subindustry")],
            by.x="TICK1", by.y="ticker", all.x=TRUE)
DT <- merge(DT, dtfSNP[, .SD, .SDcols=c("ticker", "sector", "subindustry")],
            by.x="TICK2", by.y="ticker", all.x=TRUE, suffixes = c("_t1", "_t2"))
DT[, same_sector := ifelse(sector_t1 == sector_t2, sector_t1, 'ADefault')]

DTscatter <- dcast(DT[grepl("2016", vtype)], TICK1 + TICK2 + sector_t1 + same_sector ~ vtype, value.var = "value")
DTscatter <- DTscatter[order(same_sector)]
p <- ggplot(DTscatter, aes(x=corr_2015_2016, y=dtw_2015_2016)) +
  geom_point(size=0.5, alpha=0.2, color=bpalette[match(DTscatter$same_sector, names(bpalette))]) +
  facet_wrap(vars(sector_t1), nrow=3, ncol=4) +
  geom_vline(xintercept = 0, lty=3) + 
  labs(x="Corr", y="DTW", title="Corr-vs-DTW")

print(p)

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
    scale_color_manual(values=c("gray"="Group", "red"="Group Avg")) +
    scale_alpha_identity() + 
    facet_wrap(~ grps)
  print(p)
  
  
  
}

dtf <- fread("data/std_scale_2015_2016.csv")
dtf <- Filter(function(x) all(!is.na(x)), dtf)
X <- read.table("data/association_results/dtw_2015_2016", header=T, row.names=1)
ii <- which(names(X) == "ENPH")
X <- X[-ii, -ii]

tree <- cluster_similarity(X, type="dtw", visualize = TRUE)

groups <- cutree(tree, k=9)
table(groups)
calculate_group_stats(groups, dtf)
