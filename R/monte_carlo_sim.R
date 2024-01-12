setwd("~/Git/dtw-in-finance/")
source("R/normality_tests.R")
library(data.table); library(rugarch)
mod1 <- readRDS("data/model_objects/DA_2007_2008/DA_2007_2008_AAPL.rds")
mod2 <- readRDS("data/model_objects/DA_2007_2008/DA_2007_2008_GE.rds")

files_ <- list.files("data/label_analysis/", pattern = "model_summary", full.names=TRUE)
dtfModSum <- rbindlist(lapply(files_, fread))
names(lstModSum) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])


files_ <- list.files("data/association_results/", pattern = "model_cor", full.names=TRUE)
lstCors <- lapply(files_, read_association_table)
names(lstCors) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])

file_ <- list.files("data/association_results/", pattern="dtw", full.names=TRUE)
lstDTW <- lapply(files_, read_association_table)
names(lstDTW) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])

melt_unique_assc_matrix <- function(nm, lst) {
  M <- lst[[nm]]
  M[upper.tri(M, diag=TRUE)] <- NA_real_
  M <- as.data.table(reshape2::melt(M))[!is.na(value)]
  setnames(M, c("Var1", "Var2"), c("tick1", "tick2"))
  M <- M[, `:=`(label=nm)]
  M
}

sample_cor_bin <- function(DT, size=10)
{
  if (DT[, .N <= size]) {
    return(DT)
  } else {
    return(DT[sample.int(.N, size)])
  }
}

dtfDTW <- rbindlist(lapply(names(lstDTW), melt_unique_assc_matrix, lst=lstDTW))
dtfDTW[, year := substring(label, 1, 4)]
dtfDTW[, label := NULL]
setnames(dtfDTW, "value", "dtw")

dtfCors <- rbindlist(lapply(names(lstCors), melt_unique_assc_matrix, lst=lstCors))
setnames(dtfCors, "value", "model_cor")
dtfCors[, year := substring(label, 1, 4)]
dtfCors[, label := NULL]
dtfCors[, summary(value)]
dtfCors[, hist(value, breaks=40)]
dtfCors[, cor_bin := cut(value, dtfCors[, hist(value, breaks=40, plot=FALSE)]$breaks)]
dtfSmp1 <- dtfCors[, sample_cor_bin(.SD, 300), by=cor_bin][order(cor_bin)]


# Extract fitted distribution from uGARCH-fit
sample_distribution <- function(obj, u, to_std_price=TRUE)
{
  if (missing(u)) {
    u <- runif(obj@model$modeldata$T)
  }
  dpars <- extract_ugarchfit_dist(obj)
  smp <- qdist(distribution=dpars$dentype, 
               p=u,
               mu=obj@fit$fitted.values,
               sigma=obj@fit$sigma,
               lambda=dpars$lambda,
               skew=dpars$skew,
               shape=dpars$shape)
  if (to_std_price) {
    smp <- cumprod(1 + smp)
  }
  return(smp)
}


resample <- function(m1, m2, N, rho, tdf)
{
  U <- BiCopSim(N, family=2, par=rho, par2=tdf)
  out <- matrix(NA_real_, nrow=N, ncol=2)
  out[, 1] <- sample_distribution(m1, u=U[, 1])
  out[, 2] <- sample_distribution(m2, u=U[, 2])
  return(out)
}


replicate_pair <- function(pair, B=100)
{
  m1 <- readRDS(pair[1])$model; m2 <- readRDS(pair[2])$model
  N1 <- m1@model$modeldata$T
  stopifnot(N1 == m2@model$modeldata$T)
  
  parfit <- BiCopEst(pit(m1), pit(m2), family=2)
  samples <- replicate(B, resample(m1, m2, N=N1, rho=parfit$par, tdf=parfit$par2), simplify=FALSE)
  dtws <- lapply(samples, calculate_similarity, dist.method="DTW", window.type=sakoeChibaWindow,
                 window.size=15)
  dtws <- sapply(dtws, `[`, 2)
  out <- data.table(t(c(summary(dtws), quantile(dtws, c(0.05, 0.95)))))
  return(out)
}




# smp1_pair_model_locations
collect <- vector("list", nrow(dtfSmp1))
for (i in 1:nrow(dtfSmp1)) {
  if (i %% 10 == 0) {
    print(sprintf("Starting %d out of %d at %s.", i, dtfSmp1[, .N], Sys.time()))
  }
  t1 <- dtfSmp1[i, tick1]
  t2 <- dtfSmp1[i, tick2]
  yr_lb <- dtfSmp1[i, strsplit(label, "_")[[1]][1]]
  m1_loc <- sprintf("data/model_objects/%s/%s_%s.rds", yr_lb, yr_lb, t1)
  m2_loc <- sprintf("data/model_objects/%s/%s_%s.rds", yr_lb, yr_lb, t2)
  mcsims <- replicate_pair(c(m1_loc, m2_loc), B=100)
  mcsims[, `:=`(tick1=t1, tick2=t2, year=yr_lb)]
  collect[[i]] <- mcsims
}
dtf <- merge(dtfSmp1, rbindlist(collect))
dtf <- merge(dtf, dtfDTW)
dtf[, year := factor(year, levels=sort(unique(year)))]
dtf[, IQR := `3rd Qu.` - `1st Qu.`]
dtf[, violation := dtw < `5%` | dtw > `95%`]
dtf[, mean(violation)] # It should be close to 10%
# TODO: 1. pull in unconditional mean and variance, industry, year variables into the data frame
#       2. regress against 

library(ggplot2)

ggplot(dtf, aes(x=cor_bin, y=IQR)) + 
  geom_boxplot()

ggplot(dtf, aes(x=year, y=value)) + 
  geom_boxplot()

ggplot(dtf, aes(model_cor, dtw)) + 
  geom_bin2d(binwidth=c(0.05, 20))

ggplot(dtf, aes(model_cor, log(dtw))) + 
  geom_bin2d(binwidth=c(0.05, 0.25)) +
  scale_fill_viridis_c() +
  geom_smooth(method="lm", formula=y ~ poly(x, 2), col="red")

summary(lm(log(dtw) ~ poly(model_cor, 2), data=dtf))

dtf[cor_bin=="(0.75,0.80]"][order(Mean)]


statSigDiffPairs <- which(dtf[, outer(`5%`, `95%`, function(x, y) x > y)], arr.ind=TRUE)

collect2 <- vector("list", nrow(statSigDiffPairs))
for (rr in seq_len(nrow(statSigDiffPairs))) {
  tmp1 <- dtf[statSigDiffPairs[rr,]]
  tmp2 <- data.table(p1_t1=tmp1[1, tick1],
                     p1_t2=tmp1[1, tick2],
                     p1_cor_bin=tmp1[1, cor_bin],
                     p2_t1=tmp1[2, tick1],
                     p2_t2=tmp1[2, tick2],
                     p2_cor_bin=tmp1[2, cor_bin])
  collect2[[rr]] <- tmp2
}

dtf2 <- rbindlist(collect2)


