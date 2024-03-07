setwd("~/Git/dtw-in-finance/")
source("R/normality_tests.R")
source("R/similarity_functions.R")
library(data.table); library(rugarch); library(VineCopula)

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


files_ <- list.files("data/label_analysis/", pattern = "model_summary", full.names=TRUE)
dtfModSum <- rbindlist(lapply(files_, fread))# , colClasses=c("label"="character")))
dtfModSum[, time_varying_mean := as.numeric(ar > 0 | ma > 0)]
setnames(dtfModSum, "label", "year")


files_ <- list.files("data/association_results/", pattern = "model_cor", full.names=TRUE)
lstMCors <- lapply(files_, read_association_table)
names(lstMCors) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])
dtfMCors <- rbindlist(lapply(names(lstMCors), melt_unique_assc_matrix, lst=lstMCors))
setnames(dtfMCors, "value", "model_cor")
dtfMCors[, year := substring(label, 1, 4)]
dtfMCors[, label := NULL]
dtfMCors[, mcor_bin := cut(model_cor, dtfMCors[, hist(model_cor, breaks=40, plot=FALSE)]$breaks)]


files_ <- list.files("data/association_results/", pattern = "unadjusted_cor", full.names=TRUE)
lstUCors <- lapply(files_, read_association_table)
names(lstUCors) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])
dtfUCors <- rbindlist(lapply(names(lstUCors), melt_unique_assc_matrix, lst=lstUCors))
setnames(dtfUCors, "value", "unadjusted_cor")
dtfUCors[, year := substring(label, 1, 4)]
dtfUCors[, label := NULL]
dtfUCors[, ucor_bin := cut(unadjusted_cor, dtfUCors[, hist(unadjusted_cor, breaks=40, plot=FALSE)]$breaks)]


files_ <- list.files("data/association_results/", pattern="dtw", full.names=TRUE)
lstDTW <- lapply(files_, read_association_table)
names(lstDTW) <- sapply(files_, function(x) strsplit(x, "//")[[1]][2])
dtfDTW <- rbindlist(lapply(names(lstDTW), melt_unique_assc_matrix, lst=lstDTW))
dtfDTW[, year := substring(label, 1, 4)]
dtfDTW[, label := NULL]
setnames(dtfDTW, "value", "dtw")

dtfDist <- Reduce(merge, list(dtfMCors, dtfDTW), dtfUCors)
dtfDist[, year := as.integer(year)]
dtfDist[, tick1 := as.character(tick1)]
dtfDist[, tick2 := as.character(tick2)]

# dtfSmp1 <- dtfCors[, sample_cor_bin(.SD, 300), by=cor_bin][order(cor_bin)]
# dtfSmp1 <- merge(dtfSmp1, dtfDTW)
# fwrite(dtfSmp1, file="data/monte_carlo/sample_20240118_1.csv")


# Extract fitted distribution from uGARCH-fit and perform monte-carlo sampling
# ------------------------------------------------
sample_distribution <- function(obj, u, to_std_price=TRUE, conditional_mu=TRUE,
                                conditional_sigma=TRUE)
{
  if (missing(u)) {
    u <- runif(obj@model$modeldata$T)
  }
  if (conditional_mu) {
    mu_ <- obj@fit$fitted.values
  } else {
    if (identical(conditional_mu, 0)) {
      mu_ <- 0
    } else {
      mu_ <- uncmean(obj)
    }
  }
  if (conditional_sigma) {
    sigma_ <- obj@fit$sigma
  } else {
    sigma_ <- sqrt(uncvariance(obj))
  }
  dpars <- extract_ugarchfit_dist(obj)
  smp <- qdist(distribution=dpars$dentype, 
               p=u,
               mu=mu_,
               sigma=sigma_,
               lambda=dpars$lambda,
               skew=dpars$skew,
               shape=dpars$shape)
  if (to_std_price) {
    smp <- cumprod(1 + smp)
  }
  return(smp)
}


resample <- function(m1, m2, N, rho, tdf, to_std_price=TRUE, conditional_mu=TRUE,
                     conditional_sigma=TRUE)
{
  U <- BiCopSim(N, family=2, par=rho, par2=tdf)
  out <- matrix(NA_real_, nrow=N, ncol=2)
  out[, 1] <- sample_distribution(m1, u=U[, 1], to_std_price=to_std_price,
                                  conditional_mu=conditional_mu, conditional_sigma=conditional_sigma)
  out[, 2] <- sample_distribution(m2, u=U[, 2], to_std_price=to_std_price,
                                  conditional_mu=conditional_mu, conditional_sigma=conditional_sigma)
  return(out)
}


replicate_pair <- function(pair, B=100, to_std_price=TRUE, conditional_mu=TRUE,
                           conditional_sigma=TRUE)
{
  m1 <- readRDS(pair[1])$model; m2 <- readRDS(pair[2])$model
  N1 <- m1@model$modeldata$T
  stopifnot(N1 == m2@model$modeldata$T)
  
  parfit <- BiCopEst(pit(m1), pit(m2), family=2)
  samples <- replicate(B, resample(m1, m2, N=N1, rho=parfit$par, tdf=parfit$par2, to_std_price=to_std_price,
                                   conditional_mu=conditional_mu, conditional_sigma=conditional_sigma),
                       simplify=FALSE)
  dtws <- lapply(samples, calculate_similarity, dist.method="DTW", window.type=sakoeChibaWindow,
                 window.size=15)
  dtws <- sapply(dtws, `[`, 2)
  # out <- data.table(t(c(summary(dtws), quantile(dtws, c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99)))))
  return(dtws)
}


replicate_study <- function(pair, B)
{
  m1_loc <- sprintf("data/model_objects/%s/%s_%s.rds", pair$year, pair$year, pair$tick1)
  m2_loc <- sprintf("data/model_objects/%s/%s_%s.rds", pair$year, pair$year, pair$tick2)
  mcsims <- replicate_pair(c(m1_loc, m2_loc), B=B, to_std_price=TRUE,
                           conditional_mu=TRUE, conditional_sigma=TRUE)
  mcsims2 <- replicate_pair(c(m1_loc, m2_loc), B=B, to_std_price=TRUE,
                            conditional_mu=FALSE, conditional_sigma=TRUE)
  mcsims3 <- replicate_pair(c(m1_loc, m2_loc), B=B, to_std_price=TRUE,
                            conditional_mu=TRUE, conditional_sigma=FALSE)
  tmpDT <- rbindlist(list(data.table(label="Cond_Mean_and_Var", value=mcsims),
                          data.table(label="Cond_Var", value=mcsims2),
                          data.table(label="Cond_Mean", value=mcsims3)))
  tmpDT[, `:=`(tick1=pair$tick1, tick2=pair$tick2, year=pair$year)]
  return(tmpDT)
}


BB <- 250
cluster <- parallel::makeCluster(3)
clusterEvalQ(cluster, library("data.table"))
clusterEvalQ(cluster, library("rugarch"))
clusterEvalQ(cluster, library("VineCopula"))
clusterEvalQ(cluster, source("R/normality_tests.R"))
clusterEvalQ(cluster, source("R/similarity_functions.R"))
clusterExport(cluster, c("sample_distribution", "resample", "replicate_pair", "replicate_study"))
batches <- split(dtfSmp1, 1:nrow(dtfSmp1) %% 400)
collect <- vector("list", length(batches))
for (b in seq_along(batches)) {
  lst_ <- split(batches[[b]], 1:nrow(batches[[b]]))
  print(sprintf("Batch %d of %d total started at %s", b, length(batches), Sys.time()))
  collect[[b]] <- rbindlist(parLapply(cl=cluster, lst_, replicate_study, B=BB))
  dtfSave <- rbindlist(collect)
  fwrite(dtfSave, file="data/monte_carlo/Sample_pairs_corr_to_dtw_20240118_1.csv")
  Sys.sleep(10)
}
stopCluster(cluster)

# dtfSmp1 <- fread("data/monte_carlo/sample_20240118_1.csv")
dtfSave <- fread("data/monte_carlo/Sample_pairs_corr_to_dtw_20240118_1.csv")
dtfSave2 <- fread("data/monte_carlo/Sample_pairs_corr_to_dtw_20240118_2.csv")
dtfSave2[, value := as.numeric(value)]
dtfSave3 <- rbindlist(list(dtfSave, dtfSave2))

dtfCI <- dtfSave3[, as.list(quantile(value, c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995))),
                 by=c("tick1", "tick2", "year", "label")]
dtfCI <- merge(dtfCI, dtfDist)


dtfCI <- merge(dtfCI, dtfModSum[, .SD, .SDcols=c("tick", "year", "time_varying_mean")],
                 by.x=c("tick1", "year"), by.y=c("tick", "year"))
setnames(dtfCI, "time_varying_mean", "t1_time_varying_mean")
dtfCI <- merge(dtfCI, dtfModSum[, .SD, .SDcols=c("tick", "year", "time_varying_mean")],
                 by.x=c("tick2", "year"), by.y=c("tick", "year"))
setnames(dtfCI, "time_varying_mean", "t2_time_varying_mean")

dtfCI[, tmv_present := t1_time_varying_mean + t2_time_varying_mean]
dtfCI[, viol_01 := (dtw < `0.5%` | dtw > `99.5%`)]
dtfCI[, viol_05 := (dtw < `2.5%` | dtw > `97.5%`)]
dtfCI[, viol_10 := (dtw < `5%` | dtw > `95%`)]

dtfCI[tmv_present > 0,
      as.list(round(c(colMeans(.SD), N=.N), 3)),
      .SDcols=c("viol_01", "viol_05", "viol_10"),
      by=c("label")][order(label)]

dtfCI[tmv_present == 0,
      as.list(round(c(colMeans(.SD), N=.N), 3)),
      .SDcols=c("viol_01", "viol_05", "viol_10"),
      by=c("label")][order(label)]

dtfCI[tmv_present == 0 & (label %in% c("Unconditional", "Cond_Mean")) ,
      as.list(round(c(colMeans(.SD), N=.N), 3)),
      .SDcols=c("viol_01", "viol_05", "viol_10")]

# ------------------------------------------------



# Plots
# ------------------------------------------------
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

# ------------------------------------------------
t1 <- "CTSH"
t2 <- "VRTX"
yr_lb <- "2012"
tblUnadj <- read_association_table(sprintf("data/association_results/%s_unadjusted_cor.tsv", yr_lb))
tblDTW <- read_association_table(sprintf("data/association_results/%s_dtw.tsv", yr_lb))
tblMod <- read_association_table(sprintf("data/association_results/%s_model_cor.tsv", yr_lb))

tblUnadj[t1, t2]
tblDTW[t1, t2]
tblMod[t1, t2]



dtf2[, plot(model_cor, log(dtw), main="model_cor", xlim=c(-0.6, 1.1))]
dtf2[, plot(rho_m, log(dtw), main="rho_m", xlim=c(-0.6, 1.1))]
# ------------------------------------------------

# Regression Stats: Unc Mean, Unc Variance, Persistance, Halflife
# ------------------------------------------------
dtf <- fread("data/Sample_pairs_corr_to_dtw.csv")
dtf[, year := as.character(year)]
collect3 <- vector("list", dtf[,.N])
for (i in 1:nrow(dtf)) {
  if (i %% 100 == 0) {
    print(sprintf("%d out of %d at %s", i, dtf[,.N], Sys.time()))
  }
  t1 <- dtf[i, tick1]; t2 <- dtf[i, tick2]; yr_lb <- dtf[i, year]
  t1_file_loc <- sprintf("data/model_objects/%s/%s_%s.rds", yr_lb, yr_lb, t1)
  t2_file_loc <- sprintf("data/model_objects/%s/%s_%s.rds", yr_lb, yr_lb, t2)
  
  mod1 <- readRDS(t1_file_loc)
  mod2 <- readRDS(t2_file_loc)
  tCop <- BiCopEst(as.numeric(pit(mod1$model)),
                   as.numeric(pit(mod2$model)),
                   family = 2, max.df = 30, se=TRUE)
  out <- data.table(tick1=t1, tick2=t2, year=as.character(yr_lb),
                    t1_unc_mean=uncmean(mod1$model), t1_unc_var=uncvariance(mod1$model),
                    t1_pers=persistence(mod1$model), t1_haflife=halflife(mod1$model),
                    t1_omega=coef(mod1$model)["omega"], t1_vmodel=mod1$model@model$modeldesc$vmodel,
                    t1_distr=mod1$model@model$modeldesc$distribution,
                    t1_n_arma=sum(mod1$model@model$modelinc[1:3]) - 1,
                    t1_n_garch=sum(mod1$model@model$modelinc[7:9]) - 3,
                    t2_unc_mean=uncmean(mod2$model), t2_unc_var=uncvariance(mod2$model),
                    t2_pers=persistence(mod2$model), t2_haflife=halflife(mod2$model),
                    t2_omega=coef(mod2$model)["omega"], t2_vmodel=mod2$model@model$modeldesc$vmodel,
                    t2_distr=mod2$model@model$modeldesc$distribution,
                    t2_n_arma=sum(mod2$model@model$modelinc[1:3]) - 1,
                    t2_n_garch=sum(mod2$model@model$modelinc[7:9]) - 3,
                    rho_m=tCop$par)
  collect3[[i]] <- out
}
dtfModStats <- rbindlist(collect3)
dtfModStats[, year := as.character(year)]

dtfSPComp <- fread("data/SandP_companies3.csv")
dtf2 <- merge(dtf, dtfModStats)
dtf2 <- merge(dtf2, dtfSPComp[, .(ticker, t1_sector=sector, t1_subindustry=subindustry)],
              by.x="tick1", by.y="ticker", all.x=TRUE)
dtf2 <- merge(dtf2, dtfSPComp[, .(ticker, t2_sector=sector, t2_subindustry=subindustry)],
              by.x="tick2", by.y="ticker", all.x=TRUE)
sub_industry_levels <- sort(unique(dtfSPComp$subindustry))
sector_levels <- sort(unique(dtfSPComp$sector))
dtf2[, t1_sector := factor(t1_sector, levels=sector_levels)]
dtf2[, t1_subindustry := factor(t1_subindustry, levels=sub_industry_levels)]
dtf2[, t1_vmodel := factor(t1_vmodel)]
dtf2[, t1_distr := factor(t1_distr)]
dtf2[, t2_sector := factor(t2_sector, levels=sector_levels)]
dtf2[, t2_subindustry := factor(t2_subindustry, levels=sub_industry_levels)]
dtf2[, t2_vmodel := factor(t2_vmodel)]
dtf2[, t2_distr := factor(t2_distr)]

dtf2[, intra_sector := as.integer(t1_sector == t2_sector)]
dtf2[, abs_mean_diff := abs(t1_unc_mean - t2_unc_mean)]
dtf2[, abs_omega_diff := abs(t1_omega - t2_omega)]
dtf2[, mean_var := apply(.SD, 1, mean), .SDcols=c("t1_unc_var", "t2_unc_var")]
dtf2[, sum_var := apply(.SD, 1, sum), .SDcols=c("t1_unc_var", "t2_unc_var")]
dtf2[, max_var := apply(.SD, 1, max), .SDcols=c("t1_unc_var", "t2_unc_var")]
dtf2[, abs_var_diff := apply(.SD, 1, function(x) abs(diff(x))), .SDcols=c("t1_unc_var", "t2_unc_var")]
dtf2[, lm_year := factor(year, levels=sort(unique(year)))]
dtf2[, pers_geom_avg := apply(.SD, 1, function(x) sqrt(prod(x))), .SDcols=c("t1_pers", "t2_pers")]
dtf2[, vmodel_pair := apply(.SD, 1, function(x) paste(sort(c(x)), collapse="_")), .SDcols=c("t1_vmodel", "t2_vmodel")]
dtf2[, distr_pair := apply(.SD, 1, function(x) paste(sort(c(x)), collapse="_")), .SDcols=c("t1_distr", "t2_distr")]
dtf2[, vmodel_pair := factor(vmodel_pair, levels=c("sGARCH_sGARCH", "gjrGARCH_sGARCH", "csGARCH_sGARCH",
                                                   "csGARCH_gjrGARCH", "csGARCH_csGARCH", "gjrGARCH_gjrGARCH"))]
dtf2[, distr_pair := factor(distr_pair, levels=c("std_std", "sstd_std", "sstd_sstd"))]
dtf2[, total_mean_par := apply(.SD, 1, sum), .SDcols=c("t1_n_arma", "t2_n_arma")]
dtf2[, total_var_par := apply(.SD, 1, sum), .SDcols=c("t1_n_garch", "t2_n_garch")]


# Degree 2 Polynomial and Intercept
lm_2i_0 <- lm(log(dtw) ~ I(rho_m) + I((rho_m)**2), data=dtf2)
lm_2i_1 <- lm(log(dtw) ~ I(rho_m) + I((rho_m)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg,
              data=dtf2)
lm_2i_2 <- lm(log(dtw) ~ I(rho_m) + I((rho_m)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg + lm_year,
              data=dtf2)
lm_2i_3 <- lm(log(dtw) ~ I(rho_m) + I((rho_m)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg + intra_sector + total_mean_par + total_var_par + lm_year + vmodel_pair + distr_pair,
              data=dtf2)

round(summary(lm_2i_0)$coefficients, 3)
round(summary(lm_2i_1)$coefficients, 3)
round(summary(lm_2i_2)$coefficients, 3)
round(summary(lm_2i_3)$coefficients, 3)

# Degree 2 Polynomial and constrained 
lm_2_0 <- lm(log(dtw) ~ -1 + I(rho_m - 1) + I((rho_m - 1)**2), data=dtf2)
lm_2_1 <- lm(log(dtw) ~ -1 + I(rho_m - 1) + I((rho_m - 1)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg,
              data=dtf2)
lm_2_2 <- lm(log(dtw) ~ -1 + I(rho_m - 1) + I((rho_m - 1)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg + lm_year,
              data=dtf2)
lm_2_3 <- lm(log(dtw) ~ -1 + I(rho_m - 1) + I((rho_m - 1)**2) + abs_mean_diff + abs_var_diff + sum_var + pers_geom_avg + intra_sector + total_mean_par + total_var_par + lm_year + vmodel_pair + distr_pair,
              data=dtf2)

round(summary(lm_2_0)$coefficients, 3)
round(summary(lm_2_1)$coefficients, 3)
round(summary(lm_2_2)$coefficients, 3)
round(summary(lm_2_3)$coefficients, 3)



dtf2[, plot(rho_m, log(dtw))]
dtf2[, lines(sort(rho_m), fitted(lm_2_0)[order(rho_m)], col='gold', type='b')]


# ------------------------------------------------