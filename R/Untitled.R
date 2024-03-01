setwd("~/Git/dtw-in-finance/")
library(data.table)
source("R/auto_marginal.R")
source("R/plots.R")
source("R/similarity_functions.R")# import read_and_melt, dynamic_time_warp from R/similarity_script.R



# Delannoy Number
function(x, y)
{
  L <- min(x, y)
  out <- vector("numeric", L + 1)
  for (k in seq(0, L)) {
    out[k + 1] <- choose(x, k) * choose(y, k) * 2**k
  }
  return(sum(out))
}



dtfStdRes1920 <- fread("data/label_analysis/DA_2019_2020_std_resids.csv")
dtfStdPrc1920 <- fread("data/label_analysis/DA_2019_2020_std_price.csv")
ticks <- sort(setdiff(names(dtfStdRes1920), c("Date")))

m_dtw <- dynamic_time_warp("A", ticks, dtfStdPrc1920[, -c("Date")])
m_cor <- correlation("A", ticks, dtfStdRes1920[, -c("Date")])
dtfM1920 <- data.table(dtw=c(m_dtw), cor=c(m_cor), tick=colnames(m_dtw))
plot(dtfM$cor, dtfM$dtw, ylim=c(0, 700), xlim=c(0, 1), main="2019-2020", ylab="dtw", xlab="cor")


dtfStdRes0708 <- fread("data/label_analysis/DA_2007_2008_std_resids.csv")
dtfStdPrc0708 <- fread("data/label_analysis/DA_2007_2008_std_price.csv")
ticks <- sort(setdiff(names(dtfStdRes0708), c("Date")))


m_dtw <- dynamic_time_warp("A", ticks, dtfStdPrc0708[, -c("Date")])
m_cor <- correlation("A", ticks, dtfStdRes0708[, -c("Date")])
dtfM0708 <- data.table(dtw=c(m_dtw), cor=c(m_cor), tick=colnames(m_dtw))
plot(dtfM0708$cor, dtfM0708$dtw, ylim=c(0, 700), xlim=c(0, 1), main="2007-2008", ylab="dtw", xlab="cor")


# Vanilla DTW vs Weighted DTW vs 
t1 <- "MSFT"; t2 <- "ADBE"
dtfStdPrc <- fread("data/label_analysis/2018_std_price.csv")
m1 <- readRDS(sprintf("data/model_objects/2018/2018_%s.rds", t1))
m2 <- readRDS(sprintf("data/model_objects/2018/2018_%s.rds", t2))

p <- dtfStdPrc[[t1]]
q <- dtfStdPrc[[t2]]
N <- dtfStdPrc[, .N]

summary(c(C))
image(D**2, y=1:505, col=grDevices::terrain.colors(100), x=1:505, 
      main=sprintf("Beta: %s", beta))
contour(D**2, x = 1:505, y = 1:505, add = TRUE)


unc_dtw <- dtw(p, q, keep.internals=TRUE, window.type="none")
scb_dtw <- dtw(p, q, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=15)

unc_wgt_dtw <- dtw(C, keep.internals=TRUE, window.type="none")
scb_wgt_dtw <- dtw(C, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=round(w * N))

image(W, col = grDevices::terrain.colors(100), x = 1:505, 
      y = 1:505)


idx_diff <- (scb_dtw$index1 - scb_dtw$index2)
summary(idx_diff)
summary(abs())



par(mfrow=c(2,1))
dtwPlotDensity(unc_dtw, main="Un-constrained DTW Algo")
dtwPlotDensity(scb_dtw, main=sprintf("Sakoe-Chiba (%s) DTW Algo", w))

dtwPlotDensity(unc_wgt_dtw, main="Un-constrained Wtd DTW Algo")
dtwPlotDensity(scb_wgt_dtw, main=sprintf("Sakoe-Chiba (%s) Wtd DTW Algo", w))

par(mfrow=c(2, 1))
dtwPlotTwoWay2(unc_dtw, main="Un-constrained DTW Algo")
dtwPlotTwoWay2(scb_dtw, main="Sakoe-Chiba (%s) DTW Algo") # xts=p, yts=q,

dtwPlotTwoWay2(unc_wgt_dtw,xts=p, yts=q, main="Un-constrained Wtd DTW Algo")
dtwPlotTwoWay2(scb_wgt_dtw, xts=p, yts=q, main="Sakoe-Chiba (%s) Wtd DTW Algo")

par(mfrow=c(2,1))
plot(sigma(m1$model))
plot(sigma(m2$model))


# Calculate dtw metrics
# ------------------------------------------------------------------------------
dtfSandP <- read_SandP_data()
dtfStdPrc1920 <- fread("data/label_analysis/DA_2019_2020_std_price.csv")
Rho1920 <- read.table("data/association_results/archive/corr_2019_2020", sep="\t",
                      header=TRUE, row.names = 1)
Rho1920 <- read.table("data/association_results/archive/corr_2019_2020", header=TRUE, row.names=1)
Rho1920 <- as.matrix(Rho1920)



ticks <- setdiff(names(dtfStdPrc1920), "Date")
w <- 0.1
N <- dtfStdPrc1920[,.N]
log_fun <- function(x, g, c0, wm=1)
{
  wm / (1 + exp(-g * (x - c0)))
}
W <- outer(1:N, 1:N, function(i, j) log_fun(abs(i - j), g=0.03, c0=N/4))

W_0 <- W_SCB <- W_WTD <- matrix(NA_real_, nrow=length(ticks), ncol=length(ticks),
                                dimnames = list(ticks, ticks))
tpairs <- combn(ticks, 2, simplify = FALSE)

for (tt in seq_along(tpairs)) {
  if (tt %% 250 == 0)
    print(sprintf("%d of %d started at %s", tt, length(tpairs), Sys.time()))
  t1 <- tpairs[[tt]][1]; t2 <- tpairs[[tt]][2]
  
  p <- dtfStdPrc1920[[t1]]
  q <- dtfStdPrc1920[[t2]]
  
  D <- outer(p, q, function(i, j) abs(i - j))
  C <- (D * W)**2
  
  dtw_0 <- dtw(p, q, keep.internals=TRUE, window.type="none")
   dtw_scb <- dtw(p, q, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=round(w * N))
  dtw_wgt <- dtw(C, keep.internals=TRUE, window.type="none")
  
  W_0[t1, t2] <- W_0[t2, t1] <- dtw_0$distance
  W_SCB[t1, t2] <- W_SCB[t2, t1] <- dtw_scb$distance
  W_WTD[t1, t2] <- W_WTD[t2, t1] <- dtw_wgt$distance
}

diag(W_0) <- 0
diag(W_SCB) <- 0
diag(W_WTD) <- 0

write.table(W_0, file="data/association_results/DA_2019_2020_dtw_unconstrained", sep="\t")
write.table(W_SCB, file="data/association_results/DA_2019_2020_dtw_sakoechiba", sep="\t")
write.table(W_WTD, file="data/association_results/DA_2019_2020_dtw_logweights", sep="\t")


# Use DTW distance to cluster stock series & test group statistics
library(dbscan)

WW <- W_SCB
diag(WW) <- NA_real_
d <- c(W_SCB[upper.tri(WW)])

knnd <- kNNdist(as.dist(WW), k=5)
plot(sort(knnd)[1:493], type="l")

cross_distance_set <- function(D, xlabels, ylabels=NULL)
{
  D <- copy(D)
  diag(D) <- NA_real_
  
  if (is.null(ylabels))
    ylabels <- xlabels

  b1 <- outer(rownames(D) %in% xlabels, colnames(D) %in% ylabels, `&`)
  return(c(D[b1]))
}


perform_dtw_cluster_analysis <- function(DTW, RHO, Price, k, plot_clusters=FALSE, method="ward")
{
  hc <-agnes(as.dist(DTW), method=method)
  sub_grp <- cutree(hc, k=k)
  
  dtfsp <- read_SandP_data()
  setnames(dtfsp, "ticker", "tick")
  dtfC <- merge(dtfsp, data.table(tick=colnames(DTW), group=sub_grp))
  
  classes <- sort(unique(sub_grp))
  NA_matrix <- matrix(NA_real_, ncol=length(classes), nrow=length(classes),
                      dimnames = list(classes, classes))
  m_rho <- NA_matrix
  m_rho_tt <- NA_matrix
  m_dtw <- NA_matrix
  m_dtw_tt <- NA_matrix
  
  for (g in classes) {
    
    tt1 <- dtfC[group==g, tick]
    gg_vals <- cross_distance_set(DTW, tt1)
    m_dtw[g, g] <- mean(gg_vals, na.rm=TRUE)
    
    rho_gg_vals <- cross_distance_set(RHO, tt1)
    m_rho[g, g] <- mean(rho_gg_vals, na.rm=TRUE)
    
    for (r in setdiff(classes, g)) {
      tt2 <- dtfC[group==r, tick]
      
      gr_vals <- cross_distance_set(DTW, tt1, tt2)
      m_dtw[g, r] <- mean(gr_vals, na.rm=TRUE)
      
      rho_gr_vals <- cross_distance_set(RHO, tt1, tt2)
      m_rho[g, r] <- mean(rho_gr_vals, na.rm=TRUE)
      
      if (length(gg_vals) > 2 & length(gr_vals) > 2) {
        ttest <- t.test(gg_vals, gr_vals, alternative="two.sided", mu=0)
        m_dtw_tt[g, r] <- ttest$p.value
        
        ttest2 <- t.test(rho_gg_vals, rho_gr_vals, alternative="two.sided", mu=0)
        m_rho_tt[g, r] <- ttest2$p.value
      }
    }
  }
  if (plot_clusters) {
    cluster_and_plot_series(D=DTW, DT=Price, k=k, method=method)
  }
  return(list(m_dtw=m_dtw, m_dtw_tt=m_dtw_tt, m_rho=m_rho, m_rho_tt=m_rho_tt,
              dtfC=dtfC))
}


dtfZ <- fread("data/label_analysis/DA_2019_2020_std_resids.csv")
dtfP <- fread("data/label_analysis/DA_2019_2020_std_price.csv")
R <- cor(dtfZ[, .SD, .SDcols=setdiff(names(dtfZ), "Date")])
W_SCB <- read.table("data/association_results/DA_2019_2020_dtw_sakoechiba",
                    header=TRUE, row.names=1)
W_SCB <- as.matrix(W_SCB)
tst <- perform_dtw_cluster_analysis(DTW=W_SCB, RHO=R, Price=dtfP, k=9, plot_clusters = T, method="ward")
plot(cbind(c(tst$m_dtw), c(tst$m_rho)), xlab="dtw", ylab="rho", main="K = 9")


round(sweep(tst$m_dtw, 2, diag(tst$m_dtw), `-`), 2)
sweep(m_rho, 2, diag(m_rho), `-`)


tt <- tst$m_dtw
tt[upper.tri(tt) | tst$m_rho < 0.4] <- NA_real_
round(tt, 2)

sweep(m_rho, 2, diag(m_rho), `-`)
m_rho_tt[upper.tri(m_rho_tt) & m_rho_tt > 0.05]

# ------------------------------------------------------------------------------


# Year-on-Year cluster analysis
# ------------------------------------------------------------------------------


# read in original price table
dtfSP <- fread("~/Git/dtw-in-finance/data/SandP_tick_history.csv")
dtfSP[, year := as.factor(year(Date))]
# break apart series into N separate periods
years <- 2010:2019
lst_years <- split(dtfSP[year %in% years], by="year")

standardize_price <- function(DT, ticks=NULL)
{
  if (is.null(ticks))
    ticks <- names(DT)
  ticks <- setdiff(ticks, "Date")
  P <- DT[, lapply(.SD, function(x) x/x[1]), .SDcols=ticks]
  return(P)
}

# Standardize prices in each period; Calculate DTW array of distances
calculate_dtw_matrix <- function(DT, ticks=NULL)
{
  N <- DT[,.N]
  # Omit any stocks with NA values
  DT <- DT[, Filter(function(x) !any(is.na(x)), .SD)]
  
  if (is.null(ticks))
    ticks <- setdiff(names(DT), "Date")
  
  P <- standardize_price(DT=DT, ticks=ticks)
  P$Date <- DT$Date
  W <- matrix(NA_real_,
              nrow=nrow(P), ncol=ncol(P),
              dimnames=list(ticks, ticks))
  
  tpairs <- combn(ticks, 2, simplify = FALSE)
  for (tt in seq_along(tpairs)) {
    if (tt %% 250 == 0)
      print(sprintf("%d of %d started at %s", tt, length(tpairs), Sys.time()))
    t1 <- tpairs[[tt]][1]; t2 <- tpairs[[tt]][2]
    p <- P[[t1]]; q <- P[[t2]]
    dtw_scb <- dtw(p, q, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=round(0.1 * N))
    W[t1, t2] <- W[t2, t1] <- dtw_scb$distance
  }
  diag(W) <- 0
  return(W)
}

lst_years <- lapply(lst_years,
                    function(x) {o <- standardize_price(x); o$Date <- x$Date; o})


results <- vector("list", length(years))
names(results) <- as.character(years)
for (y in years) {
  results[[as.character(y)]] <- calculate_dtw_matrix(dtfSP[year(Date) == y])
  saveRDS(results, file="~/Git/dtw-in-finance/data/dtw_by_year.rds")
}

# Cluster using same set of K
cluster_association_matrix <- function(W, k, method, sub_ticks=NULL)
{
  if (!is.null(sub_ticks)) {
    W <- W[sub_ticks, sub_ticks]
  }
  dtfsp <- read_SandP_data()[, .(ticker, sector, subindustry, added_date)]
  setnames(dtfsp, "ticker", "tick")
  
  hc <-agnes(as.dist(W), method=method)
  sub_grp <- cutree(hc, k=k)
  dtfC <- merge(dtfsp, data.table(tick=colnames(W), group=sub_grp))
  return(list(hc=hc, dtfC=dtfC))
  
}

tick_intersect <- Reduce(intersect, lapply(results, colnames))
lst <- lapply(results, cluster_association_matrix, k=9, method="ward", sub_ticks=tick_intersect)
for (y in names(lst)) {
  lst[[y]][[2]] <- lst[[y]][[2]][, year := as.integer(y)]
}

dtf <- rbindlist(lapply(lst, `[[`, 2))


aa <- cluster_and_plot_series(D=results[["2010"]], Price=lst_years[["2010"]],
                              k=9, return_melted = TRUE, method="ward")
bb <- cluster_and_plot_series(D=results[["2011"]], Price=lst_years[["2011"]],
                              k=9, return_melted = TRUE, method="ward")

# ------------------------------------------------------------------------------


# Troubleshoot model diagnostics for failed model specifications 
# ------------------------------------------------------------------------------
spec0708 <- read_json("data/marginal_specifications_2007_2008.json")
gvbt <- good_vs_bad_symbols(spec0708)


gvbt$fail_ticks

tick_ <- "MTB"
dtfU <- data.table(parse_json(tick_))
u <- dtfU[year(Date) %in% c(2007, 2008), diff(log(close))]


fit <- spec0708[[tick_]][1:7]
spec_mod <- ugarchspec(
  variance.model = list(model=fit$garchmod[[1]], garchOrder=c(as.integer(fit$arch[[1]]), as.integer(fit$garch[[1]]))),
  mean.model = list(armaOrder=c(fit$ar[[1]], fit$ma[[1]]), include.mean=TRUE),
  distribution.model = fit$distr[[1]]
)


fit_ <- try(ugarchfit(spec=spec_mod, data=u))
marg_tests <- marginal_tests(fit_, print=TRUE, plot=TRUE)
verify_marginal_test(marg_tests)


summary(u)
(min_u <- which.min(u))
dtfU[year(Date) %in% 2007:2008][(min_u - 5):(min_u + 5)]
u <- u[-which.min(u)]

fit <- auto_fit(u, 3, 3)
spec_mod <- ugarchspec(
  variance.model = list(model=fit$garchmod, garchOrder=c(as.integer(fit$arch), as.integer(fit$garch))),
  mean.model = list(armaOrder=c(fit$ar, fit$ma), include.mean=TRUE),
  distribution.model = fit$distr
)

fit_ <- try(ugarchfit(spec=spec_mod, data=u))
marg_tests <- marginal_tests(fit_, print=TRUE, plot=TRUE)
verify_marginal_test(marg_tests)
# ------------------------------------------------------------------------------


# Simulate GARCH models
# -------------------------------------------------
spec0708 <- read_json("data/marginal_specifications_2007_2008.json")
gvbt <- good_vs_bad_symbols(spec0708)

gvbt$fail_ticks

tickU <- "GOOGL"
tickV <- "LEG"
dtfU <- data.table(parse_json(tickU))
dtfV <- data.table(parse_json(tickV))
u <- dtfU[year(Date) %in% c(2007, 2008), diff(log(close))]
v <- dtfV[year(Date) %in% c(2007, 2008), diff(log(close))]

specU <- spec0708[[tickU]][1:7]
ugarchspecU <- ugarchspec(
  variance.model = list(model=specU$garchmod[[1]], garchOrder=c(as.integer(specU$arch[[1]]), as.integer(specU$garch[[1]]))),
  mean.model = list(armaOrder=c(specU$ar[[1]], specU$ma[[1]]), include.mean=TRUE),
  distribution.model = specU$distr[[1]]
)
fitU <- try(ugarchfit(spec=ugarchspecU, data=u))

specV <- spec0708[[tickV]][1:7]
ugarchspecV <- ugarchspec(
  variance.model = list(model=specV$garchmod[[1]], garchOrder=c(as.integer(specV$arch[[1]]), as.integer(specV$garch[[1]]))),
  mean.model = list(armaOrder=c(specV$ar[[1]], specV$ma[[1]]), include.mean=TRUE),
  distribution.model = specV$distr[[1]]
)
fitV <- try(ugarchfit(spec=ugarchspecV, data=v))


simU <- ugarchsim(fit=fitU, n.sim=400, n.start=100, m.sim=1)
simV <- ugarchsim(fit=fitV, n.sim=400, n.start=100, m.sim=1)
cor(simU@simulation$seriesSim, simV@simulation$seriesSim)


# ------------------------------------------------------------------------------



