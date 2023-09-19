setwd("~/Git/k-similar-neighbor/")
library(data.table)
source("R/auto_marginal.R")
source("R/plots.R")
source("R/similarity_functions.R")# import read_and_melt, dynamic_time_warp from R/similarity_script.R

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
t1 <- "SBUX"; t2 <- "CSCO"

m1 <- readRDS(sprintf("data/model_objects/DA_2019_2020/DA_2019_2020_%s.rds", t1))
m2 <- readRDS(sprintf("data/model_objects/DA_2019_2020/DA_2019_2020_%s.rds", t2))

p <- dtfStdPrc1920[[t1]]
q <- dtfStdPrc1920[[t2]]
N <- dtfStdPrc1920[, .N]

g <- 0.03
wm <- 1
c0 <- N/4
log_fun <- function(x, g, c0, wm=1)
{
  wm / (1 + exp(-g * (x - c0)))
}

# alpha <- 0.01
# beta <- (1 / (1 - alpha))
# time_discount <- function(x, beta)
# {
#   1 - beta**(x)
# }

w <- 0.1
beta <- 0.99
l <- 2
D <- outer(p, q, function(i, j) abs(i - j))
# W <- outer(1:length(p), 1:length(q), function(i, j) 1 - beta**(abs(i - j)))
W <- outer(1:length(p), 1:length(q), function(i, j) log_fun(abs(i - j), g=0.03, c0=N/4))
# diag(W) <- (1 - beta)
# W <- sweep(W, 2, colSums(W), `/`)
C <- (D * W)**(l)

summary(c(C))
image(D**2, y=1:505, col=grDevices::terrain.colors(100), x=1:505, 
      main=sprintf("Beta: %s", beta))
contour(D**2, x = 1:505, y = 1:505, add = TRUE)


variance_weights <- function(v, u, beta=0.95) {
  nv <- length(v); nu <- length(u)
  W <- matrix(NA_real_, nrow=nv, ncol=nu)
  for (i in seq_len(nv)) {
    for (j in seq_len(nu)) {
      pow <- abs(i - j)
      
      W[i, j] <- sum(v[i:j], u[i:j]) # These individual costs need to be discounted by time
    }
  }
  return(W)
}

unc_dtw <- dtw(p, q, keep.internals=TRUE, window.type="none")
scb_dtw <- dtw(p, q, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=round(w * N))

unc_wgt_dtw <- dtw(C, keep.internals=TRUE, window.type="none")
scb_wgt_dtw <- dtw(C, keep.internals=TRUE, window.type=sakoeChibaWindow, window.size=round(w * N))

image(W, col = grDevices::terrain.colors(100), x = 1:505, 
      y = 1:505)

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

library(dbscan)

WW <- W_SCB
diag(WW) <- NA_real_
d <- c(W_SCB[upper.tri(WW)])

knnd <- kNNdist(as.dist(WW), k=5)
plot(sort(knnd)[1:493], type="l")



# db <- dbscan(as.dist(WW), eps=20, minPts=5)
WW <- W_SCB
hc <-agnes(as.dist(WW), method="ward")
sub_grp <- cutree(hc, k=12)

dtfC <- data.table(tick=colnames(WW), group=sub_grp)
dtfC <- merge(dtfSandP, dtfC, by.x = "ticker", by.y = "tick")
setnames(dtfC, "ticker", "tick")

cross_distance_set <- function(D, xlabels, ylabels=NULL)
{
  D <- copy(D)
  diag(D) <- NA_real_
  
  if (is.null(ylabels))
    ylabels <- xlabels
  
  b1 <- outer(rownames(D) %in% xlabels, colnames(D) %in% ylabels, `&`)
  return(c(D[b1]))
}


classes <- sort(unique(sub_grp))

NA_matrix <- matrix(NA_real_, ncol=length(classes), nrow=length(classes), dimnames = list(classes, classes))
m_rho <- NA_matrix
m_rho_tt <- NA_matrix
m_dtw <- NA_matrix
m_dtw_tt <- NA_matrix

for (g in classes) {

  tt1 <- dtfC[group==g, tick]
  gg_vals <- cross_distance_set(WW, tt1)
  m_dtw[g, g] <- mean(gg_vals, na.rm=TRUE)
  
  rho_gg_vals <- cross_distance_set(Rho1920, tt1)
  m_rho[g, g] <- mean(rho_gg_vals, na.rm=TRUE)
  
  for (r in setdiff(classes, g)) {
    tt2 <- dtfC[group==r, tick]
    
    gr_vals <- cross_distance_set(WW, tt1, tt2)
    m_dtw[g, r] <- mean(gr_vals, na.rm=TRUE)
    
    rho_gr_vals <- cross_distance_set(Rho1920, tt1, tt2)
    m_rho[g, r] <- mean(rho_gr_vals, na.rm=TRUE)
    
    if (length(gg_vals) > 2 & length(gr_vals) > 2) {
      ttest <- t.test(gg_vals, gr_vals, alternative="two.sided", mu=0)
      m_dtw_tt[g, r] <- ttest$p.value
      
      ttest2 <- t.test(rho_gg_vals, rho_gr_vals, alternative="two.sided", mu=0)
      m_rho_tt[g, r] <- ttest2$p.value
    }
  }
}

plot(cbind(c(m_dtw), c(m_rho)))
sweep(m_dtw, 2, diag(m_dtw), `-`)
sweep(m_rho, 2, diag(m_rho), `-`)


round(m_rho * as.numeric(lower.tri(m_rho, diag = TRUE) & m_rho > 0.00),
      3)
tt <- m_rho
tt[upper.tri(tt) | m_rho < 0.00] <- NA_real_
round(tt, 3)

sweep(m_rho, 2, diag(m_rho), `-`)
m_rho_tt[upper.tri(m_rho_tt) & m_rho_tt > 0.05]

# Troubleshoot model diagnostics for failed model specifications 
# -------------------------------------------------
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




