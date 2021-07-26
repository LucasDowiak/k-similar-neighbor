# Test third and fourth moment for normality
jarque_bera <- function(x, k=1)
{
  d <- x - mean(x)
  n <- length(x)
  S1 <- sum(d**3) / n;  S2 <- (sum(d**2) / n)**(3/2)
  C1 <- sum(d**4) / n;  C2 <- (sum(d**2) / n)**2
  JB <- ((n - k + 1) / 6) * ((S1/S2)**2 + (C1/C2 - 3)**2 / 4)
  pvalue <- 1 - pchisq(JB, df = 2)
  return(pvalue)
}

# Lagrange multiplier test 
lm_test <- function(x, moment=1, nlags=10)
{
  d <- (x - mean(x))**moment
  m <- embed(d, nlags + 1)
  
  LM <- lsfit(x=m[,-1], y=m[,1], intercept=TRUE) # check if in intercept is needed
  rsquared <- 1 - sum(LM$residuals**2) / sum((m[,1] - mean(m[,1]))**2)
  pvalue <- pchisq(nrow(m) * rsquared, df = nlags, lower.tail = FALSE)
  structure(pvalue, name=names(x), moment=moment)
}


# provide the PIT of the residuals to the marginal models
# obj here is the @fit slot from a uGARCHfit object from the rugarch package
PIT <- function(obj)
{
  if (!inherits(obj, "uGARCHfit"))
    stop("Object must inherit from uGARCHfit model from rugarch package.")
  z <- residuals(obj, standardize=TRUE)
  shape <- coef(obj)[grep("shape", names(coef(obj)))]
  skew <- coef(obj)[grep("skew", names(coef(obj)))]
  lambda <- coef(obj)[grep("lambda", names(coef(obj)))]
  dentype <- obj@model$modeldesc$distribution
  pdist(distribution=dentype, as.numeric(z), lambda=lambda, skew=skew, shape=shape)
}


# given a model produced by rugarch, pull out the residuals, perform PIT and run
# gauntlet of marginal tests
marginal_tests <- function(obj, lm_lags=10, print=TRUE, plot=TRUE, hl_cl=0.95)
{
  isuGARCHfit <- inherits(obj, "uGARCHfit")
  if (isuGARCHfit) {
    # Standardized residuals force mean = 0 and var = 1
    # u <- PIT(obj)
    z <- residuals(obj, standardize=TRUE) # standardize refers to condition variance
    shape <- coef(obj)[grep("shape", names(coef(obj)))]
    skew <- coef(obj)[grep("skew", names(coef(obj)))]
    lambda <- coef(obj)[grep("lambda", names(coef(obj)))]
    dentype <- obj@model$modeldesc$distribution
    u <- pdist(distribution=dentype, z, lambda=lambda, skew=skew, shape=shape)
    d <- ddist(distribution=dentype, z, lambda=lambda, skew=skew, shape=shape)
  } else {
    u <- obj
  }
  
  # visual PIT transform
  if (plot) {
    on.exit(par(mfrow=c(1,1)))
    if (isuGARCHfit) {
      par(mfrow = c(1,2))
      support <- seq(-max(z), max(z), length=500)
      plot(density(z), col = "darkgoldenrod2", lwd = 2, xlab="", ylab="",
           main="Empirical vs theoretical")
      lines(support, ddist(dentype, support, lambda=lambda, skew=skew, shape=shape))
    }
    hist(u, freq = FALSE, breaks = 10, col = "darkgoldenrod2", xlab="", ylab="",
         main="Hist of PIT")
    abline(h=1, col="red", lty=2)
  }

  ####### Specification Tests
  ## Moment LM tests of PIT
  lm_results <- lapply(1:4, function(x) lm_test(u, x, nlags=lm_lags))
  lm_matrix <- matrix(unlist(lm_results),
                      dimnames=list(paste0("lm moment ", 1:4, ":"), "p-Value"))
  hl_test <- HLTest(u, lags = lm_lags, conf.level = hl_cl)  # Hong and Li Non-Parametric Density Test (2005, RFS) (rugarch package)
  # Normality Tests
  "add cramer von mises?"
  ks_test <- ks.test(u, punif)     # Kolmogorov-Smirnov test (stats package)
  sw_test <- shapiro.test(qnorm(as.vector(u))) # Shapiro Test for Normality (stats package)
  jb_test <- jarque_bera(qnorm(u), length(coef(obj)))  # Jarque-Bera test for joint normality of Skew and Kurtosis
  "berk_test <- BerkowitzTest(qnorm(u), lags = 15)"  # Obviously the Berkowitz test
  if (print) {
    cat("\n-------------------------------------------------------------------------\n")
    cat(sprintf("Lagrange Multiplier Test (Std Residuals) with %d lags\n", lm_lags))
    cat("Null: Independent Moment Condition\n\n")
    print(lm_matrix)
    cat("\n-------------------------------------------------------------------------\n")
    cat("Hong and Li non-parametric normality test of the qnorm(PIT)\n")
    print(hl_test)
    cat("\n-------------------------------------------------------------------------\n")
    cat("Shapiro-Wilks normality test of the qnorm(PIT)\n")
    cat("Null: Distribution is normally distributed\n\n")
    msg <- sprintf("W = %.5f, p-value = %.4f", sw_test[["statistic"]], sw_test[["p.value"]])
    cat(msg,"\n")
    cat("\n-------------------------------------------------------------------------\n")
    cat("Jarque Bera test of the qnorm(PIT)\n")
    cat("Null: Skew and kurtosis are jointly normal\n\n")
    print(matrix(jb_test, dimnames=list("", "p-Value")))
    cat("\n-------------------------------------------------------------------------\n")
    cat("Test that PIT are U(0,1)\n")
    print(ks_test)
  }
  return(invisible(list(
    lm_tests=lm_matrix,
    k_s=ks_test,
    shapiro_wilks=sw_test,
    jarque_bera=jb_test,
    hl_test=hl_test
  )))
}


verify_marginal_test <- function(mt, alpha=0.1)
{
  lmtests <- all(mt$lm_tests > alpha)
  kstest <- mt$k_s$p.value > alpha
  jbtest <- mt$jarque_bera > alpha
  return(all(lmtests, kstest, jbtest))
}

