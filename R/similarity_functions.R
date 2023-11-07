library(dtw)

# Wrapper to read in 
read_association_table <- function(file_)
{
  as.matrix(read.table(file_, header = T, row.names=1))
}


# Calculates the similarity of a batch of stock prices
#
#   M - matrix or data.table of pre-processed stock price series
#   
#   tickers - possible subset of stock tickers. defaults to all of M
#
#   dist.method - one of the defined similarity metrics in the proxy::dist function
#
#   ... - options to pass to the proxy::dist function
#
calculate_similarity <- function(M, tickers, dist.method=c("Euclidean", "Correlation", "DTW"), ...)
{
  dist.method <- match.arg(dist.method)
  stopifnot(inherits(M, c("matrix", "data.table")))

  # Force M to be a data.frame
  if (inherits(M, "matrix")) {
    M <- as.data.table(M)
  }
  
  if (missing(tickers)) {
    tickers <- setdiff(names(M), c("Date", "label"))
  } else {
    # check that every value in `tickers` is a name in M
    if(!all(tickers %in% names(M))) {
      missticks <- tickers[!tickers %in% names(M)]
      stop(sprintf("Some tickers are no present in data.frame M\nMissing Ticks: %s",
                   paste(missticks, collapse=", ")))
    }
  }
  
  if (dist.method=="Correlation") {
    D <- cor(M[, tickers, with=FALSE], ...)
  } else {
    D <- proxy::dist(t(as.matrix(M[, tickers, with=FALSE])),
                     t(as.matrix(M[, tickers, with=FALSE])),
                     method=dist.method,
                     ...)
  }
  return(D)
}


# Function that leverages hierarchical clustering to order similarity matrices
# into blocks
# References implementation of https://wil.yegelwel.com/cluster-correlation-matrix/
# And this one https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
#
#   X - Similarity/distance/correlation matrix
#
#   type - type of similarity
#
#   visualize - should the plot of the clustered similarity matrix be rendered
#
#   ... - additional variables to pass to stats::hclust
#
cluster_similarity <- function(X, type=c("cor", "dtw"), visualize=FALSE, ...)
{
  require(ggplot2)
  type <- match.arg(type)
  if (type == "cor") {
    midpoint <- 0
    limit <- c(-1, 1)
  } else if (type == "dtw") {
    midpoint <- 100
    limit <- c(0, 400)
  }
  # Default euclidean distance
  pairwise_dist <- proxy::dist(X)
  tree <- hclust(pairwise_dist, ...)
  
  if (visualize) {
    tt <- tree$labels[tree$order]
    XX <- X[tree$order, tree$order]
    XP <- as.data.table(XX)
    XP[, ticks := factor(tt, levels=tt)]
    XP <- melt(XP, id.vars="ticks", measure.vars=setdiff(names(X), "ticks"))
    XP[, variable := factor(as.character(variable), levels=tt)]
    XP <- XP[order(variable, ticks)]
    
    p <- ggplot(data=XP, aes(x=ticks, y=variable, fill=value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = midpoint, limit = limit, space = "Lab", 
                           name="") +
      theme_minimal() +
      theme(axis.text=element_blank())
    plot(p)
  }
  return(tree)
}

