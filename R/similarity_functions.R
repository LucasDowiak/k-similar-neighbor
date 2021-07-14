library(dtw)

# Similarity metrics -----------------------------------------------------------

# For all similarity metrics:
#
#   nm - specific ticker to set as reference
#
#   columns - subset of columns in M to calculate the sim to against the reference
#
#   M - data.table of pre-processed stock price series as columnar series
#
#   ... 
euclidean_dist <- function(nm, columns, M)
{
  j <- M[[nm]]
  D <- sweep(M[, columns, with=FALSE], 2, j, FUN="-")
  return(sqrt(rowSums(D**2)))
}


correlation <- function(nm, columns, M, method="pearson")
{
  j <- M[[nm]]
  f_ <- function(x) cor(x, j, method=method)
  out <- apply(as.matrix(M[, columns, with=FALSE]), 2, FUN=f_)
  out <- matrix(out, ncol=length(out), dim=list(nm, columns))
  return(out)
}


dynamic_time_warp <- function(nm, columns, M, ...) 
{
  j <-  matrix(M[[nm]], ncol=nrow(M))
  #f_ <- function(x) dtw::dtw(x, j, ...)$distance
  #out <- apply(as.matrix(M[, columns, with=FALSE]), 2, FUN=f_)
  #out <- matrix(out, ncol=length(out), dim=list(nm, columns))
  out <- t(proxy::dist(t(as.matrix(M[, columns, with=FALSE])), j, method="DTW"))
  row.names(out) <- nm
  return(out)
}


# Functions to run similarity metrics ------------------------------------------

# Reads in a stored similarity matrix located at the supplied file path (or creates
# a new matrix) and does some validation checks
#
#   file_ - path to stored similarity matrices
#
#   M - data.table of pre-processed stock price series as columnar series
#
check_ref_file <- function(file_, M)
{
  nms <- names(M)
  nc <- ncol(M)
  
  if (is(file_, "character")) {
    
    if (file.exists(file_)) {
      # Read in NxN similarity/distanct/correlation matrix
      OMG <- as.matrix(read.table(file_, header=TRUE, row.names = 1))
      cnms <- colnames(OMG)
      # row names must match column names
      B1 <- all(row.names(OMG) == colnames(OMG))
      # names in OMG must be subset of M
      B2 <- all(is.element(cnms, nms))
      
      if (!B1) {
        stop(sprintf("Row and column names do not match exactly in ref_file: %s", file_))
      }
      if (!B2) {
        warning(sprintf("The columns in `ref_file` are not a strict subset of the names in M: %s", file_))
      }
    } else {
      stop(sprintf("User supplied ref_file does not exists: %s", file_))
    }
  } else {
    OMG <- matrix(NA_real_, nrow=nc, ncol=nc, dimnames=list(nms, nms))
  }
  return(OMG)
}


# Performs some validation checks between the supplied pairs and data.frame
#
#   pairs - named list of lists 
#           e.g. $AAPL
#                [1] "AAPL"
#                [2] "AAPL" "ABC" "ABMD"
#
#                $ABC
#                [1] "ABC"
#                [2] "AAPL" "ABC" "ABMD"
#
#   M - data.table of pre-processed stock price series as columnar series
#
check_pairs <- function(pairs, M)
{
  nms <- setdiff(names(M), "Date")
  # if no pairs are provided, default to full NxN comparisons of the N columns in M
  if (is.null(pairs)) {
    out <- lapply(nms, function(x) list(x, nms))
    names(out) <- nms
  } else {
    if (is.null(names(pairs))) {
      stop("`pairs` needs to be a named list where the names must appear in M")
    }
    pnms <- c(names(pairs), unique(unlist(pairs)))
    B1 <- all(is.element(pnms, nms))
    if (!B1) {
      include <- intersect(pnms, nms)
      exclude <- setdiff(pnms, nms)
      out <- pairs[include]
      out <- lapply(out, function(x) {x[[2]] <- x[[2]][x[[2]] %in% include]; x})
      warning(sprintf("The values in `pairs` are not a strict subset of the names in M.\n  Including: [%s]\n  Excluding: [%s] \n",
                      paste(include, collapse=", "), paste(exclude, collapse=", ")))
    } else {
      out <- pairs
    }
  }
  return(out)
}


# Calculates the similarity of a batch of stock prices
#
#   M - data.table of pre-processed stock price series as columnar series
#   
#   pairs - named list of lists 
#           e.g. $AAPL
#                [1] "AAPL"
#                [2] "AAPL" "ABC" "ABMD"
#
#                $ABC
#                [1] "ABC"
#                [2] "AAPL" "ABC" "ABMD"
#
#   metric - one of the defined similarity metrics
#
#   ref_file - optional path to stored similarity matrices
#
#   overwrite - should the recent calculations be written back to disk
#
#   ncores - how many processors should be tasked to the job
#
#   ... - options to pass to metric
#
calculate_similarity <- function(M, metric, pairs=NULL, ref_file=NULL, overwrite=FALSE,
                                 ncores=1L, ...)
{
  if (is.null(ref_file) && overwrite)
    stop("If `overwrite` is TRUE then ref_file cannot be NULL")
  
  # Force M to be a data.table
  M <- as.data.table(M)
  
  # test for reference file
  OMG <- check_ref_file(ref_file, M)
  
  # check that every value in `pairs` is a name in M
  pairs <- check_pairs(pairs, M)
  
  # lapply over the pairs
  grab_metric <- function(x) metric(x[[1L]], x[[2L]], M=M, ...)
  if (ncores == 1L) {
    # do in sequence
    rho <- lapply(pairs, grab_metric)
    
  } else  if (ncores > 1L) {
    # do in parallel
    require(parallel)
    nc <- detectCores()
    if (ncores > nc) {
      warning(sprintf("`ncores` set to max number of threads available on this machine: %d", nc))
      ncores <- nc
    }
    cl <- makeForkCluster(ncores)
    on.exit(stopCluster(cl))
    rho <- parLapply(cl=cl, pairs, grab_metric)
  } else {
    stop("`ncores` needs to be an integer greater than or equal to 1")
  }
  
  # fill in Matrix OMG with the results
  for (nmi in names(rho)) {
    for (nmj in colnames(rho[[nmi]])) {
      OMG[nmi, nmj] <- rho[[nmi]][nmi, nmj] 
    }
  }
  
  # Should you write changes back to `ref_file`
  if (overwrite) {
    write.table(OMG, file=ref_file, sep="\t", col.names=NA)
  }
  return(list(rho, OMG))
}


# Function that leverages hierarchical clustring to  order similarity matrices
# into blocks
# References implementation of https://wil.yegelwel.com/cluster-correlation-matrix/
# And this one https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
#
#   X - Similarity matrix
#
#   type - type of similarity
#
#   visualize - should the plot of the clustered similarity matrix be rendered
#
cluster_similarity <- function(X, type=c("cor", "dtw"), visualize=FALSE)
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
  pairwise_dist <- proxy::dist(X)
  tree <- hclust(pairwise_dist)
  
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




