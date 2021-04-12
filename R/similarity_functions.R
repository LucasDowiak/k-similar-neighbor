euclidean_dist <- function(nm, columns, M)
{
  j <- M[[nm]]
  D <- sweep(M[, columns, with=FALSE], 2, j, FUN="-")
  return(sqrt(rowSums(D**2)))
}


covariance <- function(nm, columns, M, method="pearson")
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


check_ref_file <- function(file_, M)
{
  
  nms <- names(M)
  nc <- ncol(M)
  
  if (is(file_, "character")) {
    
    if (file.exists(file_)) {
      
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


calculate_similarity <- function(M, metric, pairs=NULL, ref_file=NULL, overwrite=FALSE,
                                 ncores=1L, ...)
{
  # M - data.table of data
  # pairs - optional list of binary pairs to calculate the metric for
  # metric - what similarity metric do we want
  # ref_file - optional file to read in 
  # overwrite - should the recent calculations be written back to disk
  # ncores - how many processors should be tasks to the job
  # ... - options to pass to metric
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
    nc <- detectCores()
    if (ncores > nc) {
      warning(sprintf("`ncores` set to max number of threads available on this machine: %d", nc))
      ncores <- nc
    }
    cl <- makeForkCluster(ncores)
    on.exit(stopCluster(cl))
    rho <- parLapply(cl=cl, pairs, grab_metric)
  } else {
    stop("`rho` needs to be an integer greater than 1")
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
