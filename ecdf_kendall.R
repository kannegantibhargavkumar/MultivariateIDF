### empirical kendall function

genEmpKenFun <- function(copula, sample=NULL) {
  if(is.null(sample)) 
    sample <- rCopula(1e6,copula)
  if(missing(copula)) {
    # taken from package copula function "Kn"
    stopifnot((n <- nrow(sample)) >= 1, (d <- ncol(sample)) >= 1)
    ken <- vapply(seq_len(n), function(i) sum(colSums(t(sample) < sample[i, ]) == d)/(n + 1), NA_real_)
  } else {
    ken <- pCopula(sample, copula)
  }
  
  return(ecdf(ken))
}

