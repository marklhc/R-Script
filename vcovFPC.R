vcovFPC <- function(object, popsize2 = NULL, 
                    popsize1 = NULL, KR = FALSE) {
  # Obtained finite-population-adjusted standard errors for fixed effect 
  # estimates for a fitted multilevel model
  #
  # Args:
  #   object: an R object of class merMod as resulting from lmer()
  #   popsize2: population size at level-2; if NULL, an infinite level-2
  #             population is assumed
  #   popsize1: population size at level-1; if NULL, an infinite level-1
  #             population is assumed
  #   KR: Whether Kenward-Roger approximation of standard errors should be used,
  #       which is recommended for sample number of clusters and cluster size. 
  #       Default to FALSE. 
  #
  # Returns: 
  #   The variance-covariance matrix of the fixed effect estimates, as
  #   returned by vcov()
  if (!inherits(object, "merMod")) {
    stop("Wrong input: Not a fitted model from lmer() with class merMod")
  }
  if (length(object@flist) != 1) {
    stop("Wrong input: Only models with two levels are supported")
  }
  if (is.null(popsize1) & is.null(popsize2)) {
    message("No FPC specified; return results from lme4::vcov.merMod()")
    return(vcov(object))
  }
  PR <- object@pp
  N <- unname(object@devcomp$dims["n"])
  nclus <- unname(ngrps(object))
  if (isTRUE(popsize2 > nclus)) fpc2 <- 1 - nclus / popsize2
  else {
    fpc2 <- 1
    message("No FPC needed at level-2")
  }
  if (isTRUE(popsize1 > N)) fpc1 <- 1 - N / popsize1
  else {
    fpc1 <- 1
    message("No FPC needed at level-1")
  }
  if (fpc1 == 1 & fpc2 ==1) {
    message("Return results from lme4::vcov.merMod()")
    return(vcov(object))
  }
  A <- PR$Lambdat %*% PR$Zt
  Astar <- A * sqrt(fpc2)
  X <- PR$X
  Astar_X <- Astar %*% X
  D <- Matrix::Diagonal(nrow(Astar), fpc1) + tcrossprod(Astar)
  Fisher_I <- (crossprod(X) - crossprod(solve(t(chol(D)), Astar_X))) / fpc1
  Phi <- solve(Fisher_I) * sigma(object)^2
  Phi <- as(Phi, "dpoMatrix")
  nmsX <- colnames(X)
  dimnames(Phi) <- list(nmsX, nmsX)
  if (!KR) {
    return(Phi)
  } else {
    if (!require("pbkrtest")) {
      stop("Please install the `pbkrtest` package for the use of Kenward-Roger correction!")
    } else {
      SigmaG <- pbkrtest::get_SigmaG(object, details = 0)
      vcov_kr <- pbkrtest:::vcovAdj16_internal(Phi, SigmaG, X, details = 0)
      vcov_kr <- as(Phi, "dpoMatrix")
      return(vcov_kr)
    }
  }
}
