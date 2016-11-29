#' Update Cholesky Decomposition
#'
#' Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
#' function \code{chol_update} updates L such that it corresponds to the decomposition of A + u*u'.
#'
#' @useDynLib ramcmc
#' @importFrom Rcpp evalCpp
#' @param L A lower triangular matrix. Strictly upper diagonal part is not referenced.
#' @param u A vector with with length matching with the dimensions of L.
#' @return Updated L.
#' @export
chol_update <- function(L, u) {
  check_args(L, u)
  chol_updateR(L, u)
}
#' Downdate Cholesky Decomposition
#'
#' Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
#' function \code{chol_downdate} updates L such that it corresponds to the decomposition of A - u*u'
#' (if such decomposition exists).
#'
#' @note The function does not check that the resulting matrix is positive semidefinite.
#'
#' @param L A lower triangular matrix. Strictly upper diagonal part is not referenced.
#' @param u A vector with with length matching with the dimensions of L.
#' @return Updated L.
#' @export
chol_downdate <- function(L, u) {
  check_args(L, u)
  L <- chol_downdateR(L, u)
  if(any(!is.finite(L)) || any(diag(L) < 0)) {
    stop("Resulting matrix is not positive definite.")
  }
  L
}
#' Update the Proposal of RAM Algorithm
#'
#' Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
#' function \code{adapt_L} updates L according to the RAM algorithm.
#'
#' @note The function does not check that the resulting matrix is positive semidefinite.
#'
#' @param L A lower triangular matrix corresponding to the Cholesky decomposition.
#' @param u A vector with with length matching with the dimensions of L.
#' @param current The current acceptance probability.
#' @param n Scaling parameter corresponding to the current iteration number.
#' @param target The target acceptance rate. Default is 0.234.
#' @param gamma Scaling parameter. Default is 2/3.
#' @return Updated L.
#' @references Matti Vihola (2012). "Robust adaptive Metropolis algorithm with coerced acceptance rate".
#' Statistics and Computing, 22: 997. doi:10.1007/s11222-011-9269-5
#' @export
#' @examples
#'
#' # sample from standard normal distribution
#'
#' adapt_mcmc <- function(n = 10000, sigma) {
#'   x <- numeric(n)
#'   loglik_old <- dnorm(x[1], log = TRUE)
#'   for (i in 2:n) {
#'     u <- rnorm(1, sd = sigma)
#'     prop <- x[i] + u
#'     loglik <- dnorm(prop, log = TRUE)
#'     accept_prob <- min(1, exp(loglik - loglik_old))
#'     if (runif(1) < accept_prob) {
#'       x[i] <- prop
#'       loglik_old <- loglik
#'     } else {
#'       x[i] <- x[i - 1]
#'     }
#'     if (i < n/2) {
#'       sigma <- adapt_L(sigma, u, accept_prob, i)
#'     }
#'   }
#'   list(x = x[(n/2):n], sigma = sigma)
#' }
#'
#' out <- adapt_mcmc(1e5, 2)
#' out$sigma
#' hist(out$x)
#' # acceptance rate:
#' 1 / mean(rle(out$x)$lengths)
#'
adapt_L <- function(L, u, current, n, target = 0.234, gamma = 2/3) {

  if(!is.matrix(L) && length(L) == 1) {
    L <- as.matrix(L)
  }
  check_args(L, u)
  if(any(L[upper.tri(L)] != 0)) {
    stop("L must be lower triangular matrix.")
  }
  if(n < 0) {
    stop("Argument 'n' must be non-negative.")
  }
  if(length(current) > 1 || current > 1 || current < 0) {
    stop("Argument 'current' must be a value between 0 and 1.")
  }
  if(length(target) > 1 || target > 1 || target < 0) {
    stop("Argument 'target' must be a value between 0 and 1.")
  }
  if(length(gamma) > 1 || gamma > 1 || gamma < 0) {
    stop("Argument 'gamma' must be a value between 0 and 1.")
  }
  adapt_LR(L, u, current, target, n, gamma)
}

check_args <- function(L, u){
  if (!is.matrix(L)) {
    stop("Argument 'L' must be a matrix.")
  }
  if (!is.numeric(u)) {
    stop("Argument 'u' must be a numeric vector.")
  }
  if (nrow(L) != ncol(L)) {
    stop("Argument 'L' must be square matrix.")
  }
  if (nrow(L) != length(u)) {
    stop("Length of 'u' must be equal to the number of rows in 'L'.")
  }
}