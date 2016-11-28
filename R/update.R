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
  L <- chol_downdateR(L, u)
  if(any(!is.finite(L)) || any(diag(L) < 0)) {
    stop("Resulting matrix is not positive definite.")
  }
  L
}
#' Update the Proposal of RAM Algorithm
#'
#' Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
#' function chol_downdate updates L such that it corresponds to the decomposition of A - u*u'
#' (if such decomposition exists).
#'
#' @note The function does not check that the resulting matrix is positive semidefinite.
#'
#' @param S A lower triangular matrix corresponding to the Cholesky decomposition.
#' @param u A vector with with length matching with the dimensions of L.
#' @param current The current acceptance probability.
#' @param n Scaling parameter corresponding to the current iteration number.
#' @param target The target acceptance rate. Default is 0.234.
#' @param gamma Scaling parameter. Default is 2/3.
#' @return Updated S.
#' @export
#' @examples
#'
#' # sample from standard normal distribution
#'
#' adapt_mcmc <- function(n = 10000, sigma) {
#'   x <- numeric(n)
#'   loglik_old <- dnorm(x[1])
#'   for (i in 2:n) {
#'     u <- rnorm(1, sd = sigma)
#'     prop <- x[i] + u
#'     loglik <- dnorm(prop)
#'     accept_prob <- min(1, loglik/loglik_old)
#'     if (runif(1) < accept_prob) {
#'       x[i] <- prop
#'       loglik_old <- loglik
#'     } else {
#'       x[i] <- x[i - 1]
#'     }
#'     if (i < n/2) {
#'       sigma <- update_S(sigma, u, accept_prob, i, 0.44)
#'     }
#'   }
#'   list(x = x[(n/2):n], sigma = sigma)
#' }
#'
#' out <- adapt_mcmc(1e4, 2)
#' out$sigma
#' hist(out$x)
#' # acceptance rate:
#' 1 / mean(rle(out$x)$lengths)
update_S <- function(S, u, current, n, target = 0.234, gamma = 2/3) {

  if(!is.matrix(S) && length(S) == 1) {
    S <- as.matrix(S)
  }

  if(any(S[upper.tri(S)] != 0)) {
    stop("S must be lower triangular matrix.")
  }
  adjust_SR(S, u, current, target, n, gamma)
}