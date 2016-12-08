#' Rank-one Update of Cholesky Decomposition
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
#' @examples
#'
#' L <- matrix(c(4,3,0,5), 2, 2)
#' u <- c(1, 2)
#' chol_update(L, u)
#' t(chol(L %*% t(L) + u %*% t(u)))
#'
chol_update <- function(L, u) {
  check_args(L, u)
  chol_updateR(L, u)
}
#' Rank-one Downdate of Cholesky Decomposition
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
#' Given the lower triangular matrix S obtained from the Cholesky decomposition of the shape
#' of the proposal distribution, function \code{adapt_S} updates S according to the RAM algorithm.
#'
#' @note If the downdating would result non-positive definite matrix, no adaptation is performed.
#'
#' @param S A lower triangular matrix corresponding to the Cholesky decomposition of the
#' scale of the proposal distribution.
#' @param u A vector with with length matching with the dimensions of S.
#' @param current The current acceptance probability.
#' @param n Scaling parameter corresponding to the current iteration number.
#' @param target The target acceptance rate. Default is 0.234.
#' @param gamma Scaling parameter. Default is 2/3.
#' @return If the resulting matrix is positive definite, an updated value of S.
#'  Otherwise original S is returned.
#' @references Matti Vihola (2012). "Robust adaptive Metropolis algorithm with coerced acceptance rate".
#' Statistics and Computing, 22: 997. doi:10.1007/s11222-011-9269-5
#' @export
#' @examples
#'
#' # sample from standard normal distribution
#' # use proposals from the uniform distribution on
#' # interval (-s, s), where we adapt s
#'
#' adapt_mcmc <- function(n = 10000, s) {
#'   x <- numeric(n)
#'   loglik_old <- dnorm(x[1], log = TRUE)
#'   for (i in 2:n) {
#'     u <- s * runif(1, -1, 1)
#'     prop <- x[i] + u
#'     loglik <- dnorm(prop, log = TRUE)
#'     accept_prob <- min(1, exp(loglik - loglik_old))
#'     if (runif(1) < accept_prob) {
#'       x[i] <- prop
#'       loglik_old <- loglik
#'     } else {
#'       x[i] <- x[i - 1]
#'     }
#'     # Adapt only during the burn-in
#'     if (i < n/2) {
#'       s <- adapt_S(s, u, accept_prob, i)
#'     }
#'   }
#'   list(x = x[(n/2):n], s = s)
#' }
#'
#' out <- adapt_mcmc(1e5, 2)
#' out$s
#' hist(out$x)
#' # acceptance rate:
#' 1 / mean(rle(out$x)$lengths)
#'
adapt_S <- function(S, u, current, n, target = 0.234, gamma = 2/3) {

  if(!is.matrix(S) && length(S) == 1) {
    S <- as.matrix(S)
  }
  check_args(S, u, matrix_name = "S")
  if(any(S[upper.tri(S)] != 0)) {
    stop("'S' must be lower triangular matrix.")
  }
  if(n < 0) {
    stop("Argument 'n' must be non-negative.")
  }
  if(length(current) > 1 || current > 1 || current < 0) {
    stop("Argument 'current' must be on interval [0, 1].")
  }
  if(length(target) > 1 || target >= 1 || target <= 0) {
    stop("Argument 'target' must be on interval (0, 1).")
  }
  if(length(gamma) > 1 || gamma > 1 || gamma < 0) {
    stop("Argument 'gamma' must be a value between 0 and 1.")
  }
  adapt_SR(S, u, current, target, n, gamma)
}

check_args <- function(L, u, matrix_name = "L"){

  if (!is.matrix(L)) {
    stop(paste0("Argument '", matrix_name, "' must be a matrix."))
  }
  if (!is.numeric(u)) {
    stop("Argument 'u' must be a numeric vector.")
  }
  if (nrow(L) != ncol(L)) {
    stop(paste0("Argument '", matrix_name, "' must be square matrix."))
  }
  if (nrow(L) != length(u)) {
    stop(paste0("Length of 'u' must be equal to the number of rows in '", matrix_name, "'."))
  }
}