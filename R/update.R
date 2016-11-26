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
  cholupdate(L, u)
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
  L <- choldowndate(L, u)
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
update_S <- function(S, u, current, n, target = 0.234, gamma = 2/3) {
  if(any(S[upper.tri(S)] != 0)) {
    stop("S must be lower triangular matrix.")
  }
  adjust_S_copy(S, u, current, target, n, gamma)
}