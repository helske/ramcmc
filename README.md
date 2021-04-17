[![downloads](http://cranlogs.r-pkg.org/badges/ramcmc)](http://cranlogs.r-pkg.org/badges/ramcmc)
[![cran version](http://www.r-pkg.org/badges/version/ramcmc)](http://cran.r-project.org/package=ramcmc)

# ramcmc: Building blocks for Robust Adaptive Metropolis algorithm

This small package provides key functions for the [RAM algorithm by Vihola (2012)](http://link.springer.com/article/10.1007/s11222-011-9269-5). These can be used directly from R or the corresponding header files can be linked to other R/Rcpp packages.

The package cointains three functions, two of which can be useful in more general context as well: The Cholesky update and downdate functions, and the actual update function for the covariance matrix of the proposal distribution used in RAM. For details, see package vignette.
