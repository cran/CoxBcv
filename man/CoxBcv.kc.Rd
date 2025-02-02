% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoxBcv.kc.R
\name{CoxBcv.kc}
\alias{CoxBcv.kc}
\title{Kauermann and Carroll (KC) bias-corrected sandwich variance estimator}
\usage{
CoxBcv.kc(Y, Delta, X, ID)
}
\arguments{
\item{Y}{vector of observed time-to-event data.}

\item{Delta}{vector of censoring indicators.}

\item{X}{matrix of marginal mean covariates with one column for one covariate (design matrix excluding intercept).}

\item{ID}{vector of cluster identifiers.}
}
\value{
\itemize{
  \item coef - estimate of coefficients.
  \item exp(coef) - estimate of hazard ratio.
  \item KC-var - KC bias-corrected sandwich variance estimate of coef.
}
}
\description{
Calculate the Kauermann and Carroll (KC; 2001) bias-corrected sandwich variance estimator, for marginal Cox analysis of cluster randomized trials, 
proposed by Wang et al. (under review).
}
\examples{
Y <- c(11,19,43,100,7,100,100,62,52,1,7,6)
Delta <- c(1,1,1,0,1,0,0,1,1,1,1,1)
X1 <- c(0,0,0,0,0,0,1,1,1,1,1,1)
X2 <- c(-19,6,-25,48,10,-25,15,22,17,-9,45,12)
ID <- c(1,1,2,2,3,3,4,4,5,5,6,6)

X <- X1
CoxBcv.kc(Y,Delta,X,ID)

X <- cbind(X1,X2)
CoxBcv.kc(Y,Delta,X,ID)

}
\references{
Kauermann, G., & Carroll, R. J. (2001). 
A note on the efficiency of sandwich covariance matrix estimation. 
Journal of the American Statistical Association, 96(456), 1387-1396.

Wang, X., Turner, E. L., & Li, F. 
Improving sandwich variance estimation for marginal Cox analysis of cluster randomized trials. 
Under Review.
}
