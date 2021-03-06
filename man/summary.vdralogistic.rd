\name{summary.vdralogistic}
\alias{summary.vdralogistic}
\title{
	Summary Method for Vertical Distributed Logistic Regression Models
}
\description{
  Produces a summary of a fitted vdra logistic regression model.
}
\usage{
## S3 method for class 'vdralogistic'
summary(x)
}
\arguments{
	\item{x}{a \code{vdralogistic} object.}
}
\value{
  An object of \code{summary.vdralogistic} with components:
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{converged}{logical value.  If \code{TRUE}, the regression converged.  If \code{FALSE}, it did not.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be missing values in the vector corresponding to the redudant columns model matrix.}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{tvals}{the t-values of the coefficietns.}
    \item{pvals}{the p-values of the coefficients.}
    \item{nulldev}{the null deviance of the fit.}
    \item{nulldev_df}{the degrees of freedom for the null deviance.}
    \item{resdev}{the residual deviance of the fit.}
    \item{resdev_df}{the degrees of freedome for the residual deviance.}
    \item{aic}{the AIC of the fit.}
    \item{bic}{the BIC of the fit.}
    \item{iter}{the number of iterations of the cox algorithm before convergence.}
}
\seealso{
  \code{\link{vdralogistic}}
}
