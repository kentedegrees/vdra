\name{summary.vdralogistic}
\alias{summary.vdralogistic}
\alias{summary.vdralogistic.object}
\alias{print.summary.vdralogistic}
\title{
	Summary Method for Vertical Distributed Logistic Regression Models
}
\description{
  Produces a summary of a fitted vdra logistic regression model.
}
\usage{
  ## S3 method for class 'vdralogistic'
  \method{summary}{vdralogistic}(object, ...)
}
\arguments{
	\item{object}{a \code{vdralogistic} object.}
  \item{...}{futher argumetns passed to or from other methods.}
}
\value{
  Returns an object of class \code{summary.vdralogistic}. Objects of this class have a method for the function \code{print}.  The following components must be included in \code{summary.vdralogistic} object.
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
\examples{
  summary(vdra_fit_logistic_A)
}
