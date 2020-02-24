\name{summary.vdracox}
\alias{summary.vdracox}
\title{
	Summary Method for Vertical Distributed COX Models}
\description{
  Produces a summary of a fitted vdra cox model.
}
\usage{
## S3 method for class 'vdracox'
summary(x)
}
\arguments{
	\item{x}{a \code{vdracox} object.}
}
\value{
  An object of \code{summary.vdracox} with components:
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{converged}{logical value.  If \code{TRUE}, the regression converged.  If \code{FALSE}, it did not.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be missing values in the vector corresponding to the redudant columns model matrix.}
    \item{expcoef}{a vector which represents exp(coefficients).}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{zvals}{the z-values of the coefficients.}
    \item{pvals}{the p-values of the coefficients.}
    \item{expncoef}{a vector which represents exp(-coefficients).}
    \item{lower95}{a vector of the lower bounds of the 95\% confidence interval for exp(coefficients).}
    \item{upper95}{a vector of the upper bounds of the 95\% confidence interval for exp(coefficients).}
    \item{n}{the number of observations in the data.}
    \item{nevent}{the number of events used in the fit.}
    \item{concordance}{a vector containing the number of events which are concordant, discordant, tied.risk, tied.time.  Also contains the concordance statistic and its standard error.  Calculated using the \code{survival} package, if installed.  If not installed, all values are \code{NA}. }
    \item{rsquare}{a vector containing an r-square value for the fit and its p-value.}
    \item{lrt}{a vector contaiing the likelihood ratio test statistic and its p-value.}
    \item{df}{the degrees of freedom.}
    \item{wald.test}{a vector containg the Wald test statistic and its p-value.}
    \item{score}{a vector contining the score test statistic and its p-value.}
    \item{iter}{the number of iterations of the cox algorithm before convergence.}
}
\seealso{
  \code{\link{vdracox}}
}
