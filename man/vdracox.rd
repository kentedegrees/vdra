\name{vdracox}
\alias{vdracox}
\alias{vdracox.object}
\alias{print.vdracox}
\title{
    Vertical Distributed Cox Regression Results Object
}
\description{
    This class of object is returned by the two party, three party, and K-party distributed regression analysis programs when "cox" regression is specified.  Objects of this class have methods for the functions \code{print} and \code{summary}.
}
\arguments{
    The following components must be included in a legitimate \code{vdracox} object.
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{converged}{logical value.  If \code{TRUE}, the regression converged.  If \code{FALSE}, it did not.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be \code{NA} values in the vector corresponding to the redudant columns model matrix.}
    \item{expcoef}{a vector which represents exp(coefficients).}
    \item{expncoef}{a vector which represents exp(-coefficients).}
    \item{var}{the variance matrix of the coefficients. Rows and columns corresponding to any missing coefficients are set to zero.}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{zvals}{the z-values of the coefficients.}
    \item{pvals}{the p-values of the coefficients.}
    \item{lower95}{a vector of the lower bounds of the 95\% confidence interval for exp(coefficients).}
    \item{upper95}{a vector of the upper bounds of the 95\% confidence interval for exp(coefficients).}
    \item{loglik}{a vector holding the loglikelihood and null loglikelihood.}
    \item{n}{the number of observations in the data.}
    \item{nevent}{the number of events used in the fit.}
    \item{iter}{the number of iterations of the cox algorithm before convergence.}
    \item{df}{the degrees of freedom.}
    \item{score}{a vector contining the score test statistic and its p-value.}
    \item{method}{"efron". The method used for the cox algorithm.}
    \item{lrt}{a vector contaiing the likelihood ratio test statistic and its p-value.}
    \item{rsquare}{a vector containing an r-square value for the fit and its p-value.}
    \item{wald.test}{a vector containg the Wald test statistic and its p-value.}
    \item{concordance}{a vector containing the number of events which are concordant, discordant, tied.risk, tied.time.  Also contains the concordance statistic and its standard error.  Calculated using the \code{survival} package, if installed.  If not installed, all values are \code{NA}. }
    \item{survival}{a matrix of values used to compute the surivival curve.}
    \item{strata}{a data.frame of the strata used in the computation.}
}
\seealso{
\code{\link{survfitDistributed}}, \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
