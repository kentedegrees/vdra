\name{vdralogistic}
\alias{vdralogistic}
\alias{vdralogistic.object}
\alias{print.vdralogistic}
\title{
    Vertical Distributed Logistic Regression Results Object
}
\description{
    This class of object is returned by the two party, three party, and K-party distributed regression analysis programs when "logistic" regression is specified.  Objects of this class have methods for the functions \code{print} and  \code{summary}.
}
\arguments{
    The following components must be included in a legitimate \code{vdralogistic} object.
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{converged}{logical value.  If \code{TRUE}, the regression converged.  If \code{FALSE}, it did not.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be \code{NA} values in the vector corresponding to the redudant columns model matrix.}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{tvals}{the t-values of the coefficietns.}
    \item{pvals}{the p-values of the coefficients.}
    \item{n}{the number of observations in the data.}
    \item{nulldev}{the null deviance of the fit.}
    \item{resdev}{the residual deviance of the fit.}
    \item{aic}{the AIC of the fit.}
    \item{bic}{the BIC of the fit.}
    \item{nulldev_df}{the degrees of freedom for the null deviance.}
    \item{resdev_df}{the degrees of freedome for the residual deviance.}
    \item{hoslem}{the Hosmer Lemshow Test statistics.}
    \item{ROC}{a list containing the coordinates for an ROC curve.}
    \item{iter}{the number of iterations of the cox algorithm before convergence.}
    \item{Y}{a matrix of the response.  Only returned to the party which holds the response.}
    \item{FinalFitted}{a matrix of final fitted values of the regression.  Only returned to the party which holds the response.}
}
\seealso{
\code{\link{HoslemTest}}, \code{\link{RocTest}}, \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
