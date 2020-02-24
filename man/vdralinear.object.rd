\name{vdralinear}
\alias{vdralinear}
\alias{print.vdralinear}
\title{
    Vertical Distributed Linear Regression Results Object
}
\description{
    This class of object is returned by the two party, three party, and K-party distributed regression analysis programs when "linear" regression is specified.  Objects of this class have methods for the functions \code{print} and  \code{summary}.
}
\arguments{
    The following components must be included in a legitimate \code{vdralinear} object.
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{converged}{logical value.  If \code{TRUE}, the regression converged.  If \code{FALSE}, it did not.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be missing values in the vector corresponding to the redudant columns model matrix.}
    \item{tvals}{the t-values of the coefficietns.}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{pvals}{the p-values of the coefficients.}
    \item{sse}{sum of squared errors.}
    \item{rstderr}{residual standard error.}
    \item{rsquare}{r squared.}
    \item{adjrsquare}{adjusted r squared.}
    \item{Fstat}{the F-statistic for the linear regression.}
    \item{Fpval}{the p-value of the F-statistic for the linear regression.}
    \item{df1}{The numerator degrees of freedom for the F-statistic}
    \item{df2}{The denominator degrees of freedom for the F-statistic}
    \item{n}{the number of observations in the data.}
    \item{xtx}{a matrix of the transpose of the covariates times the covarites.  Used by \code{\link{differentModel}}.}
    \item{xty}{a matrix of the transpose of the covarites times the response.  Used by \code{\link{differentModel}}.}
    \item{yty}{sum of squares of the reponse.  Used by \code{\link{differentModel}}.}
    \item{meansy}{the mean of the response.  Used by \code{\link{differentModel}}.}
    \item{means}{the mean of each covaraite.  Used by \code{\link{differentModel}}.}
}
\seealso{
    \code{\link{differentModel}}, \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
