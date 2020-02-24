\name{summary.vdralinear}
\alias{summary.vdralinear}
\title{
	Summary Method for Vertical Distributed Linear Regression Models
}
\description{
  Produces a summary of a fitted vdra linear regression model.
}
\usage{
## S3 method for class 'vdralinear'
summary(x)
}
\arguments{
	\item{x}{a \code{vdralinear} object.}
}
\value{
    An object of \code{summary.vdralinear} with components:
    \item{failed}{logical value.  If \code{FALSE}, then there was an error processing the data.  if \code{TRUE}, there were no errors.}
    \item{party}{a vector which indicates the party from which each covariate came.}
    \item{coefficients}{the vector of coefficients.  If the model is over-determined, there will be missing values in the vector corresponding to the redudant columns model matrix.}
    \item{secoef}{the vector of the standard error of the coefficients.}
    \item{tvals}{the t-values of the coefficietns.}
    \item{pvals}{the p-values of the coefficients.}
    \item{rstderr}{residual standard error.}
    \item{rsquare}{r squared.}
    \item{adjrsquare}{adjusted r squared.}
    \item{Fstat}{the F-statistic for the linear regression.}
    \item{df1}{The numerator degrees of freedom for the F-statistic}
    \item{df2}{The denominator degrees of freedom for the F-statistic}
    \item{Fpval}{the p-value of the F-statistic for the linear regression.}
}
\seealso{
  \code{\link{vdralinear}}
}
