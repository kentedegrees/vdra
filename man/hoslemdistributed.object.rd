\name{hoslemdistributed.object}
\alias{hoslemdistributed.object}
\alias{print.hoslemdistributed}
\title{
    Vertical Distributed Logistic Regression Hosmer-Lemeshow Object
}
\description{
  This class of object is returned by \code{\link{HoslemTest}}. Objects of this class have a method for the function \code{print}.
}
\arguments{
  The following component must be included in a legitimate \code{vdracox} object.
  \item{hoslem}{a vector containing three numeric quantities: the chi-square value of the test, the degrees of freedom of the test, and p-value of the test, in that order.}
}
\seealso{
  \code{\link{HoslemTest}}
}
