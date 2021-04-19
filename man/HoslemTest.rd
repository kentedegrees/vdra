\name{HoslemTest}
\alias{HoslemTest}
\alias{hoslemdistributed.object}
\alias{print.hoslemdistributed}
\title{
	Hosmer-Lemeshow Test for Vertical Distributed Logistic Regression
}
\description{
  Run the Hosmer-Lemeshow test for an object created by 2-party, 3-party, or K-party vdra logistic regression.  Only the party that holds the response may invoke this function.
}
\usage{
  HoslemTest(x, nGroups = 10)
}
\arguments{
	\item{x}{an object of type \code{\link{vdralogistic}}.}
	\item{nGroups}{the number of groups that the data will be sperated into.}
}
\value{
    Returns an object of class \code{hoslemdistributed}. Objects of this class have a method for the function \code{print}.   The following component must be included in a \code{hoslemdistributed} object.
  \item{hoslem}{a vector containing three numeric quantities: the chi-square value of the test, the degrees of freedom of the test, and p-value of the test, in that order.}
}

\examples{
  HoslemTest(vdra_fit_logistic_A)

  HoslemTest(vdra_fit_logistic_A, 20)
}
