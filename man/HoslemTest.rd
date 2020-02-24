\name{HoslemTest}
\alias{HoslemTest}
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
\examples{
  HoslemTest(vdra_fit_logistic_A)

  HoslemTest(vdra_fit_logistic_A, 20)
}
