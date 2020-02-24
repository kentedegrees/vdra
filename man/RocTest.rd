\name{RocTest}
\alias{RocTest}
\title{
	Create the ROC for Vertical Distributed Logistic Regression
}
\description{
  Generate the receiver operator curve on an object created by 2-party, 3-party, or K-party vdra logistic regression.  Only the party that holds the response may invoke this function.
}

\usage{
  RocTest(x, bins = 10)
}
\arguments{
    \item{x}{an object of type \code{\link{vdralogistic}}.}
	\item{bins}{the number of bins the data will be separated into.}
}
\examples{
  RocTest(vdra_fit_logistic_A)

  RocTest(vdra_fit_logistic_A, 20)
}
