\name{RocTest}
\alias{RocTest}
\alias{rocdistributed.object}
\alias{print.rocdistributed}
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
\value{
    Returns an object of class \code{rocdistributed}. Objects of this class have a method for the function \code{print}. The following components must be included in a \code{rocdistributed} object.
    \item{roc}{a two column matrix containing the cordinates of 1 - specifity and sensitivity.}
    \item{auc}{numeric value which is area under the curve.}
}


\examples{
  RocTest(vdra_fit_logistic_A)

  RocTest(vdra_fit_logistic_A, 40)
}
