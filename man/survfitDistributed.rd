\name{survfitDistributed}
\alias{survfitDistributed}
\title{
	Create Survival Curves for Vertical Distributed Cox Regression
}
\description{
  This function creates survival curves for a previously defined \code{\link{vdracox}} object.  The survival curve is returned as a \code{\link{survfitDistributed.object}} object.  The function also accepts a formula and the original data supplied by the party allowing exploration of other potential strata.  Both \code{formula} and \code{data} must be NULL or both must be specified.
}
\usage{
  survfitDistributed(x, formula = NULL, data = NULL)
}
\arguments{
	\item{x}{a \code{\link{vdracox}}} object.
	\item{formula}{a formula to which defines alternative strata for the survival curve}
	\item{data}{if formual is specified, this should be the data that was supplied by calling party}
}
\seealso{
  \code{\link{survfitDistributed.object}}, \code{\link{plot.survfitDistributed}}
}
\examples{
  sfit = survfitDistributed(vdra_fit_cox_A)
  print(sfit)
  plot(sfit)

  # From Data Partner 1
  sfit = survfitDistributed(vdra_fit_cox_A, Exposure, data = vdra_data[, 3:20])

  # From Data Partner 2
  sfit = survfitDistributed(vdra_fit_cox_B, Race + Sex, data = vdra_data[, 21:41])
}
