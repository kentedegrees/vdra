\name{survfitDistributed}
\alias{survfitDistributed}
\alias{survfitDistributed.object}
\alias{print.survfitDistributed}
\title{
	Create Survival Curves for Vertical Distributed Cox Regression
}
\description{
  This function creates survival curves for a previously defined \code{\link{vdracox}} object.  The function also accepts a formula and the original data supplied by the calling party allowing exploration of other potential strata.  Both \code{formula} and \code{data} must be NULL or both must be specified.
}
\usage{
  survfitDistributed(x, formula = NULL, data = NULL)
}
\arguments{
	\item{x}{an object of type \code{\link{vdracox}}.}
	\item{formula}{a formula which defines alternative strata for the survival curve.}
	\item{data}{if \code{formula} is specified, this should be the data that was supplied by the calling party.}
}
\value{
    Returns an object of class \code{survfitDistributed}. Objects of this class have methods for the functions \code{print} and \code{plot}. The following components must be included in a legitimate \code{survfitDistributed} object.
    \item{n}{the total number of subjects in each curve.}
    \item{time}{the time points at which the curve has a step.}
    \item{n.risk}{the number of subjects at risk at each time point.}
    \item{n.event}{the number of events that occour at each time point.}
    \item{n.censor}{the number of subjects who are censored at each time point.}
    \item{strata}{the number of points in each strata.}
    \item{surv}{the estimate of the survival time at each time step.}
    \item{type}{the type of censoring.  Currently, always "right".}
}

\seealso{
  \code{\link{plot.survfitDistributed}}
}
\examples{
    sfit = survfitDistributed(vdra_fit_cox_A)
    print(sfit)
    plot(sfit)

    # From Data Partner 1
    sfit = survfitDistributed(vdra_fit_cox_A, ~Exposure, data = vdra_data[, c(3:4, 5:7)])
    print(sfit)
    plot(sfit)

    # From Data Partner 2
    sfit = survfitDistributed(vdra_fit_cox_B, ~Race + Sex, data = vdra_data[, 8:11])
    print(sfit)
    plot(sfit, merge = TRUE)
}
