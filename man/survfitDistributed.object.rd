\name{survfitDistributed.object}
\alias{print.survfitDistributed}
\alias{survfitDistributed.object}
\title{
  Vertical Distributed Cox Regression Survfit Object
}
\description{
  This class of object is created by \code{\link{survfitDistributed}}.  Objects of this class have methods for the functions \code{plot} and \code{print}.  Designed to mimic \code{survfit.object} from the \code{survival} package.
}
\arguments{
    The following components must be included in a legitimate \code{vdracox} object.
    \item{n}{the total number of subjects in each curve}
    \item{time}{the time points at which the curve has a step}
    \item{n.risk}{the number of subjects at risk at each time point}
    \item{n.event}{the number of events that occour at each time point}
    \item{n.censor}{the number of subjects who are censored at each time point}
    \item{strata}{the number of points in each strata}
    \item{surv}{the estimate of the survival time at each time step}
    \item{type}{the type of censoring.  Currently, always "right"}
}
\seealso{
  \code{\link{survfitDistributed}}, \code{\link{plot.survfitDistributed}}
}
