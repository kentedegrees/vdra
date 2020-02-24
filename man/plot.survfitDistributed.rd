\name{plot.survfitDistributed}
\alias{plot.survfitDistributed}
\title{
	Plotting Survival Curves for Vertical Distributed Cox Regression
}
\description{
  Plots a survivial curve as specified by \code{survfitDistributed} object.
}
\usage{
plot.survfitDistributed(x, merge = TRUE, ...)
}
\arguments{
	\item{x}{a \code{survfitDistributed} object.}
	\item{merge}{logical.  It \code{TRUE}, plots all strata of the survival curve on one plot.  If \code{FALSE}, plots all strata in different plots.}
	\item{...}{common graphical parameters.}
}
\seealso{
  \code{\link{survfitDistributed.object}}
}
\examples{
  sfit = survfitDistributed(vdra_fit_cox_A)
  plot(sfit)
  plot(sfit, merge = FALSE)
}
