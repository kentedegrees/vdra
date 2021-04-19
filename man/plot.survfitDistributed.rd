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
	\item{...}{common graphical parameters (not fully implemented.}
}
\seealso{
  \code{\link{survfitDistributed}}
}
\examples{
  sfit = survfitDistributed(vdra_fit_cox_A)
  plot(sfit)

  # From Data Partner 1
  sfit = survfitDistributed(vdra_fit_cox_A, ~Exposure, data = vdra_data[, c(3:4, 5:7)])
  plot(sfit)
  plot(sfit, merge = FALSE)

  # From Data Partner 2
  sfit = survfitDistributed(vdra_fit_cox_B, ~Race + Sex, data = vdra_data[, 8:11])
  plot(sfit, merge = FALSE)
}
