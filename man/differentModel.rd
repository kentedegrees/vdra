\name{differentModel}
\alias{differentModel}
\title{
	Fitting Different Linear Models
}
\description{
  Models are specified symbolically.  A typical model is of the form \code{response ~ term_1 + term_2 + ... + term_k} where \code{response} and \code{term_i} are variables names used in the orgional linear model which created the object \code{x}.  The \code{response} can be the orginal respose or any of the other covariates.  Interactions are not allowed.  Not all variables in the original model have to be used.
}
\usage{
differentModel(formula, x)
}
\arguments{
    \item{formuala}{an object of class \code{"\link{formula}"}: a symbolic description of the model to be fitted. The model must be additive with no interactions.}
    \item{x}{an object of class \code{\link{vdralinear}}.}
}
\value{
    Returns an object of class \code{\link{vdralinear}}.
}
\seealso{
    \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
\examples{
    fit = differentModel(Change_BMI ~ Exposure + Age + NumRx, vdra_fit_linear_A)
    summary(fit)

    fit = differentModel(Age ~ Change_BMI + Exposure + NumRx, vdra_fit_linear_A)
    summary(fit)
}
