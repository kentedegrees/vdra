\name{vdra_fit}
\alias{vdra_fit_cox_A}
\alias{vdra_fit_cox_B}
\alias{vdra_fit_linear_A}
\alias{vdra_fit_linear_B}
\alias{vdra_fit_logistic_A}
\alias{vdra_fit_logistic_B}
\title{
    Return values from the various distributed regression functions in the \code{vdra} package.
}
\description{
    The objects vdra_fit_cox_A, vdra_fit_linear_A, and vdra_logistic_A are example fitted models that are obtained by the data partner which holds the response variable(s). The objects vdra_fit_cox_B, vdra_fit_linear_B, and vdra_logistic_B are example fitted models that are obtained by the data partner which does not hold the response variable(s).  These are provided so the user may see what the summary and print outputs look like before trying to run the full vertical distributed regression.  They also allow the user to experiement with the functions  \code{\link{differentModel}} (linear regression); \code{\link{RocTest}} and \code{\link{HoslemTest}} (Logistic Regression); and \code{\link{survfitDistributed}} and \code{\link{plot.survfitDistributed}} (Cox Regression).
}
\usage{
    vdra_fit_cox_A
    vdra_fit_cox_B
    vdra_fit_linear_A
    vdra_fit_linear_B
    vdra_fit_logistic_A
    vdra_fit_logistic_B
}
