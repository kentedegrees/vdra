#' @name vdra_fit
#' @aliases vdra_fit_cox_A vdra_fit_cox_B vdra_fit_linear_A vdra_fit_linear_B
#'   vdra_fit_logistic_A vdra_fit_logistic_B
#' @title Return values from the various distributed regression functions in the
#'   \code{vdra} package.
#' @description The objects vdra_fit_cox_A, vdra_fit_linear_A, and
#'   vdra_logistic_A are example fitted models that are obtained by the data
#'   partner which holds the response variable(s). The objects vdra_fit_cox_B,
#'   vdra_fit_linear_B, and vdra_logistic_B are example fitted models that are
#'   obtained by the data partner which does not hold the response variable(s).
#'   These are provided so the user may see what the summary and print outputs
#'   look like before trying to run the full vertical distributed regression.
#'   They also allow the user to experiement with the functions
#'   \code{\link{differentModel}} (linear regression); \code{\link{RocTest}} and
#'   \code{\link{HoslemTest}} (Logistic Regression); and
#'   \code{\link{survfitDistributed}} and \code{\link{plot.survfitDistributed}}
#'   (Cox Regression).
NULL

#' @name vdralinear
#' @aliases vdralinear vdralinear.object print.vdralinear
#' @title Vertical Distributed Linear Regression Results Object
#' @description This class of object is returned by the two party, three party,
#'   and K-party distributed regression analysis programs when "linear"
#'   regression is specified.  Objects of this class have methods for the
#'   functions \code{print} and  \code{summary}.
#' @format The following components must be included in a legitimate
#'   \code{vdralinear} object. \describe{ \item{failed}{logical value.  If
#'   \code{FALSE}, then there was an error processing the data.  if \code{TRUE},
#'   there were no errors.}
#'
#'   \item{converged}{logical value.  If \code{TRUE}, the regression converged.
#'   If \code{FALSE}, it did not.}
#'
#'   \item{party}{a vector which indicates the party from which each covariate
#'   came.}
#'
#'   \item{coefficients}{the vector of coefficients.  If the model is
#'   over-determined, there will be \code{NA} values in the vector corresponding
#'   to the redudant columns model matrix.}
#'
#'   \item{tvals}{the t-values of the coefficietns.}
#'
#'   \item{secoef}{the vector of the standard error of the coefficients.}
#'
#'   \item{pvals}{the p-values of the coefficients.}
#'
#'   \item{sse}{sum of squared errors.}
#'
#'   \item{rstderr}{residual standard error.}
#'
#'   \item{rsquare}{r squared.}
#'
#'   \item{adjrsquare}{adjusted r squared.}
#'
#'   \item{Fstat}{the F-statistic for the linear regression.}
#'
#'   \item{Fpval}{the p-value of the F-statistic for the linear regression.}
#'
#'   \item{df1}{The numerator degrees of freedom for the F-statistic.}
#'
#'   \item{df2}{The denominator degrees of freedom for the F-statistic.}
#'
#'   \item{n}{the number of observations in the data.}
#'
#'   \item{xtx}{a matrix of the transpose of the covariates times the covarites.
#'   Used by \code{\link{differentModel}}.}
#'
#'   \item{xty}{a matrix of the transpose of the covarites times the response.
#'   Used by \code{\link{differentModel}}.}
#'
#'   \item{yty}{sum of squares of the reponse.  Used by
#'   \code{\link{differentModel}}.}
#'
#'   \item{meansy}{the mean of the response.  Used by
#'   \code{\link{differentModel}}.}
#'
#'   \item{means}{the mean of each covaraite.  Used by
#'   \code{\link{differentModel}}.}
#'
#'   }
#' @seealso \code{\link{differentModel}} \code{\link{AnalysisCenter.2Party}}
#'   \code{\link{AnalysisCenter.3Party}} \code{\link{AnalysisCenter.KParty}}
NULL

#' @name vdralogistic
#' @aliases vdralogistic vdralogistic.object print.vdralogistic
#' @title Vertical Distributed Logistic Regression Results Object
#' @description This class of object is returned by the two party, three party,
#'   and K-party distributed regression analysis programs when "logistic"
#'   regression is specified.  Objects of this class have methods for the
#'   functions \code{print} and  \code{summary}.
#' @format The following components must be included in a legitimate
#'   \code{vdralogistic} object. \describe{ \item{failed}{logical value.  If
#'   \code{FALSE}, then there was an error processing the data.  if \code{TRUE},
#'   there were no errors.}
#'
#'   \item{converged}{logical value.  If \code{TRUE}, the regression converged.
#'   If \code{FALSE}, it did not.}
#'
#'   \item{party}{a vector which indicates the party from which each covariate
#'   came.}
#'
#'   \item{coefficients}{the vector of coefficients.  If the model is
#'   over-determined, there will be \code{NA} values in the vector corresponding
#'   to the redudant columns model matrix.}
#'
#'   \item{secoef}{the vector of the standard error of the coefficients.}
#'
#'   \item{tvals}{the t-values of the coefficietns.}
#'
#'   \item{pvals}{the p-values of the coefficients.}
#'
#'   \item{n}{the number of observations in the data.}
#'
#'   \item{nulldev}{the null deviance of the fit.}
#'
#'   \item{resdev}{the residual deviance of the fit.}
#'
#'   \item{aic}{the AIC of the fit.}
#'
#'   \item{bic}{the BIC of the fit.}
#'
#'   \item{nulldev_df}{the degrees of freedom for the null deviance.}
#'
#'   \item{resdev_df}{the degrees of freedome for the residual deviance.}
#'
#'   \item{hoslem}{the Hosmer Lemshow Test statistics.}
#'
#'   \item{ROC}{a list containing the coordinates for an ROC curve.}
#'
#'   \item{iter}{the number of iterations of the cox algorithm before
#'   convergence.}
#'
#'   \item{Y}{a matrix of the response.  Only returned to the party which holds
#'   the response.}
#'
#'   \item{FinalFitted}{a matrix of final fitted values of the regression.  Only
#'   returned to the party which holds the response.}
#'
#'   }
#' @seealso \code{\link{HoslemTest}} \code{\link{RocTest}}
#'   \code{\link{AnalysisCenter.2Party}} \code{\link{AnalysisCenter.3Party}},
#'   \code{\link{AnalysisCenter.KParty}}
NULL

#' @name vdracox
#' @aliases vdracox vdracox.object print.vdracox
#' @title Vertical Distributed Cox Regression Results Object
#' @description This class of object is returned by the two party, three party,
#'   and K-party distributed regression analysis programs when "cox" regression
#'   is specified.  Objects of this class have methods for the functions
#'   \code{print} and \code{summary}.
#' @format The following components must be included in a legitimate
#'   \code{vdracox} object. \describe{ \item{failed}{logical value.  If
#'   \code{FALSE}, then there was an error processing the data.  if \code{TRUE},
#'   there were no errors.}
#'
#'   \item{converged}{logical value.  If \code{TRUE}, the regression converged.
#'   If \code{FALSE}, it did not.}
#'
#'   \item{party}{a vector which indicates the party from which each covariate
#'   came.}
#'
#'   \item{coefficients}{the vector of coefficients.  If the model is
#'   over-determined, there will be \code{NA} values in the vector corresponding
#'   to the redudant columns model matrix.}
#'
#'   \item{expcoef}{a vector which represents exp(coefficients).}
#'
#'   \item{expncoef}{a vector which represents exp(-coefficients).}
#'
#'   \item{var}{the variance matrix of the coefficients. Rows and columns
#'   corresponding to any missing coefficients are set to zero.}
#'
#'   \item{secoef}{the vector of the standard error of the coefficients.}
#'
#'   \item{zvals}{the z-values of the coefficients.}
#'
#'   \item{pvals}{the p-values of the coefficients.}
#'
#'   \item{lower95}{a vector of the lower bounds of the 95\% confidence interval
#'   for exp(coefficients).}
#'
#'   \item{upper95}{a vector of the upper bounds of the 95\% confidence interval
#'   for exp(coefficients).}
#'
#'   \item{loglik}{a vector holding the loglikelihood and null loglikelihood.}
#'
#'   \item{n}{the number of observations in the data.}
#'
#'   \item{nevent}{the number of events used in the fit.}
#'
#'   \item{iter}{the number of iterations of the cox algorithm before
#'   convergence.}
#'
#'   \item{df}{the degrees of freedom.}
#'
#'   \item{score}{a vector contining the score test statistic and its p-value.}
#'
#'   \item{method}{"efron". The method used for the cox algorithm.}
#'
#'   \item{lrt}{a vector containing the likelihood ratio test statistic and its
#'   p-value.}
#'
#'   \item{rsquare}{a vector containing an r-square value for the fit and its
#'   p-value.}
#'
#'   \item{wald.test}{a vector containing the Wald test statistic and its
#'   p-value.}
#'
#'   \item{concordance}{a vector containing the number of events which are
#'   concordant, discordant, tied.risk, tied.time.  Also contains the
#'   concordance statistic and its standard error.  Calculated using the
#'   \code{survival} package, if installed.  If not installed, all values are
#'   \code{NA}. }
#'
#'   \item{survival}{a matrix of values used to compute the surivival curve.}
#'
#'   \item{strata}{a data.frame of the strata used in the computation.} }
#' @seealso \code{\link{survfitDistributed}} \code{\link{AnalysisCenter.2Party}}
#'   \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
NULL
