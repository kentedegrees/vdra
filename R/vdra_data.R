#' This is data to be included in my package
#'
#' @name vdra_data
#' @title Simulated data from a weight loss study.
#' @docType data
#' @description Simulated data based on a weight loss study performed by Harvard
#'   School of Medicine.  The original study was looking for comorbidities to
#'   predict weight loss / weight gain in obese subjects.  The various outcomes
#'   in this simulated data set don't really correlate to each other.  The
#'   purpose of this simulated data set is to illustrate the funcationality of
#'   the package, not to draw valid statistical inferences.
#' @format
#'
#' \describe{
#'
#' \item{Change_BMI}{Continuous response used for linear
#' regression.}
#'
#' \item{WtLost}{Binary response used for logistic regression.}
#'
#' \item{Time}{Used for time to event in Cox regression.  Continuous.  Ranges
#' from 1 to 459.}
#'
#' \item{Status}{Used for censoring in Cox regression.  A binary categroical
#' variable.}
#'
#' \item{Exposure}{A binary categorical variable.}
#'
#' \item{Age}{A continuous variable ranging from 3 to 80.}
#'
#' \item{ComorbidScore}{A ordinal variable with 11 variables: 0 to 10.}
#'
#' \item{NumRx}{A ordinal varible with 15 levels:  -2 to 12.}
#'
#' \item{BMI_pre}{A continuous variable ranging from 35.01 to 92.79.}
#'
#' \item{Race}{A factor with 6 levels: "Race 0" to "Race 5"..}
#'
#' \item{Sex}{A binary factor with 2 levels: M and F.}
#'
#' }
#' @keywords data vdra_data
NULL

#' vdra: A package for performing vertical distributed regression analysis.
#'
#' All the functions in this package are designed to be run in a distributed
#' network using popMedNet.
#'
#'
#' @docType package
#' @name vdra
#' @useDynLib vdra, .registration=TRUE
NULL
