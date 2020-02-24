\name{rocdistributed.object}
\alias{rocdistributed.object}
\alias{print.rocdistributed}
\title{
    Vertical Distributed Logistic Regression ROC Object
}
\description{
  This class of object is returned by \code{\link{RocTest}}. Objects of this class have a method for the function \code{plot}.
}
\arguments{
    The following components must be included in a legitimate \code{vdracox} object.
    \item{ROC}{a list which conitains two elements \code{roc} and \code{auc}}
    \item{ROC$roc}{a two column matrix containing the cordinates of 1 - specifity and sensitivity.}
    \item{ROC$auc}{numeric value which is area under the curve}
}
\seealso{
  \code{\link{RocTest}}
}
