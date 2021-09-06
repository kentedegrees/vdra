\name{distributed2party}
\alias{AnalysisCenter.2Party}
\alias{DataPartner.2Party}
\title{
	Two Party Vertical Distributed Regression Analysis
}
\description{
\code{AnalysisCenter.2Party} and \code{DataPartner.2Party} are used in conjuction
    with PopMedNet to perform linear, logistic, or cox regression on data that has
    been partitioned vertically between two data partners.  The data partner
    which holds the response variable(s) uses \code{AnalysisCener.2Party} and the
    other data partner uses \code{DataPartner.2Party}.  While both data partners
    share information with each other in order to perform the regression, data is
    kept secure and not shared, nor is any information shared that would allow
    one data partner to reconstruct part of the other data partners data.
    Final coefficients and other regression statistics are computed by the analysis
    center and shared with the other data partner.
}
\usage{
AnalysisCenter.2Party(regression = "linear", data = NULL, response = NULL,
                      strata = NULL, mask = TRUE, monitorFolder = NULL,
                      msreqid = "v_default_00_000", blocksize = 500,
                      tol = 1e-8, maxIterations = 25, sleepTime = 10,
                      maxWaitingTime = 86400, popmednet = TRUE,
                      trace = FALSE, verbose = TRUE)

DataPartner.2Party(regression = "linear", data = NULL, strata = NULL,
                   mask = TRUE, monitorFolder = NULL, sleepTime = 10,
                   maxWaitingTime = 86400, popmednet = TRUE,
                   trace = FALSE, verbose = TRUE)
}
\arguments{
\item{regression}{the model to be used to fit the data.  The default regression
    \code{"linear"} fits a least squares linear model to the data.
    Alternatively, \code{"logistic"} returns a fitted logistic model, and
    \code{"cox"} returns a fitted Cox proportional hazards model.}
\item{data}{a data.frame or matrix which contains the data to be used in the
    model.  For \code{DataPartner.2Party()}, all columns will be used as
    covariates in the regression.  For \code{AnalysisCenter.2Party()}, all
    columns, with the expection of the column specified by \code{response},
    will be used as covariates in the regression.}
\item{response}{for \code{"linear"} and \code{"logistic"} regression, the name
    of the column in \code{data} which holds the response variable.  If
    \code{reponse = NULL}, then the first column of \code{data} will be used
    as the response variable.  For \code{"cox"} regression response hold the
    name of the column which is time to event and the name of the column which
    is the event type (0 = censored, 1 = event).  If \code{response = NULL},
    then the first column of \code{data} is assumed to be the time to even and
    the second column is assumed to be the event type.}
\item{strata}{for \code{"cox"} regression only.  A \code{\link{vector}} of
    character strings identifying the names of the covariates from either party
    which will be used as strata.  Both \code{AnalysisCenter.2party} and
    \code{DataPartner.2Party} must specify the same vector of strata.}
\item{mask}{logical value: If \code{FALSE}, strata levels for the strata which
    belong to the party which specified \code{FALSE} will be identified by name.
    If \code{TRUE}, levels for the strata which belong to the party which
    specified \code{TRUE} will be put in a random order and level names will be
    changed to \code{NA}.}
\item{monitorFolder}{the folder where the directories \code{dplocal},
    \code{inputfiles}, \code{macros}, \code{msoc}, and \code{rprograms} are
    located.}
\item{msreqid}{a character string specifying the name of the \emph{Request ID}
    as specified when creating the Distributed Regresion request on PopMedNet.
    Used for logging purposes only.}
\item{blocksize}{the minimium size used to horizontally partition the data for
    data transfer between the two parties.}
\item{tol}{the tolerance used to determine convergence in \code{"logistic"}
    and \code{"cox"} regression.}
\item{maxIterations}{the maximum number of iterations to perform
    \code{"logistic"} or \code{"cox"} regression before non-convergence is
    declared.}
\item{sleepTime}{the number of seconds to wait after writing the last
    file to disk before signalling the PMN Datamart Client that files are ready
    to be transferred.}
\item{maxWaitingTime}{the number of seconds to wait to receive
    files before a transfer error is declared and the program halts execution.}
\item{popmednet}{logical value:  if \code{TRUE}, assumes that PopMednet is being
    used to transfer the files and implements PopMedNet specific routines. In
    particular, a 15 second offset terminiation of routines that execute in
    parallel is implemented.}
\item{trace}{logical value: if \code{TRUE} and \code{verbose == TRUE}, prints every function called during execution. Used for debugging.}
\item{verbose}{logical value.  If \code{TRUE}, prints out information to document the progression of the computation.}

}
\value{
Returns an object of \code{\link{class}} \code{\link{vdralinear}} for linear
    regression, \code{\link{vdralogistic}} for logistic regression, or
    \code{\link{vdracox}} for cox regression.
}
\seealso{
    \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
}
\examples{
\dontrun{
## 2 party linear regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.

fit = AnalysisCenter.2Party(regression = "linear", data = vdra_data[, c(1, 5:7)],
        response = "Change_BMI", monitorFolder = tempdir())

# Data Partner -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.

fit = DataPartner.2Party(regression = "linear", data = vdra_data[, 8:11],
        monitorFolder = tempdir())

## 2 party logistic regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.

fit = AnalysisCenter.2Party(regression = "logistic", data = vdra_data[, c(2, 5:7)],
        response = "WtLost", monitorFolder = tempdir())

# Data Partner -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.

fit = DataPartner.2Party(regression = "logistic", data = vdra_data[, 8:11],
        monitorFolder = tempdir())

## 2 party cox regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.

fit = AnalysisCenter.2Party(regression = "cox", data = vdra_data[, c(3:4, 5:7)],
        response = c("Time", "Status"), strata = c("Exposure", "Sex"),
        monitorFolder = tempdir())

# Data Partner -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.

fit = DataPartner.2Party(regression = "cox", data = vdra_data[, 8:11],
        strata = c("Exposure", "Sex"), monitorFolder = tempdir())
}
}
