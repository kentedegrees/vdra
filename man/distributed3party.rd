\name{distributed3party}
\alias{AnalysisCenter.3Party}
\alias{DataPartner1.3Party}
\alias{DataPartner2.3Party}
\title{
	Three Party Vertical Distributed Regression Analysis
}
\description{
\code{AnalysisCenter.3Party}, \code{DataPartner1.3Party} and
    \code{DataPartner2.3Party} are used in conjuction
    with PopMedNet to perform linear, logistic, or cox regression on data that has
    been partitioned vertically between two data partners.  The data partner
    which holds the response variable(s) uses \code{Datapartner1.3Party} and the
    other data partner uses \code{DataPartner2.3Party}.  Data partners are not
    allowed to communicate with each other, but share inforamtion via a
    trusted third party analysis center.  While any infomration that is shared
    with the analysis center by a data partner, with the exception of some
    summary statistics, is encrypted by the sending data parter, if the infomration
    needs to be sent on to the other data parter for futher analysis, the
    analysis center further encrypts the data.  That way, any information
    that deals directly with the raw data that moves between two data partners
    is doubly encrypted to keep both the analysis center and the other
    data partner from learning it.  Thus, no information is shared between
    the data partners or analysis center that would allow one data partner
    to reconstrut part of the other data partners data.  Final coefficients and
    other regression statistics are computed by the analysis center and shared
    with the data partners.
}
\usage{
AnalysisCenter.3Party(regression = "linear", monitorFolder = NULL,
                      msreqid = "v_default_00_000", blocksize = 500,
                      tol = 1e-8, maxIterations = 25, sleepTime = 10,
                      maxWaitingTime = 86400, popmednet = TRUE, trace = FALSE)

DataPartner1.3Party(regression = "linear", data = NULL, response = NULL,
                    strata = NULL, mask = TRUE, monitorFolder = NULL,
                    sleepTime = 10, maxWaitingTime = 86400, popmednet = TRUE,
                    trace = FALSE)

DataPartner2.3Party(regression = "linear", data = NULL, strata = NULL,
                    mask = TRUE, monitorFolder = NULL, sleepTime = 10,
                    maxWaitingTime = 86400, popmednet = TRUE, trace = FALSE)
}
\arguments{
\item{regression}{the model to be used to fit the data.  The default regression
    \code{"linear"} fits a least squares linear model to the data.
    Alternatively, \code{"logistic"} returns a fitted logistic model, and
    \code{"cox"} returns a fitted Cox proportional hazards model.}
\item{data}{a data.frame or matrix which contains the data to be used in the
    model.  For \code{DataPartner2.3Party()}, all columns will be used as
    covariates in the regression.  For \code{DataPartner1.3Party()}, all
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
    which will be used as strata.  Both \code{DataPartner1.3party} and
    \code{DataPartner2.3Party} must specify the same vector of strata.}
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
    files before a transfer error is declared and the program halts execution.
    Should be the same for both parties when \code{delayOffset = TRUE}.}
\item{popmednet}{logical value:  if \code{TRUE}, assumes that PopMednet is being
    used to transfer the files and implements PopMedNet specific routines. In
    particular, a 15 second offset between terminiation of routines that execute in
    parallel is implemented.}
\item{trace}{logical value: if code{TRUE}, prints every function call. Used for debugging.}
}
\value{
Returns an object of \code{\link{class}} \code{\link{vdralinear}} for linear
    regression, \code{\link{vdralogistic}} for logistic regression, or
    \code{\link{vdracox}} for cox regression.
}
\seealso{
    \code{\link{AnalysisCenter.2Party}}, \code{\link{AnalysisCenter.KParty}}
}
\examples{
\dontrun{
## 3 party linear regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.
fit = AnalysisCenter.3Party(regression = "linear", monitorFolder = getwd())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner1.3Party(regression = "linear", data = vdra_data[, c(1, 5:7)],
          response = "Change_BMI", monitorFolder = getwd())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner2.3Party(regression = "linear", data = vdra_data[, 8:11],
          monitorFolder = getwd())

## 3 party logistic regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.
fit = AnalysisCenter.3Party(regression = "logistic", monitorFolder = getwd())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner1.3Party(regression = "logistic", data = vdra_data[, c(2, 5:7)],
          response = "WtLost", monitorFolder = getwd())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner2.3Party(regression = "logistic", data = vdra_data[, 8:11],
          monitorFolder = getwd())

## 3 party cox regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.
fit = AnalysisCenter.3Party(regression = "cox", monitorFolder = getwd())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner1.3Party(regression = "cox", data = vdra_data[, c(3:4, 5:7)],
        response = c("Time", "Status"), strata = c("Exposure", "Sex"),
        monitorFolder = getwd())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner2.3Party(regression = "cox", data = vdra_data[, 8:11],
         strata = c("Exposure", "Sex"), monitorFolder = getwd())
}
}
