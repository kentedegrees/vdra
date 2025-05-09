\name{distributedKparty}
\alias{AnalysisCenter.KParty}
\alias{DataPartner.KParty}
\title{
	K-Party Vertical Distributed Regression Analysis
}
\description{
\code{AnalysisCenter.KParty} and \code{DataPartner.KParty} are used in conjuction
    with PopMedNet to perform linear, logistic, or cox regression on data that has
    been partitioned vertically between two or more data partners.  The data partners which
    holds the data use \code{DataPartner.KParty} while a trusted "third" party uses
    \code{AnalysisCenter.KParty}.  Data partners are
    allowed to communicate with each other and the analysis center, no information is
    shared between the data partners or analysis center that would allow one data partner
    or the analysis center
    to reconstrut part of the other data partners data.  Final coefficients and
    other regression statistics are computed by the analysis center and shared
    with the data partners.
}
\usage{
AnalysisCenter.KParty(regression = "linear", numDataPartners = NULL,
                      monitorFolder = NULL, msreqid = "v_default_00_000",
                      tol = 1e-8, maxIterations = 25, sleepTime = 10,
                      maxWaitingTime = 86400, popmednet = TRUE,
                      trace = FALSE, verbose = TRUE)

DataPartner.KParty(regression = "linear", data = NULL, response = NULL,
                   strata = NULL, mask = TRUE, numDataPartners = NULL,
                   dataPartnerID = NULL, monitorFolder = NULL,
                   sleepTime = 10, maxWaitingTime = 86400, popmednet = TRUE,
                   trace = FALSE, verbose = TRUE)
}
\arguments{
\item{regression}{the model to be used to fit the data.  The default regression
    \code{"linear"} fits a least squares linear model to the data.
    Alternatively, \code{"logistic"} returns a fitted logistic model, and
    \code{"cox"} returns a fitted Cox proportional hazards model.}
\item{data}{a data.frame or matrix which contains the data to be used in the
    model.  All columns will be used as
    covariates in the regression with the exception of the data partner which
    has \code{dataPartnerID = 1}.  For this data partner, all
    columns, with the expection of the column specified by \code{response},
    will be used as covariates in the regression.}
\item{response}{only used for data parther with \code{dataPartnerID = 1}.
    For \code{"linear"} and \code{"logistic"} regression, the name
    of the column in \code{data} which holds the response variable.  If
    \code{reponse = NULL}, then the first column of \code{data} will be used
    as the response variable.  For \code{"cox"} regression response hold the
    name of the column which is time to event and the name of the column which
    is the event type (0 = censored, 1 = event).  If \code{response = NULL},
    then the first column of \code{data} is assumed to be the time to even and
    the second column is assumed to be the event type.}
\item{strata}{for \code{"cox"} regression only.  A \code{\link{vector}} of
    character strings identifying the names of the covariates from either party
    which will be used as strata.  All data partners must specify the same
    vector of strata.}
\item{mask}{logical value: If \code{FALSE}, strata levels for the strata which
    belong to the party which specified \code{FALSE} will be identified by name.
    If \code{TRUE}, levels for the strata which belong to the party which
    specified \code{TRUE} will be put in a random order and level names will be
    changed to \code{NA}.}
\item{numDataPartners}{the number of data partners which are supplying data for
    the regression.}
\item{dataPartnerID}{a unique identifier for each data partner.  The data partner
    with the response variable(s) must have \code{dataPartnerID = 1}.  All other
    data partners must have an integer value from 2 to \code{numDataPartners}.}
\item{monitorFolder}{the folder where the directories \code{dplocal},
    \code{inputfiles}, \code{macros}, \code{msoc}, and \code{rprograms} are
    located.}
\item{msreqid}{a character string specifying the name of the \emph{Request ID}
    as specified when creating the Distributed Regresion request on PopMedNet.
    Used for logging purposes only.}
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
    Should be the same for all parties when \code{delayOffset = TRUE}.}
\item{popmednet}{logical value:  if \code{TRUE}, assumes that PopMednet is being
    used to transfer the files and implements PopMedNet specific routines. In
    particular, a 15 second offset between terminiation of routines that execute in
    parallel is implemented.}
\item{trace}{logical value: if \code{TRUE} and \code{verbose == TRUE}, prints every function call. Used for debugging.}
\item{verbose}{logical value.  If \code{TRUE}, prints out information to document the progression of the computation.}
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
fit = AnalysisCenter.KParty(regression = "linear", numDataPartners = 2,
              monitorFolder = tempdir())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "linear", data = vdra_data[, c(1, 5:7)],
          response = "Change_BMI", numDataPartners = 2, dataPartnerID = 1,
          monitorFolder = tempdir())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "linear", data = vdra_data[, 8:11],
          numDataPartners = 2, dataPartnerID = 2, monitorFolder = tempdir())

## 3 party logistic regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.
fit = AnalysisCenter.KParty(regression = "logistic", numDataPartners = 2,
              monitorFolder = tempdir())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "logistic", data = vdra_data[, c(2, 5:7)],
          response = "WtLost", numDataPartners = 2, dataPartnerID = 1,
          monitorFolder = tempdir())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "logistic", data = vdra_data[, 8:11],
          numDataPartners = 2, dataPartnerID = 2, monitorFolder = tempdir())

## 3 party cox regression

# Analysis Center -- To be run in one instance of R.
# The working directory should be the same as specified in the PopMedNet
# requset for the analysis center.
fit = AnalysisCenter.KParty(regression = "cox", numDataPartners = 2,
              monitorFolder = tempdir())

# Data Partner 1 -- To be run in second instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "cox", data = vdra_data[, c(3:4, 5:7)],
        response = c("Time", "Status"), strata = c("Exposure", "Sex"),
        numDataPartners = 2, dataPartnerID = 1, monitorFolder = tempdir())

# Data Partner 2 -- To be run in third instand of R, on perhaps a different machine.
# The working directory should be the same as specified in the PopMedNet
# request for the data partner.
fit = DataPartner.KParty(regression = "cox", data = vdra_data[, 8:11],
         strata = c("Exposure", "Sex"), numDataPartners = 2, dataPartnerID = 2,
         monitorFolder = tempdir())
}
}
