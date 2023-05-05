########################### 2 PARTY GLOBAL FUNCTIONS ###########################

#' @name distributed2party
#' @title Two Party Vertical Distributed Regression Analysis
#' @description  \code{analysis_center_2_party} and \code{DataPartner.2Party}
#'   are used in conjuction with PopMedNet to perform linear, logistic, or cox
#'   regression on data that has been partitioned vertically between two data
#'   partners.  The data partner which holds the response variable(s) uses
#'   \code{AnalysisCener.2Party} and the other data partner uses
#'   \code{DataPartner.2Party}.  While both data partners share information with
#'   each other in order to perform the regression, data is kept secure and not
#'   shared, nor is any information shared that would allow one data partner to
#'   reconstruct part of the other data partners data. Final coefficients and
#'   other regression statistics are computed by the analysis center and shared
#'   with the other data partner.
#' @param regression the model to be used to fit the data.  The default
#'   regression \code{"linear"} fits a least squares linear model to the data.
#'   Alternatively, \code{"logistic"} returns a fitted logistic model, and
#'   \code{"cox"} returns a fitted Cox proportional hazards model.
#' @param data a data.frame or matrix which contains the data to be used in the
#'   model.  For \code{DataPartner.2Party()}, all columns will be used as
#'   covariates in the regression.  For \code{analysis_center_2_party()}, all
#'   columns, with the expection of the column specified by \code{response},
#'   will be used as covariates in the regression.
#' @param response for \code{"linear"} and \code{"logistic"} regression, the
#'   name of the column in \code{data} which holds the response variable.  If
#'   \code{reponse = NULL}, then the first column of \code{data} will be used as
#'   the response variable.  For \code{"cox"} regression response hold the name
#'   of the column which is time to event and the name of the column which is
#'   the event type (0 = censored, 1 = event).  If \code{response = NULL}, then
#'   the first column of \code{data} is assumed to be the time to even and the
#'   second column is assumed to be the event type.
#' @param strata for \code{"cox"} regression only.  A \code{\link{vector}} of
#'   character strings identifying the names of the covariates from either party
#'   which will be used as strata.  Both \code{AnalysisCenter.2party} and
#'   \code{DataPartner.2Party} must specify the same vector of strata.
#' @param mask logical value: If \code{FALSE}, strata levels for the strata
#'   which belong to the party which specified \code{FALSE} will be identified
#'   by name. If \code{TRUE}, levels for the strata which belong to the party
#'   which specified \code{TRUE} will be put in a random order and level names
#'   will be changed to \code{NA}.
#' @param monitor_folder the folder where the directories \code{dplocal},
#'   \code{inputfiles}, \code{macros}, \code{msoc}, and \code{rprograms} are
#'   located.
#' @param msreqid a character string specifying the name of the \emph{Request
#'   ID} as specified when creating the Distributed Regresion request on
#'   PopMedNet. Used for logging purposes only.
#' @param blocksize the minimium size used to horizontally partition the data
#'   for data transfer between the two parties.
#' @param tol the tolerance used to determine convergence in \code{"logistic"}
#'   and \code{"cox"} regression.
#' @param max_iterations the maximum number of iterations to perform
#'   \code{"logistic"} or \code{"cox"} regression before non-convergence is
#'   declared.
#' @param sleep_time the number of seconds to wait after writing the last file
#'   to disk before signalling the PMN Datamart Client that files are ready to
#'   be transferred.
#' @param max_waiting_time the number of seconds to wait to receive files before
#'   a transfer error is declared and the program halts execution.
#' @param popmednet logical value:  if \code{TRUE}, assumes that PopMednet is
#'   being used to transfer the files and implements PopMedNet specific
#'   routines. In particular, a 15 second offset terminiation of routines that
#'   execute in parallel is implemented.
#' @param trace logical value: if \code{TRUE} and \code{verbose == TRUE}, prints
#'   every function called during execution. Used for debugging.
#' @param verbose logical value.  If \code{TRUE}, prints out information to
#'   document the progression of the computation.
#' @return Returns an object of \code{\link{class}} \code{\link{vdralinear}} for
#'   linear regression, \code{\link{vdralogistic}} for logistic regression, or
#'   \code{\link{vdracox}} for cox regression.
#'
#' @seealso \code{\link{AnalysisCenter.3Party}}
#'   \code{\link{AnalysisCenter.KParty}}
#' @examples
#' \dontrun{
#' ## 2 party linear regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- analysis_center_2_party(regression = "linear",
#'                             data = vdra_data[, c(1, 5:7)],
#'                             response = "Change_BMI",
#'                             monitor_folder = tempdir())
#'
#' # Data Partner -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.2Party(regression = "linear", data = vdra_data[, 8:11],
#'                             monitor_folder = tempdir())
#'
#' ## 2 party logistic regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- analysis_center_2_party(regression = "logistic",
#'                             data = vdra_data[, c(2, 5:7)],
#'                             response = "WtLost",
#'                             monitor_folder = tempdir())
#'
#' # Data Partner -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.2Party(regression = "logistic",
#'                          data = vdra_data[, 8:11],
#'                          monitor_folder = tempdir())
#'
#' ## 2 party cox regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- analysis_center_2_party(regression = "cox",
#'                             data = vdra_data[, c(3:4, 5:7)],
#'                             response = c("Time", "Status"),
#'                             strata = c("Exposure", "Sex"),
#'                             monitor_folder = tempdir())
#'
#' # Data Partner -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#'    fit <- DataPartner.2Party(regression = "cox",
#'                             data = vdra_data[, 8:11],
#'                             strata = c("Exposure", "Sex"),
#'                             monitor_folder = tempdir())
#' }
#' @export
analysis_center_2_party <- function(regression            = "linear",
                                    data                  = NULL,
                                    response              = NULL,
                                    strata                = NULL,
                                    mask                  = TRUE,
                                    monitor_folder         = NULL,
                                    msreqid               = "v_default_00_000",
                                    blocksize             = 500,
                                    tol                   = 1e-8,
                                    max_iterations         = 25,
                                    sleep_time             = 10,
                                    max_waiting_time        = 86400,
                                    popmednet             = TRUE,
                                    trace                 = FALSE,
                                    verbose               = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (regression == "cox") {
    stats <- party_a_process_2_cox(data, response, strata, mask, monitor_folder,
                                   msreqid, blocksize, tol, max_iterations,
                                   sleep_time, max_waiting_time, popmednet,
                                   trace, verbose)
  } else if (regression == "linear") {
    stats <- party_a_process_2_linear(data, response, monitor_folder, msreqid,
                                      blocksize, sleep_time, max_waiting_time,
                                      popmednet, trace, verbose)
  } else if (regression == "logistic") {
    stats <- party_a_process_2_logistic(data, response, monitor_folder, msreqid,
                                        blocksize, tol, max_iterations,
                                        sleep_time, max_waiting_time, popmednet,
                                        trace, verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}

#' @rdname distributed2party
#' @export
DataPartner.2Party <- function(regression          = "linear",
                               data                = NULL,
                               strata              = NULL,
                               mask                = TRUE,
                               monitor_folder       = NULL,
                               sleep_time           = 10,
                               max_waiting_time      = 86400,
                               popmednet           = TRUE,
                               trace               = FALSE,
                               verbose             = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (regression == "cox") {
    stats <- party_b_process_2_cox(data, strata, mask,
                                   monitor_folder, sleep_time, max_waiting_time,
                                   popmednet, trace, verbose)
  } else if (regression == "linear") {
    stats <- party_b_process_2_linear(data, monitor_folder,
                                      sleep_time, max_waiting_time,
                                      popmednet, trace, verbose)
  } else if (regression == "logistic") {
    stats <- party_b_process_2_logistic(data, monitor_folder,
                                        sleep_time, max_waiting_time,
                                        popmednet, trace, verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time, final = TRUE,
                        time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}

########################### 3 PARTY GLOBAL FUNCTIONS ###########################

#' @rdname distributed3party
#' @export
DataPartner1.3Party <- function(regression            = "linear",
                                data                  = NULL,
                                response              = NULL,
                                strata                = NULL,
                                mask                  = TRUE,
                                monitor_folder         = NULL,
                                sleep_time             = 10,
                                max_waiting_time        = 86400,
                                popmednet             = TRUE,
                                trace                 = FALSE,
                                verbose               = TRUE) {

  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (regression == "cox") {
    stats <- party_a_process_3_cox(data, response, strata, mask, monitor_folder,
                                   sleep_time, max_waiting_time, popmednet,
                                   trace, verbose)
  } else if (regression == "linear") {
    stats <- party_a_process_3_linear(data, response, monitor_folder,
                                      sleep_time, max_waiting_time,
                                      popmednet, trace,
                                      verbose)
  } else  if (regression == "logistic") {
    stats <- party_a_process_3_logistic(data, response, monitor_folder,
                                        sleep_time, max_waiting_time,
                                        popmednet, trace,
                                        verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}


#' @rdname distributed3party
#' @export
DataPartner2.3Party <- function(regression          = "linear",
                                data                = NULL,
                                strata              = NULL,
                                mask                = TRUE,
                                monitor_folder       = NULL,
                                sleep_time           = 10,
                                max_waiting_time      = 86400,
                                popmednet           = TRUE,
                                trace               = FALSE,
                                verbose             = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (regression == "cox") {
    stats <- party_b_process_3_cox(data, strata, mask, monitor_folder,
                                   sleep_time, max_waiting_time, popmednet,
                                   trace, verbose)
  } else if (regression == "linear") {
    stats <- party_b_process_3_linear(data, monitor_folder,
                                      sleep_time, max_waiting_time, popmednet,
                                      trace, verbose)
  } else if (regression == "logistic") {
    stats <- party_b_process_3_logistic(data, monitor_folder,
                                        sleep_time, max_waiting_time, popmednet,
                                        trace, verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}


#' @name distributed3party
#' @title Three Party Vertical Distributed Regression Analysis
#' @description \code{AnalysisCenter.3Party}, \code{DataPartner1.3Party} and
#'   \code{DataPartner2.3Party} are used in conjunction with PopMedNet to
#'   perform linear, logistic, or cox regression on data that has been
#'   partitioned vertically between two data partners.  The data partner which
#'   holds the response variable(s) uses \code{Datapartner1.3Party} and the
#'   other data partner uses \code{DataPartner2.3Party}.  Data partners are not
#'   allowed to communicate with each other, but share information via a trusted
#'   third party analysis center.  While any information that is shared with the
#'   analysis center by a data partner, with the exception of some summary
#'   statistics, is encrypted by the sending data partner, if the information
#'   needs to be sent on to the other data partner for further analysis, the
#'   analysis center further encrypts the data.  That way, any information that
#'   deals directly with the raw data that moves between two data partners is
#'   doubly encrypted to keep both the analysis center and the other data
#'   partner from learning it.  Thus, no information is shared between the data
#'   partners or analysis center that would allow one data partner to
#'   reconstruct part of the other data partners data.  Final coefficients and
#'   other regression statistics are computed by the analysis center and shared
#'   with the data partners.
#' @param regression the model to be used to fit the data.  The default
#'   regression \code{"linear"} fits a least squares linear model to the data.
#'   Alternatively, \code{"logistic"} returns a fitted logistic model, and
#'   \code{"cox"} returns a fitted Cox proportional hazards model.
#' @param data a data.frame or matrix which contains the data to be used in the
#'   model.  For \code{DataPartner2.3Party()}, all columns will be used as
#'   covariates in the regression.  For \code{DataPartner1.3Party()}, all
#'   columns, with the exception of the column specified by \code{response},
#'   will be used as covariates in the regression.
#' @param response for \code{"linear"} and \code{"logistic"} regression, the
#'   name of the column in \code{data} which holds the response variable.  If
#'   \code{reponse = NULL}, then the first column of \code{data} will be used as
#'   the response variable.  For \code{"cox"} regression response hold the name
#'   of the column which is time to event and the name of the column which is
#'   the event type (0 = censored, 1 = event).  If \code{response = NULL}, then
#'   the first column of \code{data} is assumed to be the time to even and the
#'   second column is assumed to be the event type.
#' @param strata for \code{"cox"} regression only.  A \code{\link{vector}} of
#'   character strings identifying the names of the covariates from either party
#'   which will be used as strata.  Both \code{DataPartner1_3party} and
#'   \code{DataPartner2.3Party} must specify the same vector of strata.
#' @param mask logical value: If \code{FALSE}, strata levels for the strata
#'   which belong to the party which specified \code{FALSE} will be identified
#'   by name. If \code{TRUE}, levels for the strata which belong to the party
#'   which specified \code{TRUE} will be put in a random order and level names
#'   will be changed to \code{NA}.
#' @param monitor_folder the folder where the directories \code{dplocal},
#'   \code{inputfiles}, \code{macros}, \code{msoc}, and \code{rprograms} are
#'   located.
#' @param msreqid a character string specifying the name of the \emph{Request
#'   ID} as specified when creating the Distributed Regression request on
#'   PopMedNet. Used for logging purposes only.
#' @param blocksize the minimum size used to horizontally partition the data for
#'   data transfer between the two parties.
#' @param tol the tolerance used to determine convergence in \code{"logistic"}
#'   and \code{"cox"} regression.
#' @param max_iterations the maximum number of iterations to perform
#'   \code{"logistic"} or \code{"cox"} regression before non-convergence is
#'   declared.
#' @param sleep_time the number of seconds to wait after writing the last file
#'   to disk before signalling the PMN Datamart Client that files are ready to
#'   be transferred.
#' @param max_waiting_time the number of seconds to wait to receive files before
#'   a transfer error is declared and the program halts execution. Should be the
#'   same for both parties when \code{delayOffset = TRUE}.
#' @param popmednet logical value:  if \code{TRUE}, assumes that PopMedNet is
#'   being used to transfer the files and implements PopMedNet specific
#'   routines. In particular, a 15 second offset between termination of routines
#'   that execute in parallel is implemented.
#' @param trace logical value: if \code{TRUE} and \code{verbose == TRUE}, prints
#'   every function call. Used for debugging.
#' @param verbose logical value.  If \code{TRUE}, prints out information to
#'   document the progression of the computation.
#'
#' @return Returns an object of \code{\link{class}} \code{\link{vdralinear}} for
#'   linear regression, \code{\link{vdralogistic}} for logistic regression, or
#'   \code{\link{vdracox}} for cox regression.
#' @seealso \code{\link{analysis_center_2_party}}
#'   \code{\link{AnalysisCenter.KParty}}
#'
#' @examples
#' \dontrun{
#' ## 3 party linear regression
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.3Party(regression = "linear",
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner1.3Party(regression = "linear",
#'                           data = vdra_data[, c(1, 5:7)],
#'                           response = "Change_BMI",
#'                           monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner2.3Party(regression = "linear",
#'                           data = vdra_data[, 8:11],
#'                           monitor_folder = tempdir())
#'
#' ## 3 party logistic regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.3Party(regression = "logistic",
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner1.3Party(regression = "logistic",
#'                           data = vdra_data[, c(2, 5:7)],
#'                           response = "WtLost",
#'                           monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner2.3Party(regression = "logistic",
#'                           data = vdra_data[, 8:11],
#'                           monitor_folder = tempdir())
#'
#' ## 3 party cox regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.3Party(regression = "cox",
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner1.3Party(regression = "cox",
#'                           data = vdra_data[, c(3:4, 5:7)],
#'                           response = c("Time", "Status"),
#'                           strata = c("Exposure", "Sex"),
#'                           monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner2.3Party(regression = "cox",
#'                           data = vdra_data[, 8:11],
#'                           strata = c("Exposure", "Sex"),
#'                           monitor_folder = tempdir())
#' }
#' @export
AnalysisCenter.3Party <- function(regression            = "linear",
                                  monitor_folder         = NULL,
                                  msreqid               = "v_default_00_000",
                                  blocksize             = 500,
                                  tol                   = 1e-8,
                                  max_iterations         = 25,
                                  sleep_time             = 10,
                                  max_waiting_time        = 86400,
                                  popmednet             = TRUE,
                                  trace                 = FALSE,
                                  verbose               = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (regression == "cox") {
    stats <- party_t_process_3_cox(monitor_folder, msreqid, blocksize, tol,
                                   max_iterations, sleep_time, max_waiting_time,
                                   popmednet, trace, verbose)
  } else if (regression == "linear") {
    stats <- party_t_process_3_linear(monitor_folder, msreqid, blocksize,
                                      sleep_time, max_waiting_time, popmednet,
                                      trace, verbose)
  } else if (regression == "logistic") {
    stats <- party_t_process_3_logistic(monitor_folder, msreqid, blocksize, tol,
                                        max_iterations, sleep_time,
                                        max_waiting_time, popmednet, trace,
                                        verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}

########################### K PARTY GLOBAL FUNCTIONS ###########################

#' @rdname distributedKparty
#' @export
DataPartner.KParty <- function(regression            = "linear",
                               data                  = NULL,
                               response              = NULL,
                               strata                = NULL,
                               mask                  = TRUE,
                               num_data_partners       = NULL,
                               data_partner_id         = NULL,
                               monitor_folder         = NULL,
                               sleep_time             = 10,
                               max_waiting_time        = 86400,
                               popmednet             = TRUE,
                               trace                 = FALSE,
                               verbose               = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (is.null(num_data_partners)) {
    warning("num_data_partners must be specified")
  } else if (is.null(data_partner_id)) {
    warning("data_partner_id must be specified")
  } else if (regression == "cox") {
    stats <- data_partner_k_cox(data, response, strata, mask, num_data_partners,
                                data_partner_id, monitor_folder,
                                sleep_time, max_waiting_time, popmednet, trace,
                                verbose)
  } else if (regression == "linear") {
    stats <- data_partner_k_linear(data, response, num_data_partners,
                                   data_partner_id, monitor_folder,
                                   sleep_time, max_waiting_time, popmednet,
                                   trace, verbose)
  } else  if (regression == "logistic") {
    stats <- data_partner_k_logistic(data, response, num_data_partners,
                                     data_partner_id, monitor_folder,
                                     sleep_time, max_waiting_time, popmednet,
                                     trace, verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}


#' @name distributedKparty
#' @title -Party Vertical Distributed Regression Analysis
#' @description \code{AnalysisCenter.KParty} and \code{DataPartner.KParty} are
#'   used in conjunction with PopMedNet to perform linear, logistic, or cox
#'   regression on data that has been partitioned vertically between two or more
#'   data partners.  The data partners which holds the data use
#'   \code{DataPartner.KParty} while a trusted "third" party uses
#'   \code{AnalysisCenter.KParty}.  Data partners are allowed to communicate
#'   with each other and the analysis center, no information is shared between
#'   the data partners or analysis center that would allow one data partner or
#'   the analysis center to reconstruct part of the other data partners data.
#'   Final coefficients and other regression statistics are computed by the
#'   analysis center and shared with the data partners.
#' @param regression the model to be used to fit the data.  The default
#'   regression \code{"linear"} fits a least squares linear model to the data.
#'   Alternatively, \code{"logistic"} returns a fitted logistic model, and
#'   \code{"cox"} returns a fitted Cox proportional hazards model.
#' @param data a data.frame or matrix which contains the data to be used in the
#'   model.  All columns will be used as covariates in the regression with the
#'   exception of the data partner which has \code{data_partner_id = 1}.  For
#'   this data partner, all columns, with the exception of the column specified
#'   by \code{response}, will be used as covariates in the regression.
#' @param response only used for data partner with \code{data_partner_id = 1}.
#'   For \code{"linear"} and \code{"logistic"} regression, the name of the
#'   column in \code{data} which holds the response variable.  If \code{reponse
#'   = NULL}, then the first column of \code{data} will be used as the response
#'   variable. For \code{"cox"} regression response hold the name of the column
#'   which is time to event and the name of the column which is the event type
#'   (0 = censored, 1 = event).  If \code{response = NULL}, then the first
#'   column of \code{data} is assumed to be the time to even and the second
#'   column is assumed to be the event type.
#' @param strata for \code{"cox"} regression only.  A \code{\link{vector}} of
#'   character strings identifying the names of the covariates from either party
#'   which will be used as strata.  All data partners must specify the same
#'   vector of strata.
#' @param mask logical value: If \code{FALSE}, strata levels for the strata
#'   which belong to the party which specified \code{FALSE} will be identified
#'   by name. If \code{TRUE}, levels for the strata which belong to the party
#'   which specified \code{TRUE} will be put in a random order and level names
#'   will be changed to \code{NA}.
#' @param num_data_partners the number of data partners which are supplying data
#'   for the regression.
#' @param data_partner_id a unique identifier for each data partner.  The data
#'   partner with the response variable(s) must have \code{data_partner_id = 1}.
#'   All other data partners must have an integer value from 2 to
#'   \code{num_data_partners}.
#' @param monitor_folder the folder where the directories \code{dplocal},
#'   \code{inputfiles}, \code{macros}, \code{msoc}, and \code{rprograms} are
#'   located.
#' @param msreqid a character string specifying the name of the \emph{Request
#'   ID} as specified when creating the Distributed Regresion request on
#'   PopMedNet. Used for logging purposes only.
#' @param tol the tolerance used to determine convergence in \code{"logistic"}
#'   and \code{"cox"} regression.
#' @param max_iterations the maximum number of iterations to perform
#'   \code{"logistic"} or \code{"cox"} regression before non-convergence is
#'   declared.
#' @param sleep_time the number of seconds to wait after writing the last file
#'   to disk before signalling the PMN Datamart Client that files are ready to
#'   be transferred.
#' @param max_waiting_time the number of seconds to wait to receive files before
#'   a transfer error is declared and the program halts execution. Should be the
#'   same for all parties when \code{delayOffset = TRUE}.
#' @param popmednet logical value:  if \code{TRUE}, assumes that PopMednet is
#'   being used to transfer the files and implements PopMedNet specific
#'   routines. In particular, a 15 second offset between terminiation of
#'   routines that execute in parallel is implemented.
#' @param trace logical value: if \code{TRUE} and \code{verbose == TRUE}, prints
#'   every function call. Used for debugging.
#' @param verbose logical value.  If \code{TRUE}, prints out information to
#'   document the progression of the computation.
#' @return Returns an object of \code{\link{class}} \code{\link{vdralinear}} for
#'   linear regression, \code{\link{vdralogistic}} for logistic regression, or
#'   \code{\link{vdracox}} for cox regression.
#' @seealso \code{\link{analysis_center_2_party}}
#'   \code{\link{AnalysisCenter.KParty}}
#' @examples
#' \dontrun{
#' ## 3 party linear regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.KParty(regression = "linear",
#'                             num_data_partners = 2,
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "linear",
#'                          data = vdra_data[, c(1, 5:7)],
#'                          response = "Change_BMI",
#'                          num_data_partners = 2,
#'                          data_partner_id = 1,
#'                          monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "linear",
#'                          data = vdra_data[, 8:11],
#'                          num_data_partners = 2,
#'                          data_partner_id = 2,
#'                          monitor_folder = tempdir())
#'
#' ## 3 party logistic regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.KParty(regression = "logistic",
#'                             num_data_partners = 2,
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "logistic",
#'                          data = vdra_data[, c(2, 5:7)],
#'                          response = "WtLost",
#'                          num_data_partners = 2,
#'                          data_partner_id = 1,
#'                          monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "logistic",
#'                          data = vdra_data[, 8:11],
#'                          num_data_partners = 2,
#'                          data_partner_id = 2,
#'                          monitor_folder = tempdir())
#'
#' ## 3 party cox regression
#'
#' # Analysis Center -- To be run in one instance of R.
#' # The working directory should be the same as specified in the PopMedNet
#' # requset for the analysis center.
#'
#' fit <- AnalysisCenter.KParty(regression = "cox",
#'                             num_data_partners = 2,
#'                             monitor_folder = tempdir())
#'
#' # Data Partner 1 -- To be run in second instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "cox",
#'                          data = vdra_data[, c(3:4, 5:7)],
#'                          response = c("Time", "Status"),
#'                          strata = c("Exposure", "Sex"),
#'                          num_data_partners = 2,
#'                          data_partner_id = 1,
#'                          monitor_folder = tempdir())
#'
#' # Data Partner 2 -- To be run in third instand of R, on perhaps a different
#' # machine. The working directory should be the same as specified in the
#' # PopMedNet request for the data partner.
#'
#' fit <- DataPartner.KParty(regression = "cox",
#'                          data = vdra_data[, 8:11],
#'                          strata = c("Exposure", "Sex"),
#'                          num_data_partners = 2,
#'                          data_partner_id = 2,
#'                          monitor_folder = tempdir())
#' }
#' @export
AnalysisCenter.KParty <- function(regression          = "linear",
                                  num_data_partners   = NULL,
                                  monitor_folder      = NULL,
                                  msreqid             = "v_default_00_000",
                                  tol                 = 1e-8,
                                  max_iterations      = 25,
                                  sleep_time          = 10,
                                  max_waiting_time    = 86400,
                                  popmednet           = TRUE,
                                  trace               = FALSE,
                                  verbose             = TRUE) {
  start_time <- proc.time()
  stats <- list()
  if (verbose) cat("Process started on", as.character(get_utc_time()), "UTC.\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.")
  } else if (is.null(num_data_partners)) {
    warning("num_data_partners must be specified.")
  } else if (regression == "cox") {
    stats <- analysis_center_k_cox(num_data_partners, monitor_folder, msreqid,
                                   tol, max_iterations, sleep_time,
                                   max_waiting_time, popmednet, trace, verbose)
  } else if (regression == "linear") {
    stats <- analysis_center_k_linear(num_data_partners, monitor_folder,
                                      msreqid, sleep_time, max_waiting_time,
                                      popmednet, trace, verbose)
  } else if (regression == "logistic") {
    stats <- analysis_center_k_logistic(num_data_partners, monitor_folder,
                                        msreqid, tol, max_iterations,
                                        sleep_time, max_waiting_time, popmednet,
                                        trace, verbose)
  } else {
    warning("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp <- get_elapsed_time(proc.time() - start_time,
                        final = TRUE, time_only = FALSE)
  if (verbose) cat("Process completed on",
                   as.character(get_utc_time()), "UTC.\n")
  if (verbose) cat(elp, "\n")
  return(stats)
}

############################ SHARED SETUP FUNCTIONS ############################

#' @importFrom utils file_test
create_io_location <- function(monitor_folder, folder) {
  location <- file.path(monitor_folder, folder)
  if (!dir.exists(location) && !file.exists(location)) {
    # directory does not exist, so create it.
    dir.create(location)
    return(TRUE)
  }
  if (file_test("-d", location)) {
    # directory exists.  No need to create it.
    return(TRUE)
  }
  # file with directory name exists.  cannot create directory
  return(FALSE)
}


#' @importFrom utils write.csv
check_data_format <- function(params, data) {
  if (is.data.frame(data)) {
    data <- data.frame(data)
  } else if (is.matrix(data)) {
    data <- matrix(data)
  } else {
    warning("Data is not a matrix or a data frame.")
    return(TRUE)
  }

  if (nrow(data) == 0 || ncol(data) == 0) {
    warning("The data is empty.")
    return(TRUE)
  }
  bad_value <- rep(FALSE, nrow(data))
  for (i in seq_len(ncol(data))) {
    if (is.integer(data[, i]) || is.double(data[, i]) ||
        is.single(data[, i]) || is.numeric(data[, i])) {
      bad_value <- bad_value | !is.finite(data[, i])
    } else {
      bad_value <- bad_value | is.na(data[, i])
    }
  }
  idx <- data.frame(which(bad_value))
  colnames(idx) <- "Observations with invalid entries"
  if (nrow(idx) > 0) {
    warning(
      paste0("Some observations contain invalid values: NA, NaN, or Inf. ",
             "A list of all such observations has been outputted to",
             file.path(params$write_path, "invalidEntries.csv"),
             ". Terminating program."))
    write.csv(idx, file.path(params$write_path, "invalidEntries.csv"))
    return(TRUE)
  }
  if (is.null(colnames(data))) {
    warning("Variables are not named.")
    return(TRUE)
  }
  return(FALSE)
}

check_response <- function(params, data, y_name) {
  if (is.null(y_name)) {
    warning("Response is not specified.")
    return(NULL)
  }
  if (!is.character(y_name)) {
    warning("response label is not a character string.")
    return(NULL)
  }
  y_name <- unique(y_name)
  if (params$analysis == "linear" || params$analysis == "logistic") {
    if (length(y_name) != 1) {
      warning(paste("Specify only one reponse for",
                    params$analysis, "regression."))
      return(NULL)
    }
    response_col_index <- which(colnames(data) %in% y_name)
    if (length(response_col_index) == 0) {
      warning("Response variable not found.")
      return(NULL)
    }
    if (length(response_col_index) > 1) {
      warning("Response variable appears more than once.")
      return(NULL)
    }
  }
  if (params$analysis == "cox") {
    if (length(y_name) != 2) {
      warning("Specify exactly two variables ",
              "(time and censor) for Cox regression.")
      return(NULL)
    }
    response_col_index_time   <- c(which(colnames(data) %in% y_name[1]))
    response_col_index_censor <- c(which(colnames(data) %in% y_name[2]))
    if (length(response_col_index_time) == 0) {
      warning("Time variable not found.")
      return(NULL)
    }
    if (length(response_col_index_time) > 1) {
      warning("Time variable appears more than once.")
      return(NULL)
    }
    if (length(response_col_index_censor) == 0) {
      warning("Censor variable not found.")
      return(NULL)
    }
    if (length(response_col_index_censor) > 1) {
      warning("Censor variable appears more than once.")
      return(NULL)
    }
    response_col_index <- c(response_col_index_time, response_col_index_censor)
  }
  for (i in seq_along(y_name)) {
    if (!is.numeric(data[, response_col_index[i]]) &&
        !is.integer(data[, response_col_index[i]])) {
      warning(paste(y_name[i], "is not numeric."))
      return(NULL)
    }
  }
  if (params$analysis == "logistic") {
    if (sum(!(data[, response_col_index] %in% c(0, 1))) > 0) {
      warning("Response variable is not binary. ",
              "It should only be 0's and  1's.")
      return(NULL)
    }
  }
  if (params$analysis == "cox") {
    if (sum(!(data[, response_col_index[2]] %in% c(0, 1))) > 0) {
      warning("Censoring variable is not binary. ",
              "It should only be 0's and  1's.")
      return(NULL)
    }
  }
  return(response_col_index)
}

create_model_matrix_tags <- function(data) {
  if (ncol(data) == 0) {
    return(c())
  }
  num     <- numeric(ncol(data))
  classes <- character(ncol(data))
  for (i in seq_len(ncol(data))) {
    if (is.integer(data[, i]) || is.double(data[, i]) ||
        is.single(data[, i]) || is.numeric(data[, i])) {
      num[i] <- 1
      classes[i] <- "numeric"
    } else {
      num[i] <- length(unique(data[, i])) - 1
      classes[i] <- "factor"
    }
  }
  tags <- rep(names(data), times = num)
  names(tags) <- rep(classes, times = num)
  return(tags)
}

########################### 2 PARTY SETUP FUNCTIONS ############################

prepare_params_2p <- function(analysis, party, msreqid = "v_default_00_000",
                              popmednet = TRUE, trace = FALSE, verbose = TRUE) {
  params                     <- list()
  params$party_name           <- party
  params$analysis            <- analysis
  params$msreqid             <- msreqid
  params$popmednet           <- popmednet
  params$trace               <- trace & verbose
  params$verbose             <- verbose
  params$failed              <- FALSE
  params$errorMessage        <- ""
  params$pmn_step_counter      <- 0
  params$algIterationCounter <- 0
  params$completed           <- FALSE
  params$converged           <- FALSE
  params$maxIterExceeded     <- FALSE
  params$lastIteration       <- FALSE
  params$p1                  <- 0
  params$p2                  <- 0
  params$p1_old              <- 0
  params$p2_old              <- 0
  params$stats               <- list()
  class(params$stats)        <- paste0("vdra", analysis)
  params$stats$failed        <- TRUE
  params$stats$converged     <- FALSE
  return(params)
}

########################### 3 PARTY SETUP FUNCTIONS ############################

prepare_params_3p <- function(analysis, party, msreqid = "v_default_00_000",
                              popmednet = TRUE, trace = FALSE, verbose = TRUE) {
  params                     <- list()
  params$party_name           <- party
  params$analysis            <- analysis
  params$msreqid             <- msreqid
  params$popmednet           <- popmednet
  params$trace               <- trace & verbose
  params$verbose             <- verbose
  params$failed              <- FALSE
  params$errorMessage        <- ""
  params$pmn_step_counter      <- 0
  params$algIterationCounter <- 0
  params$completed           <- FALSE
  params$converged           <- FALSE
  params$maxIterExceeded     <- FALSE
  params$lastIteration       <- FALSE
  params$p1                  <- 0
  params$p2                  <- 0
  params$p1_old              <- 0
  params$p2_old              <- 0
  params$stats               <- list()
  class(params$stats)        <- paste0("vdra", analysis)
  params$stats$failed        <- TRUE
  params$stats$converged     <- FALSE
  return(params)
}

########################### K PARTY SETUP FUNCTIONS ############################

prepare_params_kp <- function(analysis, data_partner_id, num_data_partners,
                              msreqid = "v_default_00_000", cutoff = NULL,
                              max_iterations = NULL, ac = FALSE,
                              popmednet = TRUE,
                              trace = FALSE, verbose = TRUE) {
  params                     <- list()
  params$data_partner_id       <- data_partner_id
  params$num_data_partners     <- num_data_partners
  params$analysis            <- analysis
  params$msreqid             <- msreqid
  params$popmednet           <- popmednet
  params$trace               <- trace & verbose
  params$verbose             <- verbose
  params$failed              <- FALSE
  params$errorMessage        <- ""
  params$pmn_step_counter      <- 0
  params$algIterationCounter <- 0
  params$max_iterations       <- max_iterations
  params$completed           <- FALSE
  params$converged           <- FALSE
  params$maxIterExceeded     <- FALSE
  params$lastIteration       <- FALSE
  params$cutoff              <- cutoff
  params$stats               <- list()
  class(params$stats)        <- paste0("vdra", analysis)
  params$stats$failed        <- TRUE
  params$stats$converged     <- FALSE
  if (((!is.integer(num_data_partners) && !is.numeric(num_data_partners)) ||
       num_data_partners <= 0 || is.infinite(num_data_partners) ||
       round(num_data_partners) != num_data_partners)) {
    params$failed <- TRUE
    params$errormessage <-
      paste("num_data_partners must be a positive integer,",
            "and must equal the number of data partners providing data.")
  }
  if (!params$failed) {
    if (ac) {
      if (data_partner_id != 0) {
        params$failed <- TRUE
        params$errormessage <-
          "data_partner_id for Analysis Center must be 0.\n\n"
      }
    } else {
      if (data_partner_id <= 0 || data_partner_id > num_data_partners) {
        params$failed <- TRUE
        params$errormessage <- paste0("data_partner_id must be between 1 and ",
                                      num_data_partners, " inclusive.\n\n")
      }
    }
  }
  return(params)
}

########################### PRETTY OUTPUT FUNCTIONS ############################

header <- function(params) {
  large_cox <-
    c("  ____ _____  __",
      " / ___/ _ \\ \\/ /",
      "| |  | | | \\  / ",
      "| |__| |_| /  \\ ",
      " \\____\\___/_/\\_\\")
  large_linear <-
    c(" _     ___ _   _ _____    _    ____  ",
      "| |   |_ _| \\ | | ____|  / \\  |  _ \\ ",
      "| |    | ||  \\| |  _|   / _ \\ | |_) |",
      "| |___ | || |\\  | |___ / ___ \\|  _ < ",
      "|_____|___|_| \\_|_____/_/   \\_|_| \\_\\")
  large_logistic <-
    c(" _     ___   ____ ___ ____ _____ ___ ____ ",
      "| |   / _ \\ / ___|_ _/ ___|_   _|_ _/ ___|",
      "| |  | | | | |  _ | |\\___ \\ | |  | | |    ",
      "| |__| |_| | |_| || | ___) || |  | | |___ ",
      "|_____\\___/ \\____|___|____/ |_| |___\\____|")
  large_regression <-
    c(" ____  _____ ____ ____  _____ ____ ____ ___ ___  _   _ ",
      "|  _ \\| ____/ ___|  _ \\| ____/ ___/ ___|_ _/ _ \\| \\ | |",
      "| |_) |  _|| |  _| |_) |  _| \\___ \\___ \\| | | | |  \\| |",
      "|  _ <| |__| |_| |  _ <| |___ ___) ___) | | |_| | |\\  |",
      "|_| \\_|_____\\____|_| \\_|_____|____|____|___\\___/|_| \\_|")
  small_cox <-
    c("  ___ _____  __",
      " / __/ _ \\ \\/ /",
      "| (_| (_) >  < ",
      " \\___\\___/_/\\_\\")
  small_linear <-
    c(" _    ___ _  _ ___   _   ___ ",
      "| |  |_ _| \\| | __| /_\\ | _ \\",
      "| |__ | || .` | _| / _ \\|   /",
      "|____|___|_|\\_|___/_/ \\_\\_|_\\")
  small_logistic <-
    c(" _    ___   ___ ___ ___ _____ ___ ___ ",
      "| |  / _ \\ / __|_ _/ __|_   _|_ _/ __|",
      "| |_| (_) | (_ || |\\__ \\ | |  | | (__ ",
      "|____\\___/ \\___|___|___/ |_| |___\\___|")
  small_regression <-
    c(" ___ ___ ___ ___ ___ ___ ___ ___ ___  _  _ ",
      "| _ \\ __/ __| _ \\ __/ __/ __|_ _/ _ \\| \\| |",
      "|   / _| (_ |   / _|\\__ \\__ \\| | (_) | .` |",
      "|_|_\\___\\___|_|_\\___|___/___/___\\___/|_|\\_|")
  tiny_cox <-
    c("+-+-+-+",
      "|C|O|X|")
  tiny_linear <-
    c("+-+-+-+-+-+-+",
      "|L|I|N|E|A|R|")
  tiny_logistic <-
    c("+-+-+-+-+-+-+-+-+",
      "|L|O|G|I|S|T|I|C|")
  tiny_regression <-
    c("+-+-+-+-+-+-+-+-+-+-+",
      "|R|E|G|R|E|S|S|I|O|N|",
      "+-+-+-+-+-+-+-+-+-+-+")

  width <- getOption("width")
  if (width > nchar(large_regression[1])) {
    cox        <- large_cox
    linear     <- large_linear
    logistic   <- large_logistic
    regression <- large_regression
  } else if (width  > nchar(small_regression[1])) {
    cox        <- small_cox
    linear     <- small_linear
    logistic   <- small_logistic
    regression <- small_regression
  } else {
    cox        <- tiny_cox
    linear     <- tiny_linear
    logistic   <- tiny_logistic
    regression <- tiny_regression
  }

  offset_cox        <- floor((width - nchar(cox[1])) / 2)
  offset_linear     <- floor((width - nchar(linear[1])) / 2)
  offset_logistic   <- floor((width - nchar(logistic[1])) / 2)
  offset_regression <- floor((width - nchar(regression[1])) / 2)

  space_cox        <- paste(rep(" ", offset_cox), collapse = "")
  space_linear     <- paste(rep(" ", offset_linear), collapse = "")
  space_logistic   <- paste(rep(" ", offset_logistic), collapse = "")
  space_regression <- paste(rep(" ", offset_regression), collapse = "")

  if (params$analysis == "linear") {
    if (params$verbose) cat(paste0("\r", space_linear, linear, "\n"))
    if (params$verbose) cat(paste0("\r", space_regression, regression, "\n"))
  }
  if (params$analysis == "logistic") {
    if (params$verbose) cat(paste0("\r", space_logistic, logistic, "\n"))
    if (params$verbose) cat(paste0("\r", space_regression, regression, "\n"))
  }
  if (params$analysis == "cox") {
    if (params$verbose) cat(paste0("\r", space_cox, cox, "\n"))
    if (params$verbose) cat(paste0("\r", space_regression, regression, "\n"))
  }
  if (params$verbose) cat("\n")
}


BeginningIteration <- function(params) {
  width  <- getOption("width")
  msg    <- paste("*** Beginning Iteration", params$algIterationCounter, "***")
  offset <- max(floor((width - nchar(msg)) / 2) - 1, 0)
  space  <- paste(rep(" ", offset), collapse = "")
  if (params$verbose) cat(space, msg, "\n\n")
}


EndingIteration <- function(params) {
  width  <- getOption("width")
  msg    <- paste("*** Ending Iteration", params$algIterationCounter, "***")
  offset <- floor((width - nchar(msg)) / 2) - 1
  space  <- paste(rep(" ", offset), collapse = "")
  if (params$verbose) cat(space, msg, "\n\n")
}


GetLion <- function(p) {
  lion1 <- rep("", 5)
  lion1[1] <- "    (\"`-''-/\").___..--''\"`-._\"      "
  lion1[2] <- "    `6_ 6 ) `-. ( ).`-.__.`)        "
  lion1[3] <- "    (_Y_.)' ._ ) `._ `. ``-..-'     "
  lion1[4] <- "     _..`--'_..-_/ /--'_.' ,'       "
  lion1[5] <- "    (il),-'' (li),' ((!.-'          "

  lion2 <- rep("", 8)
  lion2[1] <- "        ___  ___  _  _  _  _        "
  lion2[2] <- "       | -_>||__>|\\ |||\\ ||       "
  lion2[3] <- "       | |  ||__>| \\||| \\||       "
  lion2[4] <- "       |_|  ||__>|_\\_||_\\_|       "
  lion2[5] <- "     ___  _____  ___  _____  ___    "
  lion2[6] <- "    //__>|_   _|//_\\|_   _|||__>   "
  lion2[7] <- "    \\_\\  | |  | | |  | |  ||__>   "
  lion2[8] <- "    <__//  |_|  |_|_|  |_|  ||__>   "

  lion3 <- rep("", 4)
  lion3[1] <- "        ___   ___   _   _           "
  lion3[2] <- "        | -_> //__> | | | |         "
  lion3[3] <- "        | |   \\_\\ | |_| |         "
  lion3[4] <- "        |_|   <__// \\___//         "

  lion4 <-    "    (il),-'' (li),' ((!.-'      PSU "

  lion5 <-    "              PSU    "

  if (p >= 13) {
    nittany <- c(lion1, lion2)
  } else if (p >= 9) {
    nittany <- c(lion1, lion3)
  } else if (p >= 5) {
    nittany    <- lion1
    nittany[5] <- lion4
  } else {
    nittany <- lion5
  }

  diff <- p - length(nittany)
  top <- floor(diff / 2)
  bottom <- diff - top
  nittany <- c(rep("", top), nittany, rep("", bottom))
  return(nittany)
}


#' @importFrom utils flush.console
make_progress_bar_1 <- function(steps, message, verbose) {
  pb <- list()
  message_length    <- 18
  pb$num_steps      <- steps
  pb$num_blanks     <- 20
  pb$delimeter     <- "|"
  pb$filler        <- "#"
  pb$blank         <- "."
  pb$percent       <- 0
  pb$percentstr    <- "  0%"
  pb$prints <- 0
  message <- substr(message, 1, message_length)
  message <- paste0(message,
                    paste(rep(" ", message_length - nchar(message)),
                          collapse = ""))
  pb$header <- paste0("Processing ", message, ": ")
  to_print <- paste0(pb$header, pb$percentstr, pb$delimeter,
                     paste(rep(pb$blank, pb$num_blanks), collapse = ""),
                     pb$delimeter)
  if (verbose) cat(to_print, "\r")
  if (verbose) flush.console()
  return(pb)
}


#' @importFrom utils flush.console
make_progress_bar_2 <- function(i, pb, verbose) {
  percent <- floor(100 * i / pb$num_steps)
  if (percent == pb$percent) {
    return(pb)
  }
  pb$percent <- percent
  pb$percentstr <- paste0(paste(rep(" ", 3 - nchar(percent)), collapse = ""),
                          percent, "%")
  num_filler <- floor(pb$num_blanks * i / pb$num_steps)
  to_print <- paste0(pb$header, pb$percentstr, pb$delimeter,
                     paste(rep(pb$filler, num_filler), collapse = ""),
                     paste(rep(pb$blank, pb$num_blanks - num_filler),
                           collapse = ""),
                     pb$delimeter)

  if (i == pb$num_steps) {
    if (verbose) cat(to_print, "\n\n")
  } else {
    if (verbose) cat(to_print, "\r")
  }
  if (verbose) flush.console()
  return(pb)
}

############################### MATRIX FUNCTIONS ###############################

MultiplyDiagonalWTimesX <- function(w, x) {
  if (!is.matrix(x)) {
    x <- matrix(x, length(x), 1)
    wx1 <- matrix(NA, length(x), 1)
  } else {
    wx1 <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  }
  if (is.matrix(w)) {
    for (i in seq_len(nrow(w))) {
      wx1[i, ] <- w[i] * x[i, ]
    }
  } else {
    for (i in seq_along(w)) {
      wx1[i, ] <- w[i] * x[i, ]
    }
  }

  return(wx1)
}


#' @importFrom stats runif
find_orthonormal_vectors <- function(x, g) {
  x <- as.matrix(x)
  x <- cbind(x, runif(nrow(x)))  # Randomize z
  # Save the Random Vector Here
  n <- nrow(x)
  Q <- qr.Q(qr(x), complete = TRUE)
  Q <- Q[, (n - g + 1):n]
  return(Q)
}


#' @importFrom stats runif
random_orthonormal_matrix <- function(size) {
  return(qr.Q(qr(matrix(runif(size * size), size, size)), complete = TRUE))
}

#################### SHARED PMN COMMUNICATION FUNCTIONS ###################

#' @importFrom utils write.csv
make_csv <- function(file_nm, transfer_to_site_in, dp_cd_list, write_path) {
  dframe <- data.frame(file_nm, transfer_to_site_in, dp_cd_list)
  fp <- file.path(write_path, "file_list.csv")
  write.csv(dframe, fp, row.names = FALSE, quote = FALSE)
}


seq_zw <- function(letter = "Z_", nblocks = 1) {
  return(paste0(letter, 1:nblocks, ".rdata"))
}


Standby <- function(trigger_name, trigger_location,
                    sleep_time = 1, max_waiting_time = NULL, remove = FALSE,
                    verbose = TRUE) {

  found <- FALSE

  if (is.null(max_waiting_time)) {
    max_waiting_time <- 60 * 60 * 24
  }

  fpath <- file.path(trigger_location, trigger_name)
  start_time <- proc.time()[3]
  elapsed_time <- 0

  while (!found) {
    found <- all(file.exists(fpath))

    if (elapsed_time > max_waiting_time) {
      break
    }

    if (!found) {
      Sys.sleep(sleep_time)
      elapsed_time <- round(proc.time()[3] - start_time, 0)
    }

    if (verbose) cat("Elapsed Time:", HMS(elapsed_time), "\r")
  }
  if (verbose) cat("\n")
  if (!found) {
    stop("Exceeded maximum time waiting for files to be dropped.")
  }

  Sys.sleep(sleep_time)

  if (remove) delete_trigger(trigger_name, trigger_location)

}

# This function is never called!  (?)
CopyFile <- function(read_directory, write_directory, filename) {
  source      <- file.path(read_directory, filename)
  destination <- file.path(write_directory, filename)
  if (all(file.exists(source))) {
    file.copy(source, destination, overwrite = TRUE)
  } else {
    stop(paste0("These files do not exist:\n",
                paste0(source[!file.exists(source)], collapse = ", "), "\n"))
  }
}


make_trigger <- function(trigger_name, trigger_path, message = "Trigger File") {

  fn <- file.path(trigger_path, trigger_name)
  if (file.exists(fn)) {
    file.remove(fn)
  }

  write(message, fn)
}


delete_trigger <- function(trigger_name, trigger_path) {
  Sys.sleep(1)
  targets <- file.path(trigger_path, trigger_name)
  for (target in targets) {
    if (file.exists(target)) {
      start_time <- proc.time()[3]
      repeat {
        result <- suppressWarnings(try(file.remove(target)))
        if (result) break
        Sys.sleep(1)
        if (proc.time()[3] - start_time > 60) {
          stop(paste("Could not delete the file", target, "after 60 seconds."))
        }
      }
    }
  }
}


make_transfer_message <- function(write_path) {
  message <- "A has no covariates."
  save(message, file = file.path(write_path, "transferControl.rdata"))
}


make_error_message <- function(write_path, message = "") {
  save(message, file = file.path(write_path, "errorMessage.rdata"))
}


read_error_message <- function(read_path) {
  load(file.path(read_path, "errorMessage.rdata"))
  return(message)
}

###################### 2 PARTY PMN COMMUNICATION FUNCTIONS #####################

send_pause_quit_2p <- function(params,
                               files = c(),
                               sleep_time = 10,
                               job_failed = FALSE) {

  params <- store_log_entry_2p(params, files)
  params <- StoreTrackingTableEntry.2p(params)
  WriteLogCSV(params)
  write_log_raw(params)
  params$lastIteration <- TRUE
  files <- c(files, "stamps.rdata", "log.rdata", "file_list.csv")
  transfer <- c(rep(1, length(files) - 1), 10)
  if (params$party_name == "A") {
    if (job_failed) {
      files <- c(files, "job_fail.ok")
      params <- store_stamp_entry(params, "Job failed trigger file",
                                  "Trigger File created")
    } else {
      files <- c(files, "job_done.ok")
      params <- store_stamp_entry(params, "Job done trigger file",
                                  "Trigger File created")
    }
    transfer <- c(transfer, 10)
  }
  if (params$party_name == "A") {
    files <- c(files, "dl_track_tbl.csv")
    transfer <- c(transfer, 10)
    destination <- rep(1, length(files))
    destination[transfer == 10] <- 10
  } else {
    files <- c(files, "tr_tb_updt.rdata")
    transfer <- c(transfer, 1)
    destination <- rep(0, length(files))
    destination[transfer == 10] <- 10
  }
  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params, "Files done trigger file",
                              "Trigger File Created")
  params <- store_stamp_entry(
    params,
    "R program execution complete, output files written",
    "Tracking Table")
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (job_failed)  {
    make_trigger("job_fail.ok",  params$write_path)
  } else {
    make_trigger("job_done.ok",    params$write_path)
  }
  make_trigger("files_done.ok", params$write_path)
  return(params)
}

send_pause_continue_2p <- function(params,
                                   files = c(),
                                   sleep_time = 10,
                                   max_waiting_time = NULL,
                                   job_started = FALSE) {
  params <- store_log_entry_2p(params, files)
  params <- StoreTrackingTableEntry.2p(params)
  WriteLogCSV(params)
  write_log_raw(params)
  files <- c(files, "stamps.rdata", "log.rdata", "file_list.csv")
  transfer <- c(rep(1, length(files) - 1), 10)
  if (params$party_name == "A") {
    files <- c(files, "dl_track_tbl.csv")
    transfer <- c(transfer, 10)
    destination <- rep(1, length(files))
    destination[transfer == 10] <- 10
  } else {
    files <- c("tr_tb_updt.rdata", files)
    transfer <- c(1, transfer)
    destination <- rep(0, length(files))
    destination[transfer == 10] <- 10
  }
  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params,
                              "Files done trigger file",
                              "Trigger File created")
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  params$pmn_step_counter      <- params$pmn_step_counter + 2
  if (job_started) {
    make_trigger("job_started.ok", params$write_path)
  } else {
    make_trigger("files_done.ok", params$write_path)
  }
  if (params$party_name == "A") {
    if (params$verbose) cat("Waiting for data partner\n")
  } else {
    if (params$verbose) cat("Waiting for analysis center\n")
  }
  Standby("files_done.ok", params$read_path,
          max_waiting_time = max_waiting_time,
          verbose = params$verbose)
  if (params$verbose) cat("Resuming local processing\n\n")
  delete_trigger("files_done.ok", params$read_path)
  params <- read_log_raw_2p(params)
  params <- new_log_entry_2p(params)
  params <- read_stamps_raw_2p(params)
  params <- store_stamp_entry(params,
                              "R program execution begins",
                              "Tracking Table")
  if (params$party_name == "A") {
    params <- read_tracking_table_update_2p(params)
  }
  return(params)
}


pause_continue_2p <- function(params, max_waiting_time) {
  params <- store_log_entry_2p(params, "")
  WriteLogCSV(params)
  if (params$party_name == "A") {
    if (params$verbose) cat("Waiting for data partner\n")
  } else {
    if (params$verbose) cat("Waiting for analysis center\n")
  }
  Standby("files_done.ok", params$read_path,
          max_waiting_time = max_waiting_time,
          verbose = params$verbose)
  if (params$verbose) cat("Resuming local processing\n\n")
  delete_trigger("files_done.ok", params$read_path)
  params <- merge_log_raw_2p(params)
  params <- new_log_entry_2p(params)
  params <- MergeStampsRaw.2p(params)
  params <- read_tracking_table_update_2p(params)
  WriteLogCSV(params)
  return(params)
}

###################### 3 PARTY PMN COMMUNICATION FUNCTIONS #####################

wait_for_turn_3p <- function(params, sleep_time) {
  Sys.sleep(sleep_time)
  if ((params$party_name == "T") || (!params$popmednet)) return(NULL)

  if (params$verbose) cat("Waiting For Turn\n")
  start_time <- proc.time()[3]
  if (params$verbose) cat("Elapsed Time:", HMS(0), "\r")

  if (exists("party_offset")) {
    if (params$verbose) cat("\n\n")
    return()
  }

  party_offset <- 15

  modulus   <- 2 * party_offset
  if (params$party_name == "A") targetTime <- 0
  if (params$party_name == "B") targetTime <- party_offset


  while (as.integer(Sys.time()) %% modulus != targetTime) {
    elapsed_time <- round(proc.time()[3] - start_time, 0)
    if (params$verbose) cat("Elapsed Time:", HMS(elapsed_time), "\r")
    Sys.sleep(0.1)
  }
  if (params$verbose) cat("\n\n")
}


send_pause_quit_3p <- function(params,
                               files_a = NULL,
                               files_b = NULL,
                               files_t = NULL,
                               sleep_time = 10,
                               job_failed = FALSE,
                               wait_for_turn = FALSE) {

  params$lastIteration <- TRUE
  params$completed     <- TRUE
  files <- c(files_a, files_b, files_t, "file_list.csv")
  transfer <- c(rep(1, length(files) - 1), 10)
  destination <- c(rep(1, length(files_a)),
                   rep(2, length(files_b)),
                   rep(0, length(files_t)),
                   10)
  if (params$party != "T") {
    files <- c(files, "stamps.rdata", "log.rdata")
    transfer <- c(transfer, 1, 1)
    destination <- c(destination, 0, 0)
  }
  if (params$party == "T") {
    if (job_failed) {
      files <- c(files, "job_fail.ok")
      params <- store_stamp_entry(params,
                                  "Job failed trigger file",
                                  "Trigger File created")
    } else {
      files <- c(files, "job_done.ok")
      params <- store_stamp_entry(params,
                                  "Job done trigger file",
                                  "Trigger File created")
    }
    transfer <- c(transfer, 10)
    destination <- c(destination, 10)
  }
  params <- store_log_entry_3p(params, c(files_a, files_b, files_t))
  params <- StoreTrackingTableEntry_3p(params)
  WriteLogCSV(params)
  write_log_raw(params)

  if (params$party_name == "T") {
    WriteTrackingTableCSV(params)
    files <- c(files, "dl_track_tbl.csv")
    transfer <- c(transfer, 10)
    destination <- c(destination, 10)
  } else {
    WriteTrackingTableRaw(params)
    files <- c(files, "tr_tb_updt.rdata")
    transfer <- c(transfer, 1)
    destination <- c(destination, 0)
  }
  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params,
                              "Files done trigger file",
                              "Trigger File Created")
  params <- store_stamp_entry(
    params,
    "R program execution complete, output files written",
    "Tracking Table")
  if (wait_for_turn) {
    params <- store_stamp_entry(params,
                                "R program execution delayed",
                                "Tracking Table")
    wait_for_turn_3p(params, sleep_time)
    params <- store_stamp_entry(params,
                                "R program execution restarted",
                                "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (params$party == "T") {
    if (job_failed)  {
      make_trigger("job_fail.ok",  params$write_path)
    } else {
      make_trigger("job_done.ok",  params$write_path)
    }
  }
  make_trigger("files_done.ok", params$write_path)
  return(params)
}

send_pause_continue_3p <- function(params,
                                   files_a = NULL,
                                   files_b = NULL,
                                   files_t = NULL,
                                   from   = NULL,
                                   sleep_time = 10,
                                   max_waiting_time = 24 * 60 * 60,
                                   job_started = FALSE,
                                   wait_for_turn = FALSE) {
  params <- store_log_entry_3p(params, c(files_a, files_b, files_t))
  params <- StoreTrackingTableEntry_3p(params)
  WriteLogCSV(params)
  write_log_raw(params)

  files <- c(files_a, files_b, files_t, "file_list.csv")
  transfer <- c(rep(1, length(files) - 1), 10)
  destination <- c(rep(1, length(files_a)),
                   rep(2, length(files_b)),
                   rep(0, length(files_t)), 10)
  if (length(files) > 1) {
    WriteTrackingTableRaw(params)
  }
  if (!is.null(files_a)) {
    files <- c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer <- c(transfer, 1, 1, 1)
    destination <- c(destination, 1, 1, 1)
  }
  if (!is.null(files_b)) {
    files <- c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer <- c(transfer, 1, 1, 1)
    destination <- c(destination, 2, 2, 2)
  }
  if (!is.null(files_t)) {
    files <- c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer <- c(transfer, 1, 1, 1)
    destination <- c(destination, 0, 0, 0)
  }
  if (params$party_name == "T") {
    WriteTrackingTableCSV(params)
    files <- c(files, "dl_track_tbl.csv")
    transfer    <- c(transfer, 10)
    destination <- c(destination, 10)
  }
  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params,
                              "Files done trigger file",
                              "Trigger File created")
  if (wait_for_turn) {
    params <- store_stamp_entry(params,
                                "R program execution delayed",
                                "Tracking Table")
    wait_for_turn_3p(params, sleep_time)
    params <- store_stamp_entry(params,
                                "R program execution restarted",
                                "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (job_started) {
    make_trigger("job_started.ok", params$write_path)
  } else {
    make_trigger("files_done.ok", params$write_path)
  }
  if (length(from) == 1) {
    if (from == "T") {
      if (params$verbose) cat("Waiting for analysis center\n")
    } else if (from == "A") {
      if (params$verbose) cat("Waiting for data partner 1\n")
    } else {
      if (params$verbose) cat("Waiting for data partner 2\n")
    }
  } else if (length(from) == 2) {
    if (params$verbose) cat("Waiting for data partners\n")
  }
  Standby("files_done.ok", params$read_path[from],
          max_waiting_time = max_waiting_time,
          verbose = params$verbose)
  if (params$verbose) cat("Resuming local processing\n\n")
  delete_trigger("files_done.ok", params$read_path[from])
  params <- merge_log_raw_3p(params, from)
  params <- UpdateCounters_3p(params)
  params <- new_log_entry_3p(params)
  params <- MergeStampsRaw_3p(params, from)
  params <- store_stamp_entry(params,
                              "R program execution begins",
                              "Tracking Table")
  params <- MergeTrackingTableRAW_3p(params, from)
  return(params)
}


pause_continue_3p <- function(params, from = NULL,
                              max_waiting_time = 24 * 60 * 60) {
  params <- store_log_entry_3p(params, "")
  params <- StoreTrackingTableEntry_3p(params)
  WriteLogCSV(params)
  if (length(from) == 1) {
    if (from == "T") {
      if (params$verbose) cat("Waiting for analysis center\n")
    } else if (from == "A") {
      if (params$verbose) cat("Waiting for data partner 1\n")
    } else {
      if (params$verbose) cat("Waiting for data partner 2\n")
    }
  } else if (length(from) == 2) {
    if (params$verbose) cat("Waiting for data partners\n")
  }
  Standby("files_done.ok", params$read_path[from],
          max_waiting_time = max_waiting_time,
          verbose = params$verbose)
  if (params$verbose) cat("Resuming local processing\n\n")
  delete_trigger("files_done.ok", params$read_path[from])
  params <- merge_log_raw_3p(params, from)
  params <- UpdateCounters_3p(params)
  params <- new_log_entry_3p(params)
  params <- MergeStampsRaw_3p(params, from)
  params <- MergeTrackingTableRAW_3p(params, from)
  WriteLogCSV(params)
  return(params)
}


UpdateCounters_3p <- function(params) {
  params$pmn_step_counter <- max(params$log$history$Step) + 1
  return(params)
}

###################### K PARTY PMN COMMUNICATION FUNCTIONS #####################

wait_for_turn.kp <- function(params, sleep_time) {
  Sys.sleep(sleep_time)

  if (!params$popmednet) return(NULL)

  if (params$verbose) cat("Waiting For Turn\n")
  start_time <- proc.time()[3]
  if (params$verbose) cat("Elapsed Time:", HMS(0), "\r")

  if (exists("party_offset")) {
    if (params$verbose) cat("\n\n")
    return()
  }

  party_offset <- 15

  modulus   <- (params$num_data_partners + 1) * party_offset
  targetTime <- params$data_partner_id * party_offset

  if (params$verbose) cat("Elapsed Time:", HMS(0), "\r")
  while (as.integer(Sys.time()) %% modulus != targetTime) {
    elapsed_time <- round(proc.time()[3] - start_time, 0)
    if (params$verbose) cat("Elapsed Time:", HMS(elapsed_time), "\r")
    Sys.sleep(0.1)
  }
  if (params$verbose) cat("\n\n")
}


send_pause_quit_kp <- function(params,
                               files_ac = NULL,
                               files_dp = NULL,
                               sleep_time = 10,
                               job_failed = FALSE,
                               wait_for_turn = FALSE) {

  # Assumes that upon quitting, same thing is sent to everyone, so files_dp
  # cannot be a list

  params$lastIteration <- TRUE
  params$completed     <- TRUE
  params <- store_log_entry_kp(params, c(files_ac, files_dp))
  params <- StoreTrackingTableEntry.kp(params)
  WriteLogCSV(params)
  write_log_raw(params)

  if (params$data_partner_id != 0) {
    WriteTrackingTableRaw(params)

    files_ac <- c(files_ac, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    if (!is.null(files_dp)) {
      files_dp <- c(files_dp, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    }

    dataPartnerTarget <- 1:params$num_data_partners
    dataPartnerTarget <- dataPartnerTarget[-params$data_partner_id]
    files <- c(files_ac, rep(files_dp, length(dataPartnerTarget)),
               "file_list.csv")
    transfer <- c(rep(1, length(files) - 1), 10)
    destination <- c(rep(0, length(files_ac)),
                     rep(dataPartnerTarget, each = length(files_dp)),
                     10)
  }

  if (params$data_partner_id == 0) {
    WriteTrackingTableCSV(params)
    files       <- c("dl_track_tbl.csv", "file_list.csv")
    if (job_failed) {
      files <- c(files, "job_fail.ok")
      params <- store_stamp_entry(params,
                                  "Job failed trigger file",
                                  "Trigger File created")
    } else {
      files <- c(files, "job_done.ok")
      params <- store_stamp_entry(params,
                                  "Job done trigger file",
                                  "Trigger File created")
    }
    transfer <- c(10, 10, 10)
    destination <- c(10, 10, 10)
  }

  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params, "Files done trigger file",
                              "Trigger File Created")
  params <- store_stamp_entry(
    params,
    "R program execution complete, output files written",
    "Tracking Table")
  if (wait_for_turn) {
    params <- store_stamp_entry(params,
                                "R program execution delayed",
                                "Tracking Table")
    wait_for_turn.kp(params, sleep_time)
    params <- store_stamp_entry(params,
                                "R program execution restarted",
                                "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (params$data_partner_id == 0) {
    if (job_failed)  {
      make_trigger("job_fail.ok",  params$write_path)
    } else {
      make_trigger("job_done.ok",  params$write_path)
    }
  }
  make_trigger("files_done.ok", params$write_path)
  return(params)
}


send_pause_continue_kp <- function(params,
                                   files_ac = NULL,
                                   files_dp = NULL,
                                   from   = NULL,
                                   sleep_time = 10,
                                   max_waiting_time = 24 * 60 * 60,
                                   job_started = FALSE,
                                   wait_for_turn = FALSE) {
  if (!is.list(files_dp)) {
    params <- store_log_entry_kp(params, c(files_ac, files_dp))
    params <- StoreTrackingTableEntry.kp(params)
    WriteLogCSV(params)
    write_log_raw(params)
    if (length(files_ac) + length(files_dp) > 0) {
      WriteTrackingTableRaw(params)
    }

    if (!is.null(files_ac)) {
      files_ac <- c(files_ac, "stamps.rdata",
                    "log.rdata",
                    "tr_tb_updt.rdata")
    }
    if (!is.null(files_dp)) {
      files_dp <- c(files_dp, "stamps.rdata",
                    "log.rdata",
                    "tr_tb_updt.rdata")
    }
    dataPartnerTarget <- 1:params$num_data_partners
    if (params$data_partner_id != 0) {
      dataPartnerTarget <- dataPartnerTarget[-params$data_partner_id]
    }

    files <- c(files_ac, rep(files_dp, length(dataPartnerTarget)),
               "file_list.csv")
    transfer <- c(rep(1, length(files) - 1), 10)
    destination <- c(rep(0, length(files_ac)),
                     rep(dataPartnerTarget, each = length(files_dp)),
                     10)
  } else {
    files <- files_ac
    for (dp in 1:params$num_data_partners) {
      files <- c(files, files_dp[[dp]])
    }
    params <- store_log_entry_kp(params, files)
    params <- StoreTrackingTableEntry.kp(params)
    WriteLogCSV(params)
    write_log_raw(params)
    if (length(files) > 0) {
      WriteTrackingTableRaw(params)
    }
    if (!is.null(files_ac)) {
      files_ac <- c(files_ac, "stamps.rdata",
                    "log.rdata",
                    "tr_tb_updt.rdata")
    }
    files <- files_ac
    transfer <- rep(1, length(files))
    destination <- rep(0, length(files_ac))
    for (dp in 1:params$num_data_partners) {
      if (length(files_dp[[dp]]) > 0 && dp != params$data_partner_id) {
        files <- c(files, files_dp[[dp]], "stamps.rdata",
                   "log.rdata",
                   "tr_tb_updt.rdata")
        transfer <- c(transfer, rep(1, length(files_dp[[dp]]) + 3))
        destination <- c(destination, rep(dp, length(files_dp[[dp]]) + 3))
      }
    }
    files <- c(files, "file_list.csv")
    transfer <- c(transfer, 10)
    destination <- c(destination, 10)
  }
  if (params$data_partner_id == 0) {
    WriteTrackingTableCSV(params)
    files <- c(files, "dl_track_tbl.csv")
    transfer    <- c(transfer, 10)
    destination <- c(destination, 10)
  }
  make_csv(files, transfer, destination, params$write_path)
  params <- store_stamp_entry(params,
                              "Files done trigger file",
                              "Trigger File created")
  if (wait_for_turn) {
    params <- store_stamp_entry(params,
                                "R program execution delayed",
                                "Tracking Table")
    wait_for_turn.kp(params, sleep_time)
    params <- store_stamp_entry(params,
                                "R program execution restarted",
                                "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (job_started) {
    make_trigger("job_started.ok", params$write_path)
  } else {
    make_trigger("files_done.ok", params$write_path)
  }
  if (from == "AC") {
    if (params$verbose) cat("Waiting for analysis center\n")
    Standby("files_done.ok", params$readPathAC,
            max_waiting_time = max_waiting_time,
            verbose = params$verbose)
    delete_trigger("files_done.ok", params$readPathAC)
  } else if (from == "DP") {
    if (params$verbose) cat("Waiting for data partners\n")
    if (params$data_partner_id == 0) {
      Standby("files_done.ok", params$readPathDP,
              max_waiting_time = max_waiting_time,
              verbose = params$verbose)
      delete_trigger("files_done.ok", params$readPathDP)
    } else {
      Standby("files_done.ok",
              params$readPathDP[-params$data_partner_id],
              max_waiting_time = max_waiting_time,
              verbose = params$verbose)
      delete_trigger("files_done.ok",
                     params$readPathDP[-params$data_partner_id])
    }
  } else if (from == "DP1") {
    if (params$verbose) cat("Waiting for data partner 1\n")
    Standby("files_done.ok",
            params$readPathDP[1],
            max_waiting_time = max_waiting_time,
            verbose = params$verbose)
    delete_trigger("files_done.ok", params$readPathDP[1])
  } else if (from == "DP2") {
    if (params$verbose) cat("Waiting for data partner 2\n")
    Standby("files_done.ok",
            params$readPathDP[2],
            max_waiting_time = max_waiting_time,
            verbose = params$verbose)
    delete_trigger("files_done.ok", params$readPathDP[2])
  }
  if (params$verbose) cat("Resuming local processing\n\n")
  params <- merge_log_raw_kp(params, from)
  params <- update_counters_kp(params)
  params <- new_log_entry_kp(params)
  params <- MergeStampsRaw.kp(params, from)
  params <- store_stamp_entry(params,
                              "R program execution begins",
                              "Tracking Table")
  params <- MergeTrackingTableRAW.kp(params, from)
  return(params)
}


pause_continue_kp <- function(params, from = NULL,
                              max_waiting_time = 24 * 60 * 60) {
  params <- store_log_entry_kp(params, "")
  params <- StoreTrackingTableEntry.kp(params)
  WriteLogCSV(params)
  params <- store_stamp_entry(params,
                              "R program execution paused",
                              "Tracking Table")
  if (from == "AC") {
    if (params$verbose) cat("Waiting for analysis center\n")
    Standby("files_done.ok",
            params$readPathAC,
            max_waiting_time = max_waiting_time,
            verbose = params$verbose)
    delete_trigger("files_done.ok", params$readPathAC)
  } else {
    if (params$verbose) cat("Waiting for data partners\n")
    if (params$data_partner_id == 0) {
      Standby("files_done.ok",
              params$readPathDP,
              max_waiting_time = max_waiting_time,
              verbose = params$verbose)
      delete_trigger("files_done.ok", params$readPathDP)
    } else {
      Standby("files_done.ok",
              params$readPathDP[-params$data_partner_id],
              max_waiting_time = max_waiting_time,
              verbose = params$verbose)
      delete_trigger("files_done.ok",
                     params$readPathDP[-params$data_partner_id])
    }
  }
  if (params$verbose) cat("Resuming local processing\n\n")
  params <- merge_log_raw_kp(params, from)
  params <- update_counters_kp(params)
  params <- new_log_entry_kp(params)
  params <- MergeStampsRaw.kp(params, from)
  params <- store_stamp_entry(params,
                              "R program execution begins",
                              "Tracking Table")
  params <- MergeTrackingTableRAW.kp(params, from)
  WriteLogCSV(params)
  return(params)
}


update_counters_kp <- function(params) {
  params$pmn_step_counter <- max(params$log$history$Step) + 1
  return(params)
}

received_error_kp <- function(params, from) {
  result <- list()
  message <- ""
  if (from == "AC") {
    message_exists <- file.exists(file.path(params$readPathAC,
                                           "errorMessage.rdata"))
    if (message_exists) {
      message <- read_error_message(params$readPathAC)
    }
  } else {
    message_exists <- file.exists(file.path(params$readPathDP,
                                           "errorMessage.rdata"))
    for (id in 1:params$num_data_partners) {
      if (message_exists[id]) {
        message <- paste0(message,
                          read_error_message(params$readPathDP[id]), " ")
      }
    }
  }
  result$error <- any(message_exists)
  result$message <- message
  return(result)
}

################################ TIME FUNCTIONS ################################

get_utc_time <- function() {
  t <- Sys.time()
  attr(t, "tzone") <- "UTC"
  return(as.POSIXlt(t))
}


get_utc_offset <- function() {
  t <- Sys.time()
  return(format(t, "%z"))
}


get_utc_offset_seconds <- function() {
  t <- Sys.time()
  offset <- format(t, "%z")
  hour <- as.numeric(substr(offset, 2, 3))
  min <- as.numeric(substr(offset, 4, 5))
  pm  <- ifelse(substr(offset, 1, 1) == "-", -1, 1)
  return(pm * (hour * 3600 + min * 60))
}


convert_utc_roundtrip_time <- function(t) {
  month <- ifelse(t$mon  < 9,  paste0("0", t$mon + 1), t$mon + 1)
  day   <- ifelse(t$mday < 10, paste0("0", t$mday),    t$mday)
  hour  <- ifelse(t$hour < 10, paste0("0", t$hour),    t$hour)
  min   <- ifelse(t$min  < 10, paste0("0", t$min),     t$min)
  sec   <- ifelse(t$sec  < 10, paste0("0", t$sec),     t$sec)
  t <- paste0(t$year + 1900, "-", month, "-", day, " ", hour, ":",
              min, ":", sec)
}

get_round_trip_time <- function() {
  return(convert_utc_roundtrip_time(get_utc_time()))
}


get_elapsed_time <- function(time1, final = FALSE, time_only = FALSE) {
  etime <- floor(time1[3])
  hrs <- floor(etime / 3600)
  mins <- floor((etime %% 3600) / 60)
  secs <- etime - hrs * 3600 - mins * 60

  hr1 <- if (hrs > 9) toString(hrs) else paste0("0", toString(hrs))
  min1 <- if (mins > 9) toString(mins) else paste0("0", toString(mins))
  sec1 <- if (secs > 9) toString(secs) else paste0("0", toString(secs))
  if (final) {
    return(paste0("(Total time elapsed: ", hr1, "hr ", min1, "min ",
                  sec1, "sec)"))
  } else if (time_only) {
    return(paste0("(", hr1, "hr  ", min1, "min  ", sec1, "sec)"))
  }
  return(paste0("(Time elapsed: ", hr1, "hr ", min1, "min ", sec1, "sec)"))
}


convert_secs_to_hms <- function(secs, final = FALSE, time_only = FALSE) {
  if (length(secs) != 1) {
    secs <- 0
  }
  secs <- round(secs, digits = 0)
  hrs <- floor(secs / 3600)
  mins <- floor((secs %% 3600) / 60)
  secs <- secs - hrs * 3600 - mins * 60

  hr1 <- if (hrs > 9) toString(hrs) else paste0("0", toString(hrs))
  min1 <- if (mins > 9) toString(mins) else paste0("0", toString(mins))
  sec1 <- if (secs > 9) toString(secs) else paste0("0", toString(secs))
  if (final) {
    return(paste0("(Total time elapsed: ", hr1, ":", min1, ":", sec1, ")"))
  }
  if (time_only) {
    return(paste0("(", hr1, ":", min1, ":", sec1, ")"))
  }
  return(paste0("(Time elapsed: ", hr1, ":", min1, ":", sec1, ")"))
}

HMS <- function(t) {
  paste(paste(formatC(t %/% (60 * 60), width = 2, format = "d", flag = "0"),
              formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
              formatC(t %% 60, width = 2, format = "d", flag = "0"),
              sep = ":"))
}


############################### BLOCK FUNCTIONS ################################

get_block_size <- function(p_a, p_b) {
  # This minimium is set up based on our current encryption scheme
  # We need to guarentee that p_a + p_b + g <= blocksize
  # where g = (p_a + 1) / (p_a + p_b + 1) * blocksize
  # May change in the future.
  min_blocksize <- max(25, trunc(1 + (p_a + p_b + 1)^2 / p_b))
  return(min_blocksize)
}


create_blocks <- function(p_a, p_b, n, blocksize) {

  # Divides the matrix with ncol(n) into submatrices of approximately
  # equal size. minimum size is blocksize.

  blocks <- list()

  num_blocks       <- max(trunc(n / blocksize), 1)
  newBlocksize    <- trunc(n / num_blocks)
  numBigBlocks    <- n %% newBlocksize
  numLittleBlocks <- num_blocks - numBigBlocks
  bigBlocksize    <- newBlocksize + 1
  gLittleBlock    <- trunc(newBlocksize * (p_a + 1) / (p_a + p_b + 1))
  gBigBlock       <- trunc(bigBlocksize * (p_a + 1) / (p_a + p_b + 1))

  blocks$num_blocks       <- num_blocks
  blocks$little_blocksize <- newBlocksize
  blocks$bigBlocksize    <- bigBlocksize
  blocks$numLittleBlocks <- numLittleBlocks
  blocks$numBigBlocks    <- numBigBlocks
  blocks$gLittleBlock    <- gLittleBlock
  blocks$gBigBlock       <- gBigBlock

  blocks$stops <- integer()
  if (numBigBlocks > 0) {
    blocks$stops <- bigBlocksize * 1:numBigBlocks
  }
  if (numLittleBlocks > 0) {
    blocks$stops <- c(blocks$stops, bigBlocksize * numBigBlocks +
                        newBlocksize * 1:numLittleBlocks)
  }

  if (num_blocks == 1) {
    blocks$starts <- c(1)
  } else {
    blocks$starts <- c(1, 1 + blocks$stops)[1:num_blocks]
  }

  blocks$g <- c(rep(gBigBlock, numBigBlocks),
                rep(gLittleBlock, numLittleBlocks))

  return(blocks)
}


create_containers <- function(p_a, p_b, blocks) {
  containers <- list()

  maximum_filesize <- 25 * 1024^2

  num_blocks <- blocks$num_blocks
  little_blocksize <- blocks$little_blocksize

  little_block_g <- blocks$gLittleBlock

  little_filesize_z   <- 8 * little_blocksize * little_block_g
  # used for w, v, RW, wr, rv, Cox
  little_filesize_w   <- 8 * little_blocksize * p_b
  little_filesize_rz  <- 8 * little_blocksize^2
  # I think this is not used anymore
  little_filesize_pr  <- 8 * (p_a + 1) * p_b
  little_filesize_xr <- 8 * p_a * p_b

  num_containers_z <- ceiling(num_blocks * little_filesize_z / maximum_filesize)
  num_blocks_small_container_z <- trunc(num_blocks / num_containers_z)
  num_blocks_large_container_z <- num_blocks_small_container_z + 1
  num_large_container_z <- num_blocks %% num_containers_z
  num_small_container_z <- num_containers_z - num_large_container_z

  num_containers_w <- ceiling(num_blocks * little_filesize_w / maximum_filesize)
  num_blocks_small_container_w <- trunc(num_blocks / num_containers_w)
  num_blocks_large_container_w <- num_blocks_small_container_w + 1
  num_large_container_w <- num_blocks %% num_containers_w
  num_small_container_w <- num_containers_w - num_large_container_w

  num_containers_rz <- ceiling(num_blocks * little_filesize_rz /
                                 maximum_filesize)
  num_blocks_small_containers_rz <- trunc(num_blocks / num_containers_rz)
  num_blocks_large_container_rz <- num_blocks_small_containers_rz + 1
  num_large_container_rz <- num_blocks %% num_containers_rz
  num_small_container_rz <- num_containers_rz - num_large_container_rz

  num_containers_pr <- ceiling(num_blocks * little_filesize_pr /
                                 maximum_filesize)
  num_blocks_small_container_pr <- trunc(num_blocks / num_containers_pr)
  num_blocks_large_container_pr <- num_blocks_small_container_pr + 1
  num_large_container_pr <- num_blocks %% num_containers_pr
  num_small_container_pr <- num_containers_pr - num_large_container_pr

  num_containers_xr <- ceiling(num_blocks * little_filesize_xr /
                                 maximum_filesize)
  num_blocks_small_containers_xr <- trunc(num_blocks / num_containers_xr)
  num_blocks_large_container_xr <- num_blocks_small_containers_xr + 1
  num_large_container_xr <- num_blocks %% num_containers_xr
  num_small_container_xr <- num_containers_xr - num_large_container_xr

  if (num_large_container_z > 0) {
    file_break_z <- c(0:(num_large_container_z - 1) *
                        num_blocks_large_container_z + 1,
                      0:(num_small_container_z - 1) *
                        num_blocks_small_container_z + 1 +
                        num_large_container_z * num_blocks_large_container_z)
  } else {
    file_break_z <- c(0:(num_small_container_z - 1) *
                        num_blocks_small_container_z + 1 +
                        num_large_container_z * num_blocks_large_container_z)
  }

  if (num_large_container_w > 0) {
    filebreak_w <- c(0:(num_large_container_w - 1) *
                       num_blocks_large_container_w + 1,
                     0:(num_small_container_w - 1) *
                       num_blocks_small_container_w + 1 +
                       num_large_container_w * num_blocks_large_container_w)
  } else {
    filebreak_w <- c(0:(num_small_container_w - 1) *
                       num_blocks_small_container_w + 1 +
                       num_large_container_w * num_blocks_large_container_w)
  }

  if (num_large_container_rz > 0) {
    filebreak_rz <- c(0:(num_large_container_rz - 1) *
                        num_blocks_large_container_rz + 1,
                      0:(num_small_container_rz - 1) *
                        num_blocks_small_containers_rz + 1 +
                        num_large_container_rz * num_blocks_large_container_rz)
  } else {
    filebreak_rz <- c(0:(num_small_container_rz - 1) *
                        num_blocks_small_containers_rz + 1 +
                        num_large_container_rz * num_blocks_large_container_rz)
  }

  if (num_large_container_pr > 0) {
    filebreak_pr <- c(0:(num_large_container_pr - 1) *
                        num_blocks_large_container_pr + 1,
                      0:(num_small_container_pr - 1) *
                        num_blocks_small_container_pr + 1 +
                        num_large_container_pr * num_blocks_large_container_pr)
  } else {
    filebreak_pr <- c(0:(num_small_container_pr - 1) *
                        num_blocks_small_container_pr + 1 +
                        num_large_container_pr * num_blocks_large_container_pr)
  }

  if (num_large_container_xr > 0) {
    filebreak_xr <- c(0:(num_large_container_xr - 1) *
                        num_blocks_large_container_xr + 1,
                      0:(num_small_container_xr - 1) *
                        num_blocks_small_containers_xr + 1 +
                        num_large_container_xr * num_blocks_large_container_xr)
  } else {
    filebreak_xr <- c(0:(num_small_container_xr - 1) *
                        num_blocks_small_containers_xr + 1 +
                        num_large_container_xr * num_blocks_large_container_xr)
  }

  containers$file_break_z  <- file_break_z
  containers$filebreak_w   <- filebreak_w
  containers$filebreak_rz  <- filebreak_rz
  # I think we are not using this anymore
  containers$filebreak_pr  <- filebreak_pr
  containers$filebreak_v   <- filebreak_w
  containers$filebreak_RW  <- filebreak_w
  containers$filebreak_wr  <- filebreak_w
  containers$filebreak_rv  <- filebreak_w
  containers$filebreak_vr  <- filebreak_w
  containers$filebreak_cox <- filebreak_w
  containers$filebreak_xr  <- filebreak_xr

  return(containers)
}

########################### OUTPUT FORMAT FUNCTIONS ############################

formatPValue <- function(pvals, width = 7) {
  p <- c()
  for (x in pvals) {
    if (is.na(x)) {
      x <- format("NA", width = width, justify = "right")
    } else if (x > 1e-3) {
      x <- format(round(x, 5),  width = width, justify = "right", nsmall = 5)
    } else if (x > 2e-16) {
      x <- formatC(x, format = "e", digits = 1)
    } else {
      x <- format("<2e-16", width = width, justify = "right")
    }
    p <- c(p, x)
  }
  return(p)
}


formatStrings <- function(x, minWidth = NULL, justify = "left") {
  width <- max(max(nchar(x)), minWidth)
  x <- format(x, justify = justify, width = width)
  return(x)
}


formatStat <- function(x) {
  if (is.na(x)) {
    "NA"
  } else if (x >= 1000000) {
    formatC(x, format = "e", digits = 3)
  } else if (x >= 1000) {
    as.character(round(x, 0))
  } else if (x > 1e-3) {
    format(signif(x, 4))
  } else {
    formatC(x, format = "e", digits = 3)
  }
}


formatStatList <- function(vals) {
  # Assumes that x is non-empty set of numeric or NA and there are no NaN's
  # width = 10, justify = right => standard output, so no worries about justify
  # nor width
  notNA <- which(!is.na(vals))
  notZero <- which(vals != 0)
  keep <- intersect(notNA, notZero)
  if (length(keep) == 0) {
    f <- c()
    for (x in vals) {
      if (is.na(x)) {
        f <- c(f, "NA")
      } else {
        f <- c(f, "0")
      }
    }
    return(f)
  }
  temp <- vals[keep]  # All non-zero, non-NA
  minval <- min(abs(temp))
  maxval <- max(abs(temp))
  #Where most significant digit is located
  decmin <- floor(log10(minval))
  decmax <- floor(log10(maxval))
  if (minval >= 1) decmin <- decmin + 1
  if (maxval >= 1) decmax <- decmax + 1
  if ((decmin < -3) || (decmax > 6) || (decmax - decmin > 3)) {
    # scientific
    f <- c()
    for (x in vals) {
      if (is.na(x)) {
        f <- c(f, "NA")
      } else {
        f <- c(f, formatC(x, format = "e", digits = 3))
      }
    }
    return(f)
  } else {
    # standard
    if (decmin < 0) {
      nsmall <- 6
    } else {
      nsmall <- 6 - decmin
    }
    f <- c()
    for (x in vals) {
      if (is.na(x)) {
        f <- c(f, "NA")
      } else {
        f <- c(f, format(round(x, nsmall), scientific = FALSE, nsmall = nsmall))
      }
    }
    return(f)
  }
}

########################### SHARED STAMPS FUNCTIONS ############################

store_stamp_entry <- function(params, description = "", type = "") {
  newEntry             <- params$stamps$blank
  newEntry$Step        <- params$pmn_step_counter
  newEntry$Description <- description
  newEntry$Time        <- get_round_trip_time()
  newEntry$Type        <- type
  params$stamps$history <- rbind(params$stamps$history, newEntry)
  return(params)
}


WriteStampsRaw <- function(params) {
  stamps <- params$stamps$history
  save(stamps, file = file.path(params$write_path, "stamps.rdata"))
}


#' @importFrom utils write.csv
WriteStampsCSV <- function(params) {
  write.csv(params$stamps$history, file.path(params$write_path, "stamps.csv"),
            row.names = FALSE)
}

########################### 2 PARTY STAMPS FUNCTIONS ###########################

initialize_time_stamps_2p <- function(params) {
  stamps <- list()
  stamps$blank <- data.frame(
    Step        = params$pmn_step_counter,
    Source      = paste("Org", params$party_name, "Dist Reg"),
    Description = "R program execution begins",
    Time        = get_round_trip_time(),
    Type        = "Tracking Table")
  stamps$history <- stamps$blank
  params$stamps <- stamps
  return(params)
}


read_stamps_raw_2p <- function(params) {
  stamps <- NULL
  load(file.path(params$read_path, "stamps.rdata"))
  params$stamps$history <- stamps
  return(params)
}


MergeStampsRaw.2p <- function(params) {
  # This function will only be used in the function Pause Continue
  # When party A and party B run simultaneously, but Party A can run first
  # even if Party B starts the whole thing.  We append party B's log
  # to the end of Party A's log.
  stamps <- NULL
  load(file.path(params$read_path, "stamps.rdata"))
  params$stamps$history <- rbind(params$stamps$history, stamps)
  return(params)
}

########################### 3 PARTY STAMPS FUNCTIONS ###########################

initialize_time_stamps_3p <- function(params) {
  stamps <- list()
  stamps$blank <- data.frame(
    Step        = params$pmn_step_counter,
    Source      = paste("Org", params$party_name, "Dist Reg"),
    Description = "R program execution begins",
    Time        = get_round_trip_time(),
    Type        = "Tracking Table")
  stamps$history <- stamps$blank
  params$stamps <- stamps
  return(params)
}


MergeStampsRaw_3p <- function(params, from) {
  stamps <- NULL
  for (party in from) {
    load(file.path(params$read_path[[party]], "stamps.rdata"))
    key1 <- paste0(params$stamps$history$Step,
                   params$stamps$history$Source,
                   params$stamps$history$Description)
    key2 <- paste0(stamps$Step,
                   stamps$Source,
                   stamps$Description)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$stamps$history <- rbind(params$stamps$history, stamps)
    } else if (length(idx) < length(key2)) {
      params$stamps$history <- rbind(params$stamps$history, stamps[-idx, ])
    }
  }
  idx <- order(as.character(params$stamps$history$Time))
  params$stamps$history <- params$stamps$history[idx, ]
  return(params)
}

########################### K PARTY STAMPS FUNCTIONS ###########################

initialize_time_stamps_kp <- function(params) {
  stamps <- list()
  stamps$blank <- data.frame(Step        = params$pmn_step_counter,
                             Source      = paste0("Org dp",
                                                  params$data_partner_id,
                                                  " Dist Reg"),
                             Description = "R program execution begins",
                             Time        = get_round_trip_time(),
                             Type        = "Tracking Table")
  stamps$history <- stamps$blank
  params$stamps <- stamps
  return(params)
}


MergeStampsRaw.kp <- function(params, from) {
  stamps <- NULL
  if (from == "AC") {
    load(file.path(params$readPathAC, "stamps.rdata"))
    key1 <- paste0(params$stamps$history$Step,
                   params$stamps$history$Source,
                   params$stamps$history$Description)
    key2 <- paste0(stamps$Step,
                   stamps$Source,
                   stamps$Description)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$stamps$history <- rbind(params$stamps$history, stamps)
    } else if (length(idx) < length(key2)) {
      params$stamps$history <- rbind(params$stamps$history, stamps[-idx, ])
    }
  } else if (from == "DP1") {
    load(file.path(params$readPathDP[1], "stamps.rdata"))
    key1 <- paste0(params$stamps$history$Step,
                   params$stamps$history$Source,
                   params$stamps$history$Description)
    key2 <- paste0(stamps$Step,
                   stamps$Source,
                   stamps$Description)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$stamps$history <- rbind(params$stamps$history, stamps)
    } else if (length(idx) < length(key2)) {
      params$stamps$history <- rbind(params$stamps$history, stamps[-idx, ])
    }
  } else {
    for (id in 1:params$num_data_partners) {
      if (id == params$data_partner_id) next
      load(file.path(params$readPathDP[id], "stamps.rdata"))
      key1 <- paste0(params$stamps$history$Step,
                     params$stamps$history$Source,
                     params$stamps$history$Description)
      key2 <- paste0(stamps$Step,
                     stamps$Source,
                     stamps$Description)
      idx <- which(key2 %in% key1)
      if (length(idx) == 0) {
        params$stamps$history <- rbind(params$stamps$history, stamps)
      } else if (length(idx) < length(key2)) {
        params$stamps$history <- rbind(params$stamps$history, stamps[-idx, ])
      }
    }
  }
  idx <- order(as.character(params$stamps$history$Time))
  params$stamps$history <- params$stamps$history[idx, ]
  return(params)
}

############################# SHARED LOG FUNCTIONS #############################

add_to_log <- function(params, functionName, readTime, readSize,
                       writeTime, writeSize) {
  readTime  <- round(as.numeric(readTime),  digits = 2)
  writeTime <- round(as.numeric(writeTime), digits = 2)
  readSize  <- round(as.numeric(readSize),  digits = 0)
  writeSize <- round(as.numeric(writeSize), digits = 0)
  if (params$log$current$functions == "") {
    params$log$current$functions <- functionName
  } else {
    params$log$current$functions <- paste0(params$log$current$functions,
                                           ", ", functionName)
  }
  params$log$current$read_time  <- params$log$current$read_time + readTime
  params$log$current$read_Size  <- params$log$current$read_Size + readSize
  params$log$current$write_time <- params$log$current$write_time + writeTime
  params$log$current$write_Size <- params$log$current$write_Size + writeSize
  return(params)
}


write_log_raw <- function(params) {
  log <- params$log$history
  save(log, file = file.path(params$write_path, "log.rdata"))
}


#' @importFrom utils write.csv
WriteLogCSV <- function(params) {
  write.csv(params$log$history, file.path(params$write_path, "log.csv"),
            row.names = FALSE)
}


#' @importFrom utils write.table
write_to_log_summary <- function(c1 = "", c2 = "", c3 = "",
                                 write_path = NULL, append = TRUE) {
  if (is.numeric(c2)) {
    c2 <- round(c2, 2)
  }
  write.table(data.frame(c1, c2, c3),
              file.path(write_path, "log_summary.csv"), sep = ",",
              col.names = FALSE,
              row.names = FALSE, append = append)
}

############################# 2 PARTY LOG FUNCTIONS ############################

initialize_log_2p <- function(params) {
  log <- list()
  log$blank <- data.frame(Step             = 0,
                          iteration_alg    = 0,
                          Party            = "",
                          functions        = "",
                          wait_time        = 0,
                          start_time       = get_utc_time(),
                          end_time         = get_utc_time(),
                          read_time        = 0,
                          read_Size        = 0,
                          write_time       = 0,
                          write_Size       = 0,
                          computetation_time = 0,
                          files_sent       = "",
                          bytes_sent       = 0)
  log$current <- log$blank
  log$history <- log$blank
  params$log <- log
  return(params)
}


new_log_entry_2p <- function(params) {
  params$log$current <- params$log$blank
  params$log$current$Party         <- params$party_name
  params$log$current$start_time    <- get_utc_time()
  return(params)
}


store_log_entry_2p <- function(params, files) {
  params$log$current$Step          <- params$pmn_step_counter
  params$log$current$iteration_alg <- params$algIterationCounter
  params$log$current$Party <- params$party_name
  params$log$current$end_time <- get_utc_time()
  params$log$current$computetation_time <-
    round(as.numeric(difftime(
      params$log$current$end_time,
      params$log$current$start_time, units = "secs")) -
        params$log$current$read_time - params$log$current$write_time, 2)
  params$log$current$files_sent <- paste(files, collapse = ", ")
  params$log$current$bytes_sent <-
    sum(file.size(file.path(params$write_path, files)))
  if (is.na(params$log$current$bytes_sent)) {
    params$log$current$bytes_sent <- 0
  }
  nrows <- nrow(params$log$history)
  if (nrows >= 2) {
    params$log$current$wait_time <-
      round(as.numeric(difftime(
        params$log$current$start_time,
        max(params$log$history$end_time[which(params$log$history$Party ==
                                                params$log$current$Party)]),
        units = "secs")), 2)
  }
  if (params$log$history$Party[nrows] == "") {
    params$log$history <- params$log$current
  } else {
    params$log$history <- rbind(params$log$history, params$log$current)
  }
  nrows <- nrow(params$log$history)
  return(params)
}

read_log_raw_2p <- function(params) {
  load(file.path(params$read_path, "log.rdata"))
  params$log$history <- log
  return(params)
}


merge_log_raw_2p <- function(params) {
  # This function will only be used in the function Pause Continue
  # When party A and party B run simultaneously, but Party A can run first
  # even if Party B starts the whole thing.  We append party B's log
  # to the end of Party A's log.
  load(file.path(params$read_path, "log.rdata"))
  params$log$history <- rbind(params$log$history, log)
  return(params)
}


summarize_log_2p <- function(params) {
  write_path <- params$write_path
  log    <- params$log$history
  indexA <- which(log$Party == "A")
  indexB <- which(log$Party == "B")
  party_a_start_time <- log$start_time[indexA[1]]
  party_a_end_time   <- log$end_time[indexA[length(indexA)]]
  party_a_total_time <- round(as.numeric(difftime(
    party_a_end_time, party_a_start_time, units = "secs")), digits = 2)
  party_a_reading_time <- sum(log$read_time[indexA])
  party_a_writing_time <- sum(log$write_time[indexA])
  party_a_computing_time <- sum(log$computetation_time[indexA])
  party_a_waiting_time <- sum(log$wait_time[indexA])
  party_a_total_time_hms <-
    convert_secs_to_hms(party_a_total_time, time_only = TRUE)
  party_a_reading_time_hms <-
    convert_secs_to_hms(party_a_reading_time, time_only = TRUE)
  party_a_writing_time_hms <-
    convert_secs_to_hms(party_a_writing_time, time_only = TRUE)
  party_a_computing_time_hms <-
    convert_secs_to_hms(party_a_computing_time, time_only = TRUE)
  party_a_waiting_time_hms <-
    convert_secs_to_hms(party_a_waiting_time, time_only = TRUE)
  party_a_Bytes_read <- sum(log$read_Size[indexA])
  party_a_Bytes_written <- sum(log$write_Size[indexA])

  party_b_start_time <- log$start_time[indexB[1]]
  party_b_end_time   <- log$end_time[indexB[length(indexB)]]
  party_b_total_time <- round(as.numeric(difftime(
    party_b_end_time, party_b_start_time, units = "secs")), digits = 2)
  party_b_reading_time <- sum(log$read_time[indexB])
  party_b_writing_time <- sum(log$write_time[indexB])
  party_b_computing_time <- sum(log$computetation_time[indexB])
  party_b_waiting_time <- party_b_total_time - party_b_reading_time -
    party_b_writing_time - party_b_computing_time
  party_b_total_time_hms <-
    convert_secs_to_hms(party_b_total_time, time_only = TRUE)
  party_b_reading_time_hms <-
    convert_secs_to_hms(party_b_reading_time, time_only = TRUE)
  party_b_writing_time_hms <-
    convert_secs_to_hms(party_b_writing_time, time_only = TRUE)
  party_b_computing_time_hms <-
    convert_secs_to_hms(party_b_computing_time, time_only = TRUE)
  party_b_waiting_time_hms <-
    convert_secs_to_hms(party_b_waiting_time, time_only = TRUE)
  party_b_Bytes_read <- sum(log$read_Size[indexB])
  party_b_Bytes_written <- sum(log$write_Size[indexB])

  total_transfer_time <- 0
  if (max(log$Step) > 1) {
    for (i in 2:max(log$Step)) {
      idx1 <- which(log$Step == i - 1)
      idx2 <- which(log$Step == i)
      total_transfer_time <- total_transfer_time +
        as.numeric(difftime(min(log$start_time[idx2]),
                            max(log$end_time[idx1]), units = "secs"))
    }
  }
  total_transfer_time <- round(total_transfer_time, 2)

  total_reading_time <- sum(log$read_time)
  total_writing_time <- sum(log$write_time)
  total_computing_time <- sum(log$computetation_time)
  elapsed_computing_time <-
    party_a_total_time - total_transfer_time
  total_reading_time_hms <-
    convert_secs_to_hms(total_reading_time, time_only = TRUE)
  total_writing_time_hms <-
    convert_secs_to_hms(total_writing_time, time_only = TRUE)
  total_computing_time_hms <-
    convert_secs_to_hms(total_computing_time, time_only = TRUE)
  elapsed_computing_time_hms <-
    convert_secs_to_hms(elapsed_computing_time, time_only = TRUE)
  total_transfer_time_hms <-
    convert_secs_to_hms(total_transfer_time, time_only = TRUE)
  total_Bytes_transferred <- sum(log$bytes_sent)
  kb_per_second <- round(total_Bytes_transferred /
                           (total_transfer_time * 1024), digits = 2)
  write_to_log_summary(c1 = "Analysis",
                       c2 = params$analysis,
                       write_path = write_path, append = FALSE)
  if (!is.null(params$blocks)) {
    write_to_log_summary(c1 = "Blocksize",
                         c2 = params$blocks$little_blocksize,
                         write_path = write_path)
    write_to_log_summary(c1 = "Number of Blocks",
                         c2 = params$blocks$numLittleBlocks +
                           params$blocks$numBigBlocks,
                         write_path = write_path)
  }
  if (!is.null(params$n))   write_to_log_summary(c1 = "N",
                                                 c2 = params$n,
                                                 write_path = write_path)

  p <- max(0, params$p1_old - (params$analysis != "cox"))
  write_to_log_summary(c1 = "p_a", c2 = p, write_path = write_path)
  p <- params$p2_old
  write_to_log_summary(c1 = "p_b", c2 = p, write_path = write_path)

  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Party A Start Time",
                       c2 = party_a_start_time, write_path = write_path)
  write_to_log_summary(c1 = "Party A End Time",
                       c2 = party_a_end_time, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Run Time",
                       c2 = party_a_total_time,
                       c3 = party_a_total_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Reading Time",
                       c2 = party_a_reading_time,
                       c3 = party_a_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Bytes Read",
                       c2 = party_a_Bytes_read, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Writing Time",
                       c2 = party_a_writing_time,
                       c3 = party_a_writing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Bytes Written",
                       c2 = party_a_Bytes_written, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Computing Time",
                       c2 = party_a_computing_time,
                       c3 = party_a_computing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Waiting Time",
                       c2 = party_a_waiting_time,
                       c3 = party_a_waiting_time_hms, write_path = write_path)
  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Party B Start Time",
                       c2 = party_b_start_time, write_path = write_path)
  write_to_log_summary(c1 = "Party B End Time",
                       c2 = party_b_end_time, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Run Time",
                       c2 = party_b_total_time,
                       c3 = party_b_total_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Reading Time",
                       c2 = party_b_reading_time,
                       c3 = party_b_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Bytes Read",
                       c2 = party_b_Bytes_read, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Writing Time",
                       c2 = party_b_writing_time,
                       c3 = party_b_writing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Bytes Written",
                       c2 = party_b_Bytes_written, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Computing Time",
                       c2 = party_b_computing_time,
                       c3 = party_b_computing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Waiting Time",
                       c2 = party_b_waiting_time,
                       c3 = party_b_waiting_time_hms, write_path = write_path)
  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Total Reading Time",
                       c2 = total_reading_time,
                       c3 = total_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Total Writing Time",
                       c2 = total_writing_time,
                       c3 = total_writing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Computing Time",
                       c2 = total_computing_time,
                       c3 = total_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Elapsed Computing Time",
                       c2 = elapsed_computing_time,
                       c3 = elapsed_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Transfer Time",
                       c2 = total_transfer_time,
                       c3 = total_transfer_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Total Bytes Transferred",
                       c2 = total_Bytes_transferred, write_path = write_path)
  write_to_log_summary(c1 = "KB / Sec Transfer Rate",
                       c2 = kb_per_second, write_path = write_path)

}

############################# 3 PARTY LOG FUNCTIONS ############################

initialize_log_3p <- function(params) {
  log <- list()
  log$blank <- data.frame(Step             = 0,
                          iteration_alg    = 0,
                          Party            = "",
                          functions        = "",
                          wait_time        = 0,
                          start_time       = get_utc_time(),
                          end_time         = get_utc_time(),
                          read_time        = 0,
                          read_Size        = 0,
                          write_time       = 0,
                          write_Size       = 0,
                          computetation_time = 0,
                          files_sent       = "",
                          bytes_sent       = 0)
  log$current <- log$blank
  log$history <- log$blank
  params$log <- log
  return(params)
}


new_log_entry_3p <- function(params) {
  params$log$current <- params$log$blank
  params$log$current$Party         <- params$party_name
  params$log$current$start_time    <- get_utc_time()
  return(params)
}


store_log_entry_3p <- function(params, files) {
  params$log$current$Step          <- params$pmn_step_counter
  params$log$current$iteration_alg <- params$algIterationCounter
  params$log$current$Party <- params$party_name
  params$log$current$end_time <- get_utc_time()
  params$log$current$computetation_time <- round(as.numeric(difftime(
    params$log$current$end_time,
    params$log$current$start_time, units = "secs")) -
      params$log$current$read_time - params$log$current$write_time, 2)
  params$log$current$files_sent <- paste(files, collapse = ", ")
  params$log$current$bytes_sent <-
    sum(file.size(file.path(params$write_path, files)))
  if (is.na(params$log$current$bytes_sent)) {
    params$log$current$bytes_sent <- 0
  }
  nrows <- nrow(params$log$history)
  if (nrows >= 3) {
    params$log$current$wait_time <-
      round(as.numeric(difftime(
        params$log$current$start_time,
        max(params$log$history$end_time[which(params$log$history$Party ==
                                                params$log$current$Party)]),
        units = "secs")), 2)
  }
  if (params$log$history$Party[nrows] == "") {
    params$log$history <- params$log$current
  } else {
    params$log$history <- rbind(params$log$history, params$log$current)
  }
  return(params)
}


merge_log_raw_3p <- function(params, from) {
  # This function will only be used in the function Pause Continue
  # When party A and party B run simultaneously, but Party A can run first
  # even if Party B starts the whole thing.  We append party B's log
  # to the end of Party A's log.
  for (party in from) {
    load(file.path(params$read_path[[party]], "log.rdata"))
    key1 <- paste0(params$log$history$Step, params$log$history$Party)
    key2 <- paste0(log$Step, log$Party)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$log$history <- rbind(params$log$history, log)
    } else if (length(idx) < length(key2)) {
      params$log$history <- rbind(params$log$history, log[-idx, ])
    }
  }
  idx <- order(params$log$history$Step, params$log$history$Party)
  params$log$history <- params$log$history[idx, ]

  return(params)
}


summarize_log_3p <- function(params) {
  write_path <- params$write_path

  log    <- params$log$history
  indexA <- which(log$Party == "A")
  indexB <- which(log$Party == "B")
  indexT <- which(log$Party == "T")
  party_a_start_time <- log$start_time[indexA[1]]
  party_a_end_time   <- log$end_time[indexA[length(indexA)]]
  party_a_total_time <- round(as.numeric(difftime(
    party_a_end_time, party_a_start_time, units = "secs")), digits = 2)
  party_a_reading_time <- sum(log$read_time[indexA])
  party_a_writing_time <- sum(log$write_time[indexA])
  party_a_computing_time <- sum(log$computetation_time[indexA])
  party_a_waiting_time <- sum(log$wait_time[indexA])
  party_a_total_time_hms <-
    convert_secs_to_hms(party_a_total_time, time_only = TRUE)
  party_a_reading_time_hms <-
    convert_secs_to_hms(party_a_reading_time, time_only = TRUE)
  party_a_writing_time_hms <-
    convert_secs_to_hms(party_a_writing_time, time_only = TRUE)
  party_a_computing_time_hms <-
    convert_secs_to_hms(party_a_computing_time, time_only = TRUE)
  party_a_waiting_time_hms <-
    convert_secs_to_hms(party_a_waiting_time, time_only = TRUE)
  party_a_Bytes_read <- sum(log$read_Size[indexA])
  party_a_Bytes_written <- sum(log$write_Size[indexA])

  party_b_start_time <- log$start_time[indexB[1]]
  party_b_end_time   <- log$end_time[indexB[length(indexB)]]
  party_b_total_time <- round(as.numeric(difftime(
    party_b_end_time, party_b_start_time, units = "secs")), digits = 2)
  party_b_reading_time <- sum(log$read_time[indexB])
  party_b_writing_time <- sum(log$write_time[indexB])
  party_b_computing_time <- sum(log$computetation_time[indexB])
  party_b_waiting_time <- sum(log$wait_time[indexB])
  party_b_total_time_hms <-
    convert_secs_to_hms(party_b_total_time, time_only = TRUE)
  party_b_reading_time_hms <-
    convert_secs_to_hms(party_b_reading_time, time_only = TRUE)
  party_b_writing_time_hms <-
    convert_secs_to_hms(party_b_writing_time, time_only = TRUE)
  party_b_computing_time_hms <-
    convert_secs_to_hms(party_b_computing_time, time_only = TRUE)
  party_b_waiting_time_hms <-
    convert_secs_to_hms(party_b_waiting_time, time_only = TRUE)
  party_b_Bytes_read <- sum(log$read_Size[indexB])
  party_b_Bytes_written <- sum(log$write_Size[indexB])

  party_t_start_time <- log$start_time[indexT[1]]
  party_t_end_time   <- log$end_time[indexT[length(indexT)]]
  party_t_total_time <- round(as.numeric(difftime(
    party_t_end_time, party_t_start_time, units = "secs")), digits = 2)
  party_t_reading_time <- sum(log$read_time[indexT])
  party_t_writing_time <- sum(log$write_time[indexT])
  party_t_computing_time <- sum(log$computetation_time[indexT])
  party_t_waiting_time <- sum(log$wait_time[indexT])
  party_t_total_time_hms <-
    convert_secs_to_hms(party_t_total_time, time_only = TRUE)
  party_t_reading_time_hms <-
    convert_secs_to_hms(party_t_reading_time, time_only = TRUE)
  party_t_writing_time_hms <-
    convert_secs_to_hms(party_t_writing_time, time_only = TRUE)
  party_t_computing_time_hms <-
    convert_secs_to_hms(party_t_computing_time, time_only = TRUE)
  party_t_waiting_time_hms <-
    convert_secs_to_hms(party_t_waiting_time, time_only = TRUE)
  party_t_Bytes_read <- sum(log$read_Size[indexT])
  party_t_Bytes_written <- sum(log$write_Size[indexT])

  total_transfer_time <- 0
  if (max(log$Step) > 1) {
    for (i in 2:max(log$Step)) {
      idx1 <- which(log$Step == i - 1)
      idx2 <- which(log$Step == i)
      total_transfer_time <- total_transfer_time +
        as.numeric(difftime(min(log$start_time[idx2]),
                            max(log$end_time[idx1]), units = "secs"))
    }
  }
  total_transfer_time <- round(total_transfer_time, 2)
  elapsed_computing_time <- party_t_total_time - total_transfer_time

  total_reading_time <- sum(log$read_time)
  total_writing_time <- sum(log$write_time)
  total_computing_time <- sum(log$computetation_time)
  total_reading_time_hms <-
    convert_secs_to_hms(total_reading_time, time_only = TRUE)
  total_writing_time_hms <-
    convert_secs_to_hms(total_writing_time, time_only = TRUE)
  total_transfer_time_hms <-
    convert_secs_to_hms(total_transfer_time, time_only = TRUE)
  total_computing_time_hms <-
    convert_secs_to_hms(total_computing_time, time_only = TRUE)
  elapsed_computing_time_hms <-
    convert_secs_to_hms(elapsed_computing_time, time_only = TRUE)
  total_Bytes_transferred <- sum(log$bytes_sent)
  kb_per_second <- round(total_Bytes_transferred /
                           (total_transfer_time * 1024), digits = 2)

  write_to_log_summary(c1 = "Analysis",
                       c2 = params$analysis,
                       write_path = write_path,
                       append = FALSE)
  if (!is.null(params$blocks)) {
    write_to_log_summary(c1 = "Blocksize",
                         c2 = params$blocks$little_blocksize,
                         write_path = write_path)
    write_to_log_summary(c1 = "Number of Blocks",
                         c2 = params$blocks$numLittleBlocks +
                           params$blocks$numBigBlocks,
                         write_path = write_path)
  }
  if (!is.null(params$n))   write_to_log_summary(c1 = "N",
                                                 c2 = params$n,
                                                 write_path = write_path)

  p <- max(0, params$p1_old - (params$analysis != "cox"))
  write_to_log_summary(c1 = "p_a", c2 = p, write_path = write_path)
  p <- params$p2_old
  write_to_log_summary(c1 = "p_b", c2 = p, write_path = write_path)

  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Party A Start Time",
                       c2 = party_a_start_time, write_path = write_path)
  write_to_log_summary(c1 = "Party A End Time",
                       c2 = party_a_end_time, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Run Time",
                       c2 = party_a_total_time,
                       c3 = party_a_total_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Reading Time",
                       c2 = party_a_reading_time,
                       c3 = party_a_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Bytes Read",
                       c2 = party_a_Bytes_read, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Writing Time",
                       c2 = party_a_writing_time,
                       c3 = party_a_writing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Bytes Written",
                       c2 = party_a_Bytes_written, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Computing Time",
                       c2 = party_a_computing_time,
                       c3 = party_a_computing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party A Total Waiting Time",
                       c2 = party_a_waiting_time,
                       c3 = party_a_waiting_time_hms, write_path = write_path)
  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Party B Start Time",
                       c2 = party_b_start_time, write_path = write_path)
  write_to_log_summary(c1 = "Party B End Time",
                       c2 = party_b_end_time, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Run Time",
                       c2 = party_b_total_time,
                       c3 = party_b_total_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Reading Time",
                       c2 = party_b_reading_time,
                       c3 = party_b_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Bytes Read",
                       c2 = party_b_Bytes_read, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Writing Time",
                       c2 = party_b_writing_time,
                       c3 = party_b_writing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Bytes Written",
                       c2 = party_b_Bytes_written, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Computing Time",
                       c2 = party_b_computing_time,
                       c3 = party_b_computing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party B Total Waiting Time",
                       c2 = party_b_waiting_time,
                       c3 = party_b_waiting_time_hms, write_path = write_path)
  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Party T Start Time",
                       c2 = party_t_start_time, write_path = write_path)
  write_to_log_summary(c1 = "Party T End Time",
                       c2 = party_t_end_time, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Run Time",
                       c2 = party_t_total_time,
                       c3 = party_t_total_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Reading Time",
                       c2 = party_t_reading_time,
                       c3 = party_t_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Bytes Read",
                       c2 = party_t_Bytes_read, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Writing Time",
                       c2 = party_t_writing_time,
                       c3 = party_t_writing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Bytes Written",
                       c2 = party_t_Bytes_written, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Computing Time",
                       c2 = party_t_computing_time,
                       c3 = party_t_computing_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Party T Total Waiting Time",
                       c2 = party_t_waiting_time,
                       c3 = party_t_waiting_time_hms, write_path = write_path)
  write_to_log_summary(write_path = write_path)
  write_to_log_summary(c1 = "Total Reading Time",
                       c2 = total_reading_time,
                       c3 = total_reading_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Total Writing Time",
                       c2 = total_writing_time,
                       c3 = total_writing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Computing Time",
                       c2 = total_computing_time,
                       c3 = total_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Elapsed Computing Time",
                       c2 = elapsed_computing_time,
                       c3 = elapsed_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Transfer Time",
                       c2 = total_transfer_time,
                       c3 = total_transfer_time_hms, write_path = write_path)
  write_to_log_summary(c1 = "Total Bytes Transferred",
                       c2 = total_Bytes_transferred, write_path = write_path)
  write_to_log_summary(c1 = "KB / Sec Transfer Rate",
                       c2 = kb_per_second, write_path = write_path)

}

############################# K PARTY LOG FUNCTIONS ############################

initialize_log_kp <- function(params) {
  log <- list()
  log$blank <- data.frame(Step             = 0,
                          iteration_alg    = 0,
                          Party            = "",
                          functions        = "",
                          wait_time        = 0,
                          start_time       = get_utc_time(),
                          end_time         = get_utc_time(),
                          read_time        = 0,
                          read_Size        = 0,
                          write_time       = 0,
                          write_Size       = 0,
                          computetation_time = 0,
                          files_sent       = "",
                          bytes_sent       = 0)
  log$current <- log$blank
  log$history <- log$blank
  params$log <- log
  return(params)
}


new_log_entry_kp <- function(params) {
  params$log$current <- params$log$blank
  params$log$current$Party         <- paste0("dp", params$data_partner_id)
  params$log$current$start_time    <- get_utc_time()
  return(params)
}


store_log_entry_kp <- function(params, files) {
  params$log$current$Step          <- params$pmn_step_counter
  params$log$current$iteration_alg <- params$algIterationCounter
  params$log$current$Party <- paste0("dp", params$data_partner_id)
  params$log$current$end_time <- get_utc_time()
  params$log$current$computetation_time <- round(as.numeric(difftime(
    params$log$current$end_time,
    params$log$current$start_time, units = "secs")) -
      params$log$current$read_time - params$log$current$write_time, 2)
  params$log$current$files_sent <- paste(files, collapse = ", ")
  params$log$current$bytes_sent <-
    sum(file.size(file.path(params$write_path, files)))
  if (is.na(params$log$current$bytes_sent)) {
    params$log$current$bytes_sent <- 0
  }
  nrows <- nrow(params$log$history)
  if (nrows >= 3) {
    params$log$current$wait_time <-
      round(as.numeric(difftime(
        params$log$current$start_time,
        max(params$log$history$end_time[which(params$log$history$Party ==
                                                params$log$current$Party)]),
        units = "secs")), 2)
  }
  if (params$log$history$Party[nrows] == "") {
    params$log$history <- params$log$current
  } else {
    params$log$history <- rbind(params$log$history, params$log$current)
  }
  return(params)
}


merge_log_raw_kp <- function(params, from) {
  # This function will only be used in the function Pause Continue
  # When party A and party B run simultaneously, but Party A can run first
  # even if Party B starts the whole thing.  We append party B's log
  # to the end of Party A's log.
  if (from == "AC") {
    load(file.path(params$readPathAC, "log.rdata"))
    key1 <- paste0(params$log$history$Step, params$log$history$Party)
    key2 <- paste0(log$Step, log$Party)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$log$history <- rbind(params$log$history, log)
    } else if (length(idx) < length(key2)) {
      params$log$history <- rbind(params$log$history, log[-idx, ])
    }
  } else if (from == "DP1") {
    load(file.path(params$readPathDP[1], "log.rdata"))
    key1 <- paste0(params$log$history$Step, params$log$history$Party)
    key2 <- paste0(log$Step, log$Party)
    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$log$history <- rbind(params$log$history, log)
    } else if (length(idx) < length(key2)) {
      params$log$history <- rbind(params$log$history, log[-idx, ])
    }
  } else {
    for (id in 1:params$num_data_partners) {
      if (id == params$data_partner_id) next
      load(file.path(params$readPathDP[id], "log.rdata"))
      key1 <- paste0(params$log$history$Step, params$log$history$Party)
      key2 <- paste0(log$Step, log$Party)
      idx <- which(key2 %in% key1)
      if (length(idx) == 0) {
        params$log$history <- rbind(params$log$history, log)
      } else if (length(idx) < length(key2)) {
        params$log$history <- rbind(params$log$history, log[-idx, ])
      }
    }
  }
  idx <- order(params$log$history$Step, params$log$history$Party)
  params$log$history <- params$log$history[idx, ]
  return(params)
}


summarize_log_kp <- function(params) {
  write_path <- params$write_path
  log       <- params$log$history

  write_to_log_summary(c1 = "Analysis",
                       c2 = params$analysis,
                       write_path = write_path, append = FALSE)
  if (!is.null(params$n))   write_to_log_summary(c1 = "N",
                                                 c2 = params$n,
                                                 write_path = write_path)

  for (i in 1:params$num_data_partners) {
    if (is.null(params$pi))  {
      p <- 0
    } else {
      p <- params$pi[i] - (i == 1) * (2 + (params$analysis == "cox"))
    }
    write_to_log_summary(c1 = paste0("p", i), c2 = p, write_path = write_path)
  }

  write_to_log_summary(write_path = write_path)

  total.time.0 <- 0
  for (party in 0:params$num_data_partners) {
    party_name <- paste0("dp", party)
    index <- which(log$Party == party_name)
    if (length(index) > 0) {
      start_time <- log$start_time[index[1]]
      end_time   <- log$end_time[index[length(index)]]
      total_time <- round(as.numeric(difftime(end_time,
                                              start_time,
                                              units = "secs")), digits = 2)
      if (party == 0) {
        total.time.0 <- total_time
      }
      reading_time <- sum(log$read_time[index])
      writing_time <- sum(log$write_time[index])
      computing_time <- sum(log$computetation_time[index])
      waiting_time <- sum(log$wait_time[index])
      total_time_hms <-
        convert_secs_to_hms(total_time, time_only = TRUE)
      reading_time_hms <-
        convert_secs_to_hms(reading_time, time_only = TRUE)
      writing_time_hms <-
        convert_secs_to_hms(writing_time, time_only = TRUE)
      computing_time_hms <-
        convert_secs_to_hms(computing_time, time_only = TRUE)
      waiting_time_hms <-
        convert_secs_to_hms(waiting_time, time_only = TRUE)
      Bytes_read <- sum(log$read_Size[index])
      Bytes_written <- sum(log$write_Size[index])
      write_to_log_summary(c1 = paste(party_name, "Start Time"),
                           c2 = start_time, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "End Time"),
                           c2 = end_time, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Run Time"),
                           c2 = total_time,
                           c3 = total_time_hms, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Reading Time"),
                           c2 = reading_time,
                           c3 = reading_time_hms, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Bytes Read"),
                           c2 = Bytes_read, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Writing Time"),
                           c2 = writing_time,
                           c3 = writing_time_hms, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Bytes Written"),
                           c2 = Bytes_written, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Computing Time"),
                           c2 = computing_time,
                           c3 = computing_time_hms, write_path = write_path)
      write_to_log_summary(c1 = paste(party_name, "Total Waiting Time"),
                           c2 = waiting_time,
                           c3 = waiting_time_hms, write_path = write_path)
      write_to_log_summary(write_path = write_path)
    }
  }

  total_transfer_time <- 0
  if (max(log$Step) > 1) {
    for (i in 2:max(log$Step)) {
      idx1 <- which(log$Step == i - 1)
      idx2 <- which(log$Step == i)
      total_transfer_time <- total_transfer_time +
        as.numeric(difftime(min(log$start_time[idx2]),
                            max(log$end_time[idx1]), units = "secs"))
    }
  }
  total_transfer_time <- round(total_transfer_time, 2)
  elapsed_computing_time <- total.time.0 - total_transfer_time

  total_reading_time <- sum(log$read_time)
  total_writing_time <- sum(log$write_time)
  total_computing_time <- sum(log$computetation_time)
  total_reading_time_hms <-
    convert_secs_to_hms(total_reading_time, time_only = TRUE)
  total_writing_time_hms <-
    convert_secs_to_hms(total_writing_time, time_only = TRUE)
  total_transfer_time_hms <-
    convert_secs_to_hms(total_transfer_time, time_only = TRUE)
  total_computing_time_hms <-
    convert_secs_to_hms(total_computing_time, time_only = TRUE)
  elapsed_computing_time_hms <-
    convert_secs_to_hms(elapsed_computing_time, time_only = TRUE)
  total_Bytes_transferred <-
    sum(log$bytes_sent)
  kb_per_second <-
    round(total_Bytes_transferred / (total_transfer_time * 1024), digits = 2)

  write_to_log_summary(c1 = "Total Reading Time",
                       c2 = total_reading_time,
                       c3 = total_reading_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Writing Time",
                       c2 = total_writing_time,
                       c3 = total_writing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Computing Time",
                       c2 = total_computing_time,
                       c3 = total_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Elapsed Computing Time",
                       c2 = elapsed_computing_time,
                       c3 = elapsed_computing_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Transfer Time",
                       c2 = total_transfer_time,
                       c3 = total_transfer_time_hms,
                       write_path = write_path)
  write_to_log_summary(c1 = "Total Bytes Transferred",
                       c2 = total_Bytes_transferred,
                       write_path = write_path)
  write_to_log_summary(c1 = "KB / Sec Transfer Rate",
                       c2 = kb_per_second,
                       write_path = write_path)
}

####################### SHARED TRACKING TABLE FUNCTIONS ########################

WriteTrackingTableRaw <- function(params) {
  trackingTable <- params$trackingTable$history
  save(trackingTable, file = file.path(params$write_path, "tr_tb_updt.rdata"))
  return(params)
}


#' @importFrom utils write.csv
WriteTrackingTableCSV <- function(params) {
  write.csv(params$trackingTable$history,
            file.path(params$write_path, "dl_track_tbl.csv"),
            row.names = FALSE)
  return(params)
}

####################### 2 PARTY TRACKING TABLE FUNCTIONS #######################

initialize_tracking_table_2p <- function(params) {
  trackingTable <- list()
  trackingTable$current <-
    data.frame(DP_CD              = ifelse(params$party_name == "A", 0, 1),
               MSREQID            = params$msreqid,
               RUNID              = "dl",
               # from params$pmnIterationCounter
               ITER_NB            = 0,
               # if params$party_name is "A" then 2 else 1
               STEP_NB            = 0,
               # from log$start_time
               START_DTM          = get_utc_time(),
               # from log$end_time
               END_DTM            = get_utc_time(),
               CURR_STEP_IN       = 0,
               STEP_RETURN_CD     = 0,
               # copy errorMessage.rdata here if exists
               STEP_RETURN_MSG    = "PASS",
               # 1 = converge, 0 = no converge
               REG_CONV_IN        = 0,
               # Success or Failed when decided
               REG_CONV_MSG       = "",
               # 1 at last iteration, so right before quit
               LAST_ITER_IN       = 0,
               LAST_RUNID_IN      = 0,
               UTC_OFFSET_DISPLAY = get_utc_offset(),
               UTC_OFFSET_SEC     = get_utc_offset_seconds(),
               REGR_TYPE_CD       = params$analysis
    )
  trackingTable$history <- NA
  params$trackingTable <- trackingTable
  return(params)
}

#' @importFrom utils write.csv
StoreTrackingTableEntry.2p <- function(params) {
  params$trackingTable$current$ITER_NB <- params$pmn_step_counter
  params$trackingTable$current$START_DTM <- params$log$current$start_time
  params$trackingTable$current$END_DTM <- params$log$current$end_time
  if (file.exists(file.path(params$read_path, "errorMessage.rdata"))) {
    load(file.path(params$read_path, "errorMessage.rdata"))
    params$trackingTable$current$STEP_RETURN_MSG <- message
  } else if (file.exists(file.path(params$write_path, "errorMessage.rdata"))) {
    load(file.path(params$write_path, "errorMessage.rdata"))
    params$trackingTable$current$STEP_RETURN_MSG <- message
  }
  params$trackingTable$current$REG_CONV_IN <- ifelse(params$completed, 1, 0)
  if (params$completed) {
    params$trackingTable$current$REG_CONV_MSG <-
      ifelse(params$converged, "Success", "Failed")
  }
  params$trackingTable$current$LAST_ITER_IN <-
    ifelse(params$lastIteration, 1, 0)

  if (params$party_name == "A") {
    if (is.data.frame(params$trackingTable$history)) {
      params$trackingTable$history <- rbind(params$trackingTable$history,
                                            params$trackingTable$current)
    } else {
      params$trackingTable$history <- params$trackingTable$current
    }
    write.csv(params$trackingTable$history,
              file.path(params$write_path, "dl_track_tbl.csv"),
              row.names = FALSE)
  } else {
    trackingTableEntry <- params$trackingTable$current
    save(trackingTableEntry,
         file = file.path(params$write_path, "tr_tb_updt.rdata"))
  }
  return(params)
}

read_tracking_table_update_2p <- function(params) {
  trackingTableEntry <- NULL
  load(file.path(params$read_path, "tr_tb_updt.rdata"))
  trackingTableEntry$MSREQID <- params$msreqid
  if (is.data.frame(params$trackingTable$history)) {
    params$trackingTable$history <- rbind(params$trackingTable$history,
                                          trackingTableEntry)
  } else {
    params$trackingTable$history <- trackingTableEntry
  }
  return(params)
}

####################### 3 PARTY TRACKING TABLE FUNCTIONS #######################

initialize_tracking_table_3p <- function(params) {
  trackingTable <- list()
  trackingTable$current <-
    data.frame(DP_CD           = ifelse(params$party_name == "T", 0,
                                        ifelse(params$party_name == "A", 1, 2)),
               MSREQID         = params$msreqid,
               RUNID           = "dl",
               ITER_NB         = 0,  # params$pmnIterationCounter
               STEP_NB         = 0,
               START_DTM       = get_utc_time(), # from log$start_time
               END_DTM         = get_utc_time(), # from log$end_time
               CURR_STEP_IN    = 0,
               STEP_RETURN_CD  = 0,
               # copy errorMessage.rdata here if exists
               STEP_RETURN_MSG = "PASS",
               # 1 = converge, 0 = no converge
               REG_CONV_IN     = 0,
               # Success or Failed when decided
               REG_CONV_MSG    = "",
               # 1 at last iteration, so right before quit
               LAST_ITER_IN    = 0,
               LAST_RUNID_IN   = 0,
               UTC_OFFSET_DISPLAY = get_utc_offset(),
               UTC_OFFSET_SEC  = get_utc_offset_seconds(),
               REGR_TYPE_CD    = params$analysis
    )
  trackingTable$history <- NA
  params$trackingTable <- trackingTable
  return(params)
}

StoreTrackingTableEntry_3p <- function(params) {
  params$trackingTable$current$ITER_NB <- params$pmn_step_counter
  params$trackingTable$current$START_DTM <- params$log$current$start_time
  params$trackingTable$current$END_DTM <- params$log$current$end_time
  if (file.exists(file.path(params$write_path, "errorMessage.rdata"))) {
    load(file.path(params$write_path, "errorMessage.rdata"))
    params$trackingTable$current$STEP_RETURN_MSG <- message
  } else {
    msg <- ""
    for (party in c("A", "B", "T")) {
      if (!is.na(params$read_path[[party]]) &&
          file.exists(file.path(params$read_path[[party]],
                                "errorMessage.rdata"))) {
        load(file.path(params$read_path[[party]], "errorMessage.rdata"))
        msg <- paste0(msg, message)
      }
    }
    params$trackingTable$current$STEP_RETURN_MSG <- msg
  }
  params$trackingTable$current$REG_CONV_IN <- ifelse(params$completed, 1, 0)
  if (params$completed) {
    params$trackingTable$current$REG_CONV_MSG <-
      ifelse(params$converged, "Success", "Failed")
  }
  params$trackingTable$current$LAST_ITER_IN <-
    ifelse(params$lastIteration, 1, 0)
  if (params$pmn_step_counter == 0) {
    params$trackingTable$history <- params$trackingTable$current
  } else {
    params$trackingTable$history <- rbind(params$trackingTable$history,
                                          params$trackingTable$current)
  }
  return(params)
}


MergeTrackingTableRAW_3p <- function(params, from) {
  trackingTable <- NULL
  for (party in from) {
    load(file.path(params$read_path[[party]], "tr_tb_updt.rdata"))
    key1 <- paste0(params$trackingTable$history$ITER_NB,
                   params$trackingTable$history$DP_CD)
    key2 <- paste0(trackingTable$ITER_NB,
                   trackingTable$DP_CD)

    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable)
    } else if (length(idx) < length(key2)) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable[-idx, ])
    }
  }
  idx <- order(params$trackingTable$history$START_DTM)
  params$trackingTable$history <- params$trackingTable$history[idx, ]
  params$trackingTable$history$MSREQID <- params$msreqid
  return(params)
}

####################### K PARTY TRACKING TABLE FUNCTIONS #######################

initialize_tracking_table_kp <- function(params) {
  trackingTable <- list()
  trackingTable$current <- data.frame(
    DP_CD              = params$data_partner_id,
    MSREQID            = params$msreqid,
    RUNID              = "dl",
    # params$pmnIterationCounter
    ITER_NB            = 0,
    STEP_NB            = 0,
    # from log$start_time
    START_DTM          = get_utc_time(),
    # from log$end_time
    END_DTM            = get_utc_time(),
    CURR_STEP_IN       = 0,
    STEP_RETURN_CD     = 0,
    # copy errorMessage.rdata here if exists
    STEP_RETURN_MSG    = "PASS",
    # 1 = converge, 0 = no converge
    REG_CONV_IN        = 0,
    # Success or Failed when decided
    REG_CONV_MSG       = "",
    # 1 at last iteration, right before quit
    LAST_ITER_IN       = 0,
    LAST_RUNID_IN      = 0,
    UTC_OFFSET_DISPLAY = get_utc_offset(),
    UTC_OFFSET_SEC     = get_utc_offset_seconds(),
    REGR_TYPE_CD       = params$analysis
  )
  trackingTable$history <- NA
  params$trackingTable <- trackingTable
  return(params)
}

StoreTrackingTableEntry.kp <- function(params) {
  params$trackingTable$current$ITER_NB <- params$pmn_step_counter
  params$trackingTable$current$START_DTM <- params$log$current$start_time
  params$trackingTable$current$END_DTM <- params$log$current$end_time
  if (file.exists(file.path(params$write_path, "errorMessage.rdata"))) {
    load(file.path(params$write_path, "errorMessage.rdata"))
    params$trackingTable$current$STEP_RETURN_MSG <- message
  } else {
    msg <- ""
    for (id in 1:params$num_data_partners) {
      if (!is.na(params$readPathDP[id]) &&
          file.exists(file.path(params$readPathDP[id], "errorMessage.rdata"))) {
        load(file.path(params$readPathDP[id], "errorMessage.rdata"))
        msg <- paste0(msg, message)
      }
    }
    if (!is.na(params$readPathAC) &&
        file.exists(file.path(params$readPathAC, "errorMessage.rdata"))) {
      load(file.path(params$readPathAC, "errorMessage.rdata"))
      msg <- paste0(msg, message)
    }
    params$trackingTable$current$STEP_RETURN_MSG <- msg
  }
  params$trackingTable$current$REG_CONV_IN <- ifelse(params$completed, 1, 0)
  if (params$completed) {
    params$trackingTable$current$REG_CONV_MSG <-
      ifelse(params$converged, "Success", "Failed")
  }
  params$trackingTable$current$LAST_ITER_IN <-
    ifelse(params$lastIteration, 1, 0)
  if (params$pmn_step_counter == 0) {
    params$trackingTable$history <- params$trackingTable$current
  } else {
    params$trackingTable$history <- rbind(params$trackingTable$history,
                                          params$trackingTable$current)
  }

  return(params)
}


MergeTrackingTableRAW.kp <- function(params, from) {
  trackingTable <- NULL
  if (from == "AC") {
    load(file.path(params$readPathAC, "tr_tb_updt.rdata"))
    key1 <- paste0(params$trackingTable$history$ITER_NB,
                   params$trackingTable$history$DP_CD)
    key2 <- paste0(trackingTable$ITER_NB,
                   trackingTable$DP_CD)

    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable)
    } else if (length(idx) < length(key2)) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable[-idx, ])
    }
  } else if (from == "DP1") {
    load(file.path(params$readPathDP[1], "tr_tb_updt.rdata"))
    key1 <- paste0(params$trackingTable$history$ITER_NB,
                   params$trackingTable$history$DP_CD)
    key2 <- paste0(trackingTable$ITER_NB,
                   trackingTable$DP_CD)

    idx <- which(key2 %in% key1)
    if (length(idx) == 0) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable)
    } else if (length(idx) < length(key2)) {
      params$trackingTable$history <-
        rbind(params$trackingTable$history, trackingTable[-idx, ])
    }
  } else {
    for (id in 1:params$num_data_partners) {
      if (id == params$data_partner_id) next
      load(file.path(params$readPathDP[id], "tr_tb_updt.rdata"))
      key1 <- paste0(params$trackingTable$history$ITER_NB,
                     params$trackingTable$history$DP_CD)
      key2 <- paste0(trackingTable$ITER_NB,
                     trackingTable$DP_CD)

      idx <- which(key2 %in% key1)
      if (length(idx) == 0) {
        params$trackingTable$history <-
          rbind(params$trackingTable$history, trackingTable)
      } else if (length(idx) < length(key2)) {
        params$trackingTable$history <-
          rbind(params$trackingTable$history, trackingTable[-idx, ])
      }
    }
  }
  idx <- order(params$trackingTable$history$START_DTM)
  params$trackingTable$history <- params$trackingTable$history[idx, ]
  params$trackingTable$history$MSREQID <- params$msreqid
  return(params)
}

########################### VALID FORMULA FUNCTIONS ############################

validFormula <- function(expression) {
  # This function takes an expresion and checks that it is of the form var1 ~
  # var2 + var3 + ... varN It does not check for constants.  Constants are
  # ignored and treated if the are not there. Dupliate variables are ignored.
  # That is, as in lm(), formulas of the form var1 ~ var2 + var2 are equivalent
  # to var1 ~ var2

  # Check to make sure this is a valid expression
  if (tryCatch({
    is.expression(expression)
    FALSE
  },
  error = function(err) {
    TRUE
  })) {
    return(FALSE)
  }
  vars <- all.vars(expression)
  names <- all.names(expression)
  #Check to see if expresion only contains variables, ~, and +.  no other
  #symbols allowed.
  res1 <- all(names %in% c("~", "+", vars))
  #Check to see if expresion contains exactly one ~
  res2 <- (sum(names %in% "~") == 1)
  #Check to see if expression is of the form "variable ~ ....."
  res3 <- (names[1] == "~") & (names[2] %in% vars)
  #Check to see if the LHS variable does not occur on the RHS
  res4 <- !(names[2] %in% names[3:length(names)])
  #check to see if the LHS variable is not .
  res5 <- vars[1] != "."
  return(res1 & res2 & res3 & res4 & res5)
}

validFormula2 <- function(expression) {
  # This function takes and expression and checks that it is of the form ~ var1
  # + var2 + ... + varN Duplicate variables are ignored.  That is, ~ var1 + var1
  # is equivalent to ~ var1
  if (tryCatch({
    is.expression(expression)
    FALSE
  },
  error = function(err) {
    TRUE
  })) {
    return(FALSE)
  }
  vars <- all.vars(expression)
  names <- all.names(expression)
  # Check to see if expresion only contains variables, ~, and +.  no other
  # symbols allowed.
  res1 <- all(names %in% c("~", "+", vars))
  # Check to see if expresion contains exactly one ~
  res2 <- (sum(names %in% "~") == 1)
  # Check to see if expression has no LHS (should not)
  res3 <- length(expression) == 2
  return(res1 && res2 && res3)
}

###################### SHARED SUMMARY AND PRINT FUNCTIONS ######################

#' @export
print.vdralinear <- function(x, ...) {
  if (x$failed) {
    warning("Distributed linear regression failed.  No results to print.")
    return(invisible(NULL))
  }
  cat("Coefficients:\n")
  print(x$coefficients, digits = 4)
  return(invisible(NULL))
}

#' @title Summary Method for Vertical Distributed Linear Regression Models
#' @aliases summary.vdralinear summary.vdralinear.object
#'   print.summary.vdralinear
#' @description Produces a summary of a fitted vdra linear regression model.
#' @param object a \code{vdralinear} object.
#' @param ... further arguments passed to or from other methods.
#' @return Returns an object of class \code{summary.vdralinear}. Objects of this
#'   class have a method for the function \code{print}.  The following
#'   components must be included in \code{summary.vdralinear} object. \describe{
#'
#'   \item{failed}{logical value.  If \code{FALSE}, then there was an error
#'   processing the data.  if \code{TRUE}, there were no errors.}
#'
#'   \item{party}{a vector which indicates the party from which each covariate
#'   came.}
#'
#'   \item{coefficients}{the vector of coefficients.  If the model is
#'   over-determined, there will be missing values in the vector corresponding
#'   to the redundant columns model matrix.}
#'
#'   \item{secoef}{the vector of the standard error of the coefficients.}
#'
#'   \item{tvals}{the t-values of the coefficients.}
#'
#'   \item{pvals}{the p-values of the coefficients.}
#'
#'   \item{rstderr}{residual standard error.}
#'
#'   \item{rsquare}{r squared.}
#'
#'   \item{adjrsquare}{adjusted r squared.}
#'
#'   \item{Fstat}{the F-statistic for the linear regression.}
#'
#'   \item{df1}{the numerator degrees of freedom for the F-statistic.}
#'
#'   \item{df2}{the denominator degrees of freedom for the F-statistic.}
#'
#'   \item{Fpval}{the p-value of the F-statistic for the linear regression.}
#'
#'   }
#'
#' @seealso \code{\link{vdralinear}}
#' @examples
#' summary(vdra_fit_linear_A)
#' @export
summary.vdralinear <- function(object, ...) {
  temp <- list()
  class(temp)         <- "summary.vdralinear"
  temp$failed         <- object$failed
  if (object$failed) {
    return(temp)
  }
  temp$party          <- object$party
  temp$coefficients   <- object$coefficients
  temp$secoef         <- object$secoef
  temp$tvals          <- object$tvals
  temp$pvals          <- object$pvals
  temp$rstderr        <- object$rstderr
  temp$df2            <- object$df2
  temp$rsquare        <- object$rsquare
  temp$adjrsquare     <- object$adjrsquare
  temp$Fstat         <- object$Fstat
  temp$df1            <- object$df1
  temp$Fpval         <- object$Fpval
  return(temp)
}

#' @export
print.summary.vdralinear <- function(x, lion = FALSE, ...) {
  if (x$failed) {
    warning("Distributed linear regression failed.  No results to print.")
    return(invisible(NULL))
  }

  x$stars          <- sapply(x$pvals, function(x) {
    if (is.na(x))       ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else                " "
  })

  temp <- data.frame(formatStrings(names(x$party)),
                     formatStrings(x$party, minWidth = 5, justify = "centre"),
                     formatStatList(x$coefficients),
                     formatStatList(x$secoef),
                     formatStatList(x$tvals),
                     format.pval(x$pvals),
                     formatStrings(x$stars))
  colnames(temp) <- c("", "Party", "Estimate", "Std. Error",
                      "t value", "Pr(>|t|)", "")
  if (lion) {
    temp <- cbind(temp, GetLion(length(x$party)))
    colnames(temp)[8] <- ""
  }
  print(temp, row.names = FALSE, right = TRUE)
  cat("---",
      "\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("Residual standard error: ", formatStat(x$rstderr),
      "on", x$df2, "degrees of freedom\n")
  cat("Multiple R-squared: ", formatStat(x$rsquare),
      ", Adjusted R-squared: ", formatStat(x$adjrsquare), "\n")
  cat("F-statistic:", formatStat(x$Fstat), "on", x$df1,
      "and", x$df2, "DF, p-value:", format.pval(x$Fpval), "\n\n")
}

#' @export
print.vdralogistic <- function(x, ...) {
  if (x$failed) {
    warning("Distributed logistic regression failed.  No results to print.")
    return(invisible(NULL))
  }
  if (!x$converged) {
    warning(
      paste("Warning: Distributed logistic regression did not converge in",
            x$iter, "iterations. Reported statistics are approximate."))
  }

  cat("Coefficients:\n")
  print(x$coefficients, digits = 4)
  cat("\n")
  cat("Degrees of Freedom:", x$nulldev_df,
      "Total (i.e. Null); ", x$resdev_df, "Residual\n")
  cat("Null Deviance:    ", formatStat(x$nulldev), "\n")
  cat("Residual Deviance:", formatStat(x$resdev), "   AIC:", formatStat(x$aic))
  return(invisible(NULL))
}

#' @title Summary Method for Vertical Distributed Logistic Regression Models
#' @aliases  summary.vdralogistic summary.vdralogistic.object
#'   print.summary.vdralogistic
#' @description Produces a summary of a fitted vdra logistic regression model.
#' @param object a \code{vdralogistic} object.
#' @param ... futher argumetns passed to or from other methods.
#' @return Returns an object of class \code{summary.vdralogistic}. Objects of
#'   this class have a method for the function \code{print}.  The following
#'   components must be included in \code{summary.vdralogistic} object.
#'   \describe{
#'
#'   \item{failed}{ogical value.  If \code{FALSE}, then there was an error
#'   processing the data.  if \code{TRUE}, there were no errors.}
#'
#'   \item{converged}{ogical value.  If \code{TRUE}, the regression converged.
#'   If \code{FALSE}, it did not.}
#'
#'   \item{party}{ vector which indicates the party from which each covariate
#'   came.}
#'
#'   \item{coefficients}{he vector of coefficients.  If the model is
#'   over-determined, there will be missing values in the vector corresponding
#'   to the redudant columns model matrix.}
#'
#'   \item{secoef}{he vector of the standard error of the coefficients.}
#'
#'   \item{tvals}{he t-values of the coefficietns.}
#'
#'   \item{pvals}{he p-values of the coefficients.}
#'
#'   \item{nulldev}{he null deviance of the fit.}
#'
#'   \item{nulldev_df}{he degrees of freedom for the null deviance.}
#'
#'   \item{resdev}{he residual deviance of the fit.}
#'
#'   \item{resdev_df}{he degrees of freedome for the residual deviance.}
#'
#'   \item{aic the}{IC of the fit.}
#'
#'   \item{bic the}{IC of the fit.}
#'
#'   \item{iter }{ number of iterations of the cox algorithm before
#'   convergence.}
#'
#'   }
#' @seealso \code{\link{vdralogistic}}
#' @examples
#' summary(vdra_fit_logistic_A)
#' @export
summary.vdralogistic <- function(object, ...) {
  temp <- list()
  class(temp)         <- "summary.vdralogistic"
  temp$failed         <- object$failed
  temp$converged      <- object$converged
  if (object$failed) {
    return(temp)
  }
  temp$party          <- object$party
  temp$coefficients   <- object$coefficients
  temp$secoef         <- object$secoef
  temp$tvals          <- object$tvals
  temp$pvals          <- object$pvals
  temp$nulldev        <- object$nulldev
  temp$nulldev_df     <- object$nulldev_df
  temp$resdev         <- object$resdev
  temp$resdev_df      <- object$resdev_df
  temp$aic            <- object$aic
  temp$bic            <- object$bic
  temp$iter           <- object$iter
  return(temp)
}

#' @export
print.summary.vdralogistic <- function(x, lion = FALSE, ...) {
  arguments <- list(...)

  if (!is.na(arguments$lion) &&
      is.logical(arguments$lion)) lion <- arguments$lion

  if (x$failed) {
    warning("Distributed logistic regression failed.  No results to print.")
    return(invisible(NULL))
  }
  if (!x$converged) {
    warning(
      paste("Warning: Distributed logistic regression did not converge in",
            x$iter, "iterations. Reported statistics are approximate."))
  }
  x$stars <- sapply(x$pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  })

  temp <- data.frame(formatStrings(names(x$party)),
                     formatStrings(x$party, minWidth = 5, justify = "centre"),
                     formatStatList(x$coefficients),
                     formatStatList(x$secoef),
                     formatStatList(x$tvals),
                     format.pval(x$pvals),
                     formatStrings(x$stars))
  colnames(temp) <- c("", "Party", "Estimate",
                      "Std. Error", "t value", "Pr(>|t|)", "")
  if (lion) {
    temp <- cbind(temp, GetLion(length(x$party)))
    colnames(temp)[8] <- ""
  }
  print(temp, row.names = FALSE, right = TRUE)
  cat("---",
      "\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("(Dispersion parameter for binomial family taken to be 1)\n\n")
  cat("    Null Deviance:", formatStat(x$nulldev), " on ",
      x$nulldev_df, " degrees of freedom\n")
  cat("Residual deviance:", formatStat(x$resdev), " on ",
      x$resdev_df, " degrees of freedom\n")
  cat("AIC:", formatStat(x$aic), "\n")
  cat("BIC:", formatStat(x$bic), "\n\n")
  cat("Number of Newton-Raphson iterations:", x$iter, "\n\n")
  return(invisible(NULL))
}


#' @importFrom stats printCoefmat
#' @export
print.vdracox <- function(x, ...) {
  if (x$failed) {
    warning("Distributed Cox regression failed.  No results to print.")
    return(invisible(NULL))
  }
  if (!x$converged) {
    warning(paste("Warning: Distributed Cox regression did not converge in",
                  x$iter, "iterations. Reported statistics are approximate."))
  }

  coeftab <- data.frame(x$coefficients, x$expcoef, x$secoef, x$zvals, x$pvals)
  colnames(coeftab) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  printCoefmat(coeftab, P.values = TRUE,
               has.Pvalue = TRUE, signif.stars = FALSE)
  cat("\n")
  cat(paste0("Likelihood ratio test=", formatStat(x$lrt[1])), "on",
      x$df, paste0("df, p=", format.pval(x$lrt[2])), "\n")
  cat("n=", paste0(x$n, ","), "number of events=", x$nevent, "\n\n")
  return(invisible(NULL))
}

#' @title Summary Method for Vertical Distributed COX Models
#' @description Produces a summary of a fitted vdra Cox model.
#' @aliases summary.vdracox summary.vdracox.object print.summary.vdracox
#' @param object a \code{vdracox} object.
#' @param ... further arguments passed to or from other methods.
#' @return  Returns an object of class \code{summary.vdracox}. Objects of this
#'  class have a method for the function \code{print}.  The following components
#'  must be included in \code{summary.vdracox} object.
#' \describe{
#'  \item{failed}{logical value.  If \code{FALSE}, then there was an error
#'     processing the data.  if \code{TRUE}, there were no errors.}
#'  \item{converged}{logical value.  If \code{TRUE}, the regression converged.
#'     If \code{FALSE}, it did not.}
#'  \item{party}{a vector which indicates the party from which each covariate
#'     came.}
#'  \item{coefficients}{the vector of coefficients.  If the model is
#'     over-determined, there will be missing values in the vector corresponding
#'     to the redudant columns model matrix.}
#'  \item{expcoef}{a vector which represents exp(coefficients).}
#'  \item{secoef}{the vector of the standard error of the coefficients.}
#'  \item{zvals}{the z-values of the coefficients.}
#'  \item{pvals}{the p-values of the coefficients.}
#'  \item{expncoef}{a vector which represents exp(-coefficients).}
#'  \item{lower95}{a vector of the lower bounds of the 95\% confidence interval
#'     for exp(coefficients).}
#'  \item{upper95}{a vector of the upper bounds of the 95\% confidence interval
#'     for exp(coefficients).}
#'  \item{n}{the number of observations in the data.}
#'  \item{nevent}{the number of events used in the fit.}
#'  \item{concordance}{a vector containing the number of events which are
#'     concordant, discordant, tied.risk, tied.time.  Also contains the
#'     concordance statistic and its standard error.  Calculated using the
#'     \code{survival} package, if installed.  If not installed, all values are
#'     \code{NA}.}
#'  \item{rsquare}{a vector containing an r-square value for the fit and its
#'     p-value.}
#'  \item{lrt}{a vector containing the likelihood ratio test statistic and its
#'     p-value.}
#'  \item{df}{the degrees of freedom.}
#'  \item{wald.test}{a vector containing the Wald test statistic and its
#'     p-value.}
#'  \item{score}{a vector containing the score test statistic and its p-value.}
#'  \item{iter}{the number of iterations of the cox algorithm before
#'     convergence.}
#' }
#' @seealso \code{\link{vdracox}}
#' @examples
#' summary(vdra_fit_cox_A)
#' @export
summary.vdracox <- function(object, ...) {
  temp <- list()
  class(temp)         <- "summary.vdracox"
  temp$failed         <- object$failed
  temp$converged      <- object$converged
  if (object$failed) {
    return(temp)
  }
  temp$party          <- object$party
  temp$coefficients   <- object$coefficients
  temp$expcoef        <- object$expcoef
  temp$secoef         <- object$secoef
  temp$zval           <- object$zval
  temp$pvals          <- object$pvals
  temp$expncoef       <- object$expncoef
  temp$lower          <- object$lower
  temp$upper          <- object$upper
  temp$n              <- object$n
  temp$nevent         <- object$nevent
  temp$concordance    <- object$concordance
  temp$rsquare        <- object$rsquare
  temp$lrt            <- object$lrt
  temp$df             <- object$df
  temp$wald.test      <- object$wald.test
  temp$score          <- object$score
  temp$iter           <- object$iter
  return(temp)
}

#' @export
print.summary.vdracox <- function(x, lion = FALSE, ...) {
  arguments <- list(...)

  if (!is.na(arguments$lion) &&
      is.logical(arguments$lion)) lion <- arguments$lion

  if (x$failed) {
    warning("Distributed Cox regression failed.  No results to print.")
    return(invisible(NULL))
  }
  if (!x$converged) {
    warning(paste("Warning: Distributed Cox regression did not converge in",
                  x$iter, "iterations. Reported statistics are approximate."))
  }

  x$stars <- sapply(x$pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  })

  temp1 <- data.frame(formatStrings(names(x$party)),
                      formatStrings(x$party, minWidth = 5, justify = "centre"),
                      formatStatList(x$coefficients),
                      formatStatList(x$expcoef),
                      formatStatList(x$secoef),
                      formatStatList(x$zval),
                      format.pval(x$pvals),
                      formatStrings(x$stars))
  colnames(temp1) <- c("", "party", "   coef", "exp(coef)",
                       "se(coef)", "   z", "Pr(>|z|)", "")
  temp2 <- data.frame(formatStrings(names(x$party)),
                      formatStrings(x$party, minWidth = 5, justify = "centre"),
                      formatStatList(x$expcoef),
                      formatStatList(x$expncoef),
                      formatStatList(x$lower),
                      formatStatList(x$upper))
  colnames(temp2) <- c("", "party", "exp(coef)", "exp(-coef)",
                       "lower .95", "upper .95")
  cat("  n=", paste0(x$n, ","), "number of events=", x$nevent, "\n\n")
  print(temp1, row.names = FALSE, right = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  print(temp2, row.names = FALSE, right = TRUE)
  cat("\n")
  if (!is.na(x$concordance[5])) {
    cat("Concordance=", formatStat(x$concordance[5]),
        "(se =", formatStat(x$concordance[6]), ")\n")
  }
  cat("Likelihood ratio test=", formatStat(x$lrt[1]),
      "on",                     x$df,
      "df,", "p=",              format.pval(x$lrt[2]), "\n")
  cat("Wald test            =", formatStat(x$wald.test[1]),
      "on",                     x$df,
      "df, p=",                 format.pval(x$wald.test[2]), "\n")
  cat("Score test           =", formatStat(x$score[1]),
      "on",                     x$df,
      "df, p=",                 format.pval(x$score[2]), "\n\n")
  cat("Number of Newton-Raphson iterations:", x$iter, "\n\n")
}


############################ LINEAR DIFFERENT MODEL ############################

#' @title Fitting Different Linear Models
#' @description Models are specified symbolically.  A typical model is of the
#'   form \code{response ~ term_1 + term_2 + ... + term_k} where \code{response}
#'   and \code{term_i} are variables names used in the orgional linear model
#'   which created the object \code{x}.  The \code{response} can be the orginal
#'   respose or any of the other covariates.  Interactions are not allowed.  Not
#'   all variables in the original model have to be used.
#' @param formula an object of class \code{"\link{formula}"}: a symbolic
#'   description of the model to be fitted. The model must be additive with no
#'   interactions.
#' @param x an object of class \code{\link{vdralinear}}.
#' @return Returns an object of class \code{\link{vdralinear}}.
#' @seealso \code{\link{analysis_center_2_party}},
#'   \code{\link{AnalysisCenter.3Party}}, \code{\link{AnalysisCenter.KParty}}
#' @examples
#'  fit <- differentModel(Change_BMI ~ Exposure + Age +
#'                        NumRx, vdra_fit_linear_A)
#'  summary(fit)
#'
#'  fit <- differentModel(Age ~ Change_BMI + Exposure +
#'                        NumRx, vdra_fit_linear_A)
#'  summary(fit)
#' @importFrom  stats pf pt
#' @export
differentModel <- function(formula = NULL, x = NULL) {
  if (!inherits(x, "vdralinear")) {
    warning("This function can only be on objects of class vdralinear. ",
            "Returning original model.")
    return(invisible(x))
  }
  if (x$failed) {
    warning("Distributed linear regression failed. ",
            "Cannot compute a different model.")
    return(invisible(x))
  }
  if (max(table(names(x$party))) > 1) {
    warning("Duplicate variable names exist. ",
            "All variable names must be unique. Returning original model.")
    return(invisible(x))
  }
  if (!validFormula(formula)) {
    warning("Invalid formula, returning original model.")
    return(x)
  }

  valid_names    <- c(colnames(x$xty), colnames(x$xtx)[-1])
  responseName   <- all.vars(formula)[1]
  covariateNames <- all.vars(formula)[-1]
  variableNames  <- all.vars(formula)
  variableNames  <- variableNames[which(variableNames != ".")]

  if (!all(variableNames %in% valid_names)) {
    vars <- variableNames[which(variableNames %in% valid_names == FALSE)]
    if (length(vars) == 1) {
      warning("Variable", vars, "not found. Returning original model.")
    } else {
      temp <- c(paste0(vars[-length(vars)], ","), vars[length(vars)])
      warning(paste("Variables", temp, "not found. Returning original model."))
    }
    return(invisible(x))
  }

  if ("." %in% covariateNames) {
    covariateNames <- valid_names[which(valid_names != responseName)]
  }

  xytxy <- rbind(cbind(x$yty, t(x$xty)), cbind(x$xty, x$xtx))
  scramble <- c(2, 1, 3:ncol(xytxy))
  xytxy[scramble, scramble] <- xytxy

  # Put (intercept) first
  all_names <- c(colnames(x$xty), colnames(x$xtx))[scramble]
  colnames(xytxy) <- all_names
  rownames(xytxy) <- all_names

  response_index <- match(responseName, all_names)
  covariate_index <- c(1, match(covariateNames, all_names))

  xtx    <- xytxy[covariate_index, covariate_index]
  xty    <- matrix(xytxy[covariate_index, response_index], ncol = 1)
  yty    <- xytxy[response_index, response_index]
  means  <- c(x$means_y, x$means)[scramble][covariate_index]
  means_y <- c(x$means_y, x$means)[scramble][response_index]

  nrow <- nrow(xtx)
  indicies <- c(1)
  for (i in 2:nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }

  xtx_old   <- xtx
  xty_old   <- xty
  xtx       <- xtx[indicies, indicies]
  xty       <- matrix(xty[indicies, 1], ncol = 1)
  means_old <- means
  means     <- means[indicies]
  p         <- length(indicies)
  n         <- x$n

  invxtx <- solve(xtx)
  betas  <- drop(invxtx %*% xty)

  num_covariates <- p - 1

  sse     <- max(drop(yty - 2 * t(xty) %*% betas +
                        (t(betas) %*% xtx) %*% betas), 0)
  rstderr <- drop(sqrt(sse / (n - num_covariates - 1)))
  sst     <- drop(yty - means_y^2 * n)
  ssr     <- sst - sse
  df1     <- num_covariates
  df2     <- n - num_covariates - 1
  if (sse == 0) {
    Fstat <- Inf
  } else {
    Fstat <- (ssr / df1) / (sse / df2)
  }
  Fpval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  if (sse == 0) {
    r_sq <- 1
  } else {
    r_sq <- drop(1 - sse / sst)
  }
  adj_r_sq <- drop(1 - (n - 1) / (n - num_covariates - 1) * (1 - r_sq))
  if (rstderr == 0) {
    tvals <- rep(Inf, num_covariates + 1)
  } else {
    tvals   <- betas / (rstderr * sqrt(diag(invxtx)))
  }
  secoef  <- tvals^-1 * betas
  pvals   <- 2 * pt(abs(tvals), n - num_covariates - 1, lower.tail = FALSE)
  stars   <- matrix(sapply(pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  }))

  y <- list()
  class(y) <- "vdralinear"
  y$failed    <- x$failed
  y$converged <- x$converged

  y$party <- c(x$responseParty, x$party)[scramble][covariate_index]
  y$responseParty <- c(x$responseParty, x$party)[scramble][response_index]
  p1 <- length(covariate_index)
  y$coefficients           <- rep(NA, p1)
  y$tvals                  <- rep(NA, p1)
  y$secoef                 <- rep(NA, p1)
  y$pvals                  <- rep(NA, p1)

  y$sse                    <- sse
  y$coefficients[indicies] <- betas
  y$tvals[indicies]        <- tvals
  y$secoef[indicies]       <- secoef
  y$pvals[indicies]        <- pvals
  y$rstderr                <- rstderr
  y$rsquare                <- r_sq
  y$adjrsquare             <- adj_r_sq
  y$Fstat                  <- Fstat
  y$Fpval <- Fpval
  y$df1                    <- df1
  y$df2                    <- df2
  y$n                      <- x$n
  y$xtx                    <- xtx_old
  y$xty                    <- xty_old
  y$yty                    <- yty
  y$means_y                 <- means_y
  y$means                  <- means_old

  names_old                <- all_names[covariate_index]
  names(y$party)           <- names_old
  names(y$coefficients)    <- names_old
  names(y$secoef)          <- names_old
  names(y$tvals)           <- names_old
  names(y$pvals)           <- names_old
  colnames(y$xtx)          <- names_old
  rownames(y$xtx)          <- names_old
  colnames(y$xty)          <- responseName
  rownames(y$xty)          <- names_old
  return(invisible(y))
}

###################### LOGISTIC roc AND HOSLEM FUNCTIONS #######################


#' @importFrom  stats pchisq quantile
HoslemInternal <- function(x, data = NULL, nGroups = 10) {
  #            y:  response (vector, length n)
  #  final_fitted:  final_fitted from getFinalCoefA(...)  (vector, length n)
  #            p:  number of covariates p_a + p_b
  #      nGroups:  number of groups, specified by user
  #                or chosen automatically if unspecified.
  #                Common to choose ngroups = 10, as long as nGroups > p + 1.
  #
  # Returns vector c(chisq, df, pval)

  n <- x$n

  if (nGroups <= 0) {
    nGroups <- 10
  }

  if (nGroups > n) {
    nGroups <- n
  }

  if (is.null(data)) {
    Y <- x$Y
  } else {
    Y <- data$Y
  }
  pi_ <- exp(x$final_fitted) / (1 + exp(x$final_fitted))
  uq <- unique(quantile(pi_, probs = seq(0, 1, 1 / nGroups)))
  group_ <- cut(pi_, breaks = uq, include.lowest = TRUE)
  dd <- data.frame(y = Y[order(pi_)], pi_ = sort(pi_),
                   group = group_[order(pi_)])

  e1 <- by(dd, dd$group, function(x) sum(x$pi_))
  o1 <- by(dd, dd$group, function(x) sum(x$y))
  gn <- table(dd$group)
  e0 <- gn - e1
  o0 <- gn - o1

  testStat <- 0
  for (i in seq_along(e1)) {
    if (o0[i] == e0[i]) {
      temp1 <- 0
    } else {
      temp1 <- (o0[i] - e0[i])^2 / e0[i]
    }
    if (o1[i] == e1[i]) {
      temp2 <- 0
    } else {
      temp2 <- (o1[i] - e1[i])^2 / e1[i]
    }
    testStat <- testStat + temp1 + temp2
  }

  df <- nGroups - 2
  rtrn <- c(testStat, df, 1 - pchisq(testStat, df))
  names(rtrn) <- c("Chi-sq", "DF", "p-value")

  return(rtrn)
}

#' @export
print.hoslemdistributed <- function(x, ...) {
  cat("Hosmer and Lemeshow goodness of fit (GOF) test\n",
      "       Chi-squared:", x$hoslem[1], "with DF",
      paste0(x$hoslem[2], ","), " p-value:", x$hoslem[3], "\n")
}


#' @title Hosmer-Lemeshow Test for Vertical Distributed Logistic Regression
#' @description Run the Hosmer-Lemeshow test for an object created by 2-party,
#'   3-party, or K-party vdra logistic regression.  Only the party that holds
#'   the response may invoke this function.
#' @aliases HoslemTest hoslemdistributed.object print.hoslemdistributed
#' @param x an object of type \code{\link{vdralogistic}}.
#' @param nGroups the number of groups that the data will be sperated into.
#' @return  Returns an object of class \code{hoslemdistributed}. Objects of this
#'   class have a method for the function \code{print}.   The following
#'   component must be included in a \code{hoslemdistributed} object.
#'
#' \describe{
#'
#' \item{hoslem}{a vector containing three numeric quantities: the chi-square
#'   value of the test, the degrees of freedom of the test, and p-value of the
#'   test, in that order.}
#'
#' }
#' @examples
#'  HoslemTest(vdra_fit_logistic_A)
#'
#'  HoslemTest(vdra_fit_logistic_A, 20)
#' @export
HoslemTest <- function(x = NULL, nGroups = 10) {
  if (!inherits(x, "vdralogistic")) {
    warning("Cannot perform test on non vdralogistic object.")
    return(invisible(NULL))
  }
  if (!(x$converged)) {
    warning("Process did not converge. ",
            "Cannot perform Hosmer and Lemeshow goodness of fit test.")
    return(invisible(NULL))
  }
  if (is.null(x$Y) || is.null(x$final_fitted)) {
    warning("HoslemTest can only be invoked by ",
            "the party which holds the response.")
    return(invisible(NULL))
  } else if (is.numeric(nGroups)) {
    temp <- list()
    class(temp) <- "hoslemdistributed"
    temp$hoslem <- HoslemInternal(x, nGroups = nGroups)
    return(temp)
  }
}


roc_internal <- function(x, data = NULL, bins = 500) {
  #             y:  response vector (numeric, not factor, length n)
  #   final_fitted:  final_fitted from getFinalCoefA(...)  (vector, length n)
  #    thresholds:  how smooth the curve should be
  #
  #  Returns myRocObject (object$auc to get AUC)
  #
  #  Object size is roughly equal to a matrix with
  # (thresholds) rows and 2 columns

  if (is.null(data)) {
    Y <- x$Y
  } else {
    Y <- data$Y
  }

  if (bins < 2) bins <- 2

  positive <- sum(Y)
  negative <- length(Y) - positive
  pi_ <- exp(x$final_fitted) / (1 + exp(x$final_fitted))
  threshold <- seq(0, 1, length.out = bins)
  rtrn <- matrix(NA, bins, 2)

  oldX <- 1
  oldY <- 1
  AUC <- 0

  for (i in 1:bins) {
    newX <- 1 - sum(Y == 0 & pi_ < threshold[i]) / negative
    newY <- sum(Y & pi_ >= threshold[i]) / positive
    rtrn[i, 1] <- newX
    rtrn[i, 2] <- newY
    AUC <- AUC + oldY * (oldX - newX)
    oldX <- newX
    oldY <- newY
  }

  temp <- list()
  temp$roc <- rtrn
  temp$auc <- AUC
  return(temp)
}

#' @importFrom graphics axis lines text plot
#' @export
print.rocdistributed <- function(x, ...) {
  rtrn <- x$roc
  plot(rtrn[, 1], rtrn[, 2], xaxt = "n", yaxt = "n",
       xlim = c(-0.2, 1.2), ylim = c(0, 1),
       type = "s", ylab = "Sensitivity", xlab = "1 - Specificity", col = "blue",
       main = "roc Curve")
  axis(side = 1, at = seq(0, 1, 0.2))
  axis(side = 2, at = seq(0, 1, 0.2))
  lines(x = c(0, 1), y = c(0, 1), col = "limegreen")
  text(0.8, 0.05, paste("Area under the curve:",
                        format(x$auc, digits = 4)))
}


#' @aliases RocTest rocdistributed.object print.rocdistributed
#' @title Create the roc for Vertical Distributed Logistic Regression
#' @description Generate the receiver operator curve on an object created by
#'   2-party, 3-party, or K-party vdra logistic regression.  Only the party that
#'   holds the response may invoke this function.
#' @param x an object of type \code{\link{vdralogistic}}.
#' @param bins the number of bins the data will be separated into.
#' @return Returns an object of class \code{rocdistributed}. Objects of this
#'   class have a method for the function \code{print}. The following components
#'   must be included in a \code{rocdistributed} object.
#' \describe{
#'  \item{roc}{a two column matrix containing the coordinates of 1 - specificity
#'  and sensitivity.}
#'  \item{auc}{numeric value which is area under the curve.}
#' }
#' @examples
#' RocTest(vdra_fit_logistic_A)
#'
#' RocTest(vdra_fit_logistic_A, 40)
#' @export
RocTest <- function(x = NULL, bins = 10) {
  if (!inherits(x, "vdralogistic")) {
    warning("Cannot create roc on non vdralogistic object.")
    return(invisible(NULL))
  }
  if (!x$converged) {
    warning("Process did not converge.  Cannot generate roc.")
    return(invisible(NULL))
  }
  if (is.null(x$Y) || is.null(x$final_fitted)) {
    warning("RocTest can only be invoked by ",
            "the party which holds the response.")
    return(invisible(NULL))
  } else if (is.numeric(bins)) {
    temp <- list()
    temp <- roc_internal(x, bins = bins)
    class(temp) <- "rocdistributed"
    return(temp)
  }
}


################### COX DISPLAY SURVFIT AND STRATA FUNCTIONS ###################

#' @importFrom grDevices rgb
GetColors <- function(n) {
  color <- matrix(0, 6, 3)
  color[1, ] <- c(0.000, 0.000, 1.000) # blue
  color[2, ] <- c(0.627, 0.125, 0.941) # purple
  color[3, ] <- c(1.000, 0.000, 0.000) # red
  color[4, ] <- c(1.000, 0.647, 0.000) # orange
  color[5, ] <- c(0.000, 1.000, 0.000) # green

  if (n == 1) {
    return(rgb(0, 0, 0))
  }
  if (n <= 3) {
    cols <- c(rgb(color[1, 1], color[1, 2], color[1, 3]),
              rgb(color[2, 1], color[2, 2], color[2, 3]),
              rgb(color[3, 1], color[3, 2], color[3, 3]))
    return(cols[1:n])
  }
  cols <- c()
  for (i in 1:n) {
    idx <- 4 * (i - 1) / (n - 1) + 1 # If we add yellow back in, change 4 to 5
    idx1 <- floor(idx)
    idx2 <- ceiling(idx)
    dx   <- idx - idx1
    tcol <- color[idx1, ] + (color[idx2, ] - color[idx1, ]) * dx
    cols <- c(cols, rgb(tcol[1], tcol[2], tcol[3]))
  }
  return(cols)
}


#' @title Plotting Survival Curves for Vertical Distributed Cox Regression
#' @description Plots a survivial curve as specified by
#'   \code{survfitDistributed} object.
#' @param x A \code{survfitDistributed} object.
#' @param merge Logical.  It \code{TRUE}, plots all strata of the survival curve
#'   on one plot.  If \code{FALSE}, plots all strata in different plots.
#' @param ... Common graphical parameters (not fully implemented).
#' @return No return value.
#' @seealso \code{\link{survfitDistributed}}
#' @examples
#'  sfit <- survfitDistributed(vdra_fit_cox_A)
#'  plot(sfit)
#'
#'  # From Data Partner 1
#'  sfit <- survfitDistributed(vdra_fit_cox_A,
#'                            ~Exposure,
#'                            data = vdra_data[, c(3:4, 5:7)])
#'  plot(sfit)
#'  plot(sfit, merge = FALSE)
#'
#'  # From Data Partner 2
#'  sfit <- survfitDistributed(vdra_fit_cox_B,
#'                            ~Race + Sex,
#'                            data = vdra_data[, 8:11])
#'  plot(sfit, merge = FALSE)
#' @importFrom graphics legend lines
#' @export
plot.survfitDistributed <- function(x, merge = FALSE, ...) {
  max <- 0
  n <- length(x$strata)
  labels <- c()
  max <- max(x$time)
  labels <- names(x$strata)
  arguments <- list(...)
  arguments$x <- 1
  arguments$type <- "n"
  if (is.null(arguments$ylim)) arguments$ylim <- c(0, 1)
  if (is.null(arguments$xlim)) arguments$xlim <- c(0, max)
  if (is.null(arguments$xlab)) arguments$xlab <- "Time"
  if (is.null(arguments$ylab)) arguments$ylab <- "Percent Survival"
  if (is.null(arguments$main)) arguments$main <- "Survival Curve"

  if (merge) {
    do.call("plot", arguments)
    cols <- GetColors(n)
    start <- 1
    for (i in 1:n) {
      end <- start + x$strata[i] - 1
      lines(c(1, x$surv[start:end]) ~
              c(0, x$time[start:end]), type = "s", col = cols[i])
      start <- end + 1
    }
    legend("bottomleft", legend = labels, col = cols, lty = 1)
  } else {
    cols <- GetColors(1)
    start <- 1
    for (i in 1:n) {
      end <- start + x$strata[i] - 1
      do.call("plot", arguments)
      lines(c(1, x$surv[start:end]) ~
              c(0, x$time[start:end]), type = "s", col = cols)
      legend("bottomleft", legend = labels[i], col = cols, lty = 1)
      start <- end + 1
    }
  }
}


#' @rdname survfitDistributed
#' @export
print.survfitDistributed <- function(x, ...) {
  start <- 1
  events <- integer(length(x$strata))
  for (i in seq_along(x$strata)) {
    end <- start + x$strata[i] - 1
    events[i] <- sum(x$n.event[start:end])
    start <- end + 1
  }
  df <- data.frame(n = x$n, events = events)
  row.names(df) <- names(x$strata)
  print(df)
}


#' @rdname survfitDistributed
#' @export
survfitDistributed.stats <- function(x) {
  surv          <- list()
  surv$n        <- x$strata$end - x$strata$start + 1
  for (i in seq_len(nrow(x$strata))) {
    start <- x$strata$start[i]
    end   <- x$strata$end[i]
    idx   <- which(c(1, diff(x$survival$rank[start:end])) != 0)
    temp0 <- table(x$survival$rank[start:end], x$survival$status[start:end])
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 <- cbind(temp0, 0)
        colnames(temp0) <- c("0", "1")
      } else {
        temp0 <- cbind(0, temp0)
        colnames(temp0) <- c("0", "1")
      }
    }

    surv$time     <- c(surv$time, x$survival$rank[start:end][idx])
    surv$n.risk   <- c(surv$n.risk, rev(cumsum(rev(temp0[, 1] + temp0[, 2]))))
    surv$n.event  <- c(surv$n.event, temp0[, 2])
    surv$n.censor <- c(surv$n.censor, temp0[, 1])
    surv$strata   <- c(surv$strata, length(idx))
    surv$surv     <- c(surv$surv, x$survival$surv[start:end][idx])
  }
  names(surv$n.risk) <- NULL
  names(surv$n.event) <- NULL
  names(surv$n.censor) <- NULL
  names(surv$strata) <- x$strata$label
  surv$type     <- "right"
  class(surv) <- "survfitDistributed"
  return(invisible(surv))
}

#' @rdname survfitDistributed
#' @export
survfitDistributed.formula <- function(x, formula, data) {
  surv <- list()
  vars <- all.vars(formula)
  if ("." %in% vars) {
    warning("This function does not allow the . symbol in formulas.")
    return(invisible(NULL))
  }
  if (!all(vars %in% colnames(data))) {
    warning("Not all strata are found in the data.")
    return(invisible(NULL))
  }
  if (length(vars) == 0) {
    data <- data.frame(const__ = rep(1, length(x$survival$rank)))
  } else {
    idx    <- which(colnames(data) %in% vars)
    data   <- data[x$survival$sorted, idx, drop = FALSE]
  }
  sorted <- do.call("order", as.data.frame(cbind(data,
                                                 x$survival$rank,
                                                 x$survival$status)))
  data   <- data[sorted, , drop = FALSE]
  rank   <- x$survival$rank[sorted]
  status <- x$survival$status[sorted]
  data2  <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  legend <- list()
  colnames(data2) <- colnames(data)
  for (i in seq_len(ncol(data))) {
    levels <- levels(as.factor(data[, i]))
    legend[[colnames(data)[i]]] <- levels
    data2[, i] <- sapply(data[, i], function(x) {
      which(levels %in% x)
    })
  }
  ranks <- which(apply(abs(apply(data2, 2, diff)), 1, sum) > 0)
  ranks <- c(ranks, nrow(data2))
  start <- 1
  for (i in seq_along(ranks)) {
    end <- ranks[i]
    surv$n <- c(surv$n, end - start + 1)
    # Calculate the Kaplan Meier Curve Here per notes from 9/4/19
    rank2 <- rank[start:end]
    event2 <- status[start:end]
    temp <- table(rank2)
    m <- length(temp)
    temp0 <- table(rank2, event2)
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 <- cbind(temp0, 0)
        colnames(temp0) <- c("0", "1")
      } else {
        temp0 <- cbind(0, temp0)
        colnames(temp0) <- c("0", "1")
      }
    }
    idx <- which(temp0[, 2] > 0)
    if (temp0[nrow(temp0), 2] == 0) idx <- c(idx, nrow(temp0))
    n_fails <- temp0[idx, 2]
    start0 <- c(1, (cumsum(temp)[1:(m - 1)] + 1))[idx]
    start1 <- start0 + temp0[idx, 1]
    stop1  <- start1 + n_fails  - 1
    final <- length(rank2)
    S <- 1
    t2 <- rep(0, length(n_fails))
    S2 <- rep(0, length(n_fails))
    for (j in seq_along(n_fails)) {
      n <- final - start0[j] + 1
      d <- stop1[j] - start1[j] + 1
      S <- S * (n - d) / n
      t2[j] <- rank2[start0[j]]
      S2[j] <- S
    }
    surv$time <- c(surv$time, t2)
    surv$n.risk   <- c(surv$n.risk, rev(cumsum(rev(temp0[, 1] + temp0[, 2]))))
    surv$n.event <- c(surv$n.event, temp0[, 2])
    surv$n.censor <- c(surv$n.censor, temp0[, 1])
    surv$surv <- c(surv$surv, S2)
    surv$strata <- c(surv$strata, length(idx))
    if (length(vars) == 0) {
      names(surv$strata)[i] <- ""
    } else {
      label <- ""
      for (j in seq_len(ncol(data))) {
        temp <- colnames(data)[j]
        label <- paste0(label, temp, "=", legend[[temp]][data2[start, j]])
        if (j < ncol(data)) {
          label <- paste0(label, ", ")
        }
      }
      names(surv$strata)[i] <- label
    }
    start <- end + 1
  }
  surv$type            <- "right"
  names(surv$n.risk)   <- NULL
  names(surv$n.event)  <- NULL
  names(surv$n.censor) <- NULL
  class(surv) <- "survfitDistributed"
  return(invisible(surv))
}

#' @title Create Survival Curves for Vertical Distributed Cox Regression
#' @aliases survfitDistributed survfitDistributed.object
#'   print.survfitDistributed
#' @description This function creates survival curves for a previously defined
#'   \code{\link{vdracox}} object.  The function also accepts a formula and the
#'   original data supplied by the calling party allowing exploration of other
#'   potential strata.  Both \code{formula} and \code{data} must be NULL or both
#'   must be specified.
#' @param x an object of type \code{\link{vdracox}}.
#' @param formula a formula which defines alternative strata for the survival
#'   curve.
#' @param data if \code{formula} is specified, this should be the data that was
#'   supplied by the calling party.
#' @return Returns an object of class \code{survfitDistributed}. Objects of this
#'   class have methods for the functions \code{print} and \code{plot}. The
#'   following components must be included in a legitimate
#'   \code{survfitDistributed} object.
#'
#'   \item{n}{the total number of subjects in each curve.}
#'
#'   \item{time}{the time points at which the curve has a step.}
#'
#'   \item{n.risk}{the number of subjects at risk at each time point.}
#'
#'   \item{n.event}{the number of events that occour at each time point.}
#'
#'   \item{n.censor}{the number of subjects who are censored at each time
#'   point.}
#'
#'   \item{strata}{the number of points in each strata.}
#'
#'   \item{surv}{the estimate of the survival time at each time step.}
#'
#'   \item{type}{the type of censoring.  Currently, always "right".}
#'
#' @seealso \code{\link{plot.survfitDistributed}}
#' @examples
#'
#' sfit <- survfitDistributed(vdra_fit_cox_A)
#' print(sfit)
#' plot(sfit)
#'
#' # From Data Partner 1
#'
#' sfit <- survfitDistributed(vdra_fit_cox_A,
#'                           ~Exposure,
#'                           data = vdra_data[, c(3:4, 5:7)])
#' print(sfit)
#' plot(sfit)
#'
#' # From Data Partner 2
#'
#' sfit <- survfitDistributed(vdra_fit_cox_B,
#'                           ~Race + Sex,
#'                           data = vdra_data[, 8:11])
#' print(sfit)
#' plot(sfit, merge = TRUE)
#' @export
survfitDistributed <- function(x = NULL, formula = NULL, data = NULL) {
  if (!inherits(x, "vdracox")) {
    warning("The first parameter must be a vdracox object.")
    return(invisible(NULL))
  }
  if (is.null(data) && is.null(formula)) {
    return(survfitDistributed.stats(x))
  }
  if (!is.matrix(data) && !is.data.frame(data)) {
    warning(paste("the data must either be a matrix or a data.frame.",
                  "Please use the same data that you used for",
                  "the distributed regression."))
    return(invisible(NULL))
  }
  if (!inherits(formula, "formula") || !validFormula2(formula)) {
    warning(paste("The formula must be of the form \"~ var1 + ... + vark\"",
                  "where the variables",
                  "are found in the data. The formula can also be \"~ 1\"."))
  }
  return(survfitDistributed.formula(x, formula, data))
}
