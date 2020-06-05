########################### 2 PARTY GLOBAL FUNCTIONS ###########################

AnalysisCenter.2Party = function(regression            = "linear",
                                 data                  = NULL,
                                 response              = NULL,
                                 strata                = NULL,
                                 mask                  = TRUE,
                                 monitorFolder         = getwd(),
                                 msreqid               = "v_default_00_000",
                                 blocksize             = 500,
                                 tol                   = 1e-8,
                                 maxIterations         = 25,
                                 sleepTime             = 10,
                                 maxWaitingTime        = 24 * 60 * 60,
                                 popmednet             = TRUE,
                                 trace                 = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (regression == "cox") {
    stats = PartyAProcess2Cox(data, response, strata, mask, monitorFolder,
                              msreqid, blocksize, tol, maxIterations,
                              sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "linear") {
    stats = PartyAProcess2Linear(data, response, monitorFolder, msreqid,
                                 blocksize, sleepTime, maxWaitingTime,
                                 popmednet, trace)
  } else  if (regression == "logistic") {
    stats = PartyAProcess2Logistic(data, response, monitorFolder, msreqid,
                                   blocksize, tol, maxIterations, sleepTime,
                                   maxWaitingTime, popmednet, trace)
  }

  else {
    cat("Regression type must be \"cox\", \"linear\" or \"logistic\"\n\n")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n")
  options(digits.secs = digits.secs)
  return(stats)
}


DataPartner.2Party = function(regression          = "linear",
                              data                = NULL,
                              strata              = NULL,
                              mask                = TRUE,
                              monitorFolder       = getwd(),
                              sleepTime           = 10,
                              maxWaitingTime      = 24 * 60 * 60,
                              popmednet           = TRUE,
                              trace               = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (regression == "cox") {
    stats = PartyBProcess2Cox(data, strata, mask,
                              monitorFolder, sleepTime, maxWaitingTime,
                              popmednet, trace)
  } else if (regression == "linear") {
    stats = PartyBProcess2Linear(data, monitorFolder, sleepTime, maxWaitingTime,
                                 popmednet, trace)
  } else if (regression == "logistic") {
    stats = PartyBProcess2Logistic(data, monitorFolder, sleepTime, maxWaitingTime,
                                   popmednet, trace)
  } else {
    cat("Regression type must be \"cox\", \"linear\" or \"logistic\"\n\n")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n")
  options(digits.secs = digits.secs)
  return(stats)
}

########################### 3 PARTY GLOBAL FUNCTIONS ###########################

DataPartner1.3Party = function(regression            = "linear",
                               data                  = NULL,
                               response              = NULL,
                               strata                = NULL,
                               mask                  = TRUE,
                               monitorFolder         = getwd(),
                               sleepTime             = 10,
                               maxWaitingTime        = 24 * 60 * 60,
                               popmednet             = TRUE,
                               trace                 = FALSE) {

  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (regression == "cox") {
    stats = PartyAProcess3Cox(data, response, strata, mask, monitorFolder,
                              sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "linear") {
    stats = PartyAProcess3Linear(data, response, monitorFolder,
                                 sleepTime, maxWaitingTime, popmednet, trace)
  } else  if (regression == "logistic") {
    stats = PartyAProcess3Logistic(data, response, monitorFolder,
                                   sleepTime, maxWaitingTime, popmednet, trace)
  } else {
    cat("Regression type must be \"cox\", \"linear\" or \"logistic\"\n\n")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n")
  options(digits.secs = digits.secs)
  return(stats)
}


DataPartner2.3Party = function(regression          = "linear",
                               data                = NULL,
                               strata              = NULL,
                               mask                = TRUE,
                               monitorFolder       = getwd(),
                               sleepTime           = 10,
                               maxWaitingTime      = 24 * 60 * 60,
                               popmednet           = TRUE,
                               trace               = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (regression == "cox") {
    stats = PartyBProcess3Cox(data, strata, mask, monitorFolder,
                              sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "linear") {
    stats = PartyBProcess3Linear(data, monitorFolder,
                                 sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "logistic") {
    stats = PartyBProcess3Logistic(data, monitorFolder,
                                   sleepTime, maxWaitingTime, popmednet, trace)
  } else {
    cat("Regression type must be \"cox\", \"linear\" or \"logistic\"\n\n")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n")
  options(digits.secs = digits.secs)
  return(stats)
}


AnalysisCenter.3Party = function(regression            = "linear",
                                 monitorFolder         = getwd(),
                                 msreqid               = "v_default_00_000",
                                 blocksize             = 500,
                                 tol                   = 1e-8,
                                 maxIterations         = 25,
                                 sleepTime             = 10,
                                 maxWaitingTime        = 24 * 60 * 60,
                                 popmednet             = TRUE,
                                 trace                 = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (regression == "cox") {
    stats = PartyTProcess3Cox(monitorFolder, msreqid, blocksize, tol,
                              maxIterations, sleepTime, maxWaitingTime,
                              popmednet, trace)
  } else if (regression == "linear") {
    stats = PartyTProcess3Linear(monitorFolder, msreqid, blocksize,
                                 sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "logistic") {
    stats = PartyTProcess3Logistic(monitorFolder, msreqid, blocksize, tol,
                                   maxIterations, sleepTime, maxWaitingTime,
                                   popmednet, trace)
  } else {
    cat("Regression type must be \"cox\", \"linear\" or \"logistic\"\n\n")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n")
  options(digits.secs = digits.secs)
  return(stats)
}

########################### K PARTY GLOBAL FUNCTIONS ###########################

DataPartner.KParty = function(regression            = "linear",
                              data                  = NULL,
                              response              = NULL,
                              strata                = NULL,
                              mask                  = TRUE,
                              numDataPartners       = NULL,
                              dataPartnerID         = NULL,
                              monitorFolder         = getwd(),
                              sleepTime             = 10,
                              maxWaitingTime        = 24 * 60 * 60,
                              popmednet             = TRUE,
                              trace                 = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (is.null(numDataPartners)) {
    stop("numDataPartners must be specified\n")
  }
  if (is.null(dataPartnerID)) {
    stop("dataPartnerID must be specified")
  }

  if (regression == "cox") {
    stats = DataPartnerKCox(data, response, strata, mask, numDataPartners,
                            dataPartnerID, monitorFolder,
                            sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "linear") {
    stats = DataPartnerKLinear(data, response, numDataPartners,
                               dataPartnerID, monitorFolder,
                               sleepTime, maxWaitingTime, popmednet, trace)
  } else  if (regression == "logistic") {
    stats = DataPartnerKLogistic(data, response, numDataPartners,
                                 dataPartnerID, monitorFolder,
                                 sleepTime, maxWaitingTime, popmednet, trace)
  } else {
    stop("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n\n")
  options(digits.secs = digits.secs)
  return(stats)
}


AnalysisCenter.KParty = function(regression            = "linear",
                                 numDataPartners       = NULL,
                                 monitorFolder         = getwd(),
                                 msreqid               = "v_default_00_000",
                                 tol                   = 1e-8,
                                 maxIterations         = 25,
                                 sleepTime             = 10,
                                 maxWaitingTime        = 24 * 60 * 60,
                                 popmednet             = TRUE,
                                 trace                 = FALSE) {
  digits.secs = getOption("digits.secs")
  options(digits.secs = 2)
  startTime = proc.time()
  cat("Process started on", as.character(GetUTCTime()), "UTC.\n")
  stats = list()
  if (is.null(numDataPartners)) {
    stop("numDataPartners must be specified.")
    return(NULL)
  }

  if (regression == "cox") {
    stats = AnalysisCenterKCox(numDataPartners, monitorFolder, msreqid, tol,
                               maxIterations, sleepTime, maxWaitingTime,
                               popmednet, trace)
  } else if (regression == "linear") {
    stats = AnalysisCenterKLinear(numDataPartners, monitorFolder, msreqid,
                                  sleepTime, maxWaitingTime, popmednet, trace)
  } else if (regression == "logistic") {
    stats = AnalysisCenterKLogistic(numDataPartners, monitorFolder, msreqid, tol,
                                    maxIterations, sleepTime, maxWaitingTime,
                                    popmednet, trace)
  } else {
    stop("Regression type must be \"cox\", \"linear\" or \"logistic\"")
  }

  elp = GetElapsedTime(proc.time() - startTime, final = TRUE, timeOnly = FALSE)
  cat("Process completed on", as.character(GetUTCTime()), "UTC.\n")
  cat(elp, "\n\n")
  options(digits.secs = digits.secs)
  return(stats)
}

############################ SHARED SETUP FUNCTIONS ############################

CreateIOLocation = function(monitorFolder, folder) {
  location = file.path(monitorFolder, folder)
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

########################### 2 PARTY SETUP FUNCTIONS ############################

PrepareParams.2p = function(analysis, party, msreqid = "v_default_00_000",
                            popmednet = TRUE, trace = FALSE) {
  params                     = list()
  params$partyName           = party
  params$analysis            = analysis
  params$msreqid             = msreqid
  params$popmednet           = popmednet
  params$trace               = trace
  params$failed              = FALSE
  params$errorMessage        = ""
  params$pmnStepCounter      = 0
  params$algIterationCounter = 0
  params$completed           = FALSE
  params$converged           = FALSE
  params$maxIterExceeded     = FALSE
  params$lastIteration       = FALSE
  params$p1                  = 0
  params$p2                  = 0
  params$p1.old              = 0
  params$p2.old              = 0
  params$stats               = list()
  class(params$stats)        = paste0("vdra", analysis)
  params$stats$failed        = TRUE
  params$stats$converged     = FALSE
  return(params)
}

########################### 3 PARTY SETUP FUNCTIONS ############################

PrepareParams.3p = function(analysis, party, msreqid = "v_default_00_000",
                            popmednet = TRUE, trace = FALSE) {
  params                     = list()
  params$partyName           = party
  params$analysis            = analysis
  params$msreqid             = msreqid
  params$trace               = trace
  params$popmednet           = popmednet
  params$failed              = FALSE
  params$errorMessage        = ""
  params$pmnStepCounter      = 0
  params$algIterationCounter = 0
  params$completed           = FALSE
  params$converged           = FALSE
  params$maxIterExceeded     = FALSE
  params$lastIteration       = FALSE
  params$p1                  = 0
  params$p2                  = 0
  params$p1.old              = 0
  params$p2.old              = 0
  params$stats               = list()
  class(params$stats)        = paste0("vdra", analysis)
  params$stats$failed        = TRUE
  params$stats$converged     = FALSE
  return(params)
}

########################### K PARTY SETUP FUNCTIONS ############################

PrepareParams.kp = function(analysis, dataPartnerID, numDataPartners,
                            msreqid = "v_default_00_000", cutoff = NULL,
                            maxIterations = NULL, ac = FALSE, popmednet = TRUE,
                            trace = FALSE) {
  params                     = list()
  params$dataPartnerID       = dataPartnerID
  params$numDataPartners     = numDataPartners
  params$analysis            = analysis
  params$msreqid             = msreqid
  params$trace               = trace
  params$popmednet           = popmednet
  params$failed              = FALSE
  params$errorMessage        = ""
  params$pmnStepCounter      = 0
  params$algIterationCounter = 0
  params$maxIterations       = maxIterations
  params$completed           = FALSE
  params$converged           = FALSE
  params$maxIterExceeded     = FALSE
  params$lastIteration       = FALSE
  params$cutoff              = cutoff
  params$stats               = list()
  class(params$stats)        = paste0("vdra", analysis)
  params$stats$failed        = TRUE
  params$stats$converged     = FALSE
  if ((class(numDataPartners) != "integer" &&
      class(numDataPartners) != "numeric") ||
      numDataPartners <= 0 ||
      is.infinite(numDataPartners) ||
      round(numDataPartners) != numDataPartners) {
    params$failed = TRUE
    params$errorMessage = "Error: numDataPartners must be a positive integer, and must equal the number of data partners providing data.\n\n"
  }
  if (!params$failed) {
    if (ac) {
      if (dataPartnerID != 0) {
        params$failed = TRUE
        params$errorMessage = "Error: dataPartnerID for Analysis Center must be 0.\n\n"
      }
    } else {
      if (dataPartnerID <= 0 || dataPartnerID > numDataPartners) {
        params$failed = TRUE
        params$errorMessage = paste0("Error: dataPartnerID must be between 1 and ", numDataPartners, " inclusive.\n\n")
      }
    }
  }
  return(params)
}

########################### PRETTY OUTPUT FUNCTIONS ############################

Header = function(params) {
  large.cox = c("  ____ _____  __",
                " / ___/ _ \\ \\/ /",
                "| |  | | | \\  / ",
                "| |__| |_| /  \\ ",
                " \\____\\___/_/\\_\\")
  large.linear = c(" _     ___ _   _ _____    _    ____  ",
                   "| |   |_ _| \\ | | ____|  / \\  |  _ \\ ",
                   "| |    | ||  \\| |  _|   / _ \\ | |_) |",
                   "| |___ | || |\\  | |___ / ___ \\|  _ < ",
                   "|_____|___|_| \\_|_____/_/   \\_|_| \\_\\")
  large.logistic = c(" _     ___   ____ ___ ____ _____ ___ ____ ",
                     "| |   / _ \\ / ___|_ _/ ___|_   _|_ _/ ___|",
                     "| |  | | | | |  _ | |\\___ \\ | |  | | |    ",
                     "| |__| |_| | |_| || | ___) || |  | | |___ ",
                     "|_____\\___/ \\____|___|____/ |_| |___\\____|")
  large.regression = c(" ____  _____ ____ ____  _____ ____ ____ ___ ___  _   _ ",
                       "|  _ \\| ____/ ___|  _ \\| ____/ ___/ ___|_ _/ _ \\| \\ | |",
                       "| |_) |  _|| |  _| |_) |  _| \\___ \\___ \\| | | | |  \\| |",
                       "|  _ <| |__| |_| |  _ <| |___ ___) ___) | | |_| | |\\  |",
                       "|_| \\_|_____\\____|_| \\_|_____|____|____|___\\___/|_| \\_|")
  small.cox = c("  ___ _____  __",
                " / __/ _ \\ \\/ /",
                "| (_| (_) >  < ",
                " \\___\\___/_/\\_\\")
  small.linear = c(" _    ___ _  _ ___   _   ___ ",
                   "| |  |_ _| \\| | __| /_\\ | _ \\",
                   "| |__ | || .` | _| / _ \\|   /",
                   "|____|___|_|\\_|___/_/ \\_\\_|_\\")
  small.logistic = c(" _    ___   ___ ___ ___ _____ ___ ___ ",
                     "| |  / _ \\ / __|_ _/ __|_   _|_ _/ __|",
                     "| |_| (_) | (_ || |\\__ \\ | |  | | (__ ",
                     "|____\\___/ \\___|___|___/ |_| |___\\___|")
  small.regression = c(" ___ ___ ___ ___ ___ ___ ___ ___ ___  _  _ ",
                       "| _ \\ __/ __| _ \\ __/ __/ __|_ _/ _ \\| \\| |",
                       "|   / _| (_ |   / _|\\__ \\__ \\| | (_) | .` |",
                       "|_|_\\___\\___|_|_\\___|___/___/___\\___/|_|\\_|")
  tiny.cox = c("+-+-+-+",
               "|C|O|X|")
  tiny.linear = c("+-+-+-+-+-+-+",
                  "|L|I|N|E|A|R|")
  tiny.logistic = c("+-+-+-+-+-+-+-+-+",
                    "|L|O|G|I|S|T|I|C|")
  tiny.regression = c("+-+-+-+-+-+-+-+-+-+-+",
                      "|R|E|G|R|E|S|S|I|O|N|",
                      "+-+-+-+-+-+-+-+-+-+-+")

  width = getOption("width")
  if (width > nchar(large.regression[1])) {
    cox        = large.cox
    linear     = large.linear
    logistic   = large.logistic
    regression = large.regression
  } else if (width  > nchar(small.regression[1])) {
    cox        = small.cox
    linear     = small.linear
    logistic   = small.logistic
    regression = small.regression
  } else {
    cox        = tiny.cox
    linear     = tiny.linear
    logistic   = tiny.logistic
    regression = tiny.regression
  }

  offset.cox        = floor((width - nchar(cox[1])) / 2)
  offset.linear     = floor((width - nchar(linear[1])) / 2)
  offset.logistic   = floor((width - nchar(logistic[1])) / 2)
  offset.regression = floor((width - nchar(regression[1])) / 2)

  space.cox        = paste(rep(" ", offset.cox), collapse = "")
  space.linear     = paste(rep(" ", offset.linear), collapse = "")
  space.logistic   = paste(rep(" ", offset.logistic), collapse = "")
  space.regression = paste(rep(" ", offset.regression), collapse = "")

  if (params$analysis == "linear") {
    cat(paste0("\r", space.linear, linear, "\n"))
    cat(paste0("\r", space.regression, regression, "\n"))
  }
  if (params$analysis == "logistic") {
    cat(paste0("\r", space.logistic, logistic, "\n"))
    cat(paste0("\r", space.regression, regression, "\n"))
  }
  if (params$analysis == "cox") {
    cat(paste0("\r", space.cox, cox, "\n"))
    cat(paste0("\r", space.regression, regression, "\n"))
  }
  cat("\n")
}


BeginningIteration = function(params) {
  width  = getOption("width")
  msg    = paste("*** Beginning Iteration", params$algIterationCounter, "***")
  offset = max(floor((width - nchar(msg)) / 2) - 1, 0)
  space  = paste(rep(" ", offset), collapse = "")
  cat(space, msg, "\n\n")
}


EndingIteration = function(params) {
  width  = getOption("width")
  msg    = paste("*** Ending Iteration", params$algIterationCounter, "***")
  offset = floor((width - nchar(msg)) / 2) - 1
  space  = paste(rep(" ", offset), collapse = "")
  cat(space, msg, "\n\n")
}


GetLion = function(p) {
  lion = rep("", p)
  lion2 = rep("", 13)
  lion2[1 ] = "    (\"`-''-/\").___..--''\"`-._\"      "
  lion2[2 ] = "    `6_ 6 ) `-. ( ).`-.__.`)        "
  lion2[3 ] = "    (_Y_.)' ._ ) `._ `. ``-..-'     "
  lion2[4 ] = "     _..`--'_..-_/ /--'_.' ,'       "
  lion2[5 ] = "    (il),-'' (li),' ((!.-'          "
  lion2[6 ] = "        ___  ___  _  _  _  _        "
  lion2[7 ] = "       | -_>||__>|\\ |||\\ ||       "
  lion2[8 ] = "       | |  ||__>| \\||| \\||       "
  lion2[9 ] = "       |_|  ||__>|_\\_||_\\_|       "
  lion2[10] = "     ___  _____  ___  _____  ___    "
  lion2[11] = "    //__>|_   _|//_\\|_   _|||__>   "
  lion2[12] = "    \\_\\  | |  | | |  | |  ||__>   "
  lion2[13] = "    <__//  |_|  |_|_|  |_|  ||__>   "

  if (p >= 12) {
    if (p < 16) kk = p - 11
    else kk = ceiling(p / 2) - 7
    lion[kk:(kk + 12)] = lion2
  } else if (p >= 8) {
    kk = p - 7
    if (p >= 10) ceiling(p / 2) - 3
    lion2[6 ]  = "        ___   ___   _   _   "
    lion2[7 ]  = "       | -_> //__> | | | |  "
    lion2[8 ]  = "       | |   \\_\\ | |_| |  "
    lion2[9 ]  = "       |_|   <__// \\___//  "
    lion[kk:(kk + 8)] = lion2[1:9]
  } else if (p >= 4) {
    kk = p - 3
    if (p >= 5) kk = ceiling(p / 2) - 1
    lion2[5 ] = "    (il),-'' (li),' ((!.-'      PSU "
    lion[kk:(kk + 4)] = lion2[1:5]
  }
  return(lion)
}


MakeProgressBar1 = function(steps, message) {
  pb = list()
  messageLength    = 18
  pb$numSteps      = steps
  pb$numBlanks     = 20
  pb$delimeter     = "|"
  pb$filler        = "#"
  pb$blank         = "."
  pb$percent       = 0
  pb$percentstr    = "  0%"
  pb$prints = 0
  message = substr(message, 1, messageLength)
  message = paste0(message,
                   paste(rep(" ", messageLength - nchar(message)),
                         collapse = ""))
  pb$header = paste0("Processing ", message, ": ")
  toPrint = paste0(pb$header, pb$percentstr, pb$delimeter,
                   paste(rep(pb$blank, pb$numBlanks), collapse = ""), pb$delimeter)
  cat(toPrint, "\r")
  flush.console()
  return(pb)
}


MakeProgressBar2 = function(i, pb) {
  percent = floor(100 * i / pb$numSteps)
  if (percent == pb$percent) {
    return(pb)
  }
  pb$percent = percent
  pb$percentstr = paste0(paste(rep(" ", 3 - nchar(percent)), collapse = ""), percent, "%")
  numFiller = floor(pb$numBlanks * i / pb$numSteps)
  toPrint = paste0(pb$header, pb$percentstr, pb$delimeter,
                   paste(rep(pb$filler, numFiller), collapse = ""),
                   paste(rep(pb$blank, pb$numBlanks - numFiller), collapse = ""),
                   pb$delimeter)

  if (i == pb$numSteps) {
    cat(toPrint, "\n\n")
  } else {
    cat(toPrint, "\r")
  }
  flush.console()
  return(pb)
}

############################### MATRIX FUNCTIONS ###############################

MultiplyDiagonalWTimesX = function(w, x) {
  if (!is.matrix(x)) {
    x = matrix(x, length(x), 1)
    wx1 = matrix(NA, length(x), 1)
  } else {
    wx1 = matrix(NA, nrow = nrow(x), ncol = ncol(x))
  }
  if (is.matrix(w)) {
    for (i in 1:nrow(w)) {
      wx1[i, ] = w[i] * x[i, ]
    }
  } else {
    for (i in 1:length(w)) {
      wx1[i, ] = w[i] * x[i, ]
    }
  }

  return(wx1)
}


FindOrthogonalVectors = function(x, g) {
  x = as.matrix(x)
  x = cbind(x, runif(nrow(x)))  # Randomize Z
  # Save the Random Vector Here
  n = nrow(x)
  Q = qr.Q(qr(x), complete = TRUE)
  Q = Q[, (n - g + 1):n]
  return(Q)
}


RandomOrthonomalMatrix = function(size) {
  return(qr.Q(qr(matrix(runif(size * size), size, size)), complete = TRUE))
}

########################## SHARED PMN COMM FUNCTIONS ###########################

MakeCSV = function(file_nm, transfer_to_site_in, dp_cd_list, writePath) {
  dframe = data.frame(file_nm, transfer_to_site_in, dp_cd_list)
  fp = file.path(writePath, "file_list.csv")
  write.csv(dframe, fp, row.names = FALSE, quote = FALSE)
}


SeqZW = function(letter = "Z_", nblocks = 1) {
  return(paste0(letter, 1:nblocks, ".rdata"))
}


Standby = function(triggerName, triggerLocation,
                   sleepTime = 1, maxWaitingTime = NULL, remove = FALSE) {

  found = FALSE

  if (is.null(maxWaitingTime)) { maxWaitingTime = 60 * 60 * 24 }

  fpath = file.path(triggerLocation, triggerName)
  startTime = proc.time()[3]
  elapsedTime = 0

  while (!found) {
    found = all(file.exists(fpath))

    if (elapsedTime > maxWaitingTime) {
      break
    }

    if (!found) {
      Sys.sleep(sleepTime)
      elapsedTime = round(proc.time()[3] - startTime, 0)
    }

    cat("Elapsed Time:", HMS(elapsedTime), "\r")
  }
  cat("\n")
  if (!found) {
    stop("Exceeded maximum time waiting for files to be dropped.")
    return("exmwt")
  }

  if (remove) DeleteTrigger(triggerName, triggerLocation)

}

CopyFile = function(readDirectory, writeDirectory, filename) {
  source      = file.path(readDirectory, filename)
  destination = file.path(writeDirectory, filename)
  if (all(file.exists(source))) {
    file.copy(source, destination, overwrite = TRUE)
  } else {
    cat("These files do not exist:\n", source[!file.exists(source)], "\n")
  }
}


MakeTrigger = function(triggerName, triggerPath, message = "Trigger File") {

  fn = file.path(triggerPath, triggerName)
  if (file.exists(fn)) {
    file.remove(fn)
  }

  write(message, fn)
}


DeleteTrigger = function(triggerName, triggerPath) {
  Sys.sleep(1)
  targets = file.path(triggerPath, triggerName)
  for (target in targets) {
    if (file.exists(target)) {
      startTime = proc.time()[3]
      repeat {
        result = suppressWarnings(try(file.remove(target)))
        if (result) break
        Sys.sleep(1)
        if (proc.time()[3] - startTime > 60) {
          stop(paste("Error: Could not delete the file", target, "after 60 seconds."))
          return(FALSE)
        }
      }
    }
  }
}


MakeTransferMessage = function(writePath) {
  message = "A has no covariates."
  save(message, file = file.path(writePath, "transferControl.rdata"))
}


MakeErrorMessage = function(writePath, message = "") {
  save(message, file = file.path(writePath, "errorMessage.rdata"))
}


ReadErrorMessage = function(readPath) {
  load(file.path(readPath, "errorMessage.rdata"))
  return(message)
}

########################## 2 PARTY PMN COMM FUNCTIONS ##########################

SendPauseQuit.2p = function(params,
                         files = c(),
                         sleepTime = 10,
                         job_failed = FALSE) {

  params = StoreLogEntry.2p(params, files)
  params = StoreTrackingTableEntry.2p(params)
  WriteLogCSV(params)
  WriteLogRaw(params)
  params$lastIteration = TRUE
  files = c(files, "stamps.rdata", "log.rdata", "file_list.csv")
  transfer = c(rep(1, length(files) - 1), 10)
  if (params$partyName == "A") {
    if (job_failed) {
      files = c(files, "job_fail.ok")
      params = StoreStampsEntry(params, "Job failed trigger file", "Trigger File created")
    } else {
      files = c(files, "job_done.ok")
      params = StoreStampsEntry(params, "Job done trigger file", "Trigger File created")
    }
    transfer = c(transfer, 10)
  }
  if (params$partyName == "A") {
    files = c(files, "dl_track_tbl.csv")
    transfer = c(transfer, 10)
    destination = rep(1, length(files))
    destination[transfer == 10] = 10
  } else {
    files = c(files, "tr_tb_updt.rdata")
    transfer = c(transfer, 1)
    destination = rep(0, length(files))
    destination[transfer == 10] = 10
  }
  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File Created")
  params = StoreStampsEntry(params, "R program execution complete, output files written",
                            "Tracking Table")
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  Sys.sleep(sleepTime)
  if (job_failed)  {
    MakeTrigger("job_fail.ok",  params$writePath)
  } else {
    MakeTrigger("job_done.ok",    params$writePath)
  }
  MakeTrigger("files_done.ok", params$writePath)
  return(params)
}

SendPauseContinue.2p = function(params,
                             files = c(),
                             sleepTime = 10,
                             maxWaitingTime = NULL,
                             job_started = FALSE) {
  params = StoreLogEntry.2p(params, files)
  params = StoreTrackingTableEntry.2p(params)
  WriteLogCSV(params)
  WriteLogRaw(params)
  files = c(files, "stamps.rdata", "log.rdata", "file_list.csv")
  transfer = c(rep(1, length(files) - 1), 10)
  if (params$partyName == "A") {
    files = c(files, "dl_track_tbl.csv")
    transfer = c(transfer, 10)
    destination = rep(1, length(files))
    destination[transfer == 10] = 10
  } else {
    files = c("tr_tb_updt.rdata", files)
    transfer = c(1, transfer)
    destination = rep(0, length(files))
    destination[transfer == 10] = 10
  }
  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File created")
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  params$pmnStepCounter      = params$pmnStepCounter + 2
  Sys.sleep(sleepTime)
  if (job_started) {
    MakeTrigger("job_started.ok", params$writePath)
  } else {
    MakeTrigger("files_done.ok", params$writePath)
  }
  if (params$partyName == "A") {
  	cat("Waiting for data partner\n")
  } else {
  	cat("Waiting for analysis center\n")
  }
  tryCatch(Standby("files_done.ok", params$readPath, maxWaitingTime = maxWaitingTime))
  cat("Resuming local processing\n\n")
  DeleteTrigger("files_done.ok", params$readPath)
  params = ReadLogRaw.2p(params)
  params = NewLogEntry.2p(params)
  params = ReadStampsRaw.2p(params)
  params = StoreStampsEntry(params, "R program execution begins", "Tracking Table")
  if (params$partyName == "A") {
    params = ReadTrackingTableUpdate.2p(params)
  }
  return(params)
}


PauseContinue.2p = function(params, maxWaitingTime) {
  params = StoreLogEntry.2p(params, "")
  WriteLogCSV(params)
  if (params$partyName == "A") {
  	cat("Waiting for data partner\n")
  } else {
  	cat("Waiting for analysis center\n")
  }
  tryCatch(Standby("files_done.ok", params$readPath, maxWaitingTime = maxWaitingTime))
  cat("Resuming local processing\n\n")
  DeleteTrigger("files_done.ok", params$readPath)
  params = MergeLogRaw.2p(params)
  params = NewLogEntry.2p(params)
  params = MergeStampsRaw.2p(params)
  params = ReadTrackingTableUpdate.2p(params)
  WriteLogCSV(params)
  return(params)
}

########################## 3 PARTY PMN COMM FUNCTIONS ##########################

WaitForTurn.3p = function(params, sleepTime) {
  Sys.sleep(sleepTime)
  if ((params$partyName == "T") || (!params$popmednet)) return(NULL)

  cat("Waiting For Turn\n")
  startTime = proc.time()[3]
  cat("Elapsed Time:", HMS(0), "\r")

  if (exists("partyOffset")) {
    cat("\n\n")
    return()
  }

  partyOffset = 15

  modulus   = 2 * partyOffset
  if (params$partyName == "A") targetTime = 0
  if (params$partyName == "B") targetTime = partyOffset


  while (as.integer(Sys.time()) %% modulus != targetTime) {
    elapsedTime = round(proc.time()[3] - startTime, 0)
    cat("Elapsed Time:", HMS(elapsedTime), "\r")
    Sys.sleep(0.1)
  }
  cat("\n\n")
}


SendPauseQuit.3p = function(params,
                         filesA = NULL,
                         filesB = NULL,
                         filesT = NULL,
                         sleepTime = 10,
                         job_failed = FALSE,
                         waitForTurn = FALSE) {

  params$lastIteration = TRUE
  params$completed     = TRUE
  files = c(filesA, filesB, filesT, "file_list.csv")
  transfer = c(rep(1, length(files) - 1), 10)
  destination = c(rep(1, length(filesA)),
                  rep(2, length(filesB)),
                  rep(0, length(filesT)),
                  10)
  if (params$party != "T") {
    files = c(files, "stamps.rdata", "log.rdata")
    transfer = c(transfer, 1, 1)
    destination = c(destination, 0, 0)
  }
  if (params$party == "T") {
    if (job_failed) {
      files = c(files, "job_fail.ok")
      params = StoreStampsEntry(params, "Job failed trigger file", "Trigger File created")
    } else {
      files = c(files, "job_done.ok")
      params = StoreStampsEntry(params, "Job done trigger file", "Trigger File created")
    }
    transfer = c(transfer, 10)
    destination = c(destination, 10)
  }
  params = StoreLogEntry.3p(params, c(filesA, filesB, filesT))
  params = StoreTrackingTableEntry.3p(params)
  WriteLogCSV(params)
  WriteLogRaw(params)

  if (params$partyName == "T") {
    WriteTrackingTableCSV(params)
    files = c(files, "dl_track_tbl.csv")
    transfer = c(transfer, 10)
    destination = c(destination, 10)
  } else {
    WriteTrackingTableRaw(params)
    files = c(files, "tr_tb_updt.rdata")
    transfer = c(transfer, 1)
    destination = c(destination, 0)
  }
  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File Created")
  params = StoreStampsEntry(params, "R program execution complete, output files written",
                            "Tracking Table")
	if (waitForTurn) {
	  params = StoreStampsEntry(params, "R program execution delayed", "Tracking Table")
    WaitForTurn.3p(params, sleepTime)
    params = StoreStampsEntry(params, "R program execution restarted", "Tracking Table")
	}
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (params$party == "T") {
    if (job_failed)  {
      MakeTrigger("job_fail.ok",  params$writePath)
    } else {
      MakeTrigger("job_done.ok",  params$writePath)
    }
  }
  MakeTrigger("files_done.ok", params$writePath)
  return(params)
}

SendPauseContinue.3p = function(params,
                             filesA = NULL,
                             filesB = NULL,
                             filesT = NULL,
                             from   = NULL,
                             sleepTime = 10,
                             maxWaitingTime = 24 * 60 * 60,
                             job_started = FALSE,
                             waitForTurn = FALSE) {
  params = StoreLogEntry.3p(params, c(filesA, filesB, filesT))
  params = StoreTrackingTableEntry.3p(params)
  WriteLogCSV(params)
  WriteLogRaw(params)

  files = c(filesA, filesB, filesT, "file_list.csv")
  transfer = c(rep(1, length(files) - 1), 10)
  destination = c(rep(1, length(filesA)),
                  rep(2, length(filesB)),
                  rep(0, length(filesT)), 10)
  if (length(files) > 1) {
    WriteTrackingTableRaw(params)
  }
  if (!is.null(filesA)) {
    files = c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer = c(transfer, 1, 1, 1)
    destination = c(destination, 1, 1, 1)
  }
  if (!is.null(filesB)) {
    files = c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer = c(transfer, 1, 1, 1)
    destination = c(destination, 2, 2, 2)
  }
  if (!is.null(filesT)) {
    files = c(files, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    transfer = c(transfer, 1, 1, 1)
    destination = c(destination, 0, 0, 0)
  }
  if (params$partyName == "T") {
    WriteTrackingTableCSV(params)
    files = c(files, "dl_track_tbl.csv")
    transfer    = c(transfer, 10)
    destination = c(destination, 10)
  }
  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File created")
  if (waitForTurn) {
    params = StoreStampsEntry(params, "R program execution delayed", "Tracking Table")
    WaitForTurn.3p(params, sleepTime)
    params = StoreStampsEntry(params, "R program execution restarted", "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (job_started) {
    MakeTrigger("job_started.ok", params$writePath)
  } else {
    MakeTrigger("files_done.ok", params$writePath)
  }
  if (length(from) == 1) {
  	if (from == "T") {
  		cat("Waiting for analysis center\n")
  	} else if (from == "A") {
  		cat("Waiting for data partner 1\n")
  	} else {
  		cat("Waiting for data partner 2\n")
  	}
  } else if (length(from) == 2) {
  	cat("Waiting for data partners\n")
  }
  tryCatch(Standby("files_done.ok", params$readPath[from], maxWaitingTime = maxWaitingTime))
  cat("Resuming local processing\n\n")
  DeleteTrigger("files_done.ok", params$readPath[from])
  params = MergeLogRaw.3p(params, from)
  params = UpdateCounters.3p(params)
  params = NewLogEntry.3p(params)
  params = MergeStampsRaw.3p(params, from)
  params = StoreStampsEntry(params, "R program execution begins", "Tracking Table")
  params = MergeTrackingTableRAW.3p(params, from)
  return(params)
}


PauseContinue.3p = function(params, from = NULL, maxWaitingTime = 24 * 60 * 60) {
  params = StoreLogEntry.3p(params, "")
  params = StoreTrackingTableEntry.3p(params)
  WriteLogCSV(params)
  if (length(from) == 1) {
  	if (from == "T") {
  		cat("Waiting for analysis center\n")
  	} else if (from == "A") {
  		cat("Waiting for data partner 1\n")
  	} else {
  		cat("Waiting for data partner 2\n")
  	}
  } else if (length(from) == 2) {
  	cat("Waiting for data partners\n")
  }
  tryCatch(Standby("files_done.ok", params$readPath[from], maxWaitingTime = maxWaitingTime))
  cat("Resuming local processing\n\n")
  DeleteTrigger("files_done.ok", params$readPath[from])
  params = MergeLogRaw.3p(params, from)
  params = UpdateCounters.3p(params)
  params = NewLogEntry.3p(params)
  params = MergeStampsRaw.3p(params, from)
  params = MergeTrackingTableRAW.3p(params, from)
  WriteLogCSV(params)
  return(params)
}


UpdateCounters.3p = function(params) {
  params$pmnStepCounter = max(params$log$history$Step) + 1
  return(params)
}

########################## K PARTY PMN COMM FUNCTIONS ##########################

WaitForTurn.kp = function(params, sleepTime) {
  Sys.sleep(sleepTime)

	if (!params$popmednet) return(NULL)

  cat("Waiting For Turn\n")
  startTime = proc.time()[3]
  cat("Elapsed Time:", HMS(0), "\r")

  if (exists("partyOffset")) {
    cat("\n\n")
    return()
  }

  partyOffset = 15

  modulus   = (params$numDataPartners + 1) * partyOffset
  targetTime = params$dataPartnerID * partyOffset

  cat("Elapsed Time:", HMS(0), "\r")
  while (as.integer(Sys.time()) %% modulus != targetTime) {
    elapsedTime = round(proc.time()[3] - startTime, 0)
    cat("Elapsed Time:", HMS(elapsedTime), "\r")
    Sys.sleep(0.1)
  }
  cat("\n\n")
}


SendPauseQuit.kp = function(params,
                         filesAC = NULL,
                         filesDP = NULL,
                         sleepTime = 10,
                         job_failed = FALSE,
                         waitForTurn = FALSE) {

  # Assumes that upon quitting, same thing is sent to everyone, so filesDP cannot be a list

  params$lastIteration = TRUE
  params$completed     = TRUE
  params = StoreLogEntry.kp(params, c(filesAC, filesDP))
  params = StoreTrackingTableEntry.kp(params)
  WriteLogCSV(params)
  WriteLogRaw(params)

  if (params$dataPartnerID != 0) {
    WriteTrackingTableRaw(params)

    filesAC = c(filesAC, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    if (!is.null(filesDP)) {
      filesDP = c(filesDP, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    }

    dataPartnerTarget = 1:params$numDataPartners
    dataPartnerTarget = dataPartnerTarget[-params$dataPartnerID]
    files = c(filesAC, rep(filesDP, length(dataPartnerTarget)), "file_list.csv")
    transfer = c(rep(1, length(files) - 1), 10)
    destination = c(rep(0, length(filesAC)),
                    rep(dataPartnerTarget, each = length(filesDP)),
                    10)
  }

  if (params$dataPartnerID == 0) {
    WriteTrackingTableCSV(params)
    files       = c("dl_track_tbl.csv", "file_list.csv")
    if (job_failed) {
      files = c(files, "job_fail.ok")
      params = StoreStampsEntry(params, "Job failed trigger file", "Trigger File created")
    } else {
      files = c(files, "job_done.ok")
      params = StoreStampsEntry(params, "Job done trigger file", "Trigger File created")
    }
    transfer = c(10, 10, 10)
    destination = c(10, 10, 10)
  }

  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File Created")
  params = StoreStampsEntry(params, "R program execution complete, output files written",
                            "Tracking Table")
  if (waitForTurn) {
    params = StoreStampsEntry(params, "R program execution delayed", "Tracking Table")
    WaitForTurn.kp(params, sleepTime)
    params = StoreStampsEntry(params, "R program execution restarted", "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (params$dataPartnerID == 0) {
    if (job_failed)  {
      MakeTrigger("job_fail.ok",  params$writePath)
    } else {
      MakeTrigger("job_done.ok",  params$writePath)
    }
  }
  MakeTrigger("files_done.ok", params$writePath)
  return(params)
}


SendPauseContinue.kp = function(params,
                             filesAC = NULL,
                             filesDP = NULL,
                             from   = NULL,
                             sleepTime = 10,
                             maxWaitingTime = 24 * 60 * 60,
                             job_started = FALSE,
                             waitForTurn = FALSE) {
  if (class(filesDP) != "list") {
    params = StoreLogEntry.kp(params, c(filesAC, filesDP))
    params = StoreTrackingTableEntry.kp(params)
    WriteLogCSV(params)
    WriteLogRaw(params)
    if (length(filesAC) + length(filesDP) > 0) {
      WriteTrackingTableRaw(params)
    }

    if (!is.null(filesAC)) {
      filesAC = c(filesAC, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    }
    if (!is.null(filesDP)) {
      filesDP = c(filesDP, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    }
    dataPartnerTarget = 1:params$numDataPartners
    if (params$dataPartnerID != 0) {
      dataPartnerTarget = dataPartnerTarget[-params$dataPartnerID]
    }

    files = c(filesAC, rep(filesDP, length(dataPartnerTarget)), "file_list.csv")
    transfer = c(rep(1, length(files) - 1), 10)
    destination = c(rep(0, length(filesAC)),
                    rep(dataPartnerTarget, each = length(filesDP)),
                    10)
  } else {
    files = filesAC
    for (dp in 1:params$numDataPartners) {
      files = c(files, filesDP[[dp]])
    }
    params = StoreLogEntry.kp(params, files)
    params = StoreTrackingTableEntry.kp(params)
    WriteLogCSV(params)
    WriteLogRaw(params)
    if (length(files) > 0) {
      WriteTrackingTableRaw(params)
    }
    if (!is.null(filesAC)) {
      filesAC = c(filesAC, "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
    }
    files = filesAC
    transfer = rep(1, length(files))
    destination = rep(0, length(filesAC))
    for (dp in 1:params$numDataPartners) {
      if (length(filesDP[[dp]]) > 0 && dp != params$dataPartnerID) {
        files = c(files, filesDP[[dp]], "stamps.rdata", "log.rdata", "tr_tb_updt.rdata")
        transfer = c(transfer, rep(1, length(filesDP[[dp]]) + 3))
        destination = c(destination, rep(dp, length(filesDP[[dp]]) + 3))
      }
    }
    files = c(files, "file_list.csv")
    transfer = c(transfer, 10)
    destination = c(destination, 10)
  }
  if (params$dataPartnerID == 0) {
    WriteTrackingTableCSV(params)
    files = c(files, "dl_track_tbl.csv")
    transfer    = c(transfer, 10)
    destination = c(destination, 10)
  }
  MakeCSV(files, transfer, destination, params$writePath)
  params = StoreStampsEntry(params, "Files done trigger file", "Trigger File created")
  if (waitForTurn) {
    params = StoreStampsEntry(params, "R program execution delayed", "Tracking Table")
    WaitForTurn.kp(params, sleepTime)
    params = StoreStampsEntry(params, "R program execution restarted", "Tracking Table")
  }
  WriteStampsCSV(params)
  WriteStampsRaw(params)
  if (job_started) {
    MakeTrigger("job_started.ok", params$writePath)
  } else {
    MakeTrigger("files_done.ok", params$writePath)
  }
  if (from == "AC") {
    cat("Waiting for analysis center\n")
    tryCatch(Standby("files_done.ok", params$readPathAC, maxWaitingTime = maxWaitingTime))
    DeleteTrigger("files_done.ok", params$readPathAC)
  } else if (from == "DP") {
    cat("Waiting for data partners\n")
    if (params$dataPartnerID == 0) {
      tryCatch(Standby("files_done.ok", params$readPathDP, maxWaitingTime = maxWaitingTime))
      DeleteTrigger("files_done.ok", params$readPathDP)
    } else {
      tryCatch(Standby("files_done.ok", params$readPathDP[-params$dataPartnerID], maxWaitingTime = maxWaitingTime))
      DeleteTrigger("files_done.ok", params$readPathDP[-params$dataPartnerID])
    }
  } else if (from == "DP1") {
    cat("Waiting for data partner 1\n")
    tryCatch(Standby("files_done.ok", params$readPathDP[1], maxWaitingTime = maxWaitingTime))
    DeleteTrigger("files_done.ok", params$readPathDP[1])
  } else if (from == "DP2") {
    cat("Waiting for data partner 2\n")
    tryCatch(Standby("files_done.ok", params$readPathDP[2], maxWaitingTime = maxWaitingTime))
    DeleteTrigger("files_done.ok", params$readPathDP[2])
  } else {
    cat("I don't know this command for \"from\"\n")
    Sys.sleep(10000)
  }
  cat("Resuming local processing\n\n")
  params = MergeLogRaw.kp(params, from)
  params = UpdateCounters.kp(params)
  params = NewLogEntry.kp(params)
  params = MergeStampsRaw.kp(params, from)
  params = StoreStampsEntry(params, "R program execution begins", "Tracking Table")
  params = MergeTrackingTableRAW.kp(params, from)
  return(params)
}


PauseContinue.kp = function(params, from = NULL, maxWaitingTime = 24 * 60 * 60) {
  params = StoreLogEntry.kp(params, "")
  params = StoreTrackingTableEntry.kp(params)
  WriteLogCSV(params)
  params = StoreStampsEntry(params, "R program execution paused", "Tracking Table")
  if (from == "AC") {
    cat("Waiting for analysis center\n")
    tryCatch(Standby("files_done.ok", params$readPathAC, maxWaitingTime = maxWaitingTime))
    DeleteTrigger("files_done.ok", params$readPathAC)
  } else {
    cat("Waiting for data partners\n")
    if (params$dataPartnerID == 0) {
      tryCatch(Standby("files_done.ok", params$readPathDP, maxWaitingTime = maxWaitingTime))
      DeleteTrigger("files_done.ok", params$readPathDP)
    } else {
      tryCatch(Standby("files_done.ok", params$readPathDP[-params$dataPartnerID], maxWaitingTime = maxWaitingTime))
      DeleteTrigger("files_done.ok", params$readPathDP[-params$dataPartnerID])
    }
  }
  cat("Resuming local processing\n\n")
  params = MergeLogRaw.kp(params, from)
  params = UpdateCounters.kp(params)
  params = NewLogEntry.kp(params)
  params = MergeStampsRaw.kp(params, from)
  params = StoreStampsEntry(params, "R program execution begins", "Tracking Table")
  params = MergeTrackingTableRAW.kp(params, from)
  WriteLogCSV(params)
  return(params)
}


UpdateCounters.kp = function(params) {
  params$pmnStepCounter = max(params$log$history$Step) + 1
  return(params)
}

ReceivedError.kp = function(params, from) {
	result = list()
	message = ""
	if (from == "AC") {
		messageExists = file.exists(file.path(params$readPathAC, "errorMessage.rdata"))
		if (messageExists) {
			message = ReadErrorMessage(params$readPathAC)
		}
	} else {
		messageExists = file.exists(file.path(params$readPathDP, "errorMessage.rdata"))
		for (id in 1:params$numDataPartners) {
			if (messageExists[id]) {
				message = paste0(message, ReadErrorMessage(params$readPathDP[id]), " ")
			}
		}
	}
	result$error = any(messageExists)
	result$message = message
	return(result)
}

################################ TIME FUNCTIONS ################################

GetUTCTime = function() {
  t = Sys.time()
  attr(t, "tzone") = "UTC"
  return(as.POSIXlt(t))
}


GetUTCOffset = function() {
  t = Sys.time()
  return(format(t, "%z"))
}


GetUTCOffsetSeconds = function() {
  t = Sys.time()
  offset = format(t, "%z")
  hour = as.numeric(substr(offset, 2, 3))
  min = as.numeric(substr(offset, 4, 5))
  pm  = ifelse(substr(offset, 1, 1) == "-", -1, 1)
  return(pm * (hour * 3600 + min * 60))
}


ConvertUTCtoRoundTripTime = function(t) {
  if (t$mon  < 9)  {  month = paste0("0", t$mon + 1)  } else { month = t$mon + 1}
  if (t$mday < 10) {  day   = paste0("0", t$mday)     } else { day   = t$mday }
  if (t$hour < 10) {  hour  = paste0("0", t$hour)     } else { hour  = t$hour }
  if (t$min  < 10) {  min   = paste0("0", t$min)      } else { min   = t$min }
  if (t$sec  < 10) {  sec   = paste0("0", t$sec)      } else { sec   = t$sec }
  t = paste0(t$year + 1900, "-", month, "-", day, " ", hour, ":",
             min, ":", sec)
}


GetRoundTripTime = function() {
  return(ConvertUTCtoRoundTripTime(GetUTCTime()))
}


GetElapsedTime = function(time1, final = FALSE, timeOnly = FALSE) {
  etime = floor(time1[3])
  hrs = floor(etime / 3600)
  mins = floor((etime %% 3600) / 60)
  secs = etime - hrs * 3600 - mins * 60

  hr1 = if (hrs > 9) toString(hrs) else paste0("0", toString(hrs))
  min1 = if (mins > 9) toString(mins) else paste0("0", toString(mins))
  sec1 = if (secs > 9) toString(secs) else paste0("0", toString(secs))
  if (final) {
    return(paste0("(Total time elapsed: ", hr1, "hr ", min1, "min ",
                  sec1, "sec)"))
  } else if (timeOnly) {
    return(paste0("(", hr1, "hr  ", min1, "min  ", sec1, "sec)"))
  }
  return(paste0("(Time elapsed: ", hr1, "hr ", min1, "min ", sec1, "sec)"))
}


ConvertSecsToHMS = function(secs, final = FALSE, timeOnly = FALSE) {
  if (length(secs) != 1) {
    secs = 0
  }
  secs = round(secs, digits = 0)
  hrs = floor(secs / 3600)
  mins = floor((secs %% 3600) / 60)
  secs = secs - hrs * 3600 - mins * 60

  hr1 = if (hrs > 9) toString(hrs) else paste0("0", toString(hrs))
  min1 = if (mins > 9) toString(mins) else paste0("0", toString(mins))
  sec1 = if (secs > 9) toString(secs) else paste0("0", toString(secs))
  if (final) {
    return(paste0("(Total time elapsed: ", hr1, ":", min1, ":", sec1, ")"))
  }
  if (timeOnly) {
    return(paste0("(", hr1, ":", min1, ":", sec1, ")"))
  }
  return(paste0("(Time elapsed: ", hr1, ":", min1, ":", sec1, ")"))
}

HMS = function(t) {
  paste(paste(formatC(t %/% (60*60), width = 2, format = "d", flag = "0"),
              formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
              formatC(t %% 60, width = 2, format = "d", flag = "0"),
              sep = ":"))
}


############################### BLOCK FUNCTIONS ################################

GetBlockSize = function(pA, pB) {
  # This minimium is set up based on our current encryption scheme
  # We need to guarentee that pA + pB + g <= blocksize
  # where g = (pA + 1) / (pA + pB + 1) * blocksize
  # May change in the future.
  minBlocksize = max(25, trunc(1 + (pA + pB + 1)^2 / pB))
  return(minBlocksize)
}


CreateBlocks = function(pA, pB, n, blocksize) {

  # Divides the matrix with ncol(n) into submatricies of approximately
  # equal size. Minmum size is blocksize.

  blocks = list()

  numBlocks       = max(trunc(n / blocksize), 1)
  newBlocksize    = trunc(n / numBlocks)
  numBigBlocks    = n %% newBlocksize
  numLittleBlocks = numBlocks - numBigBlocks
  bigBlocksize    = newBlocksize + 1
  gLittleBlock    = trunc(newBlocksize * (pA + 1) / (pA + pB + 1))
  gBigBlock       = trunc(bigBlocksize * (pA + 1) / (pA + pB + 1))

  blocks$numBlocks       = numBlocks
  blocks$littleBlocksize = newBlocksize
  blocks$bigBlocksize    = bigBlocksize
  blocks$numLittleBlocks = numLittleBlocks
  blocks$numBigBlocks    = numBigBlocks
  blocks$gLittleBlock    = gLittleBlock
  blocks$gBigBlock       = gBigBlock

  blocks$stops = integer()
  if (numBigBlocks > 0) {
    blocks$stops = bigBlocksize * 1:numBigBlocks
  }
  if (numLittleBlocks > 0) {
    blocks$stops = c(blocks$stops, bigBlocksize * numBigBlocks +
                       newBlocksize * 1:numLittleBlocks)
  }

  if (numBlocks == 1) {
    blocks$starts = c(1)
  } else {
    blocks$starts = c(1, 1 + blocks$stops)[1:numBlocks]
  }

  blocks$g = c(rep(gBigBlock, numBigBlocks), rep(gLittleBlock, numLittleBlocks))

  return(blocks)
}


CreateContainers = function(pA, pB, blocks) {
  containers = list()

  maximumFilesize = 25 * 1024^2

  numBlocks = blocks$numBlocks
  littleBlocksize = blocks$littleBlocksize

  littleBlockG = blocks$gLittleBlock

  littleFilesize.Z   = 8 * littleBlocksize * littleBlockG
  littleFilesize.W   = 8 * littleBlocksize * pB # used for W, V, RW, WR, RV, Cox
  littleFilesize.RZ  = 8 * littleBlocksize^2
  littleFilesize.PR  = 8 * (pA + 1) * pB        # I think this is not used anymore
  littleFilesize.XR = 8 * pA * pB

  numContainers.Z = ceiling(numBlocks * littleFilesize.Z / maximumFilesize)
  numBlocksSmallContainer.Z = trunc(numBlocks / numContainers.Z)
  numBlocksLargeContainer.Z = numBlocksSmallContainer.Z + 1
  numLargeContainer.Z = numBlocks %% numContainers.Z
  numSmallContainer.Z = numContainers.Z - numLargeContainer.Z

  numContainers.W = ceiling(numBlocks * littleFilesize.W / maximumFilesize)
  numBlocksSmallContainer.W = trunc(numBlocks / numContainers.W)
  numBlocksLargeContainer.W = numBlocksSmallContainer.W + 1
  numLargeContainer.W = numBlocks %% numContainers.W
  numSmallContainer.W = numContainers.W - numLargeContainer.W

  numContainers.RZ = ceiling(numBlocks * littleFilesize.RZ / maximumFilesize)
  numBlocksSmallContainer.RZ = trunc(numBlocks / numContainers.RZ)
  numBlocksLargeContainer.RZ = numBlocksSmallContainer.RZ + 1
  numLargeContainer.RZ = numBlocks %% numContainers.RZ
  numSmallContainer.RZ = numContainers.RZ - numLargeContainer.RZ

  numContainers.PR = ceiling(numBlocks * littleFilesize.PR / maximumFilesize)
  numBlocksSmallContainer.PR = trunc(numBlocks / numContainers.PR)
  numBlocksLargeContainer.PR = numBlocksSmallContainer.PR + 1
  numLargeContainer.PR = numBlocks %% numContainers.PR
  numSmallContainer.PR = numContainers.PR - numLargeContainer.PR

  numContainers.XR = ceiling(numBlocks * littleFilesize.XR / maximumFilesize)
  numBlocksSmallContainer.XR = trunc(numBlocks / numContainers.XR)
  numBlocksLargeContainer.XR = numBlocksSmallContainer.XR + 1
  numLargeContainer.XR = numBlocks %% numContainers.XR
  numSmallContainer.XR = numContainers.XR - numLargeContainer.XR

  if (numLargeContainer.Z > 0) {
    filebreak.Z = c(0:(numLargeContainer.Z - 1) * numBlocksLargeContainer.Z + 1,
                    0:(numSmallContainer.Z - 1) * numBlocksSmallContainer.Z + 1 +
                      numLargeContainer.Z * numBlocksLargeContainer.Z)
  } else {
    filebreak.Z = c(0:(numSmallContainer.Z - 1) * numBlocksSmallContainer.Z + 1 +
                      numLargeContainer.Z * numBlocksLargeContainer.Z)
  }

  if (numLargeContainer.W > 0) {
    filebreak.W = c(0:(numLargeContainer.W - 1) * numBlocksLargeContainer.W + 1,
                    0:(numSmallContainer.W - 1) * numBlocksSmallContainer.W + 1 +
                      numLargeContainer.W * numBlocksLargeContainer.W)
  } else {
    filebreak.W = c(0:(numSmallContainer.W - 1) * numBlocksSmallContainer.W + 1 +
                      numLargeContainer.W * numBlocksLargeContainer.W)
  }

  if (numLargeContainer.RZ > 0) {
    filebreak.RZ = c(0:(numLargeContainer.RZ - 1) * numBlocksLargeContainer.RZ + 1,
                     0:(numSmallContainer.RZ - 1) * numBlocksSmallContainer.RZ + 1 +
                       numLargeContainer.RZ * numBlocksLargeContainer.RZ)
  } else {
    filebreak.RZ = c(0:(numSmallContainer.RZ - 1) * numBlocksSmallContainer.RZ + 1 +
                       numLargeContainer.RZ * numBlocksLargeContainer.RZ)
  }

  if (numLargeContainer.PR > 0) {
    filebreak.PR = c(0:(numLargeContainer.PR - 1) * numBlocksLargeContainer.PR + 1,
                     0:(numSmallContainer.PR - 1) * numBlocksSmallContainer.PR + 1 +
                       numLargeContainer.PR * numBlocksLargeContainer.PR)
  } else {
    filebreak.PR = c(0:(numSmallContainer.PR - 1) * numBlocksSmallContainer.PR + 1 +
                       numLargeContainer.PR * numBlocksLargeContainer.PR)
  }

  if (numLargeContainer.XR > 0) {
    filebreak.XR = c(0:(numLargeContainer.XR - 1) * numBlocksLargeContainer.XR + 1,
                     0:(numSmallContainer.XR - 1) * numBlocksSmallContainer.XR + 1 +
                       numLargeContainer.XR * numBlocksLargeContainer.XR)
  } else {
    filebreak.XR = c(0:(numSmallContainer.XR - 1) * numBlocksSmallContainer.XR + 1 +
                       numLargeContainer.XR * numBlocksLargeContainer.XR)
  }

  containers$filebreak.Z   = filebreak.Z
  containers$filebreak.W   = filebreak.W
  containers$filebreak.RZ  = filebreak.RZ
  containers$filebreak.PR  = filebreak.PR # I think we are not using this anymore
  containers$filebreak.V   = filebreak.W
  containers$filebreak.RW  = filebreak.W
  containers$filebreak.WR  = filebreak.W
  containers$filebreak.RV  = filebreak.W
  containers$filebreak.VR  = filebreak.W
  containers$filebreak.Cox = filebreak.W
  containers$filebreak.XR = filebreak.XR

  return(containers)
}

########################### OUTPUT FORMAT FUNCTIONS ############################

formatPValue = function(pvals, width = 7) {
	p = c()
	for (x in pvals) {
		if (is.na(x)) {
			x = format("NA", width = width, justify = "right")
		} else if (x > 1e-3) {
			x = format(round(x, 5),  width = width, justify = "right", nsmall = 5)
		} else if (x > 2e-16) {
			x = formatC(x, format = "e", digits = 1)
		} else {
			x = format("<2e-16", width = width, justify = "right")
		}
		p = c(p, x)
	}
	return(p)
}


formatStrings = function(x, minWidth = NULL, justify = "left") {
	width = max(max(nchar(x)), minWidth)
	x = format(x, justify = justify, width = width)
	return(x)
}


formatStat = function(x) {
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


formatStatList = function(vals) {
	# Assumes that x is non-empty set of numeric or NA and there are no NaN's
	# width = 10, justify = right => standard output, so no worries about justify nor width
	notNA = which(!is.na(vals))
	notZero = which(vals != 0)
	keep = intersect(notNA, notZero)
	if (length(keep) == 0) {
		f = c()
		for (x in vals) {
			if (is.na(x)) {
				f = c(f, "NA")
			} else {
				f = c(f, "0")
			}
		}
		return(f)
	}
	temp = vals[keep]  # All non-zero, non-NA
	minval = min(abs(temp))
	maxval = max(abs(temp))
	#Where most significant digit is located
	decmin = floor(log10(minval))
	decmax = floor(log10(maxval))
	if (minval >= 1) decmin = decmin + 1
	if (maxval >= 1) decmax = decmax + 1
	if ((decmin < -3) || (decmax > 6) || (decmax - decmin > 3)) {
		# scientific
		f = c()
		for (x in vals) {
			if (is.na(x)) {
				f = c(f, "NA")
			} else {
				f = c(f, formatC(x, format = "e", digits = 3))
			}
		}
		return(f)
	} else {
		# standard
		if (decmin < 0) {
			nsmall = 6
		} else {
			nsmall = 6 - decmin
		}
		f = c()
		for (x in vals) {
			if (is.na(x)) {
				f = c(f, "NA")
			} else {
				f = c(f, format(round(x, nsmall), scientific = FALSE, nsmall = nsmall))
			}
		}
		return(f)
	}
}

########################### SHARED STAMPS FUNCTIONS ############################

StoreStampsEntry = function(params, description = "", type = "") {
	newEntry             = params$stamps$blank
	newEntry$Step        = params$pmnStepCounter
	newEntry$Description = description
	newEntry$Time        = GetRoundTripTime()
	newEntry$Type        = type
	params$stamps$history = rbind(params$stamps$history, newEntry)
	return(params)
}


WriteStampsRaw = function(params) {
	stamps = params$stamps$history
	save(stamps, file = file.path(params$writePath, "stamps.rdata"))
}


WriteStampsCSV = function(params) {
	write.csv(params$stamps$history, file.path(params$writePath, "stamps.csv"),
						row.names = FALSE)
}

########################### 2 PARTY STAMPS FUNCTIONS ###########################

InitializeStamps.2p = function(params) {
	stamps = list()
	stamps$blank = data.frame(#Iteration   = params$pmnIterationCounter,
														Step        = params$pmnStepCounter,
														Source      = paste("Org", params$partyName, "Dist Reg"),
														Description = "R program execution begins",
														Time        = GetRoundTripTime(),
														Type        = "Tracking Table")
	stamps$history = stamps$blank
	params$stamps = stamps
	return(params)
}


ReadStampsRaw.2p = function(params) {
	load(file.path(params$readPath, "stamps.rdata"))
	params$stamps$history = stamps
	return(params)
}


MergeStampsRaw.2p = function(params) {
	# This function will only be used in the function Pause Continue
	# When party A and party B run simultaneously, but Party A can run first
	# even if Party B starts the whole thing.  We append party B's log
	# to the end of Party A's log.
	load(file.path(params$readPath, "stamps.rdata"))
	params$stamps$history = rbind(params$stamps$history, stamps)
	return(params)
}

########################### 3 PARTY STAMPS FUNCTIONS ###########################

InitializeStamps.3p = function(params) {
	stamps = list()
	stamps$blank = data.frame(#Iteration   = params$pmnIterationCounter,
														Step        = params$pmnStepCounter,
														Source      = paste("Org", params$partyName, "Dist Reg"),
														Description = "R program execution begins",
														Time        = GetRoundTripTime(),
														Type        = "Tracking Table")
	stamps$history = stamps$blank
	params$stamps = stamps
	return(params)
}


MergeStampsRaw.3p = function(params, from) {
	for (party in from) {
		load(file.path(params$readPath[[party]], "stamps.rdata"))
		key1 = paste0(params$stamps$history$Step,
									params$stamps$history$Source,
									params$stamps$history$Description)
		key2 = paste0(stamps$Step,
									stamps$Source,
									stamps$Description)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$stamps$history = rbind(params$stamps$history, stamps)
		} else if (length(idx) < length(key2)) {
			params$stamps$history = rbind(params$stamps$history, stamps[-idx, ])
		}
	}
	idx = order(as.character(params$stamps$history$Time))
	params$stamps$history = params$stamps$history[idx, ]
	return(params)
}

########################### K PARTY STAMPS FUNCTIONS ###########################

InitializeStamps.kp = function(params) {
	stamps = list()
	stamps$blank = data.frame(Step        = params$pmnStepCounter,
														Source      = paste0("Org dp", params$dataPartnerID, " Dist Reg"),
														Description = "R program execution begins",
														Time        = GetRoundTripTime(),
														Type        = "Tracking Table")
	stamps$history = stamps$blank
	params$stamps = stamps
	return(params)
}


MergeStampsRaw.kp = function(params, from) {
	if (from == "AC") {
		load(file.path(params$readPathAC, "stamps.rdata"))
		key1 = paste0(params$stamps$history$Step,
									params$stamps$history$Source,
									params$stamps$history$Description)
		key2 = paste0(stamps$Step,
									stamps$Source,
									stamps$Description)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$stamps$history = rbind(params$stamps$history, stamps)
		} else if (length(idx) < length(key2)) {
			params$stamps$history = rbind(params$stamps$history, stamps[-idx, ])
		}
	} else if (from == "DP1") {
		load(file.path(params$readPathDP[1], "stamps.rdata"))
		key1 = paste0(params$stamps$history$Step,
									params$stamps$history$Source,
									params$stamps$history$Description)
		key2 = paste0(stamps$Step,
									stamps$Source,
									stamps$Description)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$stamps$history = rbind(params$stamps$history, stamps)
		} else if (length(idx) < length(key2)) {
			params$stamps$history = rbind(params$stamps$history, stamps[-idx, ])
		}
	} else {
		for (id in 1:params$numDataPartners) {
			if (id == params$dataPartnerID) next
			load(file.path(params$readPathDP[id], "stamps.rdata"))
			key1 = paste0(params$stamps$history$Step,
										params$stamps$history$Source,
										params$stamps$history$Description)
			key2 = paste0(stamps$Step,
										stamps$Source,
										stamps$Description)
			idx = which(key2 %in% key1)
			if (length(idx) == 0) {
				params$stamps$history = rbind(params$stamps$history, stamps)
			} else if (length(idx) < length(key2)) {
				params$stamps$history = rbind(params$stamps$history, stamps[-idx, ])
			}
		}
	}
	idx = order(as.character(params$stamps$history$Time))
	params$stamps$history = params$stamps$history[idx, ]
	return(params)
}

############################# SHARED LOG FUNCTIONS #############################

AddToLog = function(params, functionName, readTime, readSize, writeTime, writeSize) {
	readTime  = round(as.numeric(readTime),  digits = 2)
	writeTime = round(as.numeric(writeTime), digits = 2)
	readSize  = round(as.numeric(readSize),  digits = 0)
	writeSize = round(as.numeric(writeSize), digits = 0)
	if (params$log$current$Functions == "") {
		params$log$current$Functions = functionName
	} else {
		params$log$current$Functions = paste0(params$log$current$Functions,
																					", ", functionName)
	}
	params$log$current$Read.Time = params$log$current$Read.Time + readTime
	params$log$current$Read.Size = params$log$current$Read.Size + readSize
	params$log$current$Write.Time = params$log$current$Write.Time + writeTime
	params$log$current$Write.Size = params$log$current$Write.Size + writeSize
	return(params)
}


WriteLogRaw = function(params) {
	log = params$log$history
	save(log, file = file.path(params$writePath, "log.rdata"))
}


WriteLogCSV = function(params) {
	write.csv(params$log$history, file.path(params$writePath, "log.csv"),
						row.names = FALSE)
}


WriteToLogSummary = function(c1 = "", c2 = "", c3 = "",
														 writePath = getwd(), append = TRUE) {
	if (is.numeric(c2)) {
		c2 = round(c2, 2)
	}
	write.table(data.frame(c1, c2, c3),
							file.path(writePath, "log_summary.csv"), sep = ",", col.names = FALSE,
							row.names = FALSE, append = append)
}

############################# 2 PARTY LOG FUNCTIONS ############################

InitializeLog.2p = function(params) {
	log = list()
	log$blank = data.frame(Step             = 0,
												 Iteration.alg    = 0,
												 Party            = "",
												 Functions        = "",
												 Wait.Time        = 0,
												 Start.Time       = GetUTCTime(),
												 End.Time         = GetUTCTime(),
												 Read.Time        = 0,
												 Read.Size        = 0,
												 Write.Time       = 0,
												 Write.Size       = 0,
												 Computation.Time = 0,
												 Files.Sent       = "",
												 Bytes.Sent       = 0)
	log$current = log$blank
	log$history = log$blank
	params$log  = log
	return(params)
}


NewLogEntry.2p = function(params) {
	params$log$current = params$log$blank
	params$log$current$Party         = params$partyName
	params$log$current$Start.Time    = GetUTCTime()
	return(params)
}


StoreLogEntry.2p = function(params, files) {
	params$log$current$Step          = params$pmnStepCounter
	params$log$current$Iteration.alg = params$iter
	params$log$current$Party = params$partyName
	params$log$current$End.Time = GetUTCTime()
	params$log$current$Computation.Time = round(as.numeric(difftime(
		params$log$current$End.Time, params$log$current$Start.Time, units = "secs")) -
			params$log$current$Read.Time - params$log$current$Write.Time, 2)
	params$log$current$Files.Sent = paste(files, collapse = ", ")
	params$log$current$Bytes.Sent = sum(file.size(file.path(params$writePath, files)))
	if (is.na(params$log$current$Bytes.Sent)) {
		params$log$current$Bytes.Sent = 0
	}
	nrows = nrow(params$log$history)
	if (nrows >= 2) {
		params$log$current$Wait.Time =
			round(as.numeric(difftime(
				params$log$current$Start.Time,
				max(params$log$history$End.Time[which(params$log$history$Party ==
																								params$log$current$Party)]),
				units = "secs")), 2)
	}
	if (params$log$history$Party[nrows] == "") {
		params$log$history = params$log$current
	} else {
		params$log$history = rbind(params$log$history, params$log$current)
	}
	nrows = nrow(params$log$history)
	return(params)
}

ReadLogRaw.2p = function(params) {
	load(file.path(params$readPath, "log.rdata"))
	params$log$history = log
	return(params)
}


MergeLogRaw.2p = function(params) {
	# This function will only be used in the function Pause Continue
	# When party A and party B run simultaneously, but Party A can run first
	# even if Party B starts the whole thing.  We append party B's log
	# to the end of Party A's log.
	load(file.path(params$readPath, "log.rdata"))
	params$log$history = rbind(params$log$history, log)
	return(params)
}


SummarizeLog.2p = function(params) {
	writePath = params$writePath
	log    = params$log$history
	indexA = which(log$Party == "A")
	indexB = which(log$Party == "B")
	Party.A.Start.Time = log$Start.Time[indexA[1]]
	Party.A.End.Time   = log$End.Time[indexA[length(indexA)]]
	Party.A.Total.Time = round(as.numeric(difftime(
		Party.A.End.Time, Party.A.Start.Time, units = "secs")), digits = 2)
	Party.A.Reading.Time = sum(log$Read.Time[indexA])
	Party.A.Writing.Time = sum(log$Write.Time[indexA])
	Party.A.Computing.Time = sum(log$Computation.Time[indexA])
	Party.A.Waiting.Time = sum(log$Wait.Time[indexA])
	Party.A.Total.Time.HMS = ConvertSecsToHMS(Party.A.Total.Time, timeOnly = TRUE)
	Party.A.Reading.Time.HMS = ConvertSecsToHMS(Party.A.Reading.Time, timeOnly = TRUE)
	Party.A.Writing.Time.HMS = ConvertSecsToHMS(Party.A.Writing.Time, timeOnly = TRUE)
	Party.A.Computing.Time.HMS = ConvertSecsToHMS(Party.A.Computing.Time, timeOnly = TRUE)
	Party.A.Waiting.Time.HMS = ConvertSecsToHMS(Party.A.Waiting.Time, timeOnly = TRUE)
	Party.A.Bytes.Read = sum(log$Read.Size[indexA])
	Party.A.Bytes.Written = sum(log$Write.Size[indexA])

	Party.B.Start.Time = log$Start.Time[indexB[1]]
	Party.B.End.Time   = log$End.Time[indexB[length(indexB)]]
	Party.B.Total.Time = round(as.numeric(difftime(
		Party.B.End.Time, Party.B.Start.Time, units = "secs")), digits = 2)
	Party.B.Reading.Time = sum(log$Read.Time[indexB])
	Party.B.Writing.Time = sum(log$Wait.Time[indexB])
	Party.B.Computing.Time = sum(log$Computation.Time[indexB])
	Party.B.Waiting.Time = Party.B.Total.Time - Party.B.Reading.Time -
		Party.B.Writing.Time - Party.B.Computing.Time
	Party.B.Total.Time.HMS = ConvertSecsToHMS(Party.B.Total.Time, timeOnly = TRUE)
	Party.B.Reading.Time.HMS = ConvertSecsToHMS(Party.B.Reading.Time, timeOnly = TRUE)
	Party.B.Writing.Time.HMS = ConvertSecsToHMS(Party.B.Writing.Time, timeOnly = TRUE)
	Party.B.Computing.Time.HMS = ConvertSecsToHMS(Party.B.Computing.Time, timeOnly = TRUE)
	Party.B.Waiting.Time.HMS = ConvertSecsToHMS(Party.B.Waiting.Time, timeOnly = TRUE)
	Party.B.Bytes.Read = sum(log$Read.Size[indexB])
	Party.B.Bytes.Written = sum(log$Write.Size[indexB])

	Total.Transfer.Time = 0
	if (max(log$Step) > 1) {
		for (i in 2:max(log$Step)) {
			idx1 = which(log$Step == i - 1)
			idx2 = which(log$Step == i)
			Total.Transfer.Time = Total.Transfer.Time +
				as.numeric(difftime(min(log$Start.Time[idx2]),
														max(log$End.Time[idx1]), units = "secs"))
		}
	}
	Total.Transfer.Time = round(Total.Transfer.Time, 2)

	Total.Reading.Time = sum(log$Read.Time)
	Total.Writing.Time = sum(log$Write.Time)
	Total.Computing.Time = sum(log$Computation.Time)
	Elapsed.Computing.Time = Party.A.Total.Time - Total.Transfer.Time
	Total.Reading.Time.HMS = ConvertSecsToHMS(Total.Reading.Time, timeOnly = TRUE)
	Total.Writing.Time.HMS = ConvertSecsToHMS(Total.Writing.Time, timeOnly = TRUE)
	Total.Computing.Time.HMS = ConvertSecsToHMS(Total.Computing.Time, timeOnly = TRUE)
	Elapsed.Computing.Time.HMS = ConvertSecsToHMS(Elapsed.Computing.Time, timeOnly = TRUE)
	Total.Transfer.Time.HMS = ConvertSecsToHMS(Total.Transfer.Time, timeOnly = TRUE)
	Total.Bytes.Transfered = sum(log$Bytes.Sent)
	KB.Per.Second = round(Total.Bytes.Transfered / (Total.Transfer.Time * 1024), digits = 2)
	WriteToLogSummary(c1 = "Analysis", c2 = params$analysis, writePath = writePath, append = FALSE)
	if (!is.null(params$blocks)) {
		WriteToLogSummary(c1 = "Blocksize", c2 = params$blocks$littleBlocksize, writePath = writePath)
		WriteToLogSummary(c1 = "Number of Blocks",
											c2 = params$blocks$numLittleBlocks + params$blocks$numBigBlocks,
											writePath = writePath)
	}
	if (!is.null(params$n))   WriteToLogSummary(c1 = "N", c2 = params$n, writePath = writePath)

	p = max(0, params$p1.old - (params$analysis != "cox"))
	WriteToLogSummary(c1 = "p1", c2 = p, writePath = writePath)
	p = params$p2.old
	WriteToLogSummary(c1 = "p2", c2 = p, writePath = writePath)

	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Party A Start Time", c2 = Party.A.Start.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party A End Time", c2 = Party.A.End.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Run Time", c2 = Party.A.Total.Time,
										c3 = Party.A.Total.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Reading Time", c2 = Party.A.Reading.Time,
										c3 = Party.A.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Bytes Read", c2 = Party.A.Bytes.Read, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Writing Time", c2 = Party.A.Writing.Time,
										c3 = Party.A.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Bytes Written", c2 = Party.A.Bytes.Written, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Computing Time", c2 = Party.A.Computing.Time,
										c3 = Party.A.Computing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Waiting Time", c2 = Party.A.Waiting.Time,
										c3 = Party.A.Waiting.Time.HMS, writePath = writePath)
	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Party B Start Time", c2 = Party.B.Start.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party B End Time", c2 = Party.B.End.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Run Time", c2 = Party.B.Total.Time,
										c3 = Party.B.Total.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Reading Time", c2 = Party.B.Reading.Time,
										c3 = Party.B.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Bytes Read", c2 = Party.B.Bytes.Read, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Writing Time", c2 = Party.B.Writing.Time,
										c3 = Party.B.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Bytes Written", c2 = Party.B.Bytes.Written, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Computing Time", c2 = Party.B.Computing.Time,
										c3 = Party.B.Computing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Waiting Time", c2 = Party.B.Waiting.Time,
										c3 = Party.B.Waiting.Time.HMS, writePath = writePath)
	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Total Reading Time", c2 = Total.Reading.Time,
										c3 = Total.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Writing Time", c2 = Total.Writing.Time,
										c3 = Total.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Computing Time", c2 = Total.Computing.Time,
										c3 = Total.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Elapsed Computing Time", c2 = Elapsed.Computing.Time,
										c3 = Elapsed.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Total Transfer Time", c2 = Total.Transfer.Time,
										c3 = Total.Transfer.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Bytes Transfered", c2 = Total.Bytes.Transfered, writePath = writePath)
	WriteToLogSummary(c1 = "KB / Sec Transfer Rate", c2 = KB.Per.Second, writePath = writePath)

}

############################# 3 PARTY LOG FUNCTIONS ############################

InitializeLog.3p = function(params) {
	log = list()
	log$blank = data.frame(Step             = 0,
												 Iteration.alg    = 0,
												 Party            = "",
												 Functions        = "",
												 Wait.Time        = 0,
												 Start.Time       = GetUTCTime(),
												 End.Time         = GetUTCTime(),
												 Read.Time        = 0,
												 Read.Size        = 0,
												 Write.Time       = 0,
												 Write.Size       = 0,
												 Computation.Time = 0,
												 Files.Sent       = "",
												 Bytes.Sent       = 0)
	log$current = log$blank
	log$history = log$blank
	params$log  = log
	return(params)
}


NewLogEntry.3p = function(params) {
	params$log$current = params$log$blank
	params$log$current$Party         = params$partyName
	params$log$current$Start.Time    = GetUTCTime()
	return(params)
}


StoreLogEntry.3p = function(params, files) {
	params$log$current$Step          = params$pmnStepCounter
	params$log$current$Iteration.alg = params$algIterationCounter
	params$log$current$Party = params$partyName
	params$log$current$End.Time = GetUTCTime()
	params$log$current$Computation.Time = round(as.numeric(difftime(
		params$log$current$End.Time, params$log$current$Start.Time, units = "secs")) -
			params$log$current$Read.Time - params$log$current$Write.Time, 2)
	params$log$current$Files.Sent = paste(files, collapse = ", ")
	params$log$current$Bytes.Sent = sum(file.size(file.path(params$writePath, files)))
	if (is.na(params$log$current$Bytes.Sent)) {
		params$log$current$Bytes.Sent = 0
	}
	nrows = nrow(params$log$history)
	if (nrows >= 3) {
		params$log$current$Wait.Time =
			round(as.numeric(difftime(
				params$log$current$Start.Time,
				max(params$log$history$End.Time[which(params$log$history$Party ==
																								params$log$current$Party)]),
				units = "secs")), 2)
	}
	if (params$log$history$Party[nrows] == "") {
		params$log$history = params$log$current
	} else {
		params$log$history = rbind(params$log$history, params$log$current)
	}
	return(params)
}


MergeLogRaw.3p = function(params, from) {
	# This function will only be used in the function Pause Continue
	# When party A and party B run simultaneously, but Party A can run first
	# even if Party B starts the whole thing.  We append party B's log
	# to the end of Party A's log.
	for (party in from) {
		load(file.path(params$readPath[[party]], "log.rdata"))
		key1 = paste0(params$log$history$Step, params$log$history$Party)
		key2 = paste0(log$Step, log$Party)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$log$history = rbind(params$log$history, log)
		} else if (length(idx) < length(key2)) {
			params$log$history = rbind(params$log$history, log[-idx, ])
		}
	}
	idx = order(params$log$history$Step, params$log$history$Party)
	params$log$history = params$log$history[idx, ]

	return(params)
}


SummarizeLog.3p = function(params) {
	writePath = params$writePath

	log    = params$log$history
	indexA = which(log$Party == "A")
	indexB = which(log$Party == "B")
	indexT = which(log$Party == "T")
	Party.A.Start.Time = log$Start.Time[indexA[1]]
	Party.A.End.Time   = log$End.Time[indexA[length(indexA)]]
	Party.A.Total.Time = round(as.numeric(difftime(
		Party.A.End.Time, Party.A.Start.Time, units = "secs")), digits = 2)
	Party.A.Reading.Time = sum(log$Read.Time[indexA])
	Party.A.Writing.Time = sum(log$Write.Time[indexA])
	Party.A.Computing.Time = sum(log$Computation.Time[indexA])
	Party.A.Waiting.Time = sum(log$Wait.Time[indexA])
	Party.A.Total.Time.HMS = ConvertSecsToHMS(Party.A.Total.Time, timeOnly = TRUE)
	Party.A.Reading.Time.HMS = ConvertSecsToHMS(Party.A.Reading.Time, timeOnly = TRUE)
	Party.A.Writing.Time.HMS = ConvertSecsToHMS(Party.A.Writing.Time, timeOnly = TRUE)
	Party.A.Computing.Time.HMS = ConvertSecsToHMS(Party.A.Computing.Time, timeOnly = TRUE)
	Party.A.Waiting.Time.HMS = ConvertSecsToHMS(Party.A.Waiting.Time, timeOnly = TRUE)
	Party.A.Bytes.Read = sum(log$Read.Size[indexA])
	Party.A.Bytes.Written = sum(log$Write.Size[indexA])

	Party.B.Start.Time = log$Start.Time[indexB[1]]
	Party.B.End.Time   = log$End.Time[indexB[length(indexB)]]
	Party.B.Total.Time = round(as.numeric(difftime(
		Party.B.End.Time, Party.B.Start.Time, units = "secs")), digits = 2)
	Party.B.Reading.Time = sum(log$Read.Time[indexB])
	Party.B.Writing.Time = sum(log$Write.Time[indexB])
	Party.B.Computing.Time = sum(log$Computation.Time[indexB])
	Party.B.Waiting.Time = sum(log$Wait.Time[indexB])
	Party.B.Total.Time.HMS = ConvertSecsToHMS(Party.B.Total.Time, timeOnly = TRUE)
	Party.B.Reading.Time.HMS = ConvertSecsToHMS(Party.B.Reading.Time, timeOnly = TRUE)
	Party.B.Writing.Time.HMS = ConvertSecsToHMS(Party.B.Writing.Time, timeOnly = TRUE)
	Party.B.Computing.Time.HMS = ConvertSecsToHMS(Party.B.Computing.Time, timeOnly = TRUE)
	Party.B.Waiting.Time.HMS = ConvertSecsToHMS(Party.B.Waiting.Time, timeOnly = TRUE)
	Party.B.Bytes.Read = sum(log$Read.Size[indexB])
	Party.B.Bytes.Written = sum(log$Write.Size[indexB])

	Party.T.Start.Time = log$Start.Time[indexT[1]]
	Party.T.End.Time   = log$End.Time[indexT[length(indexT)]]
	Party.T.Total.Time = round(as.numeric(difftime(
		Party.T.End.Time, Party.T.Start.Time, units = "secs")), digits = 2)
	Party.T.Reading.Time = sum(log$Read.Time[indexT])
	Party.T.Writing.Time = sum(log$Write.Time[indexT])
	Party.T.Computing.Time = sum(log$Computation.Time[indexT])
	Party.T.Waiting.Time = sum(log$Wait.Time[indexT])
	Party.T.Total.Time.HMS = ConvertSecsToHMS(Party.T.Total.Time, timeOnly = TRUE)
	Party.T.Reading.Time.HMS = ConvertSecsToHMS(Party.T.Reading.Time, timeOnly = TRUE)
	Party.T.Writing.Time.HMS = ConvertSecsToHMS(Party.T.Writing.Time, timeOnly = TRUE)
	Party.T.Computing.Time.HMS = ConvertSecsToHMS(Party.T.Computing.Time, timeOnly = TRUE)
	Party.T.Waiting.Time.HMS = ConvertSecsToHMS(Party.T.Waiting.Time, timeOnly = TRUE)
	Party.T.Bytes.Read = sum(log$Read.Size[indexT])
	Party.T.Bytes.Written = sum(log$Write.Size[indexT])

	Total.Transfer.Time = 0
	if (max(log$Step) > 1) {
		for (i in 2:max(log$Step)) {
			idx1 = which(log$Step == i - 1)
			idx2 = which(log$Step == i)
			Total.Transfer.Time = Total.Transfer.Time +
				as.numeric(difftime(min(log$Start.Time[idx2]),
														max(log$End.Time[idx1]), units = "secs"))
		}
	}
	Total.Transfer.Time = round(Total.Transfer.Time, 2)
	Elapsed.Computing.Time = Party.T.Total.Time - Total.Transfer.Time

	Total.Reading.Time = sum(log$Read.Time)
	Total.Writing.Time = sum(log$Write.Time)
	Total.Computing.Time = sum(log$Computation.Time)
	Total.Reading.Time.HMS = ConvertSecsToHMS(Total.Reading.Time, timeOnly = TRUE)
	Total.Writing.Time.HMS = ConvertSecsToHMS(Total.Writing.Time, timeOnly = TRUE)
	Total.Transfer.Time.HMS = ConvertSecsToHMS(Total.Transfer.Time, timeOnly = TRUE)
	Total.Computing.Time.HMS = ConvertSecsToHMS(Total.Computing.Time, timeOnly = TRUE)
	Elapsed.Computing.Time.HMS = ConvertSecsToHMS(Elapsed.Computing.Time, timeOnly = TRUE)
	Total.Bytes.Transfered = sum(log$Bytes.Sent)
	KB.Per.Second = round(Total.Bytes.Transfered / (Total.Transfer.Time * 1024), digits = 2)

	WriteToLogSummary(c1 = "Analysis", c2 = params$analysis, writePath = writePath, append = FALSE)
	if (!is.null(params$blocks)) {
		WriteToLogSummary(c1 = "Blocksize", c2 = params$blocks$littleBlocksize, writePath = writePath)
		WriteToLogSummary(c1 = "Number of Blocks",
											c2 = params$blocks$numLittleBlocks + params$blocks$numBigBlocks,
											writePath = writePath)
	}
	if (!is.null(params$n))   WriteToLogSummary(c1 = "N", c2 = params$n, writePath = writePath)

	p = max(0, params$p1.old - (params$analysis != "cox"))
	WriteToLogSummary(c1 = "pA", c2 = p, writePath = writePath)
	p = params$p2.old
	WriteToLogSummary(c1 = "pB", c2 = p, writePath = writePath)

	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Party A Start Time", c2 = Party.A.Start.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party A End Time", c2 = Party.A.End.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Run Time", c2 = Party.A.Total.Time,
										c3 = Party.A.Total.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Reading Time", c2 = Party.A.Reading.Time,
										c3 = Party.A.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Bytes Read", c2 = Party.A.Bytes.Read, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Writing Time", c2 = Party.A.Writing.Time,
										c3 = Party.A.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Bytes Written", c2 = Party.A.Bytes.Written, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Computing Time", c2 = Party.A.Computing.Time,
										c3 = Party.A.Computing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party A Total Waiting Time", c2 = Party.A.Waiting.Time,
										c3 = Party.A.Waiting.Time.HMS, writePath = writePath)
	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Party B Start Time", c2 = Party.B.Start.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party B End Time", c2 = Party.B.End.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Run Time", c2 = Party.B.Total.Time,
										c3 = Party.B.Total.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Reading Time", c2 = Party.B.Reading.Time,
										c3 = Party.B.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Bytes Read", c2 = Party.B.Bytes.Read, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Writing Time", c2 = Party.B.Writing.Time,
										c3 = Party.B.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Bytes Written", c2 = Party.B.Bytes.Written, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Computing Time", c2 = Party.B.Computing.Time,
										c3 = Party.B.Computing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party B Total Waiting Time", c2 = Party.B.Waiting.Time,
										c3 = Party.B.Waiting.Time.HMS, writePath = writePath)
	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Party T Start Time", c2 = Party.T.Start.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party T End Time", c2 = Party.T.End.Time, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Run Time", c2 = Party.T.Total.Time,
										c3 = Party.T.Total.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Reading Time", c2 = Party.T.Reading.Time,
										c3 = Party.T.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Bytes Read", c2 = Party.T.Bytes.Read, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Writing Time", c2 = Party.T.Writing.Time,
										c3 = Party.T.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Bytes Written", c2 = Party.T.Bytes.Written, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Computing Time", c2 = Party.T.Computing.Time,
										c3 = Party.T.Computing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Party T Total Waiting Time", c2 = Party.T.Waiting.Time,
										c3 = Party.T.Waiting.Time.HMS, writePath = writePath)
	WriteToLogSummary(writePath = writePath)
	WriteToLogSummary(c1 = "Total Reading Time", c2 = Total.Reading.Time,
										c3 = Total.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Writing Time", c2 = Total.Writing.Time,
										c3 = Total.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Computing Time", c2 = Total.Computing.Time,
										c3 = Total.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Elapsed Computing Time", c2 = Elapsed.Computing.Time,
										c3 = Elapsed.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Total Transfer Time", c2 = Total.Transfer.Time,
										c3 = Total.Transfer.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Bytes Transfered", c2 = Total.Bytes.Transfered, writePath = writePath)
	WriteToLogSummary(c1 = "KB / Sec Transfer Rate", c2 = KB.Per.Second, writePath = writePath)

}

############################# K PARTY LOG FUNCTIONS ############################

InitializeLog.kp = function(params) {
  log = list()
	log$blank = data.frame(Step             = 0,
												 Iteration.alg    = 0,
												 Party            = "",
												 Functions        = "",
												 Wait.Time        = 0,
												 Start.Time       = GetUTCTime(),
												 End.Time         = GetUTCTime(),
												 Read.Time        = 0,
												 Read.Size        = 0,
												 Write.Time       = 0,
												 Write.Size       = 0,
												 Computation.Time = 0,
												 Files.Sent       = "",
												 Bytes.Sent       = 0)
	log$current = log$blank
	log$history = log$blank
	params$log  = log
	return(params)
}


NewLogEntry.kp = function(params) {
	params$log$current = params$log$blank
	params$log$current$Party         = paste0("dp", params$dataPartnerID)
	params$log$current$Start.Time    = GetUTCTime()
	return(params)
}


StoreLogEntry.kp = function(params, files) {
	params$log$current$Step          = params$pmnStepCounter
	params$log$current$Iteration.alg = params$algIterationCounter
	params$log$current$Party = paste0("dp", params$dataPartnerID)
	params$log$current$End.Time = GetUTCTime()
	params$log$current$Computation.Time = round(as.numeric(difftime(
		params$log$current$End.Time, params$log$current$Start.Time, units = "secs")) -
			params$log$current$Read.Time - params$log$current$Write.Time, 2)
	params$log$current$Files.Sent = paste(files, collapse = ", ")
	params$log$current$Bytes.Sent = sum(file.size(file.path(params$writePath, files)))
	if (is.na(params$log$current$Bytes.Sent)) {
		params$log$current$Bytes.Sent = 0
	}
	nrows = nrow(params$log$history)
	if (nrows >= 3) {
		params$log$current$Wait.Time =
			round(as.numeric(difftime(
				params$log$current$Start.Time,
				max(params$log$history$End.Time[which(params$log$history$Party ==
																								params$log$current$Party)]),
				units = "secs")), 2)
	}
	if (params$log$history$Party[nrows] == "") {
		params$log$history = params$log$current
	} else {
		params$log$history = rbind(params$log$history, params$log$current)
	}
	return(params)
}


MergeLogRaw.kp = function(params, from) {
	# This function will only be used in the function Pause Continue
	# When party A and party B run simultaneously, but Party A can run first
	# even if Party B starts the whole thing.  We append party B's log
	# to the end of Party A's log.
	if (from == "AC") {
		load(file.path(params$readPathAC, "log.rdata"))
		key1 = paste0(params$log$history$Step, params$log$history$Party)
		key2 = paste0(log$Step, log$Party)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$log$history = rbind(params$log$history, log)
		} else if (length(idx) < length(key2)) {
			params$log$history = rbind(params$log$history, log[-idx, ])
		}
	} else if (from == "DP1") {
		load(file.path(params$readPathDP[1], "log.rdata"))
		key1 = paste0(params$log$history$Step, params$log$history$Party)
		key2 = paste0(log$Step, log$Party)
		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$log$history = rbind(params$log$history, log)
		} else if (length(idx) < length(key2)) {
			params$log$history = rbind(params$log$history, log[-idx, ])
		}
	} else {
		for (id in 1:params$numDataPartners) {
			if (id == params$dataPartnerID) next
			load(file.path(params$readPathDP[id], "log.rdata"))
			key1 = paste0(params$log$history$Step, params$log$history$Party)
			key2 = paste0(log$Step, log$Party)
			idx = which(key2 %in% key1)
			if (length(idx) == 0) {
				params$log$history = rbind(params$log$history, log)
			} else if (length(idx) < length(key2)) {
				params$log$history = rbind(params$log$history, log[-idx, ])
			}
		}
	}
	idx = order(params$log$history$Step, params$log$history$Party)
	params$log$history = params$log$history[idx, ]
	return(params)
}


SummarizeLog.kp = function(params) {
	writePath = params$writePath
	log       = params$log$history

	WriteToLogSummary(c1 = "Analysis", c2 = params$analysis, writePath = writePath, append = FALSE)
	if (!is.null(params$n))   WriteToLogSummary(c1 = "N", c2 = params$n, writePath = writePath)

	for (i in 1:params$numDataPartners) {
		if (is.null(params$pi))  {
			p = 0
		} else {
			p = params$pi[i] - (i == 1) * (2 - (params$analysis == "cox"))
		}
		WriteToLogSummary(c1 = paste0("p", i), c2 = p, writePath = writePath)
	}

	WriteToLogSummary(writePath = writePath)

	total.time.0 = 0
	for (party in 0:params$numDataPartners) {
		partyName = paste0("dp", party)
		index = which(log$Party == partyName)
		if (length(index) > 0) {
		  Start.Time = log$Start.Time[index[1]]
		  End.Time   = log$End.Time[index[length(index)]]
		  Total.Time = round(as.numeric(difftime(End.Time, Start.Time, units = "secs")), digits = 2)
		  if (party == 0) { total.time.0 = Total.Time }
		  Reading.Time = sum(log$Read.Time[index])
		  Writing.Time = sum(log$Write.Time[index])
		  Computing.Time = sum(log$Computation.Time[index])
		  Waiting.Time = sum(log$Wait.Time[index])
		  Total.Time.HMS = ConvertSecsToHMS(Total.Time, timeOnly = TRUE)
		  Reading.Time.HMS = ConvertSecsToHMS(Reading.Time, timeOnly = TRUE)
		  Writing.Time.HMS = ConvertSecsToHMS(Writing.Time, timeOnly = TRUE)
		  Computing.Time.HMS = ConvertSecsToHMS(Computing.Time, timeOnly = TRUE)
		  Waiting.Time.HMS = ConvertSecsToHMS(Waiting.Time, timeOnly = TRUE)
		  Bytes.Read = sum(log$Read.Size[index])
		  Bytes.Written = sum(log$Write.Size[index])
		  WriteToLogSummary(c1 = paste(partyName, "Start Time"), c2 = Start.Time, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "End Time"), c2 = End.Time, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Run Time"), c2 = Total.Time,
		                    c3 = Total.Time.HMS, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Reading Time"), c2 = Reading.Time,
		                    c3 = Reading.Time.HMS, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Bytes Read"), c2 = Bytes.Read, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Writing Time"), c2 = Writing.Time,
		                    c3 = Writing.Time.HMS, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Bytes Written"), c2 = Bytes.Written, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Computing Time"), c2 = Computing.Time,
		                    c3 = Computing.Time.HMS, writePath = writePath)
		  WriteToLogSummary(c1 = paste(partyName, "Total Waiting Time"), c2 = Waiting.Time,
		                    c3 = Waiting.Time.HMS, writePath = writePath)
		  WriteToLogSummary(writePath = writePath)
		}
	}

	Total.Transfer.Time = 0
	if (max(log$Step) > 1) {
		for (i in 2:max(log$Step)) {
			idx1 = which(log$Step == i - 1)
			idx2 = which(log$Step == i)
			Total.Transfer.Time = Total.Transfer.Time +
				as.numeric(difftime(min(log$Start.Time[idx2]),
														max(log$End.Time[idx1]), units = "secs"))
		}
	}
	Total.Transfer.Time = round(Total.Transfer.Time, 2)
	Elapsed.Computing.Time = total.time.0 - Total.Transfer.Time

	Total.Reading.Time = sum(log$Read.Time)
	Total.Writing.Time = sum(log$Write.Time)
	Total.Computing.Time = sum(log$Computation.Time)
	Total.Reading.Time.HMS = ConvertSecsToHMS(Total.Reading.Time, timeOnly = TRUE)
	Total.Writing.Time.HMS = ConvertSecsToHMS(Total.Writing.Time, timeOnly = TRUE)
	Total.Transfer.Time.HMS = ConvertSecsToHMS(Total.Transfer.Time, timeOnly = TRUE)
	Total.Computing.Time.HMS = ConvertSecsToHMS(Total.Computing.Time, timeOnly = TRUE)
	Elapsed.Computing.Time.HMS = ConvertSecsToHMS(Elapsed.Computing.Time, timeOnly = TRUE)
	Total.Bytes.Transfered = sum(log$Bytes.Sent)
	KB.Per.Second = round(Total.Bytes.Transfered / (Total.Transfer.Time * 1024), digits = 2)

	WriteToLogSummary(c1 = "Total Reading Time", c2 = Total.Reading.Time,
										c3 = Total.Reading.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Writing Time", c2 = Total.Writing.Time,
										c3 = Total.Writing.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Computing Time", c2 = Total.Computing.Time,
										c3 = Total.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Elapsed Computing Time", c2 = Elapsed.Computing.Time,
										c3 = Elapsed.Computing.Time.HMS,  writePath = writePath)
	WriteToLogSummary(c1 = "Total Transfer Time", c2 = Total.Transfer.Time,
										c3 = Total.Transfer.Time.HMS, writePath = writePath)
	WriteToLogSummary(c1 = "Total Bytes Transfered", c2 = Total.Bytes.Transfered, writePath = writePath)
	WriteToLogSummary(c1 = "KB / Sec Transfer Rate", c2 = KB.Per.Second, writePath = writePath)
}

####################### SHARED TRACKING TABLE FUNCTIONS ########################

WriteTrackingTableRaw = function(params) {
	trackingTable = params$trackingTable$history
	save(trackingTable, file = file.path(params$writePath, "tr_tb_updt.rdata"))
	return(params)
}


WriteTrackingTableCSV = function(params) {
	write.csv(params$trackingTable$history, file.path(params$writePath, "dl_track_tbl.csv"),
						row.names = FALSE)
	return(params)
}

####################### 2 PARTY TRACKING TABLE FUNCTIONS #######################

InitializeTrackingTable.2p = function(params) {
	trackingTable = list()
	trackingTable$current = data.frame(DP_CD              = ifelse(params$partyName == "A", 0, 1),
																		 MSREQID            = params$msreqid,
																		 RUNID              = "dl",
																		 ITER_NB            = 0,  # params$pmnIterationCounter
																		 STEP_NB            = 0, # ifelse(params$partyName == "A", 2, 1),
																		 START_DTM          = GetUTCTime(), # from log$Start.Time
																		 END_DTM            = GetUTCTime(), # from log$End.Time
																		 CURR_STEP_IN       = 0,
																		 STEP_RETURN_CD     = 0,
																		 STEP_RETURN_MSG    = "PASS", # copy errorMessage.rdata here if exists
																		 REG_CONV_IN        = 0,  # 1 = converge, 0 = no converge
																		 REG_CONV_MSG       = "", # Success or Failed when decided
																		 LAST_ITER_IN       = 0,  # 1 at last iteration, so right before quit
																		 LAST_RUNID_IN      = 0,
																		 UTC_OFFSET_DISPLAY = GetUTCOffset(),
																		 UTC_OFFSET_SEC     = GetUTCOffsetSeconds(),
																		 REGR_TYPE_CD       = params$analysis
	)
	trackingTable$history = NA
	params$trackingTable = trackingTable
	return(params)
}

StoreTrackingTableEntry.2p = function(params) {
	params$trackingTable$current$ITER_NB = params$pmnStepCounter
	params$trackingTable$current$START_DTM = params$log$current$Start.Time
	params$trackingTable$current$END_DTM = params$log$current$End.Time
	if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
		load(file.path(params$readPath, "errorMessage.rdata"))
		params$trackingTable$current$STEP_RETURN_MSG = message
	} else if (file.exists(file.path(params$writePath, "errorMessage.rdata"))) {
		load(file.path(params$writePath, "errorMessage.rdata"))
		params$trackingTable$current$STEP_RETURN_MSG = message
	}
	params$trackingTable$current$REG_CONV_IN = ifelse(params$completed, 1, 0)
	if (params$completed) {
		params$trackingTable$current$REG_CONV_MSG = ifelse(params$converged, "Success", "Failed")
	}
	params$trackingTable$current$LAST_ITER_IN = ifelse(params$lastIteration, 1, 0)

	if (params$partyName == "A") {
		if (class(params$trackingTable$history) == "data.frame") {
			params$trackingTable$history = rbind(params$trackingTable$history,
																					 params$trackingTable$current)
		} else {
			params$trackingTable$history = params$trackingTable$current
		}
		write.csv(params$trackingTable$history, file.path(params$writePath, "dl_track_tbl.csv"),
							row.names = FALSE)
	} else {
		trackingTableEntry = params$trackingTable$current
		save(trackingTableEntry, file = file.path(params$writePath, "tr_tb_updt.rdata"))
	}
	return(params)
}

ReadTrackingTableUpdate.2p = function(params) {
	load(file.path(params$readPath, "tr_tb_updt.rdata"))
	trackingTableEntry$MSREQID = params$msreqid
	if (class(params$trackingTable$history) == "data.frame") {
		params$trackingTable$history = rbind(params$trackingTable$history,
																				 trackingTableEntry)
	} else {
		params$trackingTable$history = trackingTableEntry
	}
	return(params)
}

####################### 3 PARTY TRACKING TABLE FUNCTIONS #######################

InitializeTrackingTable.3p = function(params) {
	trackingTable = list()
	trackingTable$current = data.frame(DP_CD              = ifelse(params$partyName == "T", 0,
																																 ifelse(params$partyName == "A", 1, 2)),
																		 MSREQID            = params$msreqid,
																		 RUNID              = "dl",
																		 ITER_NB            = 0,  # params$pmnIterationCounter
																		 STEP_NB            = 0,
																		 START_DTM          = GetUTCTime(), # from log$Start.Time
																		 END_DTM            = GetUTCTime(), # from log$End.Time
																		 CURR_STEP_IN       = 0,
																		 STEP_RETURN_CD     = 0,
																		 STEP_RETURN_MSG    = "PASS", # copy errorMessage.rdata here if exists
																		 REG_CONV_IN        = 0,  # 1 = converge, 0 = no converge
																		 REG_CONV_MSG       = "", # Success or Failed when decided
																		 LAST_ITER_IN       = 0,  # 1 at last iteration, so right before quit
																		 LAST_RUNID_IN      = 0,
																		 UTC_OFFSET_DISPLAY = GetUTCOffset(),
																		 UTC_OFFSET_SEC     = GetUTCOffsetSeconds(),
																		 REGR_TYPE_CD       = params$analysis
	)
	trackingTable$history = NA
	params$trackingTable = trackingTable
	return(params)
}

StoreTrackingTableEntry.3p = function(params) {
	params$trackingTable$current$ITER_NB = params$pmnStepCounter
	params$trackingTable$current$START_DTM = params$log$current$Start.Time
	params$trackingTable$current$END_DTM = params$log$current$End.Time
	if (file.exists(file.path(params$writePath, "errorMessage.rdata"))) {
		load(file.path(params$writePath, "errorMessage.rdata"))
		params$trackingTable$current$STEP_RETURN_MSG = message
	} else {
		msg = ""
		for (party in c("A", "B", "T")) {
			if (!is.na(params$readPath[[party]]) &&
					file.exists(file.path(params$readPath[[party]], "errorMessage.rdata"))) {
				load(file.path(params$readPath[[party]], "errorMessage.rdata"))
				msg = paste0(msg, message)
			}
		}
		params$trackingTable$current$STEP_RETURN_MSG = msg
	}
	params$trackingTable$current$REG_CONV_IN = ifelse(params$completed, 1, 0)
	if (params$completed) {
		params$trackingTable$current$REG_CONV_MSG = ifelse(params$converged, "Success", "Failed")
	}
	params$trackingTable$current$LAST_ITER_IN = ifelse(params$lastIteration, 1, 0)
	if (params$pmnStepCounter == 0) {
		params$trackingTable$history = params$trackingTable$current
	} else {
		params$trackingTable$history = rbind(params$trackingTable$history,
																				 params$trackingTable$current)
	}
	return(params)
}


MergeTrackingTableRAW.3p = function(params, from) {
	for (party in from) {
		load(file.path(params$readPath[[party]], "tr_tb_updt.rdata"))
		key1 = paste0(params$trackingTable$history$ITER_NB,
									params$trackingTable$history$DP_CD)
		key2 = paste0(trackingTable$ITER_NB,
									trackingTable$DP_CD)

		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable)
		} else if (length(idx) < length(key2)) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable[-idx, ])
		}
	}
	idx = order(params$trackingTable$history$START_DTM)
	params$trackingTable$history = params$trackingTable$history[idx, ]
	params$trackingTable$history$MSREQID = params$msreqid
	return(params)
}

####################### K PARTY TRACKING TABLE FUNCTIONS #######################

InitializeTrackingTable.kp = function(params) {
	trackingTable = list()
	trackingTable$current = data.frame(DP_CD              = params$dataPartnerID,
																		 MSREQID            = params$msreqid,
																		 RUNID              = "dl",
																		 ITER_NB            = 0,  # params$pmnIterationCounter
																		 STEP_NB            = 0,
																		 START_DTM          = GetUTCTime(), # from log$Start.Time
																		 END_DTM            = GetUTCTime(), # from log$End.Time
																		 CURR_STEP_IN       = 0,
																		 STEP_RETURN_CD     = 0,
																		 STEP_RETURN_MSG    = "PASS", # copy errorMessage.rdata here if exists
																		 REG_CONV_IN        = 0,  # 1 = converge, 0 = no converge
																		 REG_CONV_MSG       = "", # Success or Failed when decided
																		 LAST_ITER_IN       = 0,  # 1 at last iteration, so right before quit
																		 LAST_RUNID_IN      = 0,
																		 UTC_OFFSET_DISPLAY = GetUTCOffset(),
																		 UTC_OFFSET_SEC     = GetUTCOffsetSeconds(),
																		 REGR_TYPE_CD       = params$analysis
	)
	trackingTable$history = NA
	params$trackingTable = trackingTable
	return(params)
}

StoreTrackingTableEntry.kp = function(params) {
	params$trackingTable$current$ITER_NB = params$pmnStepCounter
	params$trackingTable$current$START_DTM = params$log$current$Start.Time
	params$trackingTable$current$END_DTM = params$log$current$End.Time
	if (file.exists(file.path(params$writePath, "errorMessage.rdata"))) {
		load(file.path(params$writePath, "errorMessage.rdata"))
		params$trackingTable$current$STEP_RETURN_MSG = message
	} else {
		msg = ""
		for (id in 1:params$numDataPartners) {
			if (!is.na(params$readPathDP[id]) &&
					file.exists(file.path(params$readPathDP[id], "errorMessage.rdata"))) {
				load(file.path(params$readPathDP[id], "errorMessage.rdata"))
				msg = paste0(msg, message)
			}
		}
		if (!is.na(params$readPathAC) &&
				file.exists(file.path(params$readPathAC, "errorMessage.rdata"))) {
			load(file.path(params$readPathAC, "errorMessage.rdata"))
			msg = paste0(msg, message)
		}
		params$trackingTable$current$STEP_RETURN_MSG = msg
	}
	params$trackingTable$current$REG_CONV_IN = ifelse(params$completed, 1, 0)
	if (params$completed) {
		params$trackingTable$current$REG_CONV_MSG = ifelse(params$converged, "Success", "Failed")
	}
	params$trackingTable$current$LAST_ITER_IN = ifelse(params$lastIteration, 1, 0)
	if (params$pmnStepCounter == 0) {
		params$trackingTable$history = params$trackingTable$current
	} else {
		params$trackingTable$history = rbind(params$trackingTable$history,
																				 params$trackingTable$current)
	}

	return(params)
}


MergeTrackingTableRAW.kp = function(params, from) {
	if (from == "AC") {
		load(file.path(params$readPathAC, "tr_tb_updt.rdata"))
		key1 = paste0(params$trackingTable$history$ITER_NB,
									params$trackingTable$history$DP_CD)
		key2 = paste0(trackingTable$ITER_NB,
									trackingTable$DP_CD)

		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable)
		} else if (length(idx) < length(key2)) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable[-idx, ])
		}
	} else if (from == "DP1") {
		load(file.path(params$readPathDP[1], "tr_tb_updt.rdata"))
		key1 = paste0(params$trackingTable$history$ITER_NB,
									params$trackingTable$history$DP_CD)
		key2 = paste0(trackingTable$ITER_NB,
									trackingTable$DP_CD)

		idx = which(key2 %in% key1)
		if (length(idx) == 0) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable)
		} else if (length(idx) < length(key2)) {
			params$trackingTable$history =
				rbind(params$trackingTable$history, trackingTable[-idx, ])
		}
	} else {
		for (id in 1:params$numDataPartners) {
			if (id == params$dataPartnerID) next
			load(file.path(params$readPathDP[id], "tr_tb_updt.rdata"))
			key1 = paste0(params$trackingTable$history$ITER_NB,
										params$trackingTable$history$DP_CD)
			key2 = paste0(trackingTable$ITER_NB,
										trackingTable$DP_CD)

			idx = which(key2 %in% key1)
			if (length(idx) == 0) {
				params$trackingTable$history =
					rbind(params$trackingTable$history, trackingTable)
			} else if (length(idx) < length(key2)) {
				params$trackingTable$history =
					rbind(params$trackingTable$history, trackingTable[-idx, ])
			}
		}
	}
	idx = order(params$trackingTable$history$START_DTM)
	params$trackingTable$history = params$trackingTable$history[idx, ]
	params$trackingTable$history$MSREQID = params$msreqid
	return(params)
}

########################### VALID FORMULA FUNCTIONS ############################

validFormula = function(expression) {
	# This function takes an expresion and checks that it is of the form var1 ~ var2 + var3 + ... varN
	# It does not check for constants.  Constants are ignored and treated if the are not there.
	# Dupliate variables are ignored.  That is, as in lm(), formulas of the form
	# var1 ~ var2 + var2 are equivalent to var1 ~ var2

	# Check to make sure this is a valid expression
	if (tryCatch({is.expression(expression); FALSE},
							 error = function(err) { TRUE })) {
		return(FALSE)
	}
	vars = all.vars(expression)
	names = all.names(expression)
	#Check to see if expresion only contains variables, ~, and +.  no other symbols allowed.
	res1 = all(names %in% c("~", "+", vars))
	#Check to see if expresion contains exactly one ~
	res2 = (sum(names %in% "~") == 1)
	#Check to see if expression is of the form "variable ~ ....."
	res3 = (names[1] == "~") & (names[2] %in% vars)
	#Check to see if the LHS variable does not occur on the RHS
	res4 = !(names[2] %in% names[3:length(names)])
	#check to see if the LHS variable is not .
	res5 = vars[1] != "."
	return(res1 & res2 & res3 & res4 & res5)
}

validFormula2 = function(expression) {
	# This function takes and expression and checks that it is of the form ~ var1 + var2 + ... + varN
	# Duplicate variables are ignored.  That is, ~ var1 + var1 is equivalent to ~ var1
	if (tryCatch({is.expression(expression); FALSE},
							 error = function(err) { TRUE })) {
		return(FALSE)
	}
	vars = all.vars(expression)
	names = all.names(expression)
	# Check to see if expresion only contains variables, ~, and +.  no other symbols allowed.
	res1 = all(names %in% c("~", "+", vars))
	# Check to see if expresion contains exactly one ~
	res2 = (sum(names %in% "~") == 1)
	# Check to see if expression has no LHS (should not)
	res3 = length(expression) == 2
	return(res1 && res2 && res3)
}

###################### SHARED SUMMARY AND PRINT FUNCTIONS ######################

print.vdralinear = function(x) {
  if (x$failed) {
    cat("Two party distributed linear regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }
  cat("Coefficients:\n")
  print(x$coefficients, digits = 4)
  return(invisible(NULL))
}


summary.vdralinear = function(x) {
  temp = list()
  class(temp)         = "summary.vdralinear"
  temp$failed         = x$failed
  if (x$failed) {
    return(temp)
  }
  temp$party          = x$party
  temp$coefficients   = x$coefficients
  temp$secoef         = x$secoef
  temp$tvals          = x$tvals
  temp$pvals          = x$pvals
  temp$rstderr        = x$rstderr
  temp$df2            = x$df2
  temp$rsquare        = x$rsquare
  temp$adjrsquare     = x$adjrsquare
  temp$Fstat          = x$Fstat
  temp$df1            = x$df1
  temp$Fpval          = x$Fpval
  return(temp)
}


print.summary.vdralinear = function(x, lion = FALSE) {
  if (x$failed) {
    cat("Two party distributed linear regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }

  x$stars          = sapply(x$pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  })

  temp = data.frame(formatStrings(names(x$party)),
                    formatStrings(x$party, minWidth = 5, justify = "centre"),
                    formatStatList(x$coefficients),
                    formatStatList(x$secoef),
                    formatStatList(x$tvals),
                    format.pval(x$pvals),
                    formatStrings(x$stars))
  colnames(temp) = c("", "Party", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
  if (lion) {
    temp = cbind(temp, GetLion(length(x$party)))
    colnames(temp)[8] = ""
  }
  print(temp, row.names = FALSE, right = TRUE)
  cat("---", "\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("Residual standard error: ", formatStat(x$rstderr), "on", x$df2, "degrees of freedom\n")
  cat("Multiple R-squared: ", formatStat(x$rsquare), 	", Adjusted R-squared: ", formatStat(x$adjrsquare), "\n")
  cat("F-statistic:", formatStat(x$Fstat), "on", x$df1, "and", x$df2, "DF, p-value:", format.pval(x$Fpval), "\n\n")
}


print.vdralogistic = function(x) {
  if (x$failed) {
    cat("Two party distributed logistic regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }
  if (!x$converged) {
    cat("Warning: Two party distributed logistic regression did not converge in",
        x$iter, "iterations. Reported statistics are approximate.\n\n")
  }

  cat("Coefficients:\n")
  print(x$coefficients, digits = 4)
  cat("\n")
  cat("Degrees of Freedom:", x$nulldev_df, "Total (i.e. Null); ", x$resdev_df, "Residual\n")
  cat("Null Deviance:    ", formatStat(x$nulldev), "\n")
  cat("Residual Deviance:", formatStat(x$resdev), "   AIC:", formatStat(x$aic))
  return(invisible(NULL))
}


summary.vdralogistic = function(x, lion = FALSE) {
  temp = list()
  class(temp)         = "summary.vdralogistic"
  temp$failed         = x$failed
  temp$converged      = x$converged
  if (x$failed) {
    return(temp)
  }
  temp$party          = x$party
  temp$coefficients   = x$coefficients
  temp$secoef         = x$secoef
  temp$tvals          = x$tvals
  temp$pvals          = x$pvals
  temp$nulldev        = x$nulldev
  temp$nulldev_df     = x$nulldev_df
  temp$resdev         = x$resdev
  temp$resdev_df      = x$resdev_df
  temp$aic            = x$aic
  temp$bic            = x$bic
  temp$iter           = x$iter
  return(temp)
}


print.summary.vdralogistic = function(x, lion = FALSE) {
  if (x$failed) {
    cat("Two party distributed logistic regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }
  if (!x$converged) {
    cat("Warning: Two party distributed logistic regression did not converge in",
        x$iter, "iterations. Reported statistics are approximate.\n\n")
  }
  x$stars = sapply(x$pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  })

  temp = data.frame(formatStrings(names(x$party)),
                    formatStrings(x$party, minWidth = 5, justify = "centre"),
                    formatStatList(x$coefficients),
                    formatStatList(x$secoef),
                    formatStatList(x$tvals),
                    format.pval(x$pvals),
                    formatStrings(x$stars))
  colnames(temp) = c("", "Party", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
  if (lion) {
    temp = cbind(temp, GetLion(length(x$party)))
    colnames(temp)[8] = ""
  }
  print(temp, row.names = FALSE, right = TRUE)
  cat("---", "\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("(Dispertion parameter for binomial family taken to be 1)\n\n")
  cat("    Null Deviance:", formatStat(x$nulldev), " on ", x$nulldev_df, " degrees of freedom\n")
  cat("Residual deviance:", formatStat(x$resdev), " on ", x$resdev_df, " degrees of freedom\n")
  cat("AIC:", formatStat(x$aic), "\n")
  cat("BIC:", formatStat(x$bic), "\n\n")
  cat("Number of Newton-Raphson iterations:", x$iter, "\n\n")
  return(invisible(NULL))
}


print.vdracox = function(x) {
  if (x$failed) {
    cat("Two party cox logistic regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }
  if (!x$converged) {
    cat("Warning: Two party distributed cox regression did not converge in",
        x$iter, "iterations. Reported statistics are approximate.\n\n")
  }

  coeftab = data.frame(x$coefficients, x$expcoef, x$secoef, x$zvals, x$pvals)
  colnames(coeftab) = c("coef", "exp(coef)", "se(coef)", "z", "p")
  printCoefmat(coeftab, P.values = TRUE, has.Pvalue=TRUE, signif.stars = FALSE)
  cat("\n")
  cat(paste0("Likelihood ratio test=", formatStat(x$lrt[1])), "on",
      x$df, paste0("df, p=", format.pval(x$lrt[2])), "\n")
  cat("n=", paste0(x$n, ","), "number of events=", x$nevent, "\n\n")
  return(invisible(NULL))
}


summary.vdracox = function(x) {
  temp = list()
  class(temp)         = "summary.vdracox"
  temp$failed         = x$failed
  temp$converged      = x$converged
  if (x$failed) {
    return(temp)
  }
  temp$party          = x$party
  temp$coefficients   = x$coefficients
  temp$expcoef        = x$expcoef
  temp$secoef         = x$secoef
  temp$zval           = x$zval
  temp$pvals          = x$pvals
  temp$expncoef       = x$expncoef
  temp$lower          = x$lower
  temp$upper          = x$upper
  temp$n              = x$n
  temp$nevent         = x$nevent
  temp$concordance    = x$concordance
  temp$rsquare        = x$rsquare
  temp$lrt            = x$lrt
  temp$df             = x$df
  temp$wald.test      = x$wald.test
  temp$score          = x$score
  temp$iter           = x$iter
  return(temp)
}


print.summary.vdracox = function(x, lion = FALSE) {
  if (x$failed) {
    cat("Two party cox logistic regression failed.  No results to print.\n\n")
    return(invisible(NULL))
  }
  if (!x$converged) {
    cat("Warning: Two party distributed cox regression did not converge in",
        x$iter, "iterations. Reported statistics are approximate.\n\n")
  }

  x$stars = sapply(x$pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  })

  temp1 = data.frame(formatStrings(names(x$party)),
                     formatStrings(x$party, minWidth = 5, justify = "centre"),
                     formatStatList(x$coefficients),
                     formatStatList(x$expcoef),
                     formatStatList(x$secoef),
                     formatStatList(x$zval),
                     format.pval(x$pvals),
                     formatStrings(x$stars))
  colnames(temp1) = c("", "party", "   coef", "exp(coef)", "se(coef)", "   z", "Pr(>|z|)", "")
  temp2 = data.frame(formatStrings(names(x$party)),
                     formatStrings(x$party, minWidth = 5, justify = "centre"),
                     formatStatList(x$expcoef),
                     formatStatList(x$expncoef),
                     formatStatList(x$lower),
                     formatStatList(x$upper))
  colnames(temp2) = c("", "party", "exp(coef)", "exp(-coef)", "lower .95", "upper .95")
  cat("  n=", paste0(x$n, ","), "number of events=", x$nevent, "\n\n")
  print(temp1, row.names = FALSE, right = TRUE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  print(temp2, row.names = FALSE, right = TRUE)
  cat("\n")
  if (!is.na(x$concordance[5])) {
    cat("Concordance=", formatStat(x$concordance[5]), "(se =",formatStat(x$concordance[6]), ")\n")
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


differentModel = function(formula = NULL, x = NULL) {
  if (class(x) != "vdralinear") {
    cat("Error: This function can only be on objects of class vdralinear. Returning original model.\n\n")
    return(invisible(x))
  }
  if (x$failed) {
    cat("Error: linear regression failed.  Cannot compute a different model.\n\n")
    return(invisible(x))
  }
  if (max(table(names(x$party))) > 1) {
    cat("Error: Duplicate variable names exist.  All variable names must be unique. Returning original model.\n\n")
    return(invisible(x))
  }
  if (!validFormula(formula)) {
    cat("Invalid formula, returning original model.\n\n")
    return(x)
  }

  valid_names    = c(colnames(x$xty), colnames(x$xtx)[-1])
  responseName   = all.vars(formula)[1]
  covariateNames = all.vars(formula)[-1]
  variableNames  = all.vars(formula)
  variableNames  = variableNames[which(variableNames != ".")]

  if (!all(variableNames %in% valid_names)) {
    vars = variableNames[which(variableNames %in% valid_names == FALSE)]
    if (length(vars) == 1) {
      cat("Variable", vars, "not found. Returning original model.\n\n")
    } else {
      temp = c(paste0(vars[-length(vars)], ","), vars[length(vars)])
      cat("Variables", temp, "not found. Returning original model.\n\n")
    }
    return(invisible(x))
  }

  if ("." %in% covariateNames) {
    covariateNames = valid_names[which(valid_names != responseName)]
  }

  xytxy = rbind(cbind(x$yty, t(x$xty)), cbind(x$xty, x$xtx))
  scramble = c(2, 1, 3:ncol(xytxy))
  xytxy[scramble, scramble] = xytxy

  all_names = c(colnames(x$xty), colnames(x$xtx))[scramble] # Put (intercept) first
  colnames(xytxy) = all_names
  rownames(xytxy) = all_names

  responseIndex = match(responseName, all_names)
  covariateIndex = c(1, match(covariateNames, all_names))

  xtx    = xytxy[covariateIndex, covariateIndex]
  xty    = matrix(xytxy[covariateIndex, responseIndex], ncol = 1)
  yty    = xytxy[responseIndex, responseIndex]
  means  = c(x$meansy, x$means)[scramble][covariateIndex]
  meansy = c(x$meansy, x$means)[scramble][responseIndex]

  nrow = nrow(xtx)
  indicies = c(1)
  for (i in 2:nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(xtx[tempIndicies, tempIndicies]) > 10 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  xtx.old   = xtx
  xty.old   = xty
  xtx       = xtx[indicies, indicies]
  xty       = matrix(xty[indicies, 1], ncol = 1)
  means.old = means
  means     = means[indicies]
  p         = length(indicies)
  n         = x$n

  invxtx = solve(xtx)
  betas  = drop(invxtx %*% xty)

  numCovariates = p - 1

  sse     = max(drop(yty - 2 * t(xty) %*% betas + (t(betas) %*% xtx) %*% betas), 0)
  rstderr = drop(sqrt(sse / (n - numCovariates - 1)))
  sst     = drop(yty - meansy^2 * n)
  ssr     = sst - sse
  df1     = numCovariates
  df2     = n - numCovariates - 1
  if (sse == 0) {
    Fstat = Inf
  } else {
    Fstat   = (ssr / df1) / (sse / df2)
  }
  Fpval   = pf(Fstat, df1, df2, lower.tail = FALSE)
  if (sse == 0) {
    Rsq = 1
  } else {
    Rsq     = drop(1 - sse / sst)
  }
  adjRsq  = drop(1 - (n - 1) / (n - numCovariates - 1) * (1 - Rsq))
  if (rstderr == 0) {
    tvals = rep(Inf, numCovariates + 1)
  } else {
    tvals   = betas / (rstderr * sqrt(diag(invxtx)))
  }
  secoef  = tvals^-1 * betas
  pvals   = 2 * pt(abs(tvals), n - numCovariates - 1, lower.tail = FALSE)
  stars   = matrix(sapply(pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  }))

  y = list()
  class(y) = "vdralinear"
  y$failed    = x$failed
  y$converged = x$converged

  y$party = c(x$responseParty, x$party)[scramble][covariateIndex]
  y$responseParty = c(x$responseParty, x$party)[scramble][responseIndex]
  p1 = length(covariateIndex)
  y$coefficients           = rep(NA, p1)
  y$tvals                  = rep(NA, p1)
  y$secoef                 = rep(NA, p1)
  y$pvals                  = rep(NA, p1)

  y$sse                    = sse
  y$coefficients[indicies] = betas
  y$tvals[indicies]        = tvals
  y$secoef[indicies]       = secoef
  y$pvals[indicies]        = pvals
  y$rstderr                = rstderr
  y$rsquare                = Rsq
  y$adjrsquare             = adjRsq
  y$Fstat                  = Fstat
  y$Fpval                  = Fpval
  y$df1                    = df1
  y$df2                    = df2
  y$n                      = x$n
  y$xtx                    = xtx.old
  y$xty                    = xty.old
  y$yty                    = yty
  y$meansy                 = meansy
  y$means                  = means.old

  names.old                = all_names[covariateIndex]
  names(y$party)           = names.old
  names(y$coefficients)    = names.old
  names(y$secoef)          = names.old
  names(y$tvals)           = names.old
  names(y$pvals)           = names.old
  colnames(y$xtx)          = names.old
  rownames(y$xtx)          = names.old
  colnames(y$xty)          = responseName
  rownames(y$xty)          = names.old
  return(invisible(y))
}

###################### LOGISTIC ROC AND HOSLEM FUNCTIONS #######################


HoslemInternal = function(x, data = NULL, nGroups = 10){
  #            y:  response (vector, length n)
  #  finalFitted:  finalFitted from getFinalCoefA(...)  (vector, length n)
  #            p:  number of covariates pA + pB
  #      nGroups:  number of groups, specified by user
  #                or chosen automatically if unspecified.
  #                Common to choose ngroups = 10, as long as nGroups > p + 1.
  #
  # Returns vector c(chisq, df, pval)

  n = x$n

  if (nGroups <= 0) {
    nGroups = 10
  }

  if (nGroups > n) {
    nGroups = n
  }

  if (is.null(data)) {
    Y = x$Y
  } else {
    Y = data$Y
  }
  pi_ = exp(x$FinalFitted) / (1 + exp(x$FinalFitted))
  uq = unique(quantile(pi_, probs = seq(0, 1, 1 / nGroups)))
  group_ = cut(pi_, breaks = uq, include.lowest = TRUE)
  dd = data.frame(y = Y[order(pi_)], pi_ = sort(pi_),
                  group = group_[order(pi_)])

  e1 = by(dd, dd$group, function(x) sum(x$pi_))
  o1 = by(dd, dd$group, function(x) sum(x$y))
  gn = table(dd$group)
  e0 = gn - e1
  o0 = gn - o1

  testStat = 0
  for (i in 1:length(e1)) {
    if (o0[i] == e0[i]) {
      temp1 = 0
    } else {
      temp1 = (o0[i] - e0[i])^2 / e0[i]
    }
    if (o1[i] == e1[i]) {
      temp2 = 0
    } else {
      temp2 = (o1[i] - e1[i])^2 / e1[i]
    }
    testStat = testStat + temp1 + temp2
  }

  df = nGroups - 2
  rtrn = c(testStat, df, 1 - pchisq(testStat, df))
  names(rtrn) = c("Chi-sq", "DF", "p-value")

  return(rtrn)
}


print.hoslemdistributed = function(x) {
  # if (!x$converged) {
  #   cat("Warning: Process did not converge.\n")
  #   cat("         Cannot perform Hosmer and Lemeshow goodness of fit test.\n")
  #   return(invisible(NULL))
  # }
  cat("Hosmer and Lemeshow goodness of fit (GOF) test\n",
      "       Chi-squared:", x$hoslem[1], "with DF",
      paste0(x$hoslem[2],","), " p-value:", x$hoslem[3], "\n")
}


HoslemTest = function(x = NULL, nGroups = 10) {
  if (class(x) != "vdralogistic") {
    cat("Warning: Cannot perform test on non vdralogistic object.\n")
    return(invisible(NULL))
  }
  if (!(x$converged)) {
    cat("Warning: Process did not converge.\n")
    cat("         Cannot perform Hosmer and Lemeshow goodness of fit test.\n")
    return(invisible(NULL))
  }
  if (is.null(x$Y) || is.null(x$FinalFitted)) {
    cat("HoslemTest can only be invoked by the party which holds the response.\n")
    return(invisible(NULL))
  } else if (is.numeric(nGroups)) {
    temp = list()
    class(temp) = "hoslemdistributed"
    temp$hoslem = HoslemInternal(x, nGroups = nGroups)
    return(temp)
  }
}


RocInternal = function(x, data = NULL, bins = 500){
  #             y:  response vector (numeric, not factor, length n)
  #   finalFitted:  final_fitted from getFinalCoefA(...)  (vector, length n)
  #    thresholds:  how smooth the curve should be
  #
  #  Returns myRocObject (object$auc to get AUC)
  #  Object size is roughly equal to a matrix with (thresholds)rows and 2 columns

  if (is.null(data)) {
    Y = x$Y
  } else {
    Y = data$Y
  }

  if (bins < 2) bins = 2

  positive = sum(Y)
  negative = length(Y) - positive
  pi_ = exp(x$FinalFitted) / (1 + exp(x$FinalFitted))
  threshold = seq(0, 1, length.out = bins)
  rtrn = matrix(NA, bins, 2)

  oldX = 1
  oldY = 1
  AUC = 0

  for (i in 1:bins) {
    newX = 1 - sum(Y == 0 & pi_ < threshold[i]) / negative
    newY = sum(Y & pi_ >= threshold[i]) / positive
    rtrn[i, 1] = newX
    rtrn[i, 2] = newY
    AUC = AUC + oldY * (oldX - newX)
    oldX = newX
    oldY = newY
  }

  temp = list()
  temp$roc = rtrn
  temp$auc = AUC
  return(temp)
}


print.rocdistributed = function(x) {
  # if (!x$converged) {
  #   cat("Warning: Process did not converge.  Cannot generate ROC.\n")
  #   return(invisible(NULL))
  # }
  rtrn = x$ROC$roc
  plot(rtrn[, 1], rtrn[, 2], xaxt = "n", yaxt = "n",
       xlim = c(-0.2, 1.2), ylim = c(0,1 ),
       type = "s", ylab = "Sensitivity", xlab = "1 - Specificity", col = "blue",
       main = "ROC Curve")
  axis(side = 1, at = seq(0,1, 0.2))
  axis(side = 2, at = seq(0,1, 0.2))
  lines(x = c(0, 1), y = c(0, 1), col = "limegreen")
  text(0.8, 0.05, paste("Area under the curve:",
                        format(x$ROC$auc, digits = 4)))
}


RocTest = function(x = NULL, bins = 10) {
  if (class(x) != "vdralogistic") {
    cat("Warning: Cannot create ROC on non vdralogistic object.\n")
    return(invisible(NULL))
  }
  if (!x$converged) {
    cat("Warning: Process did not converge.  Cannot generate ROC.\n")
    return(invisible(NULL))
  }
  if (is.null(x$Y) || is.null(x$FinalFitted)) {
    cat("RocTest can only be invoked by the party which holds the response.\n")
    return(invisible(NULL))
  } else if (is.numeric(bins)) {
    temp = list()
    class(temp) = "rocdistributed"
    temp$ROC = RocInternal(x, bins = bins)
    # temp$singularMatrix = x$singularMatrix
    return(temp)
  }
}


################### COX DISPLAY SURVFIT AND STRATA FUNCTIONS ###################


GetColors = function(n) {
  color = matrix(0, 6, 3)
  color[1, ] = c(0.000, 0.000, 1.000) # blue
  color[2, ] = c(0.627, 0.125, 0.941) # purple
  color[3, ] = c(1.000, 0.000, 0.000) # red
  color[4, ] = c(1.000, 0.647, 0.000) # orange
  color[5, ] = c(0.000, 1.000, 0.000) # green

  if (n == 1) {
    return(rgb(0, 0, 0))
  }
  if (n <= 3) {
    cols = c(rgb(color[1, 1], color[1, 2], color[1, 3]),
             rgb(color[2, 1], color[2, 2], color[2, 3]),
             rgb(color[3, 1], color[3, 2], color[3, 3]))
    return(cols[1:n])
  }
  cols = c()
  for (i in 1:n) {
    idx = 4 * (i - 1) / (n - 1) + 1 # If we add yellow back in, change 4 to 5
    idx1 = floor(idx)
    idx2 = ceiling(idx)
    dx   = idx - idx1
    tcol = color[idx1, ] + (color[idx2, ] - color[idx1, ]) * dx
    cols = c(cols, rgb(tcol[1], tcol[2], tcol[3]))
  }
  return(cols)
}


plot.survfitDistributed = function(surv, merge = TRUE, ...) {
  max = 0
  n = length(surv$strata)
  labels = c()
  max = max(surv$time)
  labels = names(surv$strata)
  arguments = list(...)
  arguments$x = 1
  arguments$type = "n"
  if (is.null(arguments$ylim)) arguments$ylim = c(0, 1)
  if (is.null(arguments$xlim)) arguments$xlim = c(0, max)
  if (is.null(arguments$xlab)) arguments$xlab = "Time"
  if (is.null(arguments$ylab)) arguments$ylab = "Percent Survival"
  if (is.null(arguments$main)) arguments$main = "Survival Curve"

  if (merge) {
    do.call("plot", arguments)
    cols = GetColors(n)
    start = 1
    for (i in 1:n) {
      end = start + surv$strata[i] - 1
      lines(c(1, surv$surv[start:end]) ~ c(0, surv$time[start:end]), type = "s", col = cols[i])
      start = end + 1
    }
    legend("bottomleft", legend = labels, col = cols, lty = 1)
  } else {
    cols = GetColors(1)
    start = 1
    for (i in 1:n) {
      end = start + surv$strata[i] - 1
      do.call("plot", arguments)
      lines(c(1, surv$surv[start:end]) ~ c(0, surv$time[start:end]), type = "s", col = cols)
      legend("bottomleft", legend = labels[i], col = cols, lty = 1)
      start = end + 1
    }
  }
}


print.survfitDistributed = function(x) {
  start = 1
  events = integer(length(x$strata))
  for (i in 1:length(x$strata)) {
    end = start + x$strata[i] - 1
    events[i] = sum(x$n.event[start:end])
    start = end + 1
  }
  df = data.frame(n = x$n, events = events)
  row.names(df) = names(x$strata)
  print(df)
}


survfitDistributed.stats = function(x) {
  surv          = list()
  surv$n        = x$strata$end - x$strata$start + 1
  for (i in 1:nrow(x$strata)) {
    start = x$strata$start[i]
    end   = x$strata$end[i]
    idx   = which(c(1, diff(x$survival$rank[start:end])) != 0)
    temp0 = table(x$survival$rank[start:end], x$survival$status[start:end])
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 = cbind(temp0, 0)
        colnames(temp0) = c("0", "1")
      } else {
        temp0 = cbind(0, temp0)
        colnames(temp0) = c("0", "1")
      }
    }

    surv$time     = c(surv$time, x$survival$rank[start:end][idx])
    surv$n.risk   = c(surv$n.risk, rev(cumsum(rev(temp0[, 1] + temp0[, 2]))))
    surv$n.event  = c(surv$n.event, temp0[, 2])
    surv$n.censor = c(surv$n.censor, temp0[, 1])
    surv$strata   = c(surv$strata, length(idx))
    surv$surv     = c(surv$surv, x$survival$surv[start:end][idx])
  }
  names(surv$n.risk) = NULL
  names(surv$n.event) = NULL
  names(surv$n.censor) = NULL
  names(surv$strata) = x$strata$label
  surv$type     = "right"
  class(surv) = "survfitDistributed"
  return(invisible(surv))
}


survfitDistributed.formula = function(x, formula, data) {
  surv = list()
  vars = all.vars(formula)
  if ("." %in% vars) {
    cat("This function does not allow the . symbol in formulas.\n")
    return(invisible(NULL))
  }
  if (!all(vars %in% colnames(data))) {
    cat("Not all strata are found in the data.\n")
    return(invisible(NULL))
  }
  if (length(vars) == 0) {
    data = data.frame(const__ = rep(1, length(x$survival$rank)))
  } else {
    idx    = which(colnames(data) %in% vars)
    data   = data[x$survival$sorted, idx, drop = FALSE]
  }
  sorted = do.call("order", as.data.frame(cbind(data, x$survival$rank, x$survival$status)))
  data   = data[sorted, , drop = FALSE]
  rank   = x$survival$rank[sorted]
  status = x$survival$status[sorted]
  data2  = matrix(0, nrow = nrow(data), ncol = ncol(data))
  legend = list()
  colnames(data2) = colnames(data)
  for (i in 1:ncol(data)) {
    levels = levels(as.factor(data[, i]))
    legend[[colnames(data)[i]]] = levels
    data2[, i] = sapply(data[, i], function(x) { which(levels %in% x)})
  }
  ranks = which(apply(abs(apply(data2, 2, diff)), 1, sum) > 0)
  ranks = c(ranks, nrow(data2))
  start = 1
  for (i in 1:length(ranks)) {
    end = ranks[i]
    surv$n = c(surv$n, end - start + 1)
    # Calculate the Kaplan Meier Curve Here per notes from 9/4/19
    rank2 = rank[start:end]
    event2 = status[start:end]
    temp = table(rank2)
    M = length(temp)
    temp0 = table(rank2, event2)
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 = cbind(temp0, 0)
        colnames(temp0) = c("0", "1")
      } else {
        temp0 = cbind(0, temp0)
        colnames(temp0) = c("0", "1")
      }
    }
    idx = which(temp0[, 2] > 0)
    if (temp0[nrow(temp0), 2] == 0) idx = c(idx, nrow(temp0))
    nfails = temp0[idx, 2]
    start0 = c(1, (cumsum(temp)[1:(M - 1)] + 1))[idx]
    start1 = start0 + temp0[idx, 1]
    stop1  = start1 + nfails  - 1
    final = length(rank2)
    S = 1
    t2 = rep(0, length(nfails))
    S2 = rep(0, length(nfails))
    for (j in 1:length(nfails)) {
      n = final - start0[j] + 1
      d = stop1[j] - start1[j] + 1
      S = S * (n - d) / n
      t2[j] = rank2[start0[j]]
      S2[j] = S
    }
    surv$time = c(surv$time, t2)
    surv$n.risk   = c(surv$n.risk, rev(cumsum(rev(temp0[, 1] + temp0[, 2]))))
    surv$n.event = c(surv$n.event, temp0[, 2])
    surv$n.censor = c(surv$n.censor, temp0[, 1])
    surv$surv = c(surv$surv, S2)
    surv$strata = c(surv$strata, length(idx))
    if (length(vars) == 0) {
      names(surv$strata)[i] = ""
    } else {
      label = ""
      for (j in 1:ncol(data)) {
        temp = colnames(data)[j]
        label = paste0(label, temp, "=", legend[[temp]][data2[start, j]])
        if (j < ncol(data)) {
          label = paste0(label, ", ")
        }
      }
      names(surv$strata)[i] = label
    }
    start = end + 1
  }
  surv$type            = "right"
  names(surv$n.risk)   = NULL
  names(surv$n.event)  = NULL
  names(surv$n.censor) = NULL
  class(surv) = "survfitDistributed"
  return(invisible(surv))
}


survfitDistributed = function(x = NULL, formula = NULL, data = NULL) {
  if (class(x) != "vdracox") {
    cat("Error: the first parameter must be a vdracox object.\n")
    return(invisible(NULL))
  }
  if (is.null(data) && is.null(formula)) {
    return(survfitDistributed.stats(x))
  }
  if (class(data) != "matrix" && class(data) != "data.frame") {
    cat("Error: the data must either be a data or a data.frame.",
        "Please use the same data that you used for the distrubuted regression.\n")
    return(invisible(NULL))
  }
  if (class(formula) != "formula" || !validFormula2(formula)) {
    cat("The formula must be of the form \"~ var1 + ... + vark\" where the variables",
        "are found in the data. The formula can also \"~ 1\".\n")
  }
  return(survfitDistributed.formula(x, formula, data))
}
