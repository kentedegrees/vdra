################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

PrepareFolder.ACDP = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolder.ACDP\n\n")
  if (is.null(monitorFolder)) {
		warning("monitorFolder must be specified.  Please use the same monitorFolder as the DataMart Client.")
		params$failed = TRUE
		return(params)
	}
	if (class(monitorFolder) != "character") {
		warning("monitorFolder directory is not valid.  Please use the same monitorFolder as the DataMart Client.")
		params$failed = TRUE
		return(params)
	}
  while (!dir.exists(monitorFolder)) {
    Sys.sleep(1)
  }

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	if (params$dataPartnerID == 0) {
		params$writePath   = file.path(monitorFolder, "inputfiles")
	} else {
		params$writePath   = file.path(monitorFolder, "msoc")
	}
	params$readPathAC    = file.path(monitorFolder, "inputfiles")
	params$readPathDP    = file.path(monitorFolder, paste0("msoc", 1:params$numDataPartners))

	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$readPathAC, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.")
	}
	for (id in 1:params$numDataPartners) {
		if (!CreateIOLocation(monitorFolder, paste0("msoc", id))) {
			params$failed = TRUE
			params$errorMessage = paste(params$errorMessage,
																	"Could not create directory",
																	paste0(params$readPathDP[id], "."),
																	"Check the path and restart the program.")
		}
	}

	if (params$dataPartnerID != 0) {
	  Sys.sleep(1)
	  DeleteTrigger("files_done.ok", params$readPathAC)
	  for (id in 1:params$numDataPartners) {
	    DeleteTrigger("files_done.ok", params$readpathDP[id])
	  }
	}

	empty = NULL
	writeTime = proc.time()[3]
	save(empty, file = file.path(params$writePath, "empty.rdata"))
	writeSize = file.size(file.path(params$writePath, "empty.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareFolder.ACDP", 0, 0, writeTime, writeSize)
	return(params)
}


#' @importFrom stats model.matrix
PrepareDataLinLog.DP1 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinLog.DP1\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  responseIndex = CheckResponse(params, data, yname)

  if (is.null(responseIndex)) {
    workdata$failed = TRUE
    return(workdata)
  }
  covariateIndex = setdiff(1:ncol(data), responseIndex)
  workdata$tags = CreateModelMatrixTags(data[, covariateIndex, drop = FALSE])
  workdata$tags = c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] = "numeric"
  X = model.matrix(~ ., data[, c(responseIndex, covariateIndex), drop = FALSE])
  rownames(X) = NULL
  covariateIndex = setdiff(1:ncol(X), 2)
  workdata$X = X[, c(2, covariateIndex), drop = FALSE]

  workdata$n        = nrow(workdata$X)
  workdata$colmin   = apply(workdata$X, 2, min)
  workdata$colmax   = apply(workdata$X, 2, max)
  workdata$colsum   = apply(workdata$X, 2, sum)
  workdata$colrange = workdata$colmax - workdata$colmin
  for (i in 1:ncol(workdata$X)) {
    if (workdata$colmin[i] == workdata$colmax[i]) {
      workdata$colmin[i] = 0
      workdata$colrange[i] = 1
    }
    workdata$X[, i] = (workdata$X[, i] - workdata$colmin[i]) / workdata$colrange[i]
  }

  return(workdata)
}

#' @importFrom stats model.matrix
PrepareDataLinLog.DPk = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinLog.DPk\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  workdata$tags = CreateModelMatrixTags(data)
  # if (ncol(data) < 2 | !("numeric" %in% names(workdata$tags))) {
  #   warning("The data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.")
  #   workdata$failed = TRUE
  #   return(workdata)
  # }
  workdata$X = model.matrix(~ ., data)
  rownames(workdata$X) = NULL
  workdata$X = workdata$X[, -1, drop = FALSE]

  workdata$n        = nrow(workdata$X)
  workdata$colmin   = apply(workdata$X, 2, min)
  workdata$colmax   = apply(workdata$X, 2, max)
  workdata$colsum   = apply(workdata$X, 2, sum)
  workdata$colrange = workdata$colmax - workdata$colmin
  for (i in 1:ncol(workdata$X)) {
    if (workdata$colmin[i] == workdata$colmax[i]) {
      workdata$colmin[i] = 0
      workdata$colrange[i] = 1
    }
    workdata$X[, i] = (workdata$X[, i] - workdata$colmin[i]) / workdata$colrange[i]
  }

  return(workdata)
}

SendBasicInfo.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendbasicInfo.DP\n\n")
  n = data$n
	params$n = n
	analysis = params$analysis
	dataPartnerID = params$dataPartnerID
	writeTime = proc.time()[3]
	save(analysis, n, dataPartnerID, file = file.path(params$writePath, "n_analysis.rdata"))
	writeSize = file.size(file.path(params$writePath, "n_analysis.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendBasicInfo.DP", 0, 0, writeTime, writeSize)
	return(params)
}

CheckAgreement.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckAgreement.AC\n\n")
  readTime = 0
	readSize = 0
	analysisAll = rep("", params$numDataPartners)
	nAll        = rep(0, params$numDataPartners)
	nDataPartnerID = rep(0, params$numDataPartners)
	message1    = NULL
	message2    = NULL
	n           = NULL
	analysis    = NULL
	dataPartnerID = NULL
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "n_analysis.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "n_analysis.rdata"))
		readTime = readTime + proc.time()[3]
		analysisAll[id] = analysis
		nAll[id]        = n
		nDataPartnerID[id] = dataPartnerID
	}

	if (any(params$analysis != analysisAll)) {
		params$failed = TRUE
		message1 = "Different regressions have been specified.\n"
		message1 = paste(message1, "Analysis center specified", params$analysis, "regression.\n")
		for (id in 1:params$numDataPartners) {
			message1 = paste(message1, "Data partner", id, "specified", analysisAll[id], "regression.\n")
		}
	}

	if (min(nAll) < max(nAll)) {
		params$failed = TRUE
		message2 = "Data partners provided different numbers of observations.\n"
		for (id in 1:params$numDataPartners) {
			message2 = paste(message2, "Data partner", id, "has", nAll[id], "observations.\n")
		}
	}

	message3error = FALSE
	message3 = ""
	for (i in 1:params$numDataPartners) {
	  if (i != nDataPartnerID[i]) {
	    message3error = TRUE
	    params$failed = TRUE
	    message3 = paste0(message3, "Data Partner ", i, " reports its ID as ", nDataPartnerID[i], "\n")
	  }
	}
	if (message3error) {
	  message3 = paste0("Check PopMedNet DataMart setup.\n", message3)
	}

	if (params$failed) {
		params$errorMessage = paste0(message1, message2, message3)
	}

	params = AddToLog(params, "CheckAgreement.AC", readTime, readSize, 0, 0)
	return(params)
}


#' @importFrom stats runif
PrepareParamsLinear.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.DP\n\n")
  params$n          = nrow(data$X)
	params$p          = ncol(data$X)
	temp = as.numeric(Sys.time())
	set.seed((temp - trunc(temp)) * .Machine$integer.max)
	params$seed       = floor(runif(1) * .Machine$integer.max)
	params$scaler     = 1 + runif(1)

	p = params$p
	seed = params$seed
	scaler = params$scaler

	writeTime = proc.time()[3]
	save(p, scaler, seed, file = file.path(params$writePath, "p_scaler_seed.rdata"))
	writeSize = file.size(file.path(params$writePath, "p_scaler_seed.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareParamsLinear.DP", 0, 0, writeTime, writeSize)
	return(params)
}


#' @importFrom stats rnorm
PrepareSharesLinear.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareSharesLinear.DP\n\n")
  readTime = 0
	readSize = 0
	p = seed = scaler = NULL

	set.seed(params$seed, kind = "Mersenne-Twister")
	halfshare = matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)

	products = rep(list(list()), params$numDataPartners)

	params$ps = c()
	params$scalers = c()
	params$seeds = c()

	for (id in 1:params$numDataPartners) {
		if (id == params$dataPartnerID) {
			products[[id]] = t(data$X) %*% data$X
			params$ps      = c(params$ps, params$p)
			params$scalers = c(params$scalers, params$scaler)
			params$seeds   = c(params$seeds, params$seed)
			next
		}
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
		readTime = readTime + proc.time()[3]
		params$ps      = c(params$ps, p)
		params$scalers = c(params$scalers, scaler)
		params$seeds   = c(params$seeds, seed)

		set.seed(seed, kind = "Mersenne-Twister")
		halfShare2 = matrix(rnorm(params$n * p, sd = 20), nrow = params$n, ncol = p)

		if (id < params$dataPartnerID) {
		  products[[id]] = t(halfShare2) %*% (data$X - scaler / (scaler + params$scaler) * halfshare)
		}

		if (id > params$dataPartnerID) {
			products[[id]] = t(data$X - scaler / (scaler + params$scaler) * halfshare) %*% halfShare2
		}
	}

	halfshare = data$X - halfshare
	colmin    = data$colmin
	colrange  = data$colrange
	colsum    = data$colsum
	colnames  = colnames(data$X)
	tags      = data$tags

	writeTime = proc.time()[3]
	save(products, file = file.path(params$writePath, "products.rdata"))
	save(halfshare, file = file.path(params$writePath, "halfshare.rdata"))
	save(colmin, colrange, colsum, colnames, tags, file = file.path(params$writePath, "colstats.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("products.rdata",
																													"halfshare.rdata",
																													"colstats.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareSharesLinear.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetProductsLinear.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLinear.AC\n\n")
  readTime = 0
	readSize = 0
	p = 0
	n = 0

	allproducts  = rep(list(list()), params$numDataPartners)
	allhalfshare = rep(list(list()), params$numDataPartners)
	alltags      = rep(list(list()), params$numDataPartners)
	products  = NULL
	halfshare = NULL
	tags      = NULL
	allcolmin = allcolrange = allcolsum = allcolnames = NULL
	colmin = colrange = colsum = colnames = NULL
	party = NULL
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "products.rdata"))
		load(file.path(params$readPathDP[id], "halfshare.rdata"))
		load(file.path(params$readPathDP[id], "colstats.rdata"))
		readSize = readSize + sum(file.size(file.path(params$readPathDP[id],
																									c("products.rdata",
																										"halfshare.rdata",
																										"colstats.rdata"))))
		readTime = readTime + proc.time()[3]

		allproducts[[id]]  = products
		allhalfshare[[id]] = halfshare
		alltags[[id]]      = tags
		allcolmin          = c(allcolmin, colmin)
		allcolrange        = c(allcolrange, colrange)
		allcolsum          = c(allcolsum, colsum)
		allcolnames        = c(allcolnames, colnames)
		party              = c(party, rep(paste0("dp", id), length(colnames)))
		p = p + ncol(halfshare)
		if (id == 1) n = nrow(halfshare)
	}

	M = matrix(0, p, p)
	colnames(M) = allcolnames
	rownames(M) = allcolnames
	offset1 = 1
	params$pi = rep(0, params$numDataPartners)
	for (id1 in 1:params$numDataPartners) {
		p1 = ncol(allhalfshare[[id1]])
		params$pi[id1] = p1
		offset2 = offset1
		for (id2 in id1:params$numDataPartners) {
			p2 = ncol(allhalfshare[[id2]])
			if (id1 == id2) {
				M[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = allproducts[[id1]][[id2]]
			} else {
				temp = allproducts[[id1]][[id2]] + allproducts[[id2]][[id1]] +
				     	 t(allhalfshare[[id1]]) %*% allhalfshare[[id2]]
				M[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = temp
				M[offset2:(offset2 + p2 - 1), offset1:(offset1 + p1 - 1)] = t(temp)
			}
			offset2 = offset2 + p2
		}
		offset1 = offset1 + p1
	}

	M = diag(allcolrange) %*% M %*% diag(allcolrange) +
		outer(allcolmin, allcolsum) + outer(allcolsum, allcolmin) -
		n * outer(allcolmin, allcolmin)

	params$xtx          = M[2:p, 2:p, drop = FALSE]
	params$xty          = M[2:p, 1, drop = FALSE]
	params$yty          = M[1, 1]
	params$meansy       = allcolsum[1] / n
	params$means        = allcolsum[-1] / n
	params$n            = n
	params$p            = p
	params$colnames     = allcolnames[-1]
	params$party        = party[-1]
	params$converged    = TRUE
  params$tags         = alltags

	params = AddToLog(params, "GetProductsLinear.AC", readTime, readSize, 0, 0)
	return(params)
}


#' @importFrom  stats pf pt
ComputeResultsLinear.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLinear.AC\n\n")
  stats           = params$stats
	stats$converged = params$converged
	n        = params$n
	yty      = params$yty
	xty      = params$xty
	xtx      = params$xtx
	meansy   = params$meansy

	# First we de-standardize.

	nrow = nrow(xtx)
	indicies = c(1)
	for (i in 2:nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}

	tags = params$tags
	min = 1
	for (id in 1:params$numDataPartners) {
	  max = min + params$pi[id] - 1
	  if (id == 1) {
	    max = max - 1
	  }
	  idx = indicies[which(min <= indicies & indicies <= max)] - min + 1
	  temp = tags[[id]]
	  temp = temp[idx]
	  tags[[id]] = temp
	  min = max + 1
	}

	params$errorMessage = ""
	numeric_found = FALSE
	for (id in 2:params$numDataPartners) {
	  if (length(unique(tags[[id]])) == 0) {
	    params$failed = TRUE
	    params$errorMessage = paste0(params$errorMessage,
	                                 paste("After removing colinear covariates, Data Partner",
	                                       id, "has no covariates."))
	  } else {
	    numeric_found = numeric_found | "numeric" %in% names(tags[[id]])
	  }
	}
	if (!numeric_found) {
	  params$failed = TRUE
	  params$errorMessage = paste0(params$errorMessage,
	                               paste("After removing colinear covariates, no Data Partner > DP1 has a numeric covariate."))
	}

	stats$failed    = params$failed

	p             = length(indicies)
	p1            = ncol(xtx)
	xtx.old       = xtx
	xty.old       = xty
	xtx           = xtx[indicies, indicies, drop = FALSE]
	xty           = xty[indicies, , drop = FALSE]

	invxtx = solve(xtx)
	betas  = drop(invxtx %*% xty)

	numCovariates = p - 1

	#   # If true sse is approximately 0, random variations could cause this
	#   # calculation to be less than 0
	#   # If calculated sse is less than 0, we set it equal to 0.
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
	stats$party                  = params$party
	stats$responseParty          = "dp1"
	stats$coefficients           = rep(NA, p1)
	stats$tvals                  = rep(NA, p1)
	stats$secoef                 = rep(NA, p1)
	stats$pvals                  = rep(NA, p1)

	stats$sse                    = sse
	stats$coefficients[indicies] = betas
	stats$tvals[indicies]        = tvals
	stats$secoef[indicies]       = secoef
	stats$pvals[indicies]        = pvals
	stats$rstderr                = rstderr
	stats$rsquare                = Rsq
	stats$adjrsquare             = adjRsq
	stats$Fstat                  = Fstat
	stats$Fpval                  = Fpval
	stats$df1                    = df1
	stats$df2                    = df2
	stats$n                      = params$n
	stats$xtx                    = xtx.old
	stats$xty                    = xty.old
	stats$yty                    = yty
	stats$meansy                 = meansy
	stats$means                  = params$means

	names(stats$party)           = params$colnames
	names(stats$coefficients)    = params$colnames
	names(stats$secoef)          = params$colnames
	names(stats$tvals)           = params$colnames
	names(stats$pvals)           = params$colnames
	colnames(stats$xtx)          = params$colnames
	rownames(stats$xtx)          = params$colnames
	colnames(stats$xty)          = colnames(params$xty)
	rownames(stats$xty)          = params$colnames

	class(stats) = "vdralinear"

	params$stats = stats
	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeResultsLinear.AC", 0, 0, writeTime, writeSize)
	return(params)
}


GetResultsLinear.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLinear.DP\n\n")
  params$converged = TRUE
	stats = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "stats.rdata"))
	readSize = file.size(file.path(params$readPathAC, "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats

	params = AddToLog(params, "GetResultsLinear.DP", readTime, readSize, 0, 0)
	return(params)
}


############################## PARENT FUNCTIONS ###############################


DataPartnerKLinear = function(data,
															yname           = NULL,
															numDataPartners = NULL,
															dataPartnerID   = NULL,
															monitorFolder   = NULL,
															sleepTime       = 10,
															maxWaitingTime  = 24 * 60 * 60,
															popmednet      = TRUE,
															trace          = FALSE,
															verbose        = TRUE) {

	params = PrepareParams.kp("linear", dataPartnerID, numDataPartners, ac = FALSE,
	                          popmednet = popmednet, trace = trace, verbose = verbose)
	if (params$failed) {
	  warning(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		warning(params$errorMessage)
		return(invisible(NULL))
	}

	if (dataPartnerID == 1) {
	  data = PrepareDataLinLog.DP1(params, data, yname)
	  params = AddToLog(params, "PrepareDataLinLog.DP1", 0, 0, 0, 0)
	} else {
	  data = PrepareDataLinLog.DPk(params, data)
	  params = AddToLog(params, "PrepareDataLinLog.DP2", 0, 0, 0, 0)
	}


	if (data$failed) {
		params$errorMessage = paste("Error processing data for data partner", params$dataPartnerID, "\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		params$errorMessage = ReadErrorMessage(params$readPathAC)
		warning(params$errorMessage)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = SendBasicInfo.DP(params, data)
	files = "n_analysis.rdata"
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	possibleError = ReceivedError.kp(params, from = "AC")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = PrepareParamsLinear.DP(params, data)
	files = "p_scaler_seed.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	params = PrepareSharesLinear.DP(params, data)
	files = c("products.rdata", "halfshare.rdata", "colstats.rdata")
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	possibleError = ReceivedError.kp(params, from = "AC")
	if (possibleError$error) {
	  params$errorMessage = possibleError$message
	  warning(possibleError$message)
	  params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
	  return(params$stats)
	} else {
	  params = GetResultsLinear.DP(params)
	  params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
	  return(params$stats)
	}
}


AnalysisCenterKLinear = function(numDataPartners = NULL,
 	                               monitorFolder   = NULL,
																 msreqid         = "v_default_0_000",
																 sleepTime       = 10,
																 maxWaitingTime  = 24 * 60 * 60,
																 popmednet       = TRUE,
																 trace           = FALSE,
																 verbose         = TRUE) {
	params = PrepareParams.kp("linear", 0, numDataPartners, msreqid, ac = TRUE,
	                          popmednet = popmednet, trace = trace, verbose = verbose)
	if (params$failed) {
	  warning(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		warning(params$errorMessage)
		return(invisible(NULL))
	}

	params = PauseContinue.kp(params, from = "DP", maxWaitingTime = maxWaitingTime)

	possibleError = ReceivedError.kp(params, from = "DP")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		MakeErrorMessage(params$writePath, possibleError$message)
		files = "errorMessage.rdata"
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	params = CheckAgreement.AC(params)

	if (params$failed) {
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		warning(params$errorMessage)
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	files = "empty.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetProductsLinear.AC(params)
	params = ComputeResultsLinear.AC(params)

	if (params$failed) {
	  MakeErrorMessage(params$writePath, params$errorMessage)
	  files = "errorMessage.rdata"
	  warning(params$errorMessage)
	  params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
	                                sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
	  params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
	  SummarizeLog.kp(params)
	  return(params$stats)
	} else {
	  files = "stats.rdata"
	  params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
	                                sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
	  params = SendPauseQuit.kp(params, sleepTime = sleepTime)
	  SummarizeLog.kp(params)
	  return(params$stats)
	}
}
