################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

PrepareFolderLogistic.A2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.A2\n\n")
	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "inputfiles")
	params$readPath      = file.path(monitorFolder, "msoc1")

	if (is.null(monitorFolder)) {
	  cat("monitorFolder must be specified.  Please use the same monitorFolder as the DataMart Client.\n")
	  params$failed = TRUE
	  return(params)
	}
	if (class(monitorFolder) != "character") {
	  cat("monitorFolder directory is not valid.  Please use the same monitorFolder as the DataMart Client.\n")
	  params$failed = TRUE
	  return(params)
	}
	while (!dir.exists(monitorFolder)) {
	  Sys.sleep(1)
	}
	params$errorMessage = NULL
	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc1")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$readPath, "."),
																"Check the path and restart the program.\n\n")
	}

	params = AddToLog(params, "PrepareDataLogistic.A23, PrepareFolderLogistic.A2", 0, 0, 0, 0)
	return(params)
}


PrepareFolderLogistic.B2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.B2\n\n")

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "msoc")
	params$readPath      = file.path(monitorFolder, "inputfiles")

	if (is.null(monitorFolder)) {
	  cat("monitorFolder must be specified.  Please use the same monitorFolder as the DataMart Client.\n")
	  params$failed = TRUE
	  return(params)
	}
	if (class(monitorFolder) != "character") {
	  cat("monitorFolder directory is not valid.  Please use the same monitorFolder as the DataMart Client.\n")
	  params$failed = TRUE
	  return(params)
	}
	while (!dir.exists(monitorFolder)) {
	  Sys.sleep(1)
	}
	params$errorMessage = NULL
	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: could not create directory",
																paste0(params$readPath, "."),
																"Check the path and restart the program.\n\n")
	}

	Sys.sleep(1)
	DeleteTrigger("files_done.ok", params$readPath)

	params = AddToLog(params, "PrepareDataLogisitc.B23, PrepareFolderLogistic.B2", 0, 0, 0, 0)

	return(params)
}


PrepareDataLogistic.A23 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.A23\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  responseIndex = CheckResponse(params, data, yname)

  if (is.null(responseIndex)) {
    workdata$failed == TRUE
    return(workdata)
  }
  covariateIndex = setdiff(1:ncol(data), responseIndex)
  workdata$tags = CreateModelMatrixTags(data[, covariateIndex, drop = FALSE])
  workdata$tags = c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] = "numeric"

  X = model.matrix(~ ., data[, c(responseIndex, covariateIndex), drop = FALSE])
  rownames(X) = NULL
  covariateIndex = setdiff(1:ncol(X), 2)

  means = apply(X, 2, mean)
  sd    = apply(X, 2, sd)
  sd    = sapply(sd, function(x) { ifelse(x > 0, x, 1)})
  workdata$Y      = X[, 2, drop = FALSE]
  workdata$X      = X[, covariateIndex, drop = FALSE]
  # workdata$meansy = means[2]
  # workdata$sdy    = sd[2]
  workdata$means  = means[covariateIndex]
  workdata$sd     = sd[covariateIndex]
  workdata$yty    = t(workdata$Y) %*% workdata$Y

  # workdata$Y      = (workdata$Y - workdata$meansy) / workdata$sdy

  if (ncol(workdata$X) >= 2) {
    for (i in 2:ncol(workdata$X)) {
      workdata$X[, i] = (workdata$X[, i] - workdata$means[i]) / workdata$sd[i]
    }
  }

  return(workdata)
}

PrepareDataLogistic.B23 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.B23\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  workdata$tags = CreateModelMatrixTags(data)
  if (ncol(data) < 2 | !("numeric" %in% names(workdata$tags))) {
    cat("Error: The data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.\n")
    workdata$failed = TRUE
    return(workdata)
  }
  workdata$X = model.matrix(~ ., data)
  rownames(workdata$X) = NULL
  workdata$X = workdata$X[, -1, drop = FALSE]
  workdata$means = apply(workdata$X, 2, mean)
  workdata$sd    = apply(workdata$X, 2, sd)
  workdata$sd    = sapply(workdata$sd, function(x) { ifelse(x > 0, x, 1)})

  for (i in 1:ncol(workdata$X)) {
    workdata$X[, i] = (workdata$X[, i] - workdata$means[i]) / workdata$sd[i]
  }


  return(workdata)
}

PrepareParamsLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLogistic.B2\n\n")
  params$failed         = FALSE
  params$converged      = FALSE
  params$halted         = FALSE

  params$n             = nrow(data$X)
  params$numEvents     = 0
  params$p1            = 0
  params$p2            = ncol(data$X)
  params$p             = params$p1 + params$p2
  params$p1.old        = 0
  params$p2.old        = params$p2
  params$Acolnames     = c("")
  params$Bcolnames     = colnames(data$X)
  params$yname         = ""
  params$Acolnames.old = c("")
  params$Bcolnames.old = c("")
  params$cutoff        = 1
  params$maxIterations = 1

  params$meansA        = 0
  params$sdA           = 0
  params$meansB        = data$means
  params$sdB           = data$sd
  params$yty           = 0

  pb          = list()
  pb$p2       = params$p2
  pb$n        = params$n
  pb$means    = data$means
  pb$sd       = data$sd
  pb$analysis = params$analysis
  pb$Bcolnames = params$Bcolnames
  pb$tags      = data$tags

  writeTime = proc.time()[3]
  save(pb, file = file.path(params$writePath, "pb.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pb.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsLogistic.B2", 0, 0, writeTime, writeSize)
  return(params)
}


PrepareParamsLogistic.A2 = function(params, data, cutoff = 0.01, maxIterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLogistic.A2\n\n")

  params$converged       = FALSE
  params$halted          = FALSE
  params$pmnStepCounter  = 1
  pb                     = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "pb.rdata")) # load pb, Bcolnames
  readSize = sum(file.size(file.path(params$readPath, "pb.rdata")))
  readTime = proc.time()[3] - readTime

  if (params$analysis != pb$analysis) {
    params$errorMessage =
      paste("Party A is running", params$analysis, "regression and Party B is running", pb$analysis, "regression.")
    cat("Error:", params$errorMessage, "\n\n")
    params$failed = TRUE
    return(params)
  }

  params$n  = nrow(data$X)
  if (pb$n != params$n) {
    params$errorMessage =
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    cat("Error:", params$errorMessage, "\n\n")
    params$failed = TRUE
  }

  params$p1 = ncol(data$X)
  params$p2 = pb$p2
  params$p  = params$p1 + params$p2
  params$p1.old = params$p1
  params$p2.old = params$p2

  params$Acolnames = colnames(data$X)
  params$Bcolnames = pb$Bcolnames
  params$yname     = colnames(data$Y)
  params$Acolnames.old = c("")
  params$Bcolnames.old = c("")
  params$Atags         = data$tags
  params$Btags         = pb$tags

  if (cutoff <= 0) cutoff = 0.01
  if (cutoff >= 1) cutoff = 0.05
  params$cutoff           = cutoff

  if (maxIterations < 1) maxIterations = 1
  params$maxIterations = maxIterations

  params$meansA = data$means
  params$sdA    = data$sd
  params$meansB = pb$means
  params$sdB    = pb$sd
  params$yty    = data$yty

  pa               = list()
  pa$p1            = params$p1
  pa$means         = data$means
  pa$sd            = data$sd
  pa$yty           = data$yty
  pa$yname         = data$yname
  pa$cutoff        = params$cutoff
  pa$maxIterations = params$maxIterations
  pa$Acolnames     = params$Acolnames
  pa$tags          = data$tags

  writeTime = proc.time()[3]
  save(pa, file = file.path(params$writePath, "pa.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pa.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsLogistic.A2", readTime, readSize,
                    writeTime, writeSize)

  return(params)
}


PrepareBlocksLogistic.A2 = function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLogistic.A2\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  = params$n
  p1 = params$p1
  p2 = params$p2

  minimumBlocksize = GetBlockSize(p1, p2)
  if (n < minimumBlocksize) {
    maxACovariates = trunc(sqrt(p2 * n) - p2 - 1)

    params$errorMessage =
      paste("Warning, the minimum secure blocksize of", minimumBlocksize,
            "is larger than the number of observations", paste0(n, "."))
    cat("Error: The minimum secure blocksize of", minimumBlocksize,
        "is larger than the number of observations", paste0(n, ".\n"))
    cat("Your options are:\n")
    cat("Increase the number of observations to at least",
        paste0(minimumBlocksize, ".\n"))
    cat("Decrease the number of A covariates to", maxACovariates, "or less.\n")
    b = n - 2 * p1 - 2
    discrim = b^2 - 4 * (p1 + 1)^2
    if (discrim >= 0) {
      minBCovariates = trunc(1 + (b - sqrt(discrim)) / 2)
      maxBCovariates = trunc((b + sqrt(discrim)) / 2)
      cat("Set the number of B covariates to be between", minBCovariates, "and",
          paste0(maxBCovariates, ".\n"))
    }
    cat("Terminating the program.\n\n")
    params$failed = TRUE
    params = AddToLog(params, "PrepareBlocksLogistic.A2", 0, 0, 0, 0)
    return(params)
  }

  if (is.null(blocksize)) {
    blocksize = minimumBlocksize
  }
  if (blocksize < minimumBlocksize) {
    cat("Warning:  Block size of", blocksize,
        "is too small. Proceeding with minimum blocksize of",
        paste0(minimumBlocksize, ".\n\n"))
    blocksize = minimumBlocksize
  } else if (n < blocksize) {
    cat("Warning: Block size of", blocksize,
        "is larger than size of data.  Proceeding with blocksize of",
        paste0(n, ".\n\n"))
  }

  params$blocks    = CreateBlocks(p1, p2, n, blocksize)
  params$container = CreateContainers(p1, p2, params$blocks)
  writeTime = proc.time()[3]
  save(blocksize, file = file.path(params$writePath, "blocksize.rdata"))
  writeSize = file.size(file.path(params$writePath, "blocksize.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareBlocksLogistic.A2", 0, 0, writeTime, writeSize)
  return(params)
}


GetZLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLogistic.A2\n\n")
  writeTime = 0
  writeSize = 0

  numBlocks = params$blocks$numBlocks
  pbar = MakeProgressBar1(numBlocks, "Z")
  containerCt.Z = 0
  for (i in 1:numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename = paste0("cz_", containerCt.Z, ".rdata")
      toWrite = file(file.path(params$writePath, filename), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]
    Z = FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$X[strt:stp, ]), g)

    writeTime = writeTime - proc.time()[3]
    writeBin(as.vector(Z), con = toWrite, endian = "little")
    writeTime = writeTime + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
      close(toWrite)
      writeSize = writeSize + file.size(file.path(params$writePath, filename))
    }
    pbar = MakeProgressBar2(i, pbar)
  }
  params = AddToLog(params, "GetZLogistic.A2", 0, 0, writeTime, writeSize)
  return(params)
}


FinalizeParamsLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsLogistic.B2\n\n")
  pa = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "pa.rdata")) # Load pa, Acolnames
  readSize = sum(file.size(file.path(params$readPath, "pa.rdata")))
  readTime = proc.time()[3] - readTime
  params$p1            = pa$p1
  params$p             = params$p1 + params$p2
  params$meansA        = pa$means
  params$sdA           = pa$sd
  params$yty           = pa$yty
  params$yname         = pa$yname
  params$cutoff        = pa$cutoff
  params$maxIterations = pa$maxIterations
  params$Acolnames     = pa$Acolnames
  params$Atags         = pa$tags
  params$Btags         = data$tags

  params = AddToLog(params, "FinalizeParamsLogistic.B2", readTime, readSize, 0, 0)
  return(params)
}


PrepareBlocksLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLogistic.B2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  readTime = proc.time()[3]
  load(file.path(params$readPath, "blocksize.rdata")) # load blocksize
  readSize = file.size(file.path(params$readPath, "blocksize.rdata"))
  readTime = proc.time()[3] - readTime
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params = AddToLog(params, "PrepareBlocksLogistic.B2", readTime, readSize, 0, 0)
  return(params)
}


GetWLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWLogistic.B2\n\n")
  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I-Z*Z')XB")

  XBTXB = t(data$X) %*% data$X

  containerCt.Z = 0
  containerCt.W = 0

  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename1 = paste0("cz_", containerCt.Z, ".rdata")
      toRead = file(file.path(params$readPath, filename1), "rb")
      readSize = readSize + file.size(file.path(params$readPath, filename1))
    }
    if (i %in% params$container$filebreak.W) {
      containerCt.W = containerCt.W + 1
      filename2 = paste0("cw_", containerCt.W, ".rdata")
      toWrite = file(file.path(params$writePath, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1
    g1 = params$blocks$g[i]

    readTime = readTime - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n2 * g1,
                       endian = "little"), nrow = n2, ncol = g1)
    readTime = readTime + proc.time()[3]

    W = data$X[strt:stp, ] - Z %*% (t(Z) %*% data$X[strt:stp, ])

    writeTime = writeTime - proc.time()[3]
    writeBin(as.vector(W), con = toWrite, endian = "little")
    writeTime = writeTime + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.Z || i == params$blocks$numBlocks) {
      close(toRead)
    }
    if ((i + 1) %in% params$container$filebreak.W || i == params$blocks$numBlocks) {
      close(toWrite)
      writeSize = writeSize + file.size(file.path(params$writePath, filename2))
    }

    pbar = MakeProgressBar2(i, pbar)
  }

  writeTime = writeTime - proc.time()[3]
  save(XBTXB, file = file.path(params$writePath, "xbtxb.rdata"))
  writeSize = writeSize + file.size(file.path(params$writePath, "xbtxb.rdata"))
  writeTime = writeTime + proc.time()[3]

  params = AddToLog(params, "GetWLogistic.B2", readTime, readSize, writeTime, writeSize)

  return(params)
}


CheckColinearityLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.A2\n\n")
  p2 = params$p2
  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0
  XBTXB     = NULL

  readTime = readTime - proc.time()[3]
  load(file.path(params$readPath, "xbtxb.rdata")) # load XBTXB
  readSize = file.size(file.path(params$readPath, "xbtxb.rdata"))
  readTime = readTime + proc.time()[3]
  XATXA = t(data$X) %*% data$X
  XATXB = 0
  XATY  = t(data$X) %*% data$Y
  YTXB  = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'X")

  containerCt.W = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.W) {
      containerCt.W = containerCt.W + 1
      filename = paste0("cw_", containerCt.W, ".rdata")
      toRead = file(file.path(params$readPath, filename), "rb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1

    readTime = readTime - proc.time()[3]
    W = matrix(readBin(con = toRead, what = numeric(), n = n2 * p2,
                       endian = "little"), nrow = n2, ncol = p2)
    readTime = readTime + proc.time()[3]

    XATXB = XATXB + t(data$X[strt:stp, ]) %*% W
    YTXB  = YTXB  + t(data$Y[strt:stp, ]) %*% W

    if ((i + 1) %in% params$container$filebreak.W || i == params$blocks$numBlocks) {
      close(toRead)
      readSize = readSize + file.size(file.path(params$readPath, filename))
    }
    pbar = MakeProgressBar2(i, pbar)
  }

  XTX = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), XBTXB))
  XTY = rbind(XATY, t(YTXB))

  nrow = nrow(XTX)
  indicies = c(1)
  for (i in 2:nrow) {
  	tempIndicies = c(indicies, i)
  	if (rcond(XTX[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
  		indicies = c(indicies, i)
  	}
  }

  XTX = XTX[indicies, indicies]
  XTY = matrix(XTY[indicies], ncol = 1)

  Anames   = params$Acolnames
  Bnames   = params$Bcolnames
  Aindex   = which(indicies <= length(Anames))
  params$IndiciesKeep  = indicies
  params$AIndiciesKeep = indicies[Aindex]
  params$BIndiciesKeep = indicies[-Aindex] - length(Anames)

  AnamesKeep = Anames[params$AIndiciesKeep]
  BnamesKeep = Bnames[params$BIndiciesKeep]
  params$Acolnames.old = params$Acolnames
  params$Bcolnames.old = params$Bcolnames
  params$Acolnames     = AnamesKeep
  params$Bcolnames     = BnamesKeep
  params$p1.old        = params$p1
  params$p2.old        = params$p2
  params$p1            = length(AnamesKeep)
  params$p2            = length(BnamesKeep)
  params$p.old         = params$p1.old + params$p2.old
  params$p             = params$p1 + params$p2
  params$meansA        = params$meansA[params$AIndiciesKeep]
  params$meansB        = params$meansB[params$BIndiciesKeep]
  params$sdA           = params$sdA[params$AIndiciesKeep]
  params$sdB           = params$sdB[params$BIndiciesKeep]
  params$xtx           = XTX
  params$xty           = XTY

  Aindicies = params$AIndiciesKeep
  Bindicies = params$BIndiciesKeep

  writeTime = writeTime - proc.time()[3]
  save(Aindicies, Bindicies, file = file.path(params$writePath, "indicies.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "indicies.rdata")))
  writeTime = writeTime + proc.time()[3]

  tags = params$Btags[params$BIndiciesKeep]

  if (length(unique(tags)) < 2) {
    params$failed = TRUE
    params$errorMessage = "After removing colinear covariates, Party B has 1 or fewer covariates.\n\n"
  } else if (!("numeric" %in% names(tags))) {
    params$failed = TRUE
    params$errorMessage = "After removing colinear covariates, Party B has no continuous covariates.\n\n"
  }

  params = AddToLog(params, "CheckColinearityLogistic.A2", readTime, readSize, writeTime, writeSize)
  return(params)
}


UpdateDataLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.A2\n\n")
  data$X = as.matrix(data$X[, params$AIndiciesKeep, drop = FALSE])
  return(data)
}


ComputeInitialBetasLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.A2\n\n")
  # de-standardize xty
  p1     = params$p1
  p2     = params$p2
  xty    = params$xty
  xtx    = params$xtx

  betas = 4 * solve(xtx) %*% xty

  Abetas   = betas[1:p1]
  Bbetas   = betas[(p1 + 1):(p1 + p2)]
  Axty     = xty[1:p1]
  Bxty     = xty[(p1 + 1):(p1 + p2)]

  params$Axty      = Axty
  params$Bxty      = Bxty
  params$betas     = betas
  params$betasA    = Abetas
  params$betasAold = matrix(0, p1, 1)
  params$betasB    = Bbetas

  params$algIterationCounter      = 1
  params$deltabeta = Inf

  writeTime = proc.time()[3]
  save(Bbetas, Bxty, file = file.path(params$writePath, "Bbetas_xty.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "Bbetas_xty.rdata")))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "ComputeInitialBetasLogistic.A2", 0, 0, writeTime, writeSize)

  return(params)
}


UpdateParamsLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.B2\n\n")
  Aindicies = NULL
  Bindicies = NULL
  Bbetas    = NULL
  Bxty      = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "indicies.rdata")) # load Aindicies, Bindicies
  load(file.path(params$readPath, "Bbetas_xty.rdata"))     # Load Bbetas
  readSize = sum(file.size(file.path(params$readPath, c("indicies.rdata",
                                                        "Bbetas_xty.rdata"))))
  readTime = proc.time()[3] - readTime
  params$Acolnames.old = params$Acolnames
  params$Bcolnames.old = params$Bcolnames
  params$Acolnames     = params$Acolnames.old[Aindicies]
  params$Bcolnames     = params$Bcolnames.old[Bindicies]
  params$p1.old = params$p1
  params$p2.old = params$p2
  params$p1     = length(Aindicies)
  params$p2     = length(Bindicies)
  params$p.old  = params$p
  params$p      = params$p1 + params$p2
  params$BIndiciesKeep = Bindicies
  params$AIndiciesKeep = Aindicies
  params$betasB    = Bbetas
  params$betasBold = matrix(0, params$p2, 1)
  params$meansB = params$meansB[Bindicies]
  params$sdB    = params$sdB[Bindicies]
  params$Bxty   = Bxty
  params = AddToLog(params, "UpdateParamsLogistic.B2", readTime, readSize, 0, 0)
  return(params)
}


UpdateDataLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.B2\n\n")
  data$X = as.matrix(data$X[, params$BIndiciesKeep, drop = FALSE])
  return(data)
}


GetXBetaLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaLogistic.B2\n\n")
  XbetaB = data$X %*% params$betasB

  writeTime = proc.time()[3]
  save(XbetaB, file = file.path(params$writePath, "xbetab.rdata"))
  writeSize = file.size(file.path(params$writePath, "xbetab.rdata"))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "GetXBetaLogistic.B2", 0, 0, writeTime, writeSize)
  return(params)
}


GetWeightsLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWeightsLogistic.A2\n\n")
  n      = params$n
  p1     = params$p1
  XbetaB = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "xbetab.rdata"))  # Load XbetaB
  readSize = file.size(file.path(params$readPath, "xbetab.rdata"))
  readTime = proc.time()[3] - readTime

  XbetaA = data$X %*% params$betasA
  Xbeta = XbetaA + XbetaB
  pi_ = (1 + exp(-Xbeta))^(-1)
  params$pi_ = pi_

  writeTime = proc.time()[3]
  save(pi_, file = file.path(params$writePath, "pi_.rdata"))
  writeSize = file.size(file.path(params$writePath, "pi_.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "GetWeightsLogistic.A2", readTime, readSize, writeTime, writeSize)
  return(params)
}


GetVLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetVLogistic.B2\n\n")
  pi_       = NULL
  writeTime = 0
  writeSize = 0
  readTime  = proc.time()[3]
  load(file.path(params$readPath, "pi_.rdata"))
  readSize  = file.size(file.path(params$readPath, "pi_.rdata"))
  readTime  = proc.time()[3] - readTime

  params$pi_ = pi_
  W = pi_ * (1 - params$pi_)

  XBTWXB = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I - Z*Z')W*XB")

  containerCt.Z = 0
  containerCt.V = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename1 = paste0("cz_", containerCt.Z, ".rdata")
      toRead = file(file.path(params$readPath, filename1), "rb")
    }
    if (i %in% params$container$filebreak.V) {
      containerCt.V = containerCt.V + 1
      filename2 = paste0("cv_", containerCt.V, ".rdata")
      toWrite = file(file.path(params$writePath, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]

    Xblock  = data$X[strt:stp, ]
    Wblock  = W[strt:stp]
    WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

    readTime = readTime - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n * g,
                       endian = "little"), nrow = n, ncol = g)
    readTime = readTime + proc.time()[3]

    V = WXblock - Z %*% (t(Z) %*% WXblock)

    writeTime = writeTime - proc.time()[3]
    writeBin(as.vector(V), con = toWrite, endian = "little")
    writeTime = writeTime + proc.time()[3]

    XBTWXB = XBTWXB + t(Xblock) %*% WXblock
    if ((i + 1) %in% params$container$filebreak.Z || i == params$blocks$numBlocks) {
      close(toRead)
      readSize = readSize + file.size(file.path(params$readPath, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$numBlocks) {
      close(toWrite)
      writeSize = writeSize + file.size(file.path(params$writePath, filename2))
    }
    pbar = MakeProgressBar2(i, pbar)
  }
  # sums of each column in WX_B
  sumsWXB = apply(MultiplyDiagonalWTimesX(W, data$X), 2, sum)
  # This information needs to be shared in order to get the intercept term

  writeTime = writeTime - proc.time()[3]
  save(sumsWXB, XBTWXB, file = file.path(params$writePath, "sumswx_xbtwxb.rdata"))
  writeSize = writeSize + sum(file.size(c(file.path(params$writePath, "sumswx_xbtwxb.rdata"))))
  writeTime = writeTime + proc.time()[3]

  params = AddToLog(params, "GetVLogistic.B2", readTime, readSize, writeTime, writeSize)
  return(params)
}


GetIILogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetIILogistic.A2\n\n")
  p1 = params$p1
  p2 = params$p2
  sumsWXB = NULL
  XBTWXB  = NULL

  writeTime = 0
  writeSize = 0
  readTime = proc.time()[3]
  load(file.path(params$readPath, "sumswx_xbtwxb.rdata")) # load sumsWXB, XBTWXB
  readSize = sum(file.size(file.path(params$readPath, "sumswx_xbtwxb.rdata")))
  readTime = proc.time()[3] - readTime

  params$sumsWXB = sumsWXB

  IA = params$Axty - t(data$X) %*% params$pi_
  W = params$pi_ * (1 - params$pi_)
  sumsWXA = apply(MultiplyDiagonalWTimesX(W, data$X), 2, sum)[-1]
  params$sumsWXA = sumsWXA

  XATWXA = t(data$X) %*% MultiplyDiagonalWTimesX(W, data$X)

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'W*X")

  XATWXB = 0
  containerCt.V = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.V) {
      containerCt.V = containerCt.V + 1
      filename1 = paste0("cv_", containerCt.V, ".rdata")
      toRead = file(file.path(params$readPath, filename1), "rb")
      readSize = readSize + file.size(file.path(params$readPath, filename1))
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1

    readTime = readTime - proc.time()[3]
    V = matrix(readBin(con = toRead, what = numeric(),
                       n = n * p2, endian = "little"), n, p2)
    readTime = readTime + proc.time()[3]
    XATWXB = XATWXB + t(data$X[strt:stp, ]) %*% V

    pbar = MakeProgressBar2(i, pbar)
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$numBlocks) {
      close(toRead)
    }
  }

  XTWX = rbind(cbind(XATWXA, XATWXB), cbind(t(XATWXB), XBTWXB))

  params$xtwx = XTWX

  II = NULL
  tryCatch({II = solve(params$xtwx)}, # dims are 1 + p1 + p2
           error = function(err) { II = NULL }
          )
  if (is.null(II)) {
    params$singularMatrix = TRUE
    params$failed = TRUE
    params$errorMessage = "The matrix t(X)WX is singular.  This is probably due to divergence of the coefficients."
    cat("ERROR: The matrix t(X)*W*X is not invertible.\n")
    cat("       This may be due to one of two possible problems.\n")
    cat("       1. Poor random initialization of the security vector.\n")
    cat("       2. Near multicollinearity in the data\n")
    cat("SOLUTIONS: \n")
    cat("       1. Rerun the data analysis.\n")
    cat("       2. If the problem persists, check the variables for\n")
    cat("          duplicates for both parties and / or reduce the\n")
    cat("          number of variables used. Once this is done,\n")
    cat("          rerun the data analysis.\n\n")
  } else {
    params$II = II
    params$IA = IA

    a21i1 = II[(p1 + 1):(p1 + p2), 1:p1] %*% matrix(IA, p1, 1)
    a11i1 = II[1:p1, 1:p1] %*% matrix(IA, p1, 1)
    params$a11i1 = a11i1

    writeTime = proc.time()[3]
    save(a21i1, XTWX, file = file.path(params$writePath, "a21i1_xtwx.rdata"))
    writeSize = sum(file.size(file.path(params$writePath, "a21i1_xtwx.rdata")))
    writeTime = proc.time()[3] - writeTime
  }
  params = AddToLog(params, "GetIILogistic.A2", readTime, readSize, writeTime, writeSize)

  return(params)
}


GetCoefLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.B2\n\n")
  p1 = params$p1
  p2 = params$p2
  XTWX  = NULL
  a21i1 = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "a21i1_xtwx.rdata"))   # load a21i1, XTWX
  readSize = sum(file.size(file.path(params$readPath, "a21i1_xtwx.rdata")))
  readTime = proc.time()[3] - readTime

  IB = params$Bxty - t(data$X) %*% params$pi_

  II = solve(XTWX)

  params$II = II
  params$IB = IB

  a22i2 = II[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  a12i2 = II[1:p1, (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  params$a22i2 = a22i2

  params$betasBold = params$betasB
  params$betasB = params$betasB + a21i1 + a22i2

  deltabetaB = max( abs(params$betasB - params$betasBold) / (abs(params$betasB) + 0.1))

  writeTime = proc.time()[3]
  save(a12i2, deltabetaB, file = file.path(params$writePath, "a12_deltabetaB.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "a12_deltabetaB.rdata")))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "GetCoefLogistic.B2", readTime, readSize, writeTime, writeSize)
  return(params)
}


GetCoefLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.A2\n\n")
  a12i2      = NULL
  deltabetaB = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "a12_deltabetaB.rdata")) # load  a12i2, deltabetab
  readSize = sum(file.size(file.path(params$readPath, "a12_deltabetaB.rdata")))
  readTime = proc.time()[3] - readTime

  params$betasAold = params$betasA
  params$betasA = params$betasA + params$a11i1 + a12i2

  deltabeta = max(abs(params$betasA - params$betasAold) / (abs(params$betasA) + 0.1), deltabetaB)

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$algIterationCounter >= params$maxIterations) {
    params$maxIterExceeded = TRUE
  }

  writeTime = proc.time()[3]
  save(deltabeta, file = file.path(params$writePath, "deltabeta.rdata"))
  writeSize = file.size(file.path(params$writePath, "deltabeta.rdata"))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "GetCoefLogistic.A2", readTime, readSize, writeTime, writeSize)


  return(params)
}


GetConvergedStatusLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetconvergedStatusLogistic.B2\n\n")
  deltabeta = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "deltabeta.rdata")) # load deltabeta.rdata
  readSize = file.size(file.path(params$readPath, "deltabeta.rdata"))
  readTime = proc.time()[3] - readTime

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$algIterationCounter >= params$maxIterations) {
    params$maxIterExceeded = TRUE
  }

  params = AddToLog(params, "GetConvergedStatusLogistic.B2", readTime, readSize, 0, 0)
  return(params)
}


GetFinalCoefLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetFinalCoefLogistic.B2\n\n")
  betasB = params$betasB / params$sdB
  offsetB = sum(betasB * params$meansB)
  BFinalFitted = t(params$sdB * t(data$X) + params$meansB) %*% betasB
  writeTime = proc.time()[3]
  save(betasB, BFinalFitted, offsetB, file = file.path(params$writePath, "b_final.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "b_final.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "GetFinalCoefLogistic.B2", 0, 0, writeTime, writeSize)
  return(params)
}


ComputeResultsLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.A2\n\n")
  stats = params$stats
  stats$failed         = FALSE
  stats$converged      = params$converged

  n      = params$n
  p1     = params$p1
  p2     = params$p2
  sdA    = params$sdA
  sdB    = params$sdB
  meansA = params$meansA
  meansB = params$meansB
  Anames = params$Acolnames.old
  Bnames = params$Bcolnames.old
  p1.old = params$p1.old
  p2.old = params$p2.old
  p.old  = params$p.old
  indicies = params$IndiciesKeep


  betasB = NULL;
  offsetB = NULL;
  BFinalFitted = NULL;
  readTime = proc.time()[3]
  load(file.path(params$readPath, "b_final.rdata"))  # betasB, offsetB, BFinalFitted
  readSize = sum(file.size(file.path(params$readPath, "b_final.rdata")))
  readTime = proc.time()[3] - readTime
  betasA = params$betasA / sdA
  offsetA = sum(betasA[-1] * params$meansA[-1])
  betasA[1] = betasA[1] - offsetA - offsetB
  betas = c(betasA, betasB)

  AFinalFitted = t(sdA * t(data$X) + meansA) %*% betasA -
  	             t(sdA[1] * t(data$X[, 1]) + meansA[1]) %*% betasA[1]
  FinalFitted = AFinalFitted + BFinalFitted + betas[1]
  params$FinalFitted = FinalFitted

  n = params$n
  ct      = sum(data$Y)
  resdev  = -2 * (sum(data$Y * FinalFitted) - sum(log(1 + exp(FinalFitted))))
  nulldev = -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

  # If xtwx were singular, it would have been caught in GetII.A2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  cov1 = solve(params$xtwx)
  secoef = sqrt(diag(cov1)) / c(sdA, sdB)
  tmp = matrix(c(1, (-meansA / sdA)[-1], -meansB / sdB), ncol = 1)
  secoef[1] = sqrt(t(tmp) %*% cov1 %*% tmp)


  stats$party = c(rep("dp0", p1.old), rep("dp1", p2.old))
  stats$coefficients = rep(NA, p.old)
  stats$secoef = rep(NA, p.old)
  stats$tvals  = rep(NA, p.old)
  stats$pvals  = rep(NA, p.old)
  stats$n  = n
  stats$nulldev = nulldev
  stats$resdev = resdev
  stats$aic = resdev + 2 * (p1 + p2)
  stats$bic = resdev + (p1 + p2) * log(n)
  stats$nulldev_df = n - 1
  stats$resdev_df = n - (p1 + p2)
  stats$coefficients[indicies] = betas
  stats$secoef[indicies] = secoef
  tvals = betas / secoef
  pvals = 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[indicies] = tvals
  stats$pvals[indicies] = pvals

  stats$nulldev = nulldev
  stats$resdev  = resdev
  stats$hoslem  = HoslemInternal(params, data)
  stats$ROC     = RocInternal(params, data)
  stats$iter    = params$algIterationCounter - 1

  names.old = c(Anames, Bnames)
  names(stats$coefficients) = names.old
  names(stats$party) = names.old
  names(stats$secoef) = names.old
  names(stats$tvals) = names.old
  names(stats$pvals) = names.old

  writeTime = proc.time()[3]
  save(stats, file = file.path(params$writePath, "stats.rdata"))
  writeSize = file.size(file.path(params$writePath, "stats.rdata"))
  writeTime = proc.time()[3] - writeTime

  stats$Y           = data$Y # For Hoslem and ROC
  stats$FinalFitted = FinalFitted
  params$stats      = stats

  params = AddToLog(params, "ComputeResultsLogistic.B2", readTime, readSize, writeTime, writeSize)
  return(params)
}


GetResultsLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.B2\n\n")
  stats = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "stats.rdata"))
  readSize = file.size(file.path(params$readPath, "stats.rdata"))
  readTime = proc.time()[3] - readTime
  params$stats = stats
  params = AddToLog(params, "GetResultsLogistic.B2", readTime, readSize, 0, 0)
  return(params)
}



############################### PARENT FUNCTIONS ###############################


PartyAProcess2Logistic = function(data,
                                  yname                 = NULL,
																	monitorFolder         = NULL,
																	msreqid               = "v_default_00_000",
                                  blocksize             = 500,
                                  cutoff                = 1e-8,
                                  maxIterations         = 25,
                                  sleepTime             = 10,
                                  maxWaitingTime        = 24 * 60 * 60,
																	popmednet             = TRUE,
																	trace                 = FALSE) {
  params = PrepareParams.2p("logistic", "A", msreqid = msreqid,
                            popmednet = popmednet, trace = trace)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)
  Header(params)
  params   = PrepareFolderLogistic.A2(params, monitorFolder)
  if (params$failed) {
  	cat(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataLogistic.A23(params, data, yname)

  params = PauseContinue.2p(params,  maxWaitingTime)
  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$completed = TRUE
    cat("Error:", ReadErrorMessage(params$readPath), "\n\n")
    params$pmnStepCounter = 1
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
  	params$completed = TRUE
  	message = "Error in processing the data for Party A."
  	MakeErrorMessage(params$writePath, message)
  	files = c("errorMessage.rdata")
  	params$pmnStepCounter = 1
  	params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
  	params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
  	SummarizeLog.2p(params)
  	return(params$stats)
  }

  params = PrepareParamsLogistic.A2(params, data, cutoff, maxIterations)

  if (params$failed) {   # Check for failed from PrepareParamsLogistic.A2()
    params$completed = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = PrepareBlocksLogistic.A2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksLogistic.A2()
    params$completed = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = GetZLogistic.A2(params, data)

  files = c("pa.rdata", "blocksize.rdata",
            SeqZW("cz_", length(params$container$filebreak.Z)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = CheckColinearityLogistic.A2(params, data)

  if (params$failed) { # Check for CheckColinearityLogistic.A2() failed
    params$completed = TRUE
    cat("Error:", params$errorMessage, "\n")
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }
  data = UpdateDataLogistic.A2(params, data)
  params = AddToLog(params, "UpdateDataLogistic.A2", 0, 0, 0, 0)
  params = ComputeInitialBetasLogistic.A2(params, data)

  files = c("indicies.rdata", "Bbetas_xty.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params = GetWeightsLogistic.A2(params, data)
    files = c("pi_.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
    params = GetIILogistic.A2(params, data)

    if (params$failed) { # Check for failed from ComputeInverseLogistic.A2()
      params$completed = TRUE
      MakeErrorMessage(params$writePath, params$errorMessage)
      files = c("errorMessage.rdata")
      params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    files = c("a21i1_xtwx.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    params = GetCoefLogistic.A2(params, data)
    files = "deltabeta.rdata"
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$completed = TRUE

  params = ComputeResultsLogistic.A2(params, data)

  files = c("stats.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  SummarizeLog.2p(params)
  return(invisible(params$stats))
}

PartyBProcess2Logistic = function(data,
																	monitorFolder       = "v_default_00_000",
                                  sleepTime           = 10,
                                  maxWaitingTime      = 24 * 60 * 60,
																	popmednet           = TRUE,
																	trace               = FALSE) {
  params = PrepareParams.2p("logistic", "B",
                            popmednet = popmednet, trace = trace)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)
  Header(params)
  params   = PrepareFolderLogistic.B2(params, monitorFolder)
  if (params$failed) {
  	cat(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataLogistic.B23(params, data)

  if (data$failed) { # Check for Error from PrepareDataLogistic.B2()
    params$completed = TRUE
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$writePath, message)
    files = c("errorMessage.rdata")
    params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params   = PrepareParamsLogistic.B2(params, data)

  files = c("pb.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$completed = TRUE
    cat("Error:", ReadErrorMessage(params$readPath), "\n\n")
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params = FinalizeParamsLogistic.B2(params, data)
  params = PrepareBlocksLogistic.B2(params)
  params = GetWLogistic.B2(params, data)

  files = c("xbtxb.rdata", SeqZW("cw_", length(params$container$filebreak.W)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$completed = TRUE
    cat("Error:", ReadErrorMessage(params$readPath), "\n\n")
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params = UpdateParamsLogistic.B2(params)
  data = UpdateDataLogistic.B2(params, data)
  params = AddToLog(params, "UpdateDataLogistic.B2", 0, 0, 0, 0)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params = GetXBetaLogistic.B2(params, data)

    files = c("xbetab.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    params = GetVLogistic.B2(params, data)
    files = c("sumswx_xbtwxb.rdata",
              SeqZW("cv_", length(params$container$filebreak.V)))
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
      params$completed = TRUE
      cat("Error:", ReadErrorMessage(params$readPath), "\n\n")
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }

    params = GetCoefLogistic.B2(params, data)
    files = c("a12_deltabetaB.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    params = GetConvergedStatusLogistic.B2(params)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$completed = TRUE

  params = GetFinalCoefLogistic.B2(params, data)
  files = c("b_final.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = GetResultsLogistic.B2(params)
  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  return(invisible(params$stats))
}
