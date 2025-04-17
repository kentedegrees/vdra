################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

PrepareFolderLinear.A2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.A2\n\n")
	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "inputfiles")
	params$readPath      = file.path(monitorFolder, "msoc1")

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

	params$errorMessage = NULL
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
																paste0(params$writePath, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "msoc1")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$readPath, "."),
																"Check the path and restart the program.")
	}

	params = AddToLog(params, "PrepareDataLinear.A23, PrepareFolderLinear.A2", 0, 0, 0, 0)
	return(params)
}


PrepareFolderLinear.B2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.B2\n\n")

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "msoc")
	params$readPath      = file.path(monitorFolder, "inputfiles")

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

	params$errorMessage = NULL
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
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Could not create directory",
																paste0(params$readPath, "."),
																"Check the path and restart the program.")
	}

	Sys.sleep(1)
	DeleteTrigger("files_done.ok", params$readPath)

	params = AddToLog(params, "PrepareDataLinear.B23, PrepareFolderLinear.B2", 0, 0, 0, 0)
	return(params)
}

#' @importFrom stats model.matrix sd
PrepareDataLinear.A23 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.A23\n\n")

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
  means = apply(X, 2, mean)
  sd    = apply(X, 2, sd)
  sd    = sapply(sd, function(x) { ifelse(x > 0, x, 1)})
  workdata$Y      = X[, 2, drop = FALSE]
  workdata$X      = X[, covariateIndex, drop = FALSE]
  workdata$meansy = means[2]
  workdata$sdy    = sd[2]
  workdata$means  = means[covariateIndex]
  workdata$sd     = sd[covariateIndex]
  workdata$yty    = t(workdata$Y) %*% workdata$Y

  workdata$Y      = (workdata$Y - workdata$meansy) / workdata$sdy

  if (ncol(workdata$X) >= 2) {
    for (i in 2:ncol(workdata$X)) {
      workdata$X[, i] = (workdata$X[, i] - workdata$means[i]) / workdata$sd[i]
    }
  }

  return(workdata)
}

#' @importFrom stats model.matrix sd
PrepareDataLinear.B23 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.B23\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }


  data = data.frame(data) # convert to a clean data.frame

  workdata$tags = CreateModelMatrixTags(data)

  if (ncol(data) < 2 | !("numeric" %in% names(workdata$tags))) {
    warning("The data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.")
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

PrepareParamsLinear.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.B2\n\n")
  params$failed         = FALSE
  params$halted         = FALSE
  params$singularMatrix = FALSE

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
  pb$Bcolnames   = params$Bcolnames
  pb$tags        = data$tags

  writeTime = proc.time()[3]
  save(pb, file = file.path(params$writePath, "pb.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pb.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsLinear.B2", 0, 0, writeTime, writeSize)
  return(params)
}


PrepareParamsLinear.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.A2\n\n")

  params$halted          = FALSE
  params$singularMatrix  = FALSE
  params$pmnStepCounter  = 1
  pb                     = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "pb.rdata")) # load pb
  readSize = sum(file.size(file.path(params$readPath, "pb.rdata")))
  readTime = proc.time()[3] - readTime

  if (params$analysis != pb$analysis) {
    params$errorMessage =
      paste("Party A is running", params$analysis, "regression and Party B is running", pb$analysis, "regression.")
    warning(params$errorMessage)
    params$failed = TRUE
    return(params)
  }

  params$n  = nrow(data$X)
  if (pb$n != params$n) {
    params$errorMessage =
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    warning(params$errorMessage)
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
  params$Atags     = data$tags
  params$Btags     = pb$tags

  params$meansA = data$means
  params$sdA    = data$sd
  params$meansB = pb$means
  params$sdB    = pb$sd
  params$yty    = data$yty
  params$meansy = data$meansy
  params$sdy    = data$sdy

  pa        = list()
  pa$p1     = params$p1
  pa$means  = data$means
  pa$sd     = data$sd
  pa$yty    = data$yty
  pa$yname  = data$yname
  pa$Acolnames = params$Acolnames
  writeTime = proc.time()[3]
  save(pa, file = file.path(params$writePath, "pa.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pa.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsLinear.A2", readTime, readSize,
                    writeTime, writeSize)

  return(params)
}


PrepareBlocksLinear.A2 = function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.A2\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  = params$n
  p1 = params$p1
  p2 = params$p2

  minimumBlocksize = GetBlockSize(p1, p2)
  if (n < minimumBlocksize) {
    maxACovariates = trunc(sqrt(p2 * n) - p2 - 1)

    params$errorMessage =
      paste("The minimum secure blocksize of", minimumBlocksize,
            "is larger than the number of observations", paste0(n, ".\n"),
            "Your options are:\n",
            "Increase the number of observations to at least",
            paste0(minimumBlocksize, ".\n"),
            "Decrease the number of A covariates to", maxACovariates, "or less.")

    b = n - 2 * p1 - 2
    discrim = b^2 - 4 * (p1 + 1)^2
    if (discrim >= 0) {
      minBCovariates = trunc(1 + (b - sqrt(discrim)) / 2)
      maxBCovariates = trunc((b + sqrt(discrim)) / 2)
      params$errorMessage =
        paste0(params$errorMessage,
               "\nSet the number of B covariates to be between ", minBCovariates, "and",
               paste0(maxBCovariates, "."))
    }
    warning(params$errorMessage)
    params$failed = TRUE
    params = AddToLog(params, "PrepareBlocksCox.A2", 0, 0, 0, 0)
    return(params)
  }

  if (is.null(blocksize)) {
    blocksize = minimumBlocksize
  }
  if (blocksize < minimumBlocksize) {
    message(paste("Block size of", blocksize,
                  "is too small. Proceeding with minimum blocksize of",
                  paste0(minimumBlocksize, ".")))
    blocksize = minimumBlocksize
  } else if (n < blocksize) {
    message(paste("Block size of", blocksize,
                  "is larger than size of data.  Proceeding with blocksize of",
                  paste0(n, ".")))
  }

  params$blocks    = CreateBlocks(p1, p2, n, blocksize)
  params$container = CreateContainers(p1, p2, params$blocks)
  writeTime = proc.time()[3]
  save(blocksize, file = file.path(params$writePath, "blocksize.rdata"))
  writeSize = file.size(file.path(params$writePath, "blocksize.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareBlocksLinear.A2", 0, 0, writeTime, writeSize)
  return(params)
}


GetZLinear.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLinear.A2\n\n")
  writeTime = 0
  writeSize = 0

  numBlocks = params$blocks$numBlocks
  pbar = MakeProgressBar1(numBlocks, "Z", params$verbose)
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
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params = AddToLog(params, "GetZLinear.A2", 0, 0, writeTime, writeSize)
  return(params)
}


FinalizeParamsLinear.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsLinear.B2\n\n")
  pa = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "pa.rdata")) # read pa
  readSize = sum(file.size(file.path(params$readPath, "pa.rdata")))
  readTime = proc.time()[3] - readTime
  params$p1     = pa$p1
  params$p1.old = params$p1
  params$p      = params$p1 + params$p2
  params$meansA = pa$means
  params$sdA    = pa$sd
  params$yty    = pa$yty
  params$yname  = pa$yname

  params$Acolnames = pa$Acolnames
  params = AddToLog(params, "FinalizeParamsLinear.B2", readTime, readSize, 0, 0)
  return(params)
}


PrepareBlocksLinear.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.B2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  readTime = proc.time()[3]
  load(file.path(params$readPath, "blocksize.rdata")) # load blocksize
  readSize = file.size(file.path(params$readPath, "blocksize.rdata"))
  readTime = proc.time()[3] - readTime
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params = AddToLog(params, "PrepareBlocksLinear.B2", readTime, readSize, 0, 0)
  return(params)
}


GetWLinear.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWLinear.B2\n\n")
  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I-Z*Z')X", params$verbose)

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

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  writeTime = writeTime - proc.time()[3]
  save(XBTXB, file = file.path(params$writePath, "xbtxb.rdata"))
  writeSize = writeSize + file.size(file.path(params$writePath, "xbtxb.rdata"))
  writeTime = writeTime + proc.time()[3]

  params = AddToLog(params, "GetWLinear.B2", readTime, readSize, writeTime, writeSize)

  return(params)
}


GetProductsLinear.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLinear.A2\n\n")
  n  = params$n
  p1 = params$p1
  p2 = params$p2
  XBTXB = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "xbtxb.rdata"))
  readSize = file.size(file.path(params$readPath, "xbtxb.rdata"))
  readTime = proc.time()[3] - readTime

  XATXA = t(data$X) %*% data$X
  XATY  = t(data$X) %*% data$Y
  YTXB  = 0
  XATXB = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'X", params$verbose)

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
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  XTX = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), XBTXB))
  XTY = rbind(XATY, t(YTXB))

  # lasso: x is standardized but needs to be divided by sqrt(n - 1),
  # y is standardized
  XTXLasso = XTX / (n - 1)
  XTYLasso = params$sdy * XTY / sqrt(n - 1)

  params$xtx = XTX
  params$xty = XTY
  params$xtxLasso = XTXLasso
  params$xtyLasso = XTYLasso

  params$converged = TRUE

  params = AddToLog(params, "GetProductsLinear.A2", readTime, readSize, 0, 0)
  return(params)
}

#' @importFrom  stats pf pt
ComputeResultsLinear.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLinear.A2\n\n")
  stats    = params$stats
  stats$converged = params$converged
  stats$failed    = FALSE
  Anames   = params$Acolnames
  Bnames   = params$Bcolnames
  n        = params$n
  yty      = params$yty
  xty      = params$xty
  xtx      = params$xtx
  sdy      = params$sdy
  sdA      = params$sdA
  sdB      = params$sdB
  meansy   = params$meansy
  meansA   = params$meansA
  meansB   = params$meansB

  # First we de-standardize.
  xtx = diag(c(sdA, sdB)) %*% xtx %*% diag(c(sdA, sdB))
  offset = matrix(c(meansA, meansB), ncol = 1) %*%
    matrix(c(meansA, meansB), nrow = 1) * n
  offset[1, 1] = 0
  xtx = xtx + offset

  xty = diag(c(sdA, sdB)) %*% xty * sdy
  offset = n * meansy * matrix(c(meansA, meansB), ncol = 1)
  xty = xty + offset

  # Now, we check for colinearity
  nrow = nrow(xtx)
  indicies = c(1)
  for (i in 2:nrow) {
  	tempIndicies = c(indicies, i)
  	if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
  		indicies = c(indicies, i)
  	}
  }

  Aindex        = which(indicies <= length(Anames))
  AIndiciesKeep = indicies[Aindex]
  BIndiciesKeep = indicies[-Aindex] - length(Anames)

  names.old     = c(Anames, Bnames)
  p             = length(indicies)
  xtx.old       = xtx
  xty.old       = xty
  xtx           = xtx[indicies, indicies, drop = FALSE]
  xty           = matrix(xty[indicies, ], ncol = 1)

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
  stats$party                  = c(rep("dp0", length(Anames)),
                                   rep("dp1", length(Bnames)))
  stats$responseParty          = "dp0"
  stats$coefficients           = rep(NA, params$p)
  stats$tvals                  = rep(NA, params$p)
  stats$secoef                 = rep(NA, params$p)
  stats$pvals                  = rep(NA, params$p)

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
  stats$means                  = c(meansA, meansB)

  names(stats$party)           = names.old
  names(stats$coefficients)    = names.old
  names(stats$secoef)          = names.old
  names(stats$tvals)           = names.old
  names(stats$pvals)           = names.old

  colnames(stats$xtx)          = names.old
  rownames(stats$xtx)          = names.old
  colnames(stats$xty)          = colnames(params$xty)
  rownames(stats$xty)          = names.old

  params$stats = stats
  writeTime = proc.time()[3]
  save(stats, file = file.path(params$writePath, "stats.rdata"))
  writeSize = file.size(file.path(params$writePath, "stats.rdata"))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "ComputeResultsLinear.A2", 0, 0, writeTime, writeSize)
  return(params)
}


GetResultsLinear.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLinear.B2\n\n")
  params$converged = TRUE
  stats = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "stats.rdata"))
  readSize = file.size(file.path(params$readPath, "stats.rdata"))
  readTime = proc.time()[3] - readTime
  params$stats = stats
  params = AddToLog(params, "GetResultsLinear.B2", readTime, readSize, 0, 0)
  return(params)
}





############################## PARENT FUNCTIONS ###############################

PartyAProcess2Linear = function(data,
                                yname                 = NULL,
																monitorFolder         = NULL,
																msreqid               = "v_default_00_0000",
                                blocksize             = NULL,
																sleepTime             = 10,
                                maxWaitingTime        = 24 * 60 * 60,
																popmednet             = TRUE,
																trace                 = FALSE,
																verbose               = TRUE) {

  params = PrepareParams.2p("linear", "A", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)
  Header(params)

  params   = PrepareFolderLinear.A2(params, monitorFolder)
  if (params$failed) {
  	warning(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataLinear.A23(params, data, yname)

  params = PauseContinue.2p(params,  maxWaitingTime)
  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$readPath))
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

  params = PrepareParamsLinear.A2(params, data)

  if (params$failed) {   # Check for failed from PrepareParamsLinear.A2()
    params$completed = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = PrepareBlocksLinear.A2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksCox.A2()
    params$completed = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = GetZLinear.A2(params, data)

  files = c("pa.rdata", "blocksize.rdata",
            SeqZW("cz_", length(params$container$filebreak.Z)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params$completed = TRUE
  params = GetProductsLinear.A2(params, data)
  params = ComputeResultsLinear.A2(params, data)
  files = c("stats.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  SummarizeLog.2p(params)
  return(params$stats)
}

PartyBProcess2Linear = function(data,
																monitorFolder       = NULL,
																sleepTime           = 10,
																maxWaitingTime      = 24 * 60 * 60,
																popmednet           = TRUE,
																trace               = FALSE,
																verbose             = TRUE) {
  params = PrepareParams.2p("linear", "B",
                            popmednet = popmednet, trace = trace,
                            verbose = verbose)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)

  Header(params)
  params   = PrepareFolderLinear.B2(params, monitorFolder)
  if (params$failed) {
  	warning(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataLinear.B23(params, data)

  if (data$failed) { # Check for Error from PrepareDataCox.B2()
    params$completed = TRUE
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$writePath, message)
    files = c("errorMessage.rdata")
    params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params   = PrepareParamsLinear.B2(params, data)

  files = c("pb.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$readPath))
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params = FinalizeParamsLinear.B2(params, data)
  params = PrepareBlocksLinear.B2(params)
  params = GetWLinear.B2(params, data)

  files = c("xbtxb.rdata", SeqZW("cw_", length(params$container$filebreak.W)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = GetResultsLinear.B2(params)
  params$completed = TRUE

  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  return(params$stats)
}
