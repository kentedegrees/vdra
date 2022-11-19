#################### DISTRIBUTED COX REGRESSION FUNCTIONS ####################

PrepareFolderCox.A2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderCox.A2\n\n")
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

	params = AddToLog(params, "PrepareDataCox.A23, PrepareFolderCox.A2", 0, 0, 0, 0)
	return(params)
}


PrepareFolderCox.B2 = function(params, monitorFolder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderCox.B2\n\n")

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

	params = AddToLog(params, "PrepareDataCox.B23, PrepareFolderCox.B2", 0, 0, 0, 0)

	return(params)
}


extractStrata = function(params, data, stratas, mask) {
  strata = list()
  strata$failed = FALSE
  if (!is.null(stratas)) {
    if (!("character" %in% class(stratas))) {
      warning("Strata is not a valid variable name(s).")
      strata$failed = TRUE
      return(strata)
    }
    if (length(stratas) > 0) {
      idx = stratas %in% colnames(data)
      if(!is.null(params$partyName) && params$partyName == "A")  {
        strata$strataFromA = stratas[idx]
        strata$strataFromB = stratas[!idx]
      } else {
        strata$strataFromB = stratas[idx]
        strata$strataFromA = stratas[!idx]
      }
      if (!is.null(params$dataPartnerID) && params$dataPartnerID == "1") {
        strata$strataFromMe     = stratas[idx]
        strata$strataFromOthers = stratas[!idx]
      } else {
        strata$strataFromMe     = stratas[idx]
        strata$strataFromOthers = stratas[!idx]
      }
      strata$strataIndex = which(colnames(data) %in% stratas)
      if (length(strata$strataIndex) > 0) {
        strata$X = data[, strata$strataIndex, drop = FALSE]
        strata$legend = list()
        for (i in 1:ncol(strata$X)) {
          levels = unique(strata$X[, i])
          strata$legend[[colnames(strata$X)[i]]] = levels
          if (mask) {
            levels = sample(levels)
            strata$legend[[colnames(strata$X)[i]]] = rep("NA", length(levels))
          }
          strata$X[, i] = sapply(strata$X[, i], function(x) { which(levels %in% x)})
        }
      }
    }
  }
  return(strata)
}


#' @importFrom stats model.matrix
PrepareDataCox.23 = function(params, data, yname, strata, mask) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataCox.23\n\n")

  workdata = list()
  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  workdata$strata = extractStrata(params, data, strata, mask)
  if (workdata$strata$failed) {
    workdata$failed = TRUE
    return(workdata)
  }
  strataIndex = workdata$strata$strataIndex
  responseIndex = numeric()

  if (params$partyName == "A") {
    responseIndex = CheckResponse(params, data, yname)

    if (is.null(responseIndex)) {
      workdata$failed = TRUE
      return(workdata)
    }

    workdata$survival        = list()
    workdata$survival$rank   = data[, responseIndex[1]]
    workdata$survival$status = data[, responseIndex[2]]
    if (length(intersect(strataIndex, responseIndex)) > 0) {
      warning("Response and strata share a variable.")
      workdata$failed = TRUE
      return(workdata)
    }
  }

  covariateIndex = setdiff(1:ncol(data), union(strataIndex, responseIndex))

  if (length(covariateIndex) == 0) {
    if (params$partyName == "A") {
      workdata$X = matrix(0, nrow = nrow(data), ncol = 0)
    } else {
      warning("After removing strata, data is empty.  Party B must supply at least one non-strata covariate.")
      workdata$failed = TRUE
      return(workdata)
    }
  } else {
    workdata$tags = CreateModelMatrixTags(data[, covariateIndex, drop = FALSE])
    if (params$partyName == "B" & (ncol(data) < 2 | !("numeric" %in% names(workdata$tags)))) {
      warning("The data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.")
      workdata$failed = TRUE
      return(workdata)
    }
    workdata$X = scale(model.matrix(~ ., data[, covariateIndex, drop = FALSE]),
                       center = TRUE, scale = FALSE)
    workdata$X = workdata$X[, -1, drop = FALSE]
  }

  return(workdata)
}


PrepareParamsCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.B2\n\n")
  params$failed         = FALSE
  params$converged      = FALSE
  params$halted         = FALSE
  params$singularMatrix = FALSE

  params$n = nrow(data$X)
  params$numEvents = 0
  params$p1 = 0
  params$p2 = ncol(data$X)
  params$p  = params$p1 + params$p2
  params$p1.old = params$p1
  params$p2.old = params$p2
  params$Acolnames = c("")
  params$Bcolnames = colnames(data$X)
  params$Acolnames.old = c("")
  params$Bcolnames.old = c("")
  params$cutoff        = 1
  params$maxIterations = 1

  params$survivalInstalled = requireNamespace("survival", quietly = TRUE)
  if (params$survivalInstalled & !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  pb           = list()
  pb$p2        = params$p2
  pb$n         = params$n
  pb$analysis  = params$analysis
  pb$Bcolnames = params$Bcolnames
  pb$strataB   = list()
  pb$strataB$strataFromA = data$strata$strataFromA
  pb$strataB$strataFromB = data$strata$strataFromB
  pb$tags      = data$tags

  writeTime = proc.time()[3]
  save(pb, file = file.path(params$writePath, "pb.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pb.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsCox.B2", 0, 0, writeTime, writeSize)
  return(params)
}

CheckStrataCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckStrataCox.A2\n\n")
  pb = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "pb.rdata")) # load pb
  readSize = sum(file.size(file.path(params$readPath, "pb.rdata")))
  readTime = proc.time()[3] - readTime

  if (length(pb$strataB$strataFromA) == length(data$strata$strataFromA) &&
      length(pb$strataB$strataFromB) == length(data$strata$strataFromB) &&
      ifelse(length(pb$strataB$strataFromA) == 0, TRUE,
             order(pb$strataB$strataFromA) == order(data$strata$strataFromA)) &&
      ifelse(length(pb$strataB$strataFromB) == 0, TRUE,
             order(pb$strataB$strataFromB) == order(data$strata$strataFromB))) {
    params$getStrataFromB = length(data$strata$strataFromB) > 0
  } else {
    params$getStrataFromB = FALSE
    AcapB = intersect(data$strata$strataFromA, pb$strataB$strataFromB)
    BcapA = intersect(data$strata$strataFromB, pb$strataB$strataFromA)
    if (length(AcapB) > 0) {
      params$errorMessage =
        paste("Party A and Party B have", length(AcapB), "variable(s) with the same name which are used in the strata.",
              "These variable(s) are <", paste0(AcapB, collapse = ", "), ">.",
              "Make sure the variables from each party have distinct names.")
    }
    else if (length(BcapA) > 0) {
      params$errorMessage =
        paste("Party A and Party B have specified", length(BcapA), "variable(s) for the strata which are not found in the data.",
              "These variable(s) are <", paste0(BcapA, collapse = ", "), ">.",
              "Check the spelling of the variables names and / or remove them from the strata.")
    } else {
      params$errorMessage =
        paste("Party A and Party B have specified different strata.",
              "Verify that both parties specify the same strata.")
    }
    warning(params$errorMessage)
    params$failed = TRUE
  }
  params = AddToLog(params, "CheckStrata.A2", readTime, readSize, 0, 0)
  return(params)
}


SendStrataCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataCox.B2\n\n")
  strataTemp = data$strata
  writeTime = proc.time()[3]
  save(strataTemp, file = file.path(params$writePath, "strata.rdata"))
  writeSize = file.size(file.path(params$writePath, "strata.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "SendStrataCox.B2", 0, 0, writeTime, writeSize)
  return(params)
}


PrepareStrataCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareStrataCox.A2\n\n")
  # I will need to update both data$survival and params$logs, etc in this function.
  readTime = 0
  readSize = 0
  writeTime = 0
  writeSize = 0
  strataTemp = NULL

  if (params$getStrataFromB) {
    readTime = proc.time()[3]
    load(file = file.path(params$readPath, "strata.rdata"))
    readSize = file.size(file.path(params$readPath, "strata.rdata"))
    readTime = proc.time()[3] - readTime
  }
  if (length(data$strata$strataFromA) == 0 && length(data$strata$strataFromB) == 0) {
    strataTemp$X = data.frame(const__ = rep(1, params$n))
    strataTemp$legend = FALSE
  } else if (length(data$strata$strataFromA) > 0 && length(data$strata$strataFromB) == 0) {
    strataTemp$X = data$strata$X
    strataTemp$legend = data$strata$legend
  } else if (length(data$strata$strataFromA) > 0 && length(data$strata$strataFromB) > 0) {
    strataTemp$X = cbind(data$strata$X, strataTemp$X)
    strataTemp$legend = c(data$strata$legend, strataTemp$legend)
  }

  sorted = do.call("order", cbind(strataTemp$X, data$survival$rank, data$survival$status))
  strataTemp$X = strataTemp$X[sorted, , drop = FALSE]
  data$survival$rank   = data$survival$rank[sorted]
  data$survival$status = data$survival$status[sorted]
  data$X = data$X[sorted, , drop = FALSE]
  data$survival$sorted = sorted
  ranks = which(apply(abs(apply(strataTemp$X, 2, diff)), 1, sum) > 0)
  ranks = c(ranks, nrow(strataTemp$X))
  names(ranks) = NULL
  strata = rep(list(list()), length(ranks))
  if (length(ranks) == 1 && colnames(strataTemp$X)[1] == "const__") {
    strata[[1]]$start = 1
    strata[[1]]$end   = as.integer(nrow(data$X))
    strata[[1]]$label = ""
  } else {
    start = 1
    for (i in  1:length(ranks)) {
      strata[[i]]$start = start
      strata[[i]]$end   = as.integer(ranks[i])
      label = ""
      for (j in 1:ncol(strataTemp$X)) {
        temp = colnames(strataTemp$X)[j]
        label = paste0(label, temp, "=", strataTemp$legend[[temp]][strataTemp$X[start, j]])
        if (j < ncol(strataTemp$X)) {
          label = paste0(label, ", ")
        }
      }
      strata[[i]]$label = label
      start = as.numeric(ranks[i]) + 1
    }
  }

  emptyStrata = c()
  data$X = cbind(matrix(0, nrow = nrow(data$X), ncol = length(strata)), data$X)
  for (i in 1:length(strata)) {
    idx = strata[[i]]$start:strata[[i]]$end
    data$X[idx, i] = 1
    temp  = table(data$survival$rank[idx])
    M = length(temp)   # number of unique observed times, including where no one fails
    # Count the number of 0's and 1's for each observed time
    temp0 = table(data$survival$rank[idx], data$survival$status[idx])
    # Check if there are all 1's or all 0's .  If so, add them into the table.
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 = cbind(temp0, 0)
        colnames(temp0) = c("0", "1")
      } else {
        temp0 = cbind(0, temp0)
        colnames(temp0) = c("0", "1")
      }
    }
    strata[[i]]$J = as.integer(sum(temp0[, 2] > 0))         # number of distinct failure times
    if (strata[[i]]$J == 0) {
      emptyStrata = c(emptyStrata, i)
    }
    # The number of failures at each rank which has a failure
    strata[[i]]$nfails = as.numeric(temp0[which(temp0[, 2] > 0), 2])
    # The first index of the ranks for which the number of failures is > 0.
    strata[[i]]$start0 = c(1, (cumsum(temp)[1:(M - 1)] + 1))[which(temp0[, 2] > 0)]
    # The first index of a failure for each rank which has a failure
    strata[[i]]$start1 = strata[[i]]$start0 + temp0[which(temp0[, 2] > 0), 1]
    # The last index of a failure for each rank which has a failure
    strata[[i]]$stop1  = as.numeric(strata[[i]]$start1 + strata[[i]]$nfails  - 1)
    strata[[i]]$start0 = as.numeric(strata[[i]]$start0 + strata[[i]]$start - 1)
    strata[[i]]$start1 = as.numeric(strata[[i]]$start1 + strata[[i]]$start - 1)
    strata[[i]]$stop1  = as.numeric(strata[[i]]$stop1  + strata[[i]]$start - 1)
  }

  data$survival$strata = strata

  writeTime = proc.time()[3]
  survival = data$survival
  save(survival, file = file.path(params$writePath, "survival.rdata"))
  writeSize = file.size(file.path(params$writePath, "survival.rdata"))
  writeTime = proc.time()[3] - writeTime
  data$readTime = readTime
  data$readSize = readSize
  data$writeTime = writeTime
  data$writeSize = writeSize

  return(data)
}


PrepareParamsCox.A2 = function(params, data, cutoff = 0.01, maxIterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.A2\n\n")
  params$converged       = FALSE
  params$halted          = FALSE
  params$singularMatrix  = FALSE
  params$pmnStepCounter  = 1
  pb = NULL

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
  params$numEvents = sum(data$survival$status)
  if (pb$n != params$n) {
    params$errorMessage =
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    warning(params$errorMessage)
    params$failed = TRUE
    return(params)
  }

  params$p1 = ncol(data$X)
  params$p2 = pb$p2
  params$p  = params$p1 + params$p2
  params$p1.old = params$p1
  params$p2.old = params$p2

  params$Acolnames = colnames(data$X)
  params$Bcolnames = pb$Bcolnames
  params$Acolnames.old = c("")
  params$Bcolnames.old = c("")

  params$Atags         = data$tags
  params$Btags         = pb$tags

  if (cutoff <= 0) cutoff = 0.01
  if (cutoff >= 1) cutoff = 0.05
  params$cutoff        = cutoff

  if (maxIterations < 1) maxIterations = 1
  params$maxIterations = maxIterations

  params$survivalInstalled = requireNamespace("survival", quietly = TRUE)
  if (params$survivalInstalled & !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  # if (requireNamespace("survival", quietly = TRUE)) {
  #   params$survivalInstalled = TRUE
  # } else {
  #   params$survivalInstalled = FALSE
  # }

  pa = list()
  pa$p1 = params$p1
  pa$cutoff = params$cutoff
  pa$maxIterations = params$maxIterations
  pa$Acolnames = params$Acolnames
  pa$tags = data$tags

  writeTime = proc.time()[3]
  save(pa, file = file.path(params$writePath, "pa.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "pa.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsCox.A2", readTime, readSize,
                    writeTime, writeSize)
  return(params)
}


PrepareBlocksCox.A2 = function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksCox.A2\n\n")
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
  params = AddToLog(params, "PrepareBlocksCox.A2", 0, 0, writeTime, writeSize)
  return(params)
}


GetZCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZCox.A2\n\n")
  writeTime = 0
  writeSize = 0

  numBlocks = params$blocks$numBlocks
  pbar = MakeProgressBar1(numBlocks, "Z Files", params$verbose)
  containerCt.Z = 0
  for (i in 1:numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename = paste0("cz_", containerCt.Z, ".rdata")
      toWrite = file(file.path(params$writePath, filename), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    g = params$blocks$g[i]
    Z = FindOrthogonalVectors(data$X[strt:stp, ], g)
    writeTime = writeTime - proc.time()[3]
    writeBin(as.vector(Z), con = toWrite, endian = "little")
    writeTime = writeTime + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
      close(toWrite)
      writeSize = writeSize + file.size(file.path(params$writePath, filename))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params = AddToLog(params, "GetZCox.A2", 0, 0, writeTime, writeSize)
  return(params)
}


SortDataCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SortDataCox.B2\n\n")
  survival = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "survival.rdata"))
  readSize = file.size(file.path(params$readPath, "survival.rdata"))
  readTime = proc.time()[3] - readTime
  data$X = data$X[survival$sorted, , drop = FALSE]
  data$survival = survival
  data$readTime = readTime
  data$readSize = readSize
  return(data)
}


FinalizeParamsCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsCox.B2\n\n")
  pa = NULL
  readTime = proc.time()[3]
  load(file.path(params$readPath, "pa.rdata")) # read pa
  readSize = sum(file.size(file.path(params$readPath, "pa.rdata")))
  readTime             = proc.time()[3] - readTime
  params$numEvents     = sum(data$survival$status)
  params$p1            = pa$p1
  params$cutoff        = pa$cutoff
  params$maxIterations = pa$maxIterations
  params$p             = params$p1 + params$p2
  params$Acolnames     = pa$Acolnames
  params = AddToLog(params, "FinalizeParamsCox.B2", readTime, readSize, 0, 0)
  return(params)
}


PrepareBlocksCox.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksCox.B2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  readTime = proc.time()[3]
  load(file.path(params$readPath, "blocksize.rdata")) # load blocksize
  readSize = file.size(file.path(params$readPath, "blocksize.rdata"))
  readTime = proc.time()[3] - readTime
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params = AddToLog(params, "PrepareBlocksCox.B2", readTime, readSize, 0, 0)
  return(params)
}


GetWCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWcox.B2\n\n")
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

  params = AddToLog(params, "GetWCox.B2", readTime, readSize, writeTime, writeSize)

  return(params)
}


CheckColinearityCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityCox.A2\n\n")
  p2 = params$p2
  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0
  XBTXB = NULL

  readTime = readTime - proc.time()[3]
  load(file.path(params$readPath, "xbtxb.rdata")) # load XBTXB
  readSize = file.size(file.path(params$readPath, "xbtxb.rdata"))
  readTime = readTime + proc.time()[3]
  XATXA = t(data$X) %*% data$X
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

    if ((i + 1) %in% params$container$filebreak.W || i == params$blocks$numBlocks) {
      close(toRead)
      readSize = readSize + file.size(file.path(params$readPath, filename))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  XTX = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), XBTXB))

  nrow = nrow(XTX)
  indicies = 1:length(data$survival$strata)
  for (i in (1 + length(data$survival$strata)):nrow) {
  	tempIndicies = c(indicies, i)
  	if (rcond(XTX[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
  		indicies = c(indicies, i)
  	}
  }

  numStrata = length(data$survival$strata)
  Anames = params$Acolnames
  Bnames = params$Bcolnames
  indicies = indicies[-(1:length(data$survival$strata))] - numStrata # Get rid of the strata indicators
  Aindex = which(indicies <= length(Anames))
  Bindex = which(indicies > length(Anames))
  params$AIndiciesKeep = indicies[Aindex]
  params$BIndiciesKeep = indicies[Bindex] - length(Anames)
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

  Aindicies = params$AIndiciesKeep
  Bindicies = params$BIndiciesKeep


  writeTime = writeTime - proc.time()[3]
  save(Aindicies, Bindicies, file = file.path(params$writePath, "indicies.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "indicies.rdata")))
  writeTime = writeTime + proc.time()[3]
  tags = params$Btags[Bindicies]

  if (length(unique(tags)) == 0) {
    params$failed = TRUE
    params$errorMessage = "After removing colinear covariates, Party B has no covariates."
  }
  # if (length(unique(tags)) < 2) {
  #   params$failed = TRUE
  #   params$errorMessage = "After removing colinear covariates, Party B has 1 or fewer covariates.\n\n"
  # } else if (!("numeric" %in% names(tags))) {
  #   params$failed = TRUE
  #   params$errorMessage = "After removing colinear covariates, Party B has no continuous covariates.\n\n"
  # }

#
#   if (params$p2 == 0) {
#     params$failed = TRUE
#     params$errorMessage = "All of party B's covariates are either linear or are colinear with Party A's covariates."
#     warning("All of Party B's covariates are either linear or are colinear with Party A's covariates.")
#   }
  params = AddToLog(params, "CheckColinearityCox.A2", readTime, readSize, writeTime, writeSize)
  return(params)
}


UpdateDataCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.A2\n\n")
  numStrata = length(data$survival$strata)
  data$X = as.matrix(data$X[, params$AIndiciesKeep + numStrata, drop = FALSE])
  return(data)
}


PrepareLoopCox.A2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareLoopCox.A2\n\n")
  params$betasA    = matrix(0, params$p1, 1)
  params$betasAold = matrix(0, params$p1, 1)
  params$algIterationCounter      = 1
  params$deltabeta = Inf
  params = AddToLog(params, "PrepareLoopCox.A2", 0, 0, 0, 0)
  return(params)
}


ComputeLogLikelihoodCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeLogLikelihoodCox.A2\n\n")
  n  = params$n
  p1 = params$p1
  readTime = 0
  readSize = 0
  XB.betasB = NULL

  if (params$algIterationCounter == 1) {
    params$X.betas.old = matrix(0, n, 1)
    X.betas = matrix(0, n, 1)
    params$loglikelihood.old = -Inf
  } else {
    readTime = -proc.time()[3]
    load(file.path(params$readPath, "XB_betasB.rdata")) # load XB.betasB
    readSize = file.size(file.path(params$readPath, "XB_betasB.rdata"))
    readTime = readTime + proc.time()[3]
    params$X.betas.old = params$X.betas
    X.betas = params$XA.betasA + XB.betasB
    params$loglikelihood.old = params$loglikelihood
  }

  numEvents = sum(data$survival$status)

  stepSize = 1
  w = exp(X.betas)

  while (max(w) == Inf) {
    X.betas = (X.betas + params$X.betas.old) * 0.5
    stepSize = stepSize * 0.5
    w = exp(X.betas)
  }

  computeLoglikelihood = TRUE

  while (computeLoglikelihood) {
    stepCounter = 0
    pbar = MakeProgressBar1(numEvents, "Loglikelihood", params$verbose)
    loglikelihood = 0
    for (i in 1:length(data$survival$strata)) {                    ##!
      if (data$survival$strata[[i]]$J > 0) {                       ##!
        for (j in 1:data$survival$strata[[i]]$J) {                 ##!
          nj = data$survival$strata[[i]]$nfails[j]                 ##!
          yIndex = data$survival$strata[[i]]$start0[j]:data$survival$strata[[i]]$end      ##!
          zIndex = data$survival$strata[[i]]$start1[j]:data$survival$strata[[i]]$stop1[j] ##!
          Aj1 = sum(w[yIndex])
          Aj2 = sum(w[zIndex]) / nj
          loglikelihood = loglikelihood + sum(log(w[zIndex]))
          for (r in 0:(nj - 1)) {
            Ajr = Aj1 - r * Aj2
            loglikelihood = loglikelihood - log(Ajr)
          }
          stepCounter = stepCounter + nj
          pbar = MakeProgressBar2(stepCounter, pbar, params$verbose)
        }
      }
    }
    if (loglikelihood > params$loglikelihood.old || stepSize < 0.5^6) {
      computeLoglikelihood = FALSE
    } else {
      if (params$verbose) cat("Step Halving\n\n")
      X.betas = (X.betas + params$X.betas.old) * 0.5
      stepSize = stepSize * 0.5
      w = exp(X.betas)
    }
  }

  deltal = as.numeric(data$survival$status)
  deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  W.XA = matrix(0, n, p1)
  stepCounter = 0

  .Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XA,
  			as.integer(n), as.integer(p1), as.integer(numEvents),
  			as.integer(params$verbose))

  params$loglikelihood = loglikelihood
  params$deltal = deltal
  params$tXA.W.XA = t(data$X) %*% W.XA
  params$tXA.deltal = t(data$X) %*% deltal

  if (params$algIterationCounter == 1) {
    params$nullLoglikelihood = loglikelihood
    params$nullScore = params$tXA.deltal
  }
  params$X.betas = X.betas
  params$stepSize = stepSize

  writeTime = proc.time()[3]
  save(X.betas, stepSize, file = file.path(params$writePath, "X_betas_ss.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, "X_betas_ss.rdata")))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "ComputeLogLikelihoodCox.A2", readTime, readSize,
                    writeTime, writeSize)
  return(params)
}


UpdateParamsCox.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.B2\n\n")
  Aindicies = NULL
  Bindicies = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "indicies.rdata")) # load Aindicies, Bindicies
  readSize = sum(file.size(file.path(params$readPath, "indicies.rdata")))
  readTime = proc.time()[3] - readTime
  params$Acolnames.old = params$Acolnames
  params$Bcolnames.old = params$Bcolnames
  params$Acolnames     = params$Acolnames.old[Aindicies]
  params$Bcolnames     = params$Bcolnames.old[Bindicies]
  params$p1.old = params$p1
  params$p2.old = params$p2
  params$p1     = length(Aindicies)
  params$p2     = length(Bindicies)
  params$p.old  = params$p - 1
  params$p      = params$p1 + params$p2
  params$BIndiciesKeep = Bindicies
  params$AIndiciesKeep = Aindicies
  params$betasB    = matrix(0, params$p2, 1)
  params$betasBold = matrix(0, params$p2, 1)
  params = AddToLog(params, "UpdateParamsCox.B2", readTime, readSize, 0, 0)
  return(params)
}

UpdateDataCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.B2\n\n")
  data$X = as.matrix(data$X[, params$BIndiciesKeep, drop = FALSE])
  return(data)
}


ComputeLogLikelihoodCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeLogLikelihoodCox.B2\n\n")
  n   = params$n
  p2  = params$p2
  stepSize = NULL
  X.betas  = NULL

  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0

  readTime = readTime - proc.time()[3]
  load(file.path(params$readPath, "X_betas_ss.rdata")) # load X.betas, stepSize
  readSize = sum(file.size(file.path(params$readPath, "X_betas_ss.rdata")))
  readTime = readTime + proc.time()[3]

  params$stepSize = stepSize

  w = exp(X.betas)
  deltal = as.numeric(data$survival$status)
  deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  numEvents = sum(data$survival$status)
  W.XB = matrix(0, n, p2)

  .Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XB,
  			as.integer(n), as.integer(p2), as.integer(numEvents),
  			as.integer(params$verbose))

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I-Z*Z')WX", params$verbose)
  containerCt.Z = 0
  containerCt.Cox = 0

  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename1 = paste0("cz_", containerCt.Z, ".rdata")
      toRead = file(file.path(params$readPath, filename1), "rb")
    }
    if (i %in% params$container$filebreak.Cox) {
      containerCt.Cox = containerCt.Cox + 1
      filename2 = paste0("cCox_", containerCt.Cox, ".rdata")
      toWrite = file(file.path(params$writePath, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1
    g = params$blocks$g[i]

    readTime = readTime - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n2 * g,
                       endian = "little"), nrow = n2, ncol = g)
    readTime = readTime + proc.time()[3]

    IZ.tZ.W.XBtemp = W.XB[strt:stp, , drop = FALSE] - Z %*% (t(Z) %*% W.XB[strt:stp, , drop = FALSE])

    writeTime = writeTime - proc.time()[3]
    writeBin(as.vector(IZ.tZ.W.XBtemp), con = toWrite, endian = "little")
    writeTime = writeTime + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.Z || i == params$blocks$numBlocks) {
      close(toRead)
      readSize = readSize + file.size(file.path(params$readPath, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.Cox || i == params$blocks$numBlocks) {
      close(toWrite)
      writeSize = writeSize + file.size(file.path(params$writePath, filename2))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params$deltal = deltal
  params$tXB.W.XB   = t(data$X) %*% W.XB
  params$tXB.deltal = t(data$X) %*% deltal
  if (params$algIterationCounter == 1) {
    params$nullScore = params$tXB.deltal
  }

  tXB.W.XB = params$tXB.W.XB
  writeTime = writeTime - proc.time()[3]
  save(tXB.W.XB, file = file.path(params$writePath, "tXB_W_XB.rdata"))
  writeSize = writeSize + file.size(file.path(params$writePath, "tXB_W_XB.rdata"))
  writeTime = writeTime + proc.time()[3]
  params = AddToLog(params, "ComputeLogLikelihoodCox.B2", readTime, readSize, writeTime, writeSize)
  return(params)
}


ComputeInverseCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInverseCox.A2\n\n")
  p1 = params$p1
  p2 = params$p2
  readTime  = 0
  readSize  = 0
  writeTime = 0
  writeSize = 0
  tXB.W.XB = 0
  tXA.W.XB = 0

  readTime = readTime - proc.time()[3]
  load(file.path(params$readPath, "tXB_W_XB.rdata")) # load tXB.W.XB
  readSize = readSize + file.size(file.path(params$readPath, "tXB_W_XB.rdata"))
  readTime = readTime + proc.time()[3]

  containerCt.Cox = 0

  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Cox) {
      containerCt.Cox = containerCt.Cox + 1
      filename = paste0("cCox_", containerCt.Cox, ".rdata")
      toRead = file(file.path(params$readPath, filename), "rb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1

    readTime = readTime - proc.time()[3]
    IZ.tZ.W.XB = matrix(readBin(con = toRead, what = numeric(), n = n2 * p2,
                        endian = "little"), nrow = n2, ncol = p2)
    readTime = readTime + proc.time()[3]

    tXA.W.XB = tXA.W.XB + t(data$X[strt:stp, ]) %*% IZ.tZ.W.XB
    if ((i + 1) %in% params$container$filebreak.Cox || i == params$blocks$numBlocks) {
      close(toRead)
      readSize = readSize + file.size(file.path(params$readPath, filename))
    }
  }

  M = rbind(cbind(params$tXA.W.XA, tXA.W.XB),
            cbind(t(tXA.W.XB), tXB.W.XB))

  if (params$algIterationCounter == 1) {
    params$nullHessian = M
  }
  params$XtWX = M


 	inv = NULL
 	tryCatch({inv = solve(M)},
 					 error = function(err) { inv = NULL } )
	M = inv
	if (is.null(M)) {
		params$failed = TRUE
    params$errorMessage = "The matrix XWX is singular.  This is probably due to divergence of the coefficients."
    warning(params$errorMessage)

    betas = rep(NA, length(params$Acolnames.old))
    betas[params$AIndiciesKeep] = params$betasA
    betas = data.frame(betas)
    rownames(betas) = params$Acolnames.old
    params = AddToLog(params, "ComputeInverseCox.A2", readTime, readSize, writeTime, writeSize)
    return(params)
  }
  M3 = M[(p1 + 1):(p1 + p2), 1:p1]
  params$M = M
  params$M3.tXA.deltal = M3 %*% params$tXA.deltal
  M3.tXA.deltal = params$M3.tXA.deltal
  writeTime = writeTime - proc.time()[3]
  save(M, M3.tXA.deltal, file = file.path(params$writePath, "M.rdata"))
  writeSize = writeSize + sum(file.size(file.path(params$writePath, "M.rdata")))
  writeTime = writeTime + proc.time()[3]
  params = AddToLog(params, "ComputeInverseCox.A2", readTime, readSize, writeTime, writeSize)

  return(params)
}


ComputeBetaCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeBetaCox.B2\n\n")
  p1 = params$p1
  p2 = params$p2
  M  = 0
  M3.tXA.deltal = 0

  readTime = proc.time()[3]
  load(file.path(params$readPath, "M.rdata")) # load M, M3.tXA.deltal
  readSize = sum(file.size(file.path(params$readPath, "M.rdata")))
  readTime = proc.time()[3] - readTime

  if (params$stepSize < 1) {  # Is this in the wrong spot?
    params$betasB = params$betasBold + (params$betasB - params$betasBold) * params$stepSize
  }

  params$betasBold = params$betasB
  M4 = M[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
  M2 = M[1:p1, (p1 + 1):(p1 + p2)]

  params$M2.tXB.deltal = M2 %*% params$tXB.deltal
  params$betasB = params$betasB + M3.tXA.deltal + M4 %*% params$tXB.deltal

  params$XB.betasB = data$X %*% params$betasB

  M2.tXB.deltal = params$M2.tXB.deltal
  XB.betasB = params$XB.betasB
  writeTime = proc.time()[3]
  save(M2.tXB.deltal, file = file.path(params$writePath, "M2_tXB_deltal.rdata"))
  save(XB.betasB, file = file.path(params$writePath, "XB_betasB.rdata"))
  writeSize = sum(file.size(file.path(params$writePath, c("M2_tXB_deltal.rdata",
                                                          "XB_betasB.rdata"))))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "ComputeBetaCox.B2", readTime, readSize, writeTime, writeSize)

  return(params)
}


ComputeBetaCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeBetaCox.A2\n\n")
  p1 = params$p1
  M2.tXB.deltal = 0

  readTime = proc.time()[3]
  load(file.path(params$readPath, "M2_tXB_deltal.rdata")) # load M2.tXB.deltal
  readSize = sum(file.size(file.path(params$readPath, "M2_tXB_deltal.rdata")))
  readTime = proc.time()[3] - readTime

  if (params$stepSize < 1) { # Is this in the wrong spot?
    params$betasA = params$betasAold + (params$betasA - params$betasAold) * params$stepSize
  }

  params$betasAold = params$betasA
  params$betasA = params$betasA + params$M[1:p1, 1:p1] %*% params$tXA.deltal + M2.tXB.deltal

  converged = abs(params$loglikelihood - params$loglikelihood.old) /
              (abs(params$loglikelihood) + 0.1) < params$cutoff
  params$converged = converged

  params$XA.betasA = data$X %*% params$betasA

  if (params$algIterationCounter >= params$maxIterations) {
    params$halted = TRUE
    warning(paste("Failed to converged in", params$maxIterations, "iterations."))
  }

  writeTime = proc.time()[3]
  save(converged, file = file.path(params$writePath, "converged.rdata"))
  writeSize = file.size(file.path(params$writePath, "converged.rdata"))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "ComputeBetaCox.A2", readTime, readSize, writeTime, writeSize)
  return(params)
}


GetConvergedStatusCox.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetConvergedStatusCox.B2\n\n")
  converged = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "converged.rdata")) # load deltabeta
  readSize = file.size(file.path(params$readPath, "converged.rdata"))
  readTime = proc.time()[3] - readTime
  params$converged = converged
  if (params$algIterationCounter > params$maxIterations) {
    params$halted = TRUE
    warning(paste("Failed to converged in", params$maxIterations, "iterations."))
  }
  writeTime = 0
  writeSize = 0
  if (params$converged || params$halted) {
    betasB = params$betasB
    nullScoreB = params$nullScore
    writeTime = proc.time()[3]
    save(betasB, nullScoreB, file = file.path(params$writePath, "B_betas_ns.rdata"))
    writeSize = sum(file.size(file.path(params$writePath, "B_betas_ns.rdata")))
    writeTime = proc.time()[3] - writeTime
  }
  params = AddToLog(params, "GetConvergedStatusCox.B2", readTime, readSize, writeTime, writeSize)
  return(params)
}


SurvFitCox.A2 = function(params, survival, pred) {
  if (params$trace) cat(as.character(Sys.time()), "SurvFitCox.A2\n\n")
  surv = rep(1, length(survival$rank))
  for (i in 1:length(survival$strata)) {
    if (survival$strata[[i]]$J > 0) {
      start   = survival$strata[[i]]$start
      end     = survival$strata[[i]]$end
      risk    = exp(pred[start:end])
      dtime   = survival$rank[start:end]
      status  = survival$status[start:end]
      death   = status == 1                                  # times where death happened
      time    = sort(unique(dtime))                          # (A) get unique event times
      rcumsum = function(x) rev(cumsum(rev(x)))
      nevent  = as.vector(rowsum(as.numeric(death), dtime))  # (A) Count the number of deaths at each event time
      ndeath  = rowsum(status, dtime)                        # (A) number of deaths at each unique event time
      nrisk   = rcumsum(rowsum(risk, dtime))                 # (A) rowsum = sum of risk at each time, then reverse cum sum, sorted by time
      erisk   = rowsum(risk * death, dtime)                  # (A) risk score sums of death at each unique event time
      n       = length(nevent)
      sum1 = double(n)   # a vector of 0's, length number of unique event times
      for (i in 1:n) {
        d = ndeath[i];
        if (d == 1) {
          sum1[i] = 1 / nrisk[i];
        } else if (d > 1) {
          for (j in 0:(d - 1)) {
            sum1[i] = sum1[i] + 1 / (d * nrisk[i] - erisk[i] * j)
          }
        }
      }
      temp = exp(-cumsum(nevent * sum1))
      for (i in start:end) {
        surv[i] = temp[which(time == survival$rank[i])]
      }
    }
  }
  return(surv)
}

#' @importFrom  stats pchisq pnorm qnorm
ComputeResultsCox.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsCox.A2\n\n")
  stats = params$stats
  stats$converged = params$converged
  stats$failed    = FALSE
  readTime = 0
  readSize = 0
  betasB     = NULL
  nullScoreB = NULL
  XB.betasB  = NULL

  readTime = proc.time()[3]
  load(file.path(params$readPath, "B_betas_ns.rdata")) # load betasB, nullscoreB
  load(file.path(params$readPath, "XB_betasB.rdata")) # load XB.betasB
  readSize = sum(file.size(file.path(params$readPath, c("B_betas_ns.rdata",
                                                        "XB_betasB.rdata"))))
  readTime = proc.time()[3] - readTime
  params$betasB = betasB
  params$nullScore = rbind(params$nullScore, nullScoreB)

  names.new          = c(params$Acolnames, params$Bcolnames)
  names.old          = c(params$Acolnames.old, params$Bcolnames.old)
  idxA               = params$AIndiciesKeep
  idxB               = params$BIndiciesKeep
  idx                = c(idxA, idxB + length(params$Acolnames.old))
  stats$party        = c(rep("dp0", length(params$Acolnames.old)),
                         rep("dp1", length(params$Bcolnames.old)))
  stats$coefficients = rep(NA, length(stats$party))
  tempcoefs          = c(params$betasA, params$betasB)
  stats$coefficients[idx] = tempcoefs
  stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
  stats$expncoef     = exp(-stats$coefficients)
  tempvar            = solve(params$XtWX)
  stats$var          = matrix(0, length(names.old), length(names.old))
  stats$var[idx, idx] = tempvar
  stats$secoef       = rep(NA, length(names.old))
  stats$secoef[idx]  = sqrt(diag(tempvar))  # se(coef)

  stats$zvals        = stats$coefficients / stats$secoef  # z values
  stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE )   # pvals
  stats$stars        = matrix(sapply(stats$pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  }))
  stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  stats$loglik       = c(params$nullLoglikelihood, params$loglikelihood)
  stats$n            = params$n
  stats$nevent       = params$numEvents
  stats$iter         = params$algIterationCounter - 1
  stats$df           = params$p
  stats$score        = t(params$nullScore) %*% solve(params$nullHessian) %*% params$nullScore
  stats$score        = c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       = "efron"
  stats$lrt          = 2*(stats$loglik[2] - stats$loglik[1])
  stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      = c(1 - exp(-stats$lrt[1]/stats$n),
                         1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    = t(tempcoefs) %*% params$XtWX %*% tempcoefs
  stats$wald.test    = c(stats$wald.test,
                         1 - pchisq(stats$wald.test, stats$df))
  pred = -params$XA.betasA - XB.betasB
  if (params$survivalInstalled) {
    surv = survival::Surv(data$survival$rank, data$survival$status)
    strat = rep(0, length(surv))
    for (i in 1:length(data$survival$strata)) {
      strat[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] = i
    }
    results = survival::concordance(surv~pred + strata(strat))
    if (class(results$stats) == "matrix") {  # more than one strata
      stats$concordance = c(apply(results$count, 2, sum)[1:4], results$concordance, sqrt(results$var))
    } else {                                 # only one strata, so a numeric vector
      stats$concordance = c(results$count[1:4], results$concordance, sqrt(results$var))
    }
  } else {
    stats$concordance = c(NA, NA, NA, NA, NA, NA)
  }

  stats$survival = data.frame(
    rank   = data$survival$rank,
    status = data$survival$status,
    sorted = data$survival$sorted,
    surv   = SurvFitCox.A2(params, data$survival, pred)
  )

  stats$strata = as.data.frame(matrix(0, length(data$survival$strata), 3))
  stats$strata$label = ""
  colnames(stats$strata) = c("start", "end", "events", "label")
  for (i in 1:length(data$survival$strata)) {
    stats$strata$start[i]  = data$survival$strata[[i]]$start
    stats$strata$end[i]    = data$survival$strata[[i]]$end
    stats$strata$events[i] = sum(data$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  = data$survival$strata[[i]]$label
  }

  names(stats$party)           = names.old
  names(stats$coefficients)    = names.old
  names(stats$expcoef)         = names.old
  names(stats$expncoef)        = names.old
  rownames(stats$var)          = names.old
  colnames(stats$var)          = names.old
  names(stats$secoef)          = names.old
  names(stats$zvals)           = names.old
  names(stats$pvals)           = names.old
  names(stats$stars)           = names.old
  names(stats$lower95)         = names.old
  names(stats$upper95)         = names.old
  names(stats$loglik)          = c("loglikelihood", "null loglikelihood")
  names(stats$score)           = c("score", "p-value")
  names(stats$lrt)             = c("likelihood ratio", "p-value")
  names(stats$rsquare)         = c("r-square", "max possible")
  names(stats$wald.test)       = c("wald", "p-value")
  names(stats$concordance)     = c("concordant", "discordant", "tied.risk", "tied.time",
                                   "concordance", "stderr")

  params$stats = stats
  writeTime = proc.time()[3]
  save(stats, file = file.path(params$writePath, "stats.rdata"))
  writeSize = file.size(file.path(params$writePath, "stats.rdata"))
  writeTime = proc.time()[3] - writeTime

  params = AddToLog(params, "ComputeResultsCox.A2",
                    readTime, readSize, writeTime, writeSize)
  return(params)
}


GetResultsCox.B2 = function(params) {
  stats = NULL
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.B2\n\n")
  readTime = proc.time()[3]
  load(file.path(params$readPath, "stats.rdata"))
  readSize = file.size(file.path(params$readPath, "stats.rdata"))
  readTime = proc.time()[3] - readTime
  params$stats = stats
  params = AddToLog(params, "GetResultsCox.B2", readTime, readSize, 0, 0)
  return(params)
}


####################### REGRESSION BY B ONLY FUNCTIONS #######################

FinalizeParams2Cox.B2 = function(params, data) {
  pa = NULL
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParams2Cox.B2\n\n")
  readTime = proc.time()[3]
  load(file.path(params$readPath, "pa.rdata")) # read pa
  readSize = file.size(file.path(params$readPath, "pa.rdata"))
  readTime = proc.time()[3] - readTime
  params$numEvents = sum(data$survival$status)
  params$p1 = pa$p1
  params$cutoff = pa$cutoff
  params$maxIterations = pa$maxIterations
  params$p = params$p1 + params$p2
  params = AddToLog(params, "FinalizeParams2Cox.B2", readTime, readSize, 0, 0)
  return(params)
}


CheckColinearityCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityCox.B2\n\n")
  # Add in the strata here
  numStrata = length(data$survival$strata)
  X = cbind(matrix(0, nrow = nrow(data$X), ncol = numStrata), data$X)
  for (i in 1:numStrata) {
    X[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
  }

  XTX = t(X) %*% X
  nrow = nrow(XTX)
  indicies = 1:numStrata
  for (i in (numStrata + 1):nrow) {
  	tempIndicies = c(indicies, i)
  	if (rcond(XTX[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
  		indicies = c(indicies, i)
  	}
  }
  indicies = indicies[-(1:numStrata)] - numStrata  # Get rid of the Strata

  params$AIndiciesKeep = c()
  params$Acolnames.old = c()
  params$Acolnames     = c()
  params$p1.old        = params$p1
  params$p1            = 0

  Bnames               = params$Bcolnames
  BnamesKeep           = Bnames[indicies]
  params$BIndiciesKeep = indicies
  params$Bcolnames.old = params$Bcolnames
  params$Bcolnames     = BnamesKeep
  params$p2.old        = params$p2
  params$p2            = length(BnamesKeep)
  params$p             = params$p1 + params$p2

  if (params$p2 == 0) {
    params$failed = TRUE
    params$errorMessage = "Party A has no covariates and all of Party B's covariates are linear."
    warning(params$errorMessage)
  }
  params = AddToLog(params, "CheckColinearityCox.B2", 0, 0, 0, 0)

  return(params)
}

#' @importFrom stats as.formula
ComputeCoxFromSurvival.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCoxFromSurvival.B2\n\n")
  # We have loaded survival previously

  strata = rep(0, nrow(data$X))
  for (i in 1:length(data$survival$strata)) {
    strata[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] = i
  }

  colnames(data$X) = paste0("V", 1:ncol(data$X))
  f = paste(c("Surv(rank, status) ~ strata(strata)", paste0("V", 1:ncol(data$X))), collapse = " + ")

  error = tryCatch(
    {fit = survival::coxph(as.formula(f),
                 data = data.frame(rank = data$survival$rank,
                                   status = data$survival$status,
                                   strata = strata,
                                   data$X),
                 iter.max = params$maxIterations)},
    error = function(e) { return(TRUE)},
    warning = function(e) { return(FALSE)}
  )

  if ((class(error) == "logical" && error)) {
    params$converged = FALSE
    params$failed    = TRUE
    params$errorMessage = "Coxph in the survival package failed to converge."
    warning(params$errorMessage)
  } else {
  	params$converged = TRUE
  	if (class(error) == "logical") {
  	  fit = suppressWarnings(survival::coxph(as.formula(f),
  								data = data.frame(rank = data$survival$rank,
  																	status = data$survival$status,
  																	strata = strata,
  																	data$X),
  								iter.max = params$maxIterations))
  		params$converged = FALSE
  	}
    fit$linear.predictors = NULL
    fit$residuals = NULL
    fit$y = NULL
    fit$formula = NULL
    fit$call = NULL
    fit$assign = NULL
    fit$terms = NULL
    fit$means = NULL
    params$fit = fit
  }
  params = AddToLog(params, "ComputeCoxFromSurvival.B2", 0, 0, 0, 0)
  return(params)
}


ComputeCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCox.B2\n\n")
  n           = params$n
  p2          = params$p2
  params$algIterationCounter = 1
  X.betas.old = matrix(0, n, 1)
  X.betas     = matrix(0, n, 1)
  betasB      = matrix(0, p2, 1)
  betasBold   = betasB
  loglikelihood.old = -Inf


  while (params$algIterationCounter <= params$maxIterations && !params$converged) {
  	BeginningIteration(params)
  	loglikelihood = 0
  	stepSize = 1
  	w = exp(X.betas)
  	while (max(w) == Inf) {
  		if (params$verbose) cat("Step Halving\n\n")
  		X.betas = (X.betas + X.betas.old) * 0.5
  		stepSize = stepSize * 0.5
  		w = exp(X.betas)
  	}
  	computeLoglikelihood = TRUE
  	while (computeLoglikelihood) {
  		numEvents = sum(data$survival$status)
  		stepCounter = 0
  		pbar = MakeProgressBar1(numEvents, "Loglikelihood", params$verbose)
  		loglikelihood = 0
  		for (i in 1:length(data$survival$strata)) {                    ##!
  			if (data$survival$strata[[i]]$J > 0) {                       ##!
  				for (j in 1:data$survival$strata[[i]]$J) {                 ##!
  					nj = data$survival$strata[[i]]$nfails[j]                 ##!
  					yIndex = data$survival$strata[[i]]$start0[j]:data$survival$strata[[i]]$end      ##!
  					zIndex = data$survival$strata[[i]]$start1[j]:data$survival$strata[[i]]$stop1[j] ##!
  					Aj1 = sum(w[yIndex])
  					Aj2 = sum(w[zIndex]) / nj
  					loglikelihood = loglikelihood + sum(log(w[zIndex]))
  					for (r in 0:(nj - 1)) {
  						Ajr = Aj1 - r * Aj2
  						loglikelihood = loglikelihood - log(Ajr)
  					}
  					stepCounter = stepCounter + nj
  					pbar = MakeProgressBar2(stepCounter, pbar, params$verbose)
  				}
  			}
  		}
  		if (loglikelihood > loglikelihood.old || stepSize < 0.5^6) {
  			computeLoglikelihood = FALSE
  		} else {
  			if (params$verbose) cat("Step Halving\n\n")
  			X.betas = (X.betas + X.betas.old) * 0.5
  			stepSize = stepSize * 0.5
  			w = exp(X.betas)
  		}
  	}
  	numEvents = sum(data$survival$status)
  	deltal = as.numeric(data$survival$status)
  	deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
  	# a pass by reference with the C call.
  	W.XB = matrix(0, n, p2)

  	.Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XB,
  				as.integer(n), as.integer(p2), as.integer(numEvents),
  				as.integer(params$verbose))

    M = t(data$X) %*% W.XB

    params$XtWX = M
    if (params$algIterationCounter == 1) {
      params$nullHessian = M
    }

  	inv = NULL
  	tryCatch({inv = solve(M)},
  					 error = function(err) { inv = NULL } )
  	M = inv
		if (is.null(M)) {
      params$failed = TRUE
      params$errorMessage = "The matrix t(X)WX is singular.  This is probably due to divergence of the coefficients."
      warning(params$errorMessage)

      betas = rep(NA, length(params$Bcolnames.old))
      betas[params$BIndiciesKeep] = betasB
      betas = data.frame(betas)
      rownames(betas) = params$Bcolnames.old
      # if (params$verbose) cat("Current Parameters:\n")
      # if (params$verbose) print(betas)
      # if (params$verbose) cat("\n")
      params = AddToLog(params, "ComputeCox.B2", 0, 0, 0, 0)
      return(params)
    }
    deltaBeta = M %*% t(data$X) %*% deltal
    betasB    = betasBold + (betasB - betasBold) * stepSize
    betasBold = betasB
    betasB    = betasB + deltaBeta
    X.betas   = data$X %*% betasB

    converged = abs(loglikelihood - loglikelihood.old) /
      (abs(loglikelihood) + 0.1) < params$cutoff
    params$converged = converged

    if (params$algIterationCounter == 1) {
      params$nullScore         = t(data$X) %*% deltal
      params$nullLoglikelihood = loglikelihood
    }
    loglikelihood.old = loglikelihood
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$loglikelihood = loglikelihood
  params$betasB = betasB
  params$X.betas = X.betas
  if (!params$converged) {
    params$failed = TRUE
    params$errorMessage = "Cox failed to converge in the specified number of iterations."
    warning(params$errorMessage)
  }
  params = AddToLog(params, "ComputeCox.B2", 0, 0, 0, 0)
  return(params)
}

#' @importFrom  stats pchisq pnorm qnorm
ComputeResultsCox.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsCox.B2\n\n")
  stats = params$stats
  stats$converged = params$converged
  stats$partyName = params$partyName
  stats$failed    = FALSE

  fitExists = !is.null(params$fit)
  names.old          = c(params$Acolnames.old, params$Bcolnames.old)
  idxA               = params$AIndiciesKeep
  idxB               = params$BIndiciesKeep
  idx                = c(idxA, idxB + length(params$Acolnames.old))
  stats$party        = c(rep("dp0", length(params$Acolnames.old)),
                         rep("dp1", length(params$Bcolnames.old)))
  stats$coefficients = rep(NA, length(stats$party))
  if (fitExists) {
    stats$coefficients[idx] = params$fit$coefficients
    tempvar          = params$fit$var
  } else {
    stats$coefficients[idx] = params$betasB
    tempvar            = solve(params$XtWX)
  }
  stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
  stats$expncoef     = exp(-stats$coefficients)
  stats$var          = matrix(0, length(names.old), length(names.old))
  stats$var[idx, idx] = tempvar
  stats$secoef       = rep(NA, length(names.old))
  stats$secoef[idx]  = sqrt(diag(tempvar))  # se(coef)

  stats$zvals        = stats$coefficients / stats$secoef  # z values
  stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE )   # pvals
  stats$stars        = matrix(sapply(stats$pvals, function(x) {
    if (is.na(x)) ''
    else if (x < 0.001) '***'
    else if (x < 0.01) '**'
    else if (x < 0.05) '*'
    else if (x < 0.1) '.'
    else ' '
  }))
  stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  if (fitExists) {
    stats$loglik     = params$fit$loglik
    stats$n          = params$fit$n
    stats$nevent     = params$fit$nevent
    stats$iter       = params$fit$iter
    stats$score      = params$fit$score
    stats$wald.test  = params$fit$wald.test
    stats$concordance = params$fit$concordance[c(2, 1, 3, 4, 6, 7)]
  } else {
    stats$loglik       = c(params$nullLoglikelihood, params$loglikelihood)
    stats$n            = params$n
    stats$nevent       = params$numEvents
    stats$iter         = params$algIterationCounter - 1
    stats$score        = t(params$nullScore) %*% solve(params$nullHessian) %*%
      params$nullScore
    stats$wald.test    = t(params$betasB) %*% params$XtWX %*% params$betasB
    stats$concordance = c(NA, NA, NA, NA, NA, NA)
  }
  stats$df           = params$p
  stats$score        = c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       = "efron"
  stats$lrt          = 2*(stats$loglik[2] - stats$loglik[1])
  stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      = c(1 - exp(-stats$lrt[1]/stats$n),
                         1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    = c(stats$wald.test,
                         1 - pchisq(stats$wald.test, stats$df))

  pred = data$X %*% stats$coefficients[idx]
  stats$survival = data.frame(
    rank   = data$survival$rank,
    status = data$survival$status,
    sorted = data$survival$sorted,
    surv   = SurvFitCox.A2(params, data$survival, pred)
  )
  stats$strata = as.data.frame(matrix(0, length(data$survival$strata), 3))
  stats$strata$label = ""
  colnames(stats$strata) = c("start", "end", "events", "label")
  for (i in 1:length(data$survival$strata)) {
    stats$strata$start[i]  = data$survival$strata[[i]]$start
    stats$strata$end[i]    = data$survival$strata[[i]]$end
    stats$strata$events[i] = sum(data$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  = data$survival$strata[[i]]$label
  }

  names(stats$party)           = names.old
  names(stats$coefficients)    = names.old
  names(stats$expcoef)         = names.old
  names(stats$expncoef)        = names.old
  rownames(stats$var)          = names.old
  colnames(stats$var)          = names.old
  names(stats$secoef)          = names.old
  names(stats$zvals)           = names.old
  names(stats$pvals)           = names.old
  names(stats$stars)           = names.old
  names(stats$lower95)         = names.old
  names(stats$upper95)         = names.old
  names(stats$loglik)          = c("loglikelihood", "null loglikelihood")
  names(stats$score)           = c("score", "p-value")
  names(stats$lrt)             = c("likelihood ratio", "p-value")
  names(stats$rsquare)         = c("r-square", "max possible")
  names(stats$wald.test)       = c("wald", "p-value")
  names(stats$concordance)     = c("concordant", "discordant", "tied.risk", "tied.time",
                                   "concordance", "stderr")

  params$stats = stats
  writeTime = proc.time()[3]
  save(stats, file = file.path(params$writePath, "stats.rdata"))
  writeSize = file.size(file.path(params$writePath, "stats.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "ComputeResultsCox.B2", 0, 0, writeTime, writeSize)
  return(params)
}


GetResultsCox.A2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.A2\n\n")
  readTime = proc.time()[3]
  load(file.path(params$readPath, "stats.rdata"))
  readSize = file.size(file.path(params$readPath, "stats.rdata"))
  readTime = proc.time()[3] - readTime
  stats$partyName = params$partyName
  params$stats = stats
  params = AddToLog(params, "GetResultsCox.A2", readTime, readSize, 0, 0)
  return(params)
}


############################## PARENT FUNCTIONS ##############################


PartyAProcess2Cox = function(data,
														 yname          = NULL,
														 strata         = NULL,
														 mask           = TRUE,
														 monitorFolder  = NULL,
														 msreqid        = "v_default_00_000",
														 blocksize      = 500,
														 cutoff         = 1e-8,
														 maxIterations  = 25,
														 sleepTime      = 10,
														 maxWaitingTime = 24 * 60 * 60,
														 popmednet      = TRUE,
														 trace          = FALSE,
														 verbose        = TRUE) {
  params = PrepareParams.2p("cox", "A", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)
  Header(params)
  params   = PrepareFolderCox.A2(params, monitorFolder)
  if (params$failed) {
  	warning(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataCox.23(params, data, yname, strata, mask)

  params = PauseContinue.2p(params, maxWaitingTime)
  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$complete = TRUE
    warning(ReadErrorMessage(params$readPath))
    params$pmnStepCounter = 1
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
  	params$complete = TRUE
  	message = "Error in processing the data for Party A."
  	MakeErrorMessage(params$writePath, message)
  	files = c("errorMessage.rdata")
  	params$pmnStepCounter = 1
  	params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
  	params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
  	SummarizeLog.2p(params)
  	return(params$stats)
  }

  params = CheckStrataCox.A2(params, data)
  if (params$getStrataFromB) {
    files = c()
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
  }

  params = PrepareParamsCox.A2(params, data, cutoff, maxIterations)

  if (params$failed) {   # Check for failed from PrepareParamsCox.A2()
    params$complete = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }


  data = PrepareStrataCox.A2(params, data)
  params = AddToLog(params, "PrepareStrataCox.A2", data$readTime,
                     data$readSize, data$writeTime, data$writeSize)

  if (params$p1 == 0) { # Check for $p1 == 0 => no covariates, only strata
    MakeTransferMessage(params$writePath)
    files = c("transferControl.rdata", "pa.rdata", "survival.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
    if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
      params$complete = TRUE
      warning(ReadErrorMessage(params$readPath))
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    params = GetResultsCox.A2(params)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = PrepareBlocksCox.A2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksCox.A2()
    params$complete = TRUE
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params = GetZCox.A2(params, data)

  #This works even in the case that we send no blocks over.  Just set
  #Filebreak.Z = c()
  files = c("pa.rdata", "blocksize.rdata", "survival.rdata",
            SeqZW("cz_", length(params$container$filebreak.Z)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = CheckColinearityCox.A2(params, data)
  # if (params$verbose) cat(params$failed,"\n")



  if (params$failed) { # Check for CheckColinearityCox.A2() failed
    params$complete = TRUE
    warning(params$errorMessage)
    MakeErrorMessage(params$writePath, params$errorMessage)
    files = c("errorMessage.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (params$p1 == 0) { # No covariates left.  All colinear with Strata
    MakeTransferMessage(params$writePath)
    files = c("transferControl.rdata", "indicies.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
    if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
      params$complete = TRUE
      warning(ReadErrorMessage(params$readPath))
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    params = GetResultsCox.A2(params)
    params = SendPauseQuit.2p(params, sleepTime = sleepTime)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  data = UpdateDataCox.A2(params, data)
  params = AddToLog(params, "UpdateDataCox.A2", 0, 0, 0, 0)
  params = PrepareLoopCox.A2(params)

  while (!params$converged && !params$halted) {
    BeginningIteration(params)

    params = ComputeLogLikelihoodCox.A2(params, data)

    files = c("X_betas_ss.rdata")
    if (params$algIterationCounter == 1) {
      files = c("indicies.rdata", files)
    } else {
      files = c("converged.rdata", files)
    }

    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    params = ComputeInverseCox.A2(params, data)

    if (params$failed) { # Check for failed from ComputeInverseCox.A2()
      params$complete = TRUE
      MakeErrorMessage(params$writePath, params$errorMessage)
      files = c("errorMessage.rdata")
      params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }

    files = c("M.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    params = ComputeBetaCox.A2(params, data)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$complete = TRUE
  files = c("converged.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = ComputeResultsCox.A2(params, data)

  files = c("stats.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime = sleepTime)
  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  SummarizeLog.2p(params)
  return(params$stats)
}


PartyBProcess2Cox = function(data,
                             strata              = NULL,
                             mask                = TRUE,
														 monitorFolder       = NULL,
														 sleepTime           = 10,
                             maxWaitingTime      = 24 * 60 * 60,
														 popmednet           = TRUE,
														 trace               = FALSE,
														 verbose             = TRUE) {
  params = PrepareParams.2p("cox", "B", popmednet = popmednet, trace = trace,
                            verbose = verbose)
  params = InitializeLog.2p(params)
  params = InitializeStamps.2p(params)
  params = InitializeTrackingTable.2p(params)
  Header(params)
  params = PrepareFolderCox.B2(params, monitorFolder)
  if (params$failed) {
  	warning(params$errorMessage)
  	return(invisible(NULL))
  }
  data = PrepareDataCox.23(params, data, NULL, strata, mask)

  if (data$failed) { # Check for Error from PrepareDataCox.B2()
    params$complete = TRUE
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$writePath, message)
    files = c("errorMessage.rdata")
    params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  params = PrepareParamsCox.B2(params, data)
  files = c("pb.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
    params$complete = TRUE
    warning(ReadErrorMessage(params$readPath))
    params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
    return(params$stats)
  }

  if (!is.null(data$strata$strataFromB)) {
    params = SendStrataCox.B2(params, data)
    files = c("strata.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
      params$complete = TRUE
      warning(ReadErrorMessage(params$readPath))
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }
  }

  data = SortDataCox.B2(params, data)
  params = AddToLog(params, "SortDataCox.B2", data$readTime,
                     data$readSize, 0, 0)

  if (file.exists(file.path(params$readPath, "transferControl.rdata"))) {
    params = FinalizeParams2Cox.B2(params, data)
    params = CheckColinearityCox.B2(params, data)

    if (params$failed) {  # Happens if pB.new == 0
      params$complete = TRUE
      MakeErrorMessage(params$writePath, params$errorMessage)
      files = c("errorMessage.rdata")
      params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }
    data = UpdateDataCox.B2(params, data)
    params = AddToLog(params, "UpdateDataCox.B2", 0, 0, 0, 0)
    if (params$survivalInstalled) {
      params = ComputeCoxFromSurvival.B2(params, data)
    } else {
      params = ComputeCox.B2(params, data)
    }

    if (params$failed) {      # We could get a job_failed here from coefficient explosion
      params$complete = TRUE
      MakeErrorMessage(params$writePath, params$errorMessage)
      files = c("errorMessage.rdata")
      params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }
    params = ComputeResultsCox.B2(params, data)
    stats = params$stats
    save(stats, file = file.path(params$writePath, "stats.rdata"))
    files = c("stats.rdata")
    params = SendPauseQuit.2p(params, files, sleepTime = sleepTime)
    return(params$stats)
  }

  params = FinalizeParamsCox.B2(params, data)
  params = PrepareBlocksCox.B2(params)
  params = GetWCox.B2(params, data)

  files = c("xbtxb.rdata", SeqZW("cw_", length(params$container$filebreak.W)))
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  if (file.exists(file.path(params$readPath, "transferControl.rdata"))) {
    params = UpdateParamsCox.B2(params)
    data = UpdateDataCox.B2(params, data)
    params = AddToLog(params, "UpdateDataCox.B2", 0, 0, 0, 0)
    if (params$survivalInstalled) {
      params = ComputeCoxFromSurvival.B2(params, data)
    } else {
      params = ComputeCox.B2(params, data)
    }

    if (params$failed) {      # We could get a job_failed here from coefficient explosion
      params$complete = TRUE
      MakeErrorMessage(params$writePath, params$errorMessage)
      files = c("errorMessage.rdata")
      params = SendPauseQuit.2p(params, files, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }
    params = ComputeResultsCox.B2(params, data)
    stats = params$stats
    save(stats, file = file.path(params$writePath, "stats.rdata"))
    files = c("stats.rdata")
    params = SendPauseQuit.2p(params, files, sleepTime = sleepTime)
    return(params$stats)
  }

  params$algIterationCounter = 1
  repeat {
    if (params$algIterationCounter == 1) {
      if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
        params$complete = TRUE
        warning(ReadErrorMessage(params$readPath))
        params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
        return(params$stats)
      }
      params = UpdateParamsCox.B2(params)
      data = UpdateDataCox.B2(params, data)
      params = AddToLog(params, "UpdateDataCox.B2", 0, 0, 0, 0)
    } else {
      params = GetConvergedStatusCox.B2(params)
      if (params$converged || params$halted) {
        break
      }
    }

    BeginningIteration(params)

    params = ComputeLogLikelihoodCox.B2(params, data)

    files = c("tXB_W_XB.rdata", SeqZW("cCox_", length(params$container$filebreak.Cox)))
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    if (file.exists(file.path(params$readPath, "errorMessage.rdata"))) {
      params$complete = TRUE
      warning(ReadErrorMessage(params$readPath))
      params = SendPauseQuit.2p(params, sleepTime = sleepTime, job_failed = TRUE)
      return(params$stats)
    }

    params = ComputeBetaCox.B2(params, data)

    files = c("XB_betasB.rdata", "M2_tXB_deltal.rdata")
    params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

    EndingIteration(params)

    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$complete = TRUE

  files = c("B_betas_ns.rdata")
  params = SendPauseContinue.2p(params, files, sleepTime, maxWaitingTime)

  params = GetResultsCox.B2(params)

  params = SendPauseQuit.2p(params, sleepTime = sleepTime)
  return(params$stats)
}
