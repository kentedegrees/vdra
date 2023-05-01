################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

PrepareFolderLogistic.A2 = function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.A2\n\n")
  params$dplocalPath   = file.path(monitor_folder, "dplocal")
  params$rprogramsPath = file.path(monitor_folder, "rprograms")
  params$macrosPath    = file.path(monitor_folder, "macros")
  params$write_path     = file.path(monitor_folder, "inputfiles")
  params$read_path      = file.path(monitor_folder, "msoc1")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed = TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
    params$failed = TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }
  params$error_message = NULL
  if (!CreateIOLocation(monitor_folder, "dplocal")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dplocalPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "rprograms")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$rprogramsPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "macros")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macrosPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "inputfiles")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "msoc1")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  params <- add_to_log(params, "PrepareDataLogistic.A23, PrepareFolderLogistic.A2", 0, 0, 0, 0)
  return(params)
}


PrepareFolderLogistic.B2 = function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.B2\n\n")

  params$dplocalPath   = file.path(monitor_folder, "dplocal")
  params$rprogramsPath = file.path(monitor_folder, "rprograms")
  params$macrosPath    = file.path(monitor_folder, "macros")
  params$write_path     = file.path(monitor_folder, "msoc")
  params$read_path      = file.path(monitor_folder, "inputfiles")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed = TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
    params$failed = TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }
  params$error_message = NULL
  if (!CreateIOLocation(monitor_folder, "dplocal")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dplocalPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "rprograms")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$rprogramsPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "macros")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macrosPath, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "msoc")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "inputfiles")) {
    params$failed = TRUE
    params$error_message = paste(params$error_message,
                                "ould not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  Sys.sleep(1)
  DeleteTrigger("files_done.ok", params$read_path)

  params <- add_to_log(params, "PrepareDataLogisitc.B23, PrepareFolderLogistic.B2", 0, 0, 0, 0)

  return(params)
}


#' @importFrom stats model.matrix
PrepareDataLogistic.A23 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.A23\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  response_index = CheckResponse(params, data, yname)

  if (is.null(response_index)) {
    workdata$failed = TRUE
    return(workdata)
  }
  covariate_index = setdiff(1:ncol(data), response_index)
  workdata$tags = CreateModelMatrixTags(data[, covariate_index, drop = FALSE])
  workdata$tags = c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] = "numeric"

  X = model.matrix(~ ., data[, c(response_index, covariate_index), drop = FALSE])
  rownames(X) = NULL
  covariate_index = setdiff(1:ncol(X), 2)

  means = apply(X, 2, mean)
  sd    = apply(X, 2, sd)
  sd    = sapply(sd, function(x) {
    ifelse(x > 0, x, 1)
  })
  workdata$Y      = X[, 2, drop = FALSE]
  workdata$X      = X[, covariate_index, drop = FALSE]
  workdata$means  = means[covariate_index]
  workdata$sd     = sd[covariate_index]
  workdata$yty    = t(workdata$Y) %*% workdata$Y

  if (ncol(workdata$X) >= 2) {
    for (i in 2:ncol(workdata$X)) {
      workdata$X[, i] = (workdata$X[, i] - workdata$means[i]) / workdata$sd[i]
    }
  }

  return(workdata)
}

#' @importFrom stats model.matrix sd
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
    warning("The data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.")
    workdata$failed = TRUE
    return(workdata)
  }
  workdata$X = model.matrix(~ ., data)
  rownames(workdata$X) = NULL
  workdata$X = workdata$X[, -1, drop = FALSE]
  workdata$means = apply(workdata$X, 2, mean)
  workdata$sd    = apply(workdata$X, 2, sd)
  workdata$sd    = sapply(workdata$sd, function(x) {
    ifelse(x > 0, x, 1)
  })

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

  write_time = proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLogistic.B2", 0, 0, write_time, write_size)
  return(params)
}


PrepareParamsLogistic.A2 = function(params, data, cutoff = 0.01, maxIterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLogistic.A2\n\n")

  params$converged       = FALSE
  params$halted          = FALSE
  params$pmnStepCounter  = 1
  pb                     = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb, Bcolnames
  read_size = sum(file.size(file.path(params$read_path, "pb.rdata")))
  read_time = proc.time()[3] - read_time

  if (params$analysis != pb$analysis) {
    params$error_message =
      paste("Party A is running", params$analysis, "regression and Party B is running", pb$analysis, "regression.")
    warning(params$error_message)
    params$failed = TRUE
    return(params)
  }

  params$n  = nrow(data$X)
  if (pb$n != params$n) {
    params$error_message =
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    warning(params$error_message)
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

  write_time = proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLogistic.A2", read_time, read_size,
                    write_time, write_size)

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

    params$error_message =
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
      params$error_message =
        paste0(params$error_message,
               "\nSet the number of B covariates to be between ", minBCovariates, "and",
               paste0(maxBCovariates, "."))
    }
    warning(params$error_message)
    params$failed = TRUE
    params <- add_to_log(params, "PrepareBlocksLogistic.A2", 0, 0, 0, 0)
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
  write_time = proc.time()[3]
  save(blocksize, file = file.path(params$write_path, "blocksize.rdata"))
  write_size = file.size(file.path(params$write_path, "blocksize.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareBlocksLogistic.A2", 0, 0, write_time, write_size)
  return(params)
}


GetZLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLogistic.A2\n\n")
  write_time = 0
  write_size = 0

  numBlocks = params$blocks$numBlocks
  pbar = MakeProgressBar1(numBlocks, "Z", params$verbose)
  containerCt.Z = 0
  for (i in 1:numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename = paste0("cz_", containerCt.Z, ".rdata")
      toWrite = file(file.path(params$write_path, filename), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]
    Z = FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$X[strt:stp, ]), g)

    write_time = write_time - proc.time()[3]
    writeBin(as.vector(Z), con = toWrite, endian = "little")
    write_time = write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
      close(toWrite)
      write_size = write_size + file.size(file.path(params$write_path, filename))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "GetZLogistic.A2", 0, 0, write_time, write_size)
  return(params)
}


FinalizeParamsLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsLogistic.B2\n\n")
  pa = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # Load pa, Acolnames
  read_size = sum(file.size(file.path(params$read_path, "pa.rdata")))
  read_time = proc.time()[3] - read_time
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

  params <- add_to_log(params, "FinalizeParamsLogistic.B2", read_time, read_size, 0, 0)
  return(params)
}


PrepareBlocksLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLogistic.B2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  read_time = proc.time()[3]
  load(file.path(params$read_path, "blocksize.rdata")) # load blocksize
  read_size = file.size(file.path(params$read_path, "blocksize.rdata"))
  read_time = proc.time()[3] - read_time
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params <- add_to_log(params, "PrepareBlocksLogistic.B2", read_time, read_size, 0, 0)
  return(params)
}


GetWLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWLogistic.B2\n\n")
  read_time  = 0
  read_size  = 0
  write_time = 0
  write_size = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I-Z*Z')XB", params$verbose)

  XBTXB = t(data$X) %*% data$X

  containerCt.Z = 0
  containerCt.W = 0

  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename1 = paste0("cz_", containerCt.Z, ".rdata")
      toRead = file(file.path(params$read_path, filename1), "rb")
      read_size = read_size + file.size(file.path(params$read_path, filename1))
    }
    if (i %in% params$container$filebreak.W) {
      containerCt.W = containerCt.W + 1
      filename2 = paste0("cw_", containerCt.W, ".rdata")
      toWrite = file(file.path(params$write_path, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1
    g1 = params$blocks$g[i]

    read_time = read_time - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n2 * g1,
                       endian = "little"), nrow = n2, ncol = g1)
    read_time = read_time + proc.time()[3]

    W = data$X[strt:stp, ] - Z %*% (t(Z) %*% data$X[strt:stp, ])

    write_time = write_time - proc.time()[3]
    writeBin(as.vector(W), con = toWrite, endian = "little")
    write_time = write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.Z || i == params$blocks$numBlocks) {
      close(toRead)
    }
    if ((i + 1) %in% params$container$filebreak.W || i == params$blocks$numBlocks) {
      close(toWrite)
      write_size = write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  write_time = write_time - proc.time()[3]
  save(XBTXB, file = file.path(params$write_path, "xbtxb.rdata"))
  write_size = write_size + file.size(file.path(params$write_path, "xbtxb.rdata"))
  write_time = write_time + proc.time()[3]

  params <- add_to_log(params, "GetWLogistic.B2", read_time, read_size, write_time, write_size)

  return(params)
}


CheckColinearityLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.A2\n\n")
  p2 = params$p2
  read_time  = 0
  read_size  = 0
  write_time = 0
  write_size = 0
  XBTXB     = NULL

  read_time = read_time - proc.time()[3]
  load(file.path(params$read_path, "xbtxb.rdata")) # load XBTXB
  read_size = file.size(file.path(params$read_path, "xbtxb.rdata"))
  read_time = read_time + proc.time()[3]
  XATXA = t(data$X) %*% data$X
  XATXB = 0
  XATY  = t(data$X) %*% data$Y
  YTXB  = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'X", params$verbose)

  containerCt.W = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.W) {
      containerCt.W = containerCt.W + 1
      filename = paste0("cw_", containerCt.W, ".rdata")
      toRead = file(file.path(params$read_path, filename), "rb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n2 = stp - strt + 1

    read_time = read_time - proc.time()[3]
    W = matrix(readBin(con = toRead, what = numeric(), n = n2 * p2,
                       endian = "little"), nrow = n2, ncol = p2)
    read_time = read_time + proc.time()[3]

    XATXB = XATXB + t(data$X[strt:stp, ]) %*% W
    YTXB  = YTXB  + t(data$Y[strt:stp, ]) %*% W

    if ((i + 1) %in% params$container$filebreak.W || i == params$blocks$numBlocks) {
      close(toRead)
      read_size = read_size + file.size(file.path(params$read_path, filename))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
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

  write_time = write_time - proc.time()[3]
  save(Aindicies, Bindicies, file = file.path(params$write_path, "indicies.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "indicies.rdata")))
  write_time = write_time + proc.time()[3]

  tags = params$Btags[params$BIndiciesKeep]

  if (length(unique(tags)) < 2) {
    params$failed = TRUE
    params$error_message = "After removing colinear covariates, Party B has 1 or fewer covariates.\n\n"
  } else if (!("numeric" %in% names(tags))) {
    params$failed = TRUE
    params$error_message = "After removing colinear covariates, Party B has no continuous covariates.\n\n"
  }

  params <- add_to_log(params, "CheckColinearityLogistic.A2", read_time, read_size, write_time, write_size)
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

  write_time = proc.time()[3]
  save(Bbetas, Bxty, file = file.path(params$write_path, "Bbetas_xty.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "Bbetas_xty.rdata")))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeInitialBetasLogistic.A2", 0, 0, write_time, write_size)

  return(params)
}


UpdateParamsLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.B2\n\n")
  Aindicies = NULL
  Bindicies = NULL
  Bbetas    = NULL
  Bxty      = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "indicies.rdata")) # load Aindicies, Bindicies
  load(file.path(params$read_path, "Bbetas_xty.rdata"))     # Load Bbetas
  read_size = sum(file.size(file.path(params$read_path, c("indicies.rdata",
                                                        "Bbetas_xty.rdata"))))
  read_time = proc.time()[3] - read_time
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
  params <- add_to_log(params, "UpdateParamsLogistic.B2", read_time, read_size, 0, 0)
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

  write_time = proc.time()[3]
  save(XbetaB, file = file.path(params$write_path, "xbetab.rdata"))
  write_size = file.size(file.path(params$write_path, "xbetab.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetXBetaLogistic.B2", 0, 0, write_time, write_size)
  return(params)
}


GetWeightsLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWeightsLogistic.A2\n\n")
  n      = params$n
  p1     = params$p1
  XbetaB = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "xbetab.rdata"))  # Load XbetaB
  read_size = file.size(file.path(params$read_path, "xbetab.rdata"))
  read_time = proc.time()[3] - read_time

  XbetaA = data$X %*% params$betasA
  Xbeta = XbetaA + XbetaB
  pi_ = (1 + exp(-Xbeta))^(-1)
  params$pi_ = pi_

  write_time = proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi_.rdata"))
  write_size = file.size(file.path(params$write_path, "pi_.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "GetWeightsLogistic.A2", read_time, read_size, write_time, write_size)
  return(params)
}


GetVLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetVLogistic.B2\n\n")
  pi_       = NULL
  write_time = 0
  write_size = 0
  read_time  = proc.time()[3]
  load(file.path(params$read_path, "pi_.rdata"))
  read_size  = file.size(file.path(params$read_path, "pi_.rdata"))
  read_time  = proc.time()[3] - read_time

  params$pi_ = pi_
  W = pi_ * (1 - params$pi_)

  XBTWXB = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I - Z*Z')W*XB", params$verbose)

  containerCt.Z = 0
  containerCt.V = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.Z) {
      containerCt.Z = containerCt.Z + 1
      filename1 = paste0("cz_", containerCt.Z, ".rdata")
      toRead = file(file.path(params$read_path, filename1), "rb")
    }
    if (i %in% params$container$filebreak.V) {
      containerCt.V = containerCt.V + 1
      filename2 = paste0("cv_", containerCt.V, ".rdata")
      toWrite = file(file.path(params$write_path, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]

    Xblock  = data$X[strt:stp, ]
    Wblock  = W[strt:stp]
    WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

    read_time = read_time - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n * g,
                       endian = "little"), nrow = n, ncol = g)
    read_time = read_time + proc.time()[3]

    V = WXblock - Z %*% (t(Z) %*% WXblock)

    write_time = write_time - proc.time()[3]
    writeBin(as.vector(V), con = toWrite, endian = "little")
    write_time = write_time + proc.time()[3]

    XBTWXB = XBTWXB + t(Xblock) %*% WXblock
    if ((i + 1) %in% params$container$filebreak.Z || i == params$blocks$numBlocks) {
      close(toRead)
      read_size = read_size + file.size(file.path(params$read_path, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$numBlocks) {
      close(toWrite)
      write_size = write_size + file.size(file.path(params$write_path, filename2))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  # sums of each column in WX_B
  sumsWXB = apply(MultiplyDiagonalWTimesX(W, data$X), 2, sum)
  # This information needs to be shared in order to get the intercept term

  write_time = write_time - proc.time()[3]
  save(sumsWXB, XBTWXB, file = file.path(params$write_path, "sumswx_xbtwxb.rdata"))
  write_size = write_size + sum(file.size(c(file.path(params$write_path, "sumswx_xbtwxb.rdata"))))
  write_time = write_time + proc.time()[3]

  params <- add_to_log(params, "GetVLogistic.B2", read_time, read_size, write_time, write_size)
  return(params)
}


GetIILogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetIILogistic.A2\n\n")
  p1 = params$p1
  p2 = params$p2
  sumsWXB = NULL
  XBTWXB  = NULL

  write_time = 0
  write_size = 0
  read_time = proc.time()[3]
  load(file.path(params$read_path, "sumswx_xbtwxb.rdata")) # load sumsWXB, XBTWXB
  read_size = sum(file.size(file.path(params$read_path, "sumswx_xbtwxb.rdata")))
  read_time = proc.time()[3] - read_time

  params$sumsWXB = sumsWXB

  IA = params$Axty - t(data$X) %*% params$pi_
  W = params$pi_ * (1 - params$pi_)
  sumsWXA = apply(MultiplyDiagonalWTimesX(W, data$X), 2, sum)[-1]
  params$sumsWXA = sumsWXA

  XATWXA = t(data$X) %*% MultiplyDiagonalWTimesX(W, data$X)

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'W*X", params$verbose)

  XATWXB = 0
  containerCt.V = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.V) {
      containerCt.V = containerCt.V + 1
      filename1 = paste0("cv_", containerCt.V, ".rdata")
      toRead = file(file.path(params$read_path, filename1), "rb")
      read_size = read_size + file.size(file.path(params$read_path, filename1))
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1

    read_time = read_time - proc.time()[3]
    V = matrix(readBin(con = toRead, what = numeric(),
                       n = n * p2, endian = "little"), n, p2)
    read_time = read_time + proc.time()[3]
    XATWXB = XATWXB + t(data$X[strt:stp, ]) %*% V

    pbar = MakeProgressBar2(i, pbar, params$verbose)
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$numBlocks) {
      close(toRead)
    }
  }

  XTWX = rbind(cbind(XATWXA, XATWXB), cbind(t(XATWXB), XBTWXB))

  params$xtwx = XTWX

  II = NULL
  tryCatch({
    II = solve(params$xtwx)
  }, # dims are 1 + p1 + p2
  error = function(err) {
    II = NULL
  }
  )
  if (is.null(II)) {
    params$singularMatrix = TRUE
    params$failed = TRUE
    params$error_message =
      paste0("The matrix t(X)*W*X is not invertible.\n",
             "       This may be due to one of two possible problems.\n",
             "       1. Poor random initialization of the security vector.\n",
             "       2. Near multicollinearity in the data\n",
             "SOLUTIONS: \n",
             "       1. Rerun the data analysis.\n",
             "       2. If the problem persists, check the variables for\n",
             "          duplicates for both parties and / or reduce the\n",
             "          number of variables used. Once this is done,\n",
             "          rerun the data analysis.")
    warning(params$error_message)
  } else {
    params$II = II
    params$IA = IA

    a21i1 = II[(p1 + 1):(p1 + p2), 1:p1] %*% matrix(IA, p1, 1)
    a11i1 = II[1:p1, 1:p1] %*% matrix(IA, p1, 1)
    params$a11i1 = a11i1

    write_time = proc.time()[3]
    save(a21i1, XTWX, file = file.path(params$write_path, "a21i1_xtwx.rdata"))
    write_size = sum(file.size(file.path(params$write_path, "a21i1_xtwx.rdata")))
    write_time = proc.time()[3] - write_time
  }
  params <- add_to_log(params, "GetIILogistic.A2", read_time, read_size, write_time, write_size)

  return(params)
}


GetCoefLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.B2\n\n")
  p1 = params$p1
  p2 = params$p2
  XTWX  = NULL
  a21i1 = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "a21i1_xtwx.rdata"))   # load a21i1, XTWX
  read_size = sum(file.size(file.path(params$read_path, "a21i1_xtwx.rdata")))
  read_time = proc.time()[3] - read_time

  IB = params$Bxty - t(data$X) %*% params$pi_

  II = solve(XTWX)

  params$II = II
  params$IB = IB

  a22i2 = II[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  a12i2 = II[1:p1, (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  params$a22i2 = a22i2

  params$betasBold = params$betasB
  params$betasB = params$betasB + a21i1 + a22i2

  deltabetaB = max(abs(params$betasB - params$betasBold) / (abs(params$betasB) + 0.1))

  write_time = proc.time()[3]
  save(a12i2, deltabetaB, file = file.path(params$write_path, "a12_deltabetaB.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "a12_deltabetaB.rdata")))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetCoefLogistic.B2", read_time, read_size, write_time, write_size)
  return(params)
}


GetCoefLogistic.A2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.A2\n\n")
  a12i2      = NULL
  deltabetaB = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "a12_deltabetaB.rdata")) # load  a12i2, deltabetab
  read_size = sum(file.size(file.path(params$read_path, "a12_deltabetaB.rdata")))
  read_time = proc.time()[3] - read_time

  params$betasAold = params$betasA
  params$betasA = params$betasA + params$a11i1 + a12i2

  deltabeta = max(abs(params$betasA - params$betasAold) / (abs(params$betasA) + 0.1), deltabetaB)

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$algIterationCounter >= params$maxIterations) {
    params$maxIterExceeded = TRUE
    warning(paste("Failed to converged in", params$maxIterations, "iterations."))
  }

  write_time = proc.time()[3]
  save(deltabeta, file = file.path(params$write_path, "deltabeta.rdata"))
  write_size = file.size(file.path(params$write_path, "deltabeta.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetCoefLogistic.A2", read_time, read_size, write_time, write_size)


  return(params)
}


GetConvergedStatusLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetconvergedStatusLogistic.B2\n\n")
  deltabeta = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "deltabeta.rdata")) # load deltabeta.rdata
  read_size = file.size(file.path(params$read_path, "deltabeta.rdata"))
  read_time = proc.time()[3] - read_time

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$algIterationCounter >= params$maxIterations) {
    params$maxIterExceeded = TRUE
    warning(paste("Failed to converged in", params$maxIterations, "iterations."))
  }

  params <- add_to_log(params, "GetConvergedStatusLogistic.B2", read_time, read_size, 0, 0)
  return(params)
}


GetFinalCoefLogistic.B2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetFinalCoefLogistic.B2\n\n")
  betasB = params$betasB / params$sdB
  offsetB = sum(betasB * params$meansB)
  BFinalFitted = t(params$sdB * t(data$X) + params$meansB) %*% betasB
  write_time = proc.time()[3]
  save(betasB, BFinalFitted, offsetB, file = file.path(params$write_path, "b_final.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "b_final.rdata")))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "GetFinalCoefLogistic.B2", 0, 0, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
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


  betasB = NULL
  offsetB = NULL
  BFinalFitted = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "b_final.rdata"))  # betasB, offsetB, BFinalFitted
  read_size = sum(file.size(file.path(params$read_path, "b_final.rdata")))
  read_time = proc.time()[3] - read_time
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

  write_time = proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size = file.size(file.path(params$write_path, "stats.rdata"))
  write_time = proc.time()[3] - write_time

  stats$Y           = data$Y # For Hoslem and ROC
  stats$FinalFitted = FinalFitted
  params$stats      = stats

  params <- add_to_log(params, "ComputeResultsLogistic.B2", read_time, read_size, write_time, write_size)
  return(params)
}


GetResultsLogistic.B2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.B2\n\n")
  stats = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size = file.size(file.path(params$read_path, "stats.rdata"))
  read_time = proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsLogistic.B2", read_time, read_size, 0, 0)
  return(params)
}



############################### PARENT FUNCTIONS ###############################


PartyAProcess2Logistic = function(data,
                                  yname                 = NULL,
                                  monitor_folder         = NULL,
                                  msreqid               = "v_default_00_000",
                                  blocksize             = 500,
                                  cutoff                = 1e-8,
                                  maxIterations         = 25,
                                  sleep_time             = 10,
                                  maxWaitingTime        = 24 * 60 * 60,
                                  popmednet             = TRUE,
                                  trace                 = FALSE,
                                  verbose               = TRUE) {
  params <- PrepareParams.2p("logistic", "A", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.2p(params)
  params <- InitializeStamps.2p(params)
  params <- InitializeTrackingTable.2p(params)
  Header(params)
  params   = PrepareFolderLogistic.A2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLogistic.A23(params, data, yname)

  params <- PauseContinue.2p(params,  maxWaitingTime)
  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$read_path))
    params$pmnStepCounter = 1
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
    params$completed = TRUE
    message = "Error in processing the data for Party A."
    MakeErrorMessage(params$write_path, message)
    files = c("error_message.rdata")
    params$pmnStepCounter = 1
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- PrepareParamsLogistic.A2(params, data, cutoff, maxIterations)

  if (params$failed) {   # Check for failed from PrepareParamsLogistic.A2()
    params$completed = TRUE
    MakeErrorMessage(params$write_path, params$error_message)
    files = c("error_message.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- PrepareBlocksLogistic.A2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksLogistic.A2()
    params$completed = TRUE
    MakeErrorMessage(params$write_path, params$error_message)
    files = c("error_message.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- GetZLogistic.A2(params, data)

  files = c("pa.rdata", "blocksize.rdata",
            SeqZW("cz_", length(params$container$filebreak.Z)))
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  params <- CheckColinearityLogistic.A2(params, data)

  if (params$failed) { # Check for CheckColinearityLogistic.A2() failed
    params$completed = TRUE
    warning(params$error_message)
    MakeErrorMessage(params$write_path, params$error_message)
    files = c("error_message.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }
  data = UpdateDataLogistic.A2(params, data)
  params <- add_to_log(params, "UpdateDataLogistic.A2", 0, 0, 0, 0)
  params <- ComputeInitialBetasLogistic.A2(params, data)

  files = c("indicies.rdata", "Bbetas_xty.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetWeightsLogistic.A2(params, data)
    files = c("pi_.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)
    params <- GetIILogistic.A2(params, data)

    if (params$failed) { # Check for failed from ComputeInverseLogistic.A2()
      params$completed = TRUE
      MakeErrorMessage(params$write_path, params$error_message)
      files = c("error_message.rdata")
      params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
      params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    files = c("a21i1_xtwx.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

    params <- GetCoefLogistic.A2(params, data)
    files = "deltabeta.rdata"
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$completed = TRUE

  params <- ComputeResultsLogistic.A2(params, data)

  files = c("stats.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)
  params <- SendPauseQuit.2p(params, sleep_time = sleep_time)
  SummarizeLog.2p(params)
  return(invisible(params$stats))
}

PartyBProcess2Logistic = function(data,
                                  monitor_folder       = "v_default_00_000",
                                  sleep_time           = 10,
                                  maxWaitingTime      = 24 * 60 * 60,
                                  popmednet           = TRUE,
                                  trace               = FALSE,
                                  verbose             = TRUE) {
  params <- PrepareParams.2p("logistic", "B",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.2p(params)
  params <- InitializeStamps.2p(params)
  params <- InitializeTrackingTable.2p(params)
  Header(params)
  params   = PrepareFolderLogistic.B2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLogistic.B23(params, data)

  if (data$failed) { # Check for Error from PrepareDataLogistic.B2()
    params$completed = TRUE
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$write_path, message)
    files = c("error_message.rdata")
    params <- SendPauseQuit.2p(params, files, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params   = PrepareParamsLogistic.B2(params, data)

  files = c("pb.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$read_path))
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params <- FinalizeParamsLogistic.B2(params, data)
  params <- PrepareBlocksLogistic.B2(params)
  params <- GetWLogistic.B2(params, data)

  files = c("xbtxb.rdata", SeqZW("cw_", length(params$container$filebreak.W)))
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$read_path))
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params <- UpdateParamsLogistic.B2(params)
  data = UpdateDataLogistic.B2(params, data)
  params <- add_to_log(params, "UpdateDataLogistic.B2", 0, 0, 0, 0)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetXBetaLogistic.B2(params, data)

    files = c("xbetab.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

    params <- GetVLogistic.B2(params, data)
    files = c("sumswx_xbtwxb.rdata",
              SeqZW("cv_", length(params$container$filebreak.V)))
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$completed = TRUE
      warning(ReadErrorMessage(params$read_path))
      params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
      return(params$stats)
    }

    params <- GetCoefLogistic.B2(params, data)
    files = c("a12_deltabetaB.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

    params <- GetConvergedStatusLogistic.B2(params)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$completed = TRUE

  params <- GetFinalCoefLogistic.B2(params, data)
  files = c("b_final.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  params <- GetResultsLogistic.B2(params)
  params <- SendPauseQuit.2p(params, sleep_time = sleep_time)
  return(invisible(params$stats))
}
