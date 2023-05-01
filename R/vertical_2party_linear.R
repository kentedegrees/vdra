################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

PrepareFolderLinear.a2 = function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.a2\n\n")
  params$dp_local_path   = file.path(monitor_folder, "dplocal")
  params$r_programs_path = file.path(monitor_folder, "rprograms")
  params$macros_path    = file.path(monitor_folder, "macros")
  params$write_path     = file.path(monitor_folder, "inputfiles")
  params$read_path      = file.path(monitor_folder, "msoc1")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }

  params$error_message <- NULL
  if (!CreateIOLocation(monitor_folder, "dplocal")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dp_local_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "rprograms")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$r_programs_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "macros")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macros_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "inputfiles")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "msoc1")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  params <- add_to_log(params, "PrepareDataLinear.a23, PrepareFolderLinear.a2", 0, 0, 0, 0)
  return(params)
}


PrepareFolderLinear.b2 = function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.b2\n\n")

  params$dp_local_path   = file.path(monitor_folder, "dplocal")
  params$r_programs_path = file.path(monitor_folder, "rprograms")
  params$macros_path    = file.path(monitor_folder, "macros")
  params$write_path     = file.path(monitor_folder, "msoc")
  params$read_path      = file.path(monitor_folder, "inputfiles")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }

  params$error_message <- NULL
  if (!CreateIOLocation(monitor_folder, "dplocal")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dp_local_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "rprograms")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$r_programs_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "macros")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macros_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!CreateIOLocation(monitor_folder, "inputfiles")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  Sys.sleep(1)
  DeleteTrigger("files_done.ok", params$read_path)

  params <- add_to_log(params, "PrepareDataLinear.b23, PrepareFolderLinear.b2", 0, 0, 0, 0)
  return(params)
}

#' @importFrom stats model.matrix sd
PrepareDataLinear.a23 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.a23\n\n")

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
  workdata$meansy = means[2]
  workdata$sdy    = sd[2]
  workdata$means  = means[covariate_index]
  workdata$sd     = sd[covariate_index]
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
PrepareDataLinear.b23 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.b23\n\n")

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

PrepareParamsLinear.b2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.b2\n\n")
  params$failed         = FALSE
  params$halted         = FALSE
  params$singular_matrix = FALSE

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

  write_time = proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLinear.b2", 0, 0, write_time, write_size)
  return(params)
}


PrepareParamsLinear.a2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.a2\n\n")

  params$halted          = FALSE
  params$singular_matrix  = FALSE
  params$pmnStepCounter  = 1
  pb                     = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb
  read_size = sum(file.size(file.path(params$read_path, "pb.rdata")))
  read_time = proc.time()[3] - read_time

  if (params$analysis != pb$analysis) {
    params$error_message <-
      paste("Party A is running", params$analysis, "regression and Party B is running", pb$analysis, "regression.")
    warning(params$error_message)
    params$failed <- TRUE
    return(params)
  }

  params$n  = nrow(data$X)
  if (pb$n != params$n) {
    params$error_message <-
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    warning(params$error_message)
    params$failed <- TRUE
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
  write_time = proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size = sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLinear.a2", read_time, read_size,
                    write_time, write_size)

  return(params)
}


PrepareBlocksLinear.a2 = function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.a2\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  = params$n
  p1 = params$p1
  p2 = params$p2

  minimumBlocksize = GetBlockSize(p1, p2)
  if (n < minimumBlocksize) {
    maxACovariates = trunc(sqrt(p2 * n) - p2 - 1)

    params$error_message <-
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
      params$error_message <-
        paste0(params$error_message,
               "\nSet the number of B covariates to be between ", minBCovariates, "and",
               paste0(maxBCovariates, "."))
    }
    warning(params$error_message)
    params$failed <- TRUE
    params <- add_to_log(params, "PrepareBlocksCox.a2", 0, 0, 0, 0)
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
  params <- add_to_log(params, "PrepareBlocksLinear.a2", 0, 0, write_time, write_size)
  return(params)
}


GetZLinear.a2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLinear.a2\n\n")
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
  params <- add_to_log(params, "GetZLinear.a2", 0, 0, write_time, write_size)
  return(params)
}


FinalizeParamsLinear.b2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsLinear.b2\n\n")
  pa = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # read pa
  read_size = sum(file.size(file.path(params$read_path, "pa.rdata")))
  read_time = proc.time()[3] - read_time
  params$p1     = pa$p1
  params$p1.old = params$p1
  params$p      = params$p1 + params$p2
  params$meansA = pa$means
  params$sdA    = pa$sd
  params$yty    = pa$yty
  params$yname  = pa$yname

  params$Acolnames = pa$Acolnames
  params <- add_to_log(params, "FinalizeParamsLinear.b2", read_time, read_size, 0, 0)
  return(params)
}


PrepareBlocksLinear.b2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.b2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  read_time = proc.time()[3]
  load(file.path(params$read_path, "blocksize.rdata")) # load blocksize
  read_size = file.size(file.path(params$read_path, "blocksize.rdata"))
  read_time = proc.time()[3] - read_time
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params <- add_to_log(params, "PrepareBlocksLinear.b2", read_time, read_size, 0, 0)
  return(params)
}


GetWLinear.b2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWLinear.b2\n\n")
  read_time  = 0
  read_size  = 0
  write_time = 0
  write_size = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "(I-Z*Z')X", params$verbose)

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

  params <- add_to_log(params, "GetWLinear.b2", read_time, read_size, write_time, write_size)

  return(params)
}


GetProductsLinear.a2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLinear.a2\n\n")
  n  = params$n
  p1 = params$p1
  p2 = params$p2
  XBTXB = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path, "xbtxb.rdata"))
  read_size = file.size(file.path(params$read_path, "xbtxb.rdata"))
  read_time = proc.time()[3] - read_time

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

  # lasso: x is standardized but needs to be divided by sqrt(n - 1),
  # y is standardized
  XTXLasso = XTX / (n - 1)
  XTYLasso = params$sdy * XTY / sqrt(n - 1)

  params$xtx = XTX
  params$xty = XTY
  params$xtxLasso = XTXLasso
  params$xtyLasso = XTYLasso

  params$converged = TRUE

  params <- add_to_log(params, "GetProductsLinear.a2", read_time, read_size, 0, 0)
  return(params)
}

#' @importFrom  stats pf pt
ComputeResultsLinear.a2 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLinear.a2\n\n")
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
  write_time = proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size = file.size(file.path(params$write_path, "stats.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeResultsLinear.a2", 0, 0, write_time, write_size)
  return(params)
}


GetResultsLinear.b2 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLinear.b2\n\n")
  params$converged = TRUE
  stats = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size = file.size(file.path(params$read_path, "stats.rdata"))
  read_time = proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsLinear.b2", read_time, read_size, 0, 0)
  return(params)
}





############################## PARENT FUNCTIONS ###############################

PartyAProcess2Linear = function(data,
                                yname                 = NULL,
                                monitor_folder         = NULL,
                                msreqid               = "v_default_00_0000",
                                blocksize             = NULL,
                                sleep_time             = 10,
                                maxWaitingTime        = 24 * 60 * 60,
                                popmednet             = TRUE,
                                trace                 = FALSE,
                                verbose               = TRUE) {

  params <- PrepareParams.2p("linear", "A", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.2p(params)
  params <- InitializeStamps.2p(params)
  params <- InitializeTrackingTable.2p(params)
  Header(params)

  params   = PrepareFolderLinear.a2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLinear.a23(params, data, yname)

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

  params <- PrepareParamsLinear.a2(params, data)

  if (params$failed) {   # Check for failed from PrepareParamsLinear.a2()
    params$completed = TRUE
    MakeErrorMessage(params$write_path, params$error_message)
    files = c("error_message.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- PrepareBlocksLinear.a2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksCox.a2()
    params$completed = TRUE
    MakeErrorMessage(params$write_path, params$error_message)
    files = c("error_message.rdata")
    params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- GetZLinear.a2(params, data)

  files = c("pa.rdata", "blocksize.rdata",
            SeqZW("cz_", length(params$container$filebreak.Z)))
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  params$completed = TRUE
  params <- GetProductsLinear.a2(params, data)
  params <- ComputeResultsLinear.a2(params, data)
  files = c("stats.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time = sleep_time)
  params <- SendPauseQuit.2p(params, sleep_time = sleep_time)
  SummarizeLog.2p(params)
  return(params$stats)
}

PartyBProcess2Linear = function(data,
                                monitor_folder       = NULL,
                                sleep_time           = 10,
                                maxWaitingTime      = 24 * 60 * 60,
                                popmednet           = TRUE,
                                trace               = FALSE,
                                verbose             = TRUE) {
  params <- PrepareParams.2p("linear", "B",
                            popmednet = popmednet, trace = trace,
                            verbose = verbose)
  params <- InitializeLog.2p(params)
  params <- InitializeStamps.2p(params)
  params <- InitializeTrackingTable.2p(params)

  Header(params)
  params   = PrepareFolderLinear.b2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLinear.b23(params, data)

  if (data$failed) { # Check for Error from prepare_data_cox_b2()
    params$completed = TRUE
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$write_path, message)
    files = c("error_message.rdata")
    params <- SendPauseQuit.2p(params, files, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params   = PrepareParamsLinear.b2(params, data)

  files = c("pb.rdata")
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(ReadErrorMessage(params$read_path))
    params <- SendPauseQuit.2p(params, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params <- FinalizeParamsLinear.b2(params, data)
  params <- PrepareBlocksLinear.b2(params)
  params <- GetWLinear.b2(params, data)

  files = c("xbtxb.rdata", SeqZW("cw_", length(params$container$filebreak.W)))
  params <- SendPauseContinue.2p(params, files, sleep_time, maxWaitingTime)

  params <- GetResultsLinear.b2(params)
  params$completed = TRUE

  params <- SendPauseQuit.2p(params, sleep_time = sleep_time)
  return(params$stats)
}
