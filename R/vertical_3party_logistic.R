################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

CheckColinearityLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.t3\n\n")
  xtx = params$xtx
  xty = params$xty

  nrow = nrow(xtx)
  indicies = c(1)
  for (i in 2:nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  xtx = xtx[indicies, indicies, drop = FALSE]
  xty = xty[indicies, drop = FALSE]

  Anames               = params$colnamesA
  Bnames               = params$colnamesB
  Aindex               = which(indicies <= length(Anames))
  params$IndiciesKeep  = indicies
  params$AIndiciesKeep = indicies[Aindex]
  params$BIndiciesKeep = indicies[-Aindex] - length(Anames)
  AnamesKeep           = Anames[params$AIndiciesKeep]
  BnamesKeep           = Bnames[params$BIndiciesKeep]
  params$colnamesA.old = params$colnamesA
  params$colnamesB.old = params$colnamesB
  params$colnamesA     = AnamesKeep
  params$colnamesB     = BnamesKeep
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
  params$xtx           = xtx
  params$xty           = xty

  Aindicies = params$AIndiciesKeep
  Bindicies = params$BIndiciesKeep

  write_time = proc.time()[3]
  save(Aindicies, file = file.path(params$write_path, "Aindicies.rdata"))
  save(Bindicies, file = file.path(params$write_path, "Bindicies.rdata"))
  write_size = sum(file.size(file.path(params$write_path, c("Aindicies.rdata",
                                                          "Bindicies.rdata"))))
  write_time = proc.time()[3] - write_time

  Btags = params$Btags[params$BIndiciesKeep]
  Atags = params$Atags[params$AIndiciesKeep][-1]

  if ((length(unique(Atags)) == 1) | (length(unique(Atags)) >= 2 & !("numeric" %in% names(Atags)))) {
    params$failed <- TRUE
    params$error_message <- "A must have no covariates or at least 2 covariates at least one of which is continuous."
  } else if (length(unique(Btags)) < 2) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has 1 or fewer covariates."
  } else if (!("numeric" %in% names(Btags))) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has no continuous covariates."
  }

  params <- add_to_log(params, "CheckColinearityLogistic.t3", 0, 0, write_time, write_size)
  return(params)
}


ComputeInitialBetasLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.t3\n\n")
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
  params$converged = FALSE
  converged = FALSE
  maxIterExceeded = FALSE

  write_time = proc.time()[3]
  save(Abetas, file = file.path(params$write_path, "betasA.rdata"))
  save(p2, Axty,   file = file.path(params$write_path, "Axty.rdata"))
  save(Bbetas, file = file.path(params$write_path, "betasB.rdata"))
  save(Bxty,   file = file.path(params$write_path, "Bxty.rdata"))
  save(converged, maxIterExceeded,
       file = file.path(params$write_path, "converged.rdata"))
  write_size = sum(file.size(file.path(params$write_path, c("betasA.rdata",
                                                          "betasB.rdata",
                                                          "Axty.rdata",
                                                          "Bxty.rdata",
                                                          "converged.rdata"))))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeInitialBetasLogistic.t3", 0, 0, write_time, write_size)

  return(params)
}


UpdateParamsLogistic.a3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.a3\n\n")
  Aindicies = NULL
  Axty      = NULL
  p2        = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "Aindicies.rdata"))
  load(file.path(params$read_path[["T"]], "Axty.rdata"))
  read_size = sum(file.size(file.path(params$read_path[["T"]], c("Aindicies.rdata",
                                                               "Axty.rdata"))))

  read_time = proc.time()[3] - read_time
  params$colnamesA.old = params$colnamesA
  params$colnamesA     = params$colnamesA.old[Aindicies]
  params$p.old         = params$p
  params$p             = length(Aindicies)
  params$p2            = p2
  params$AIndiciesKeep = Aindicies
  params$means         = params$means[Aindicies]
  params$sd            = params$sd[Aindicies]
  params$Axty          = Axty
  params <- add_to_log(params, "UpdateParamsLogistic.a3, UpdateDataLogistic.a3", read_time, read_size, 0, 0)
  return(params)
}


UpdateParamsLogistic.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.b3\n\n")
  Bindicies = NULL
  Bxty      = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "Bindicies.rdata"))
  load(file.path(params$read_path[["T"]], "Bxty.rdata"))
  read_size = sum(file.size(file.path(params$read_path[["T"]], c("Bindicies.rdata",
                                                               "Bxty.rdata"))))

  read_time = proc.time()[3] - read_time
  params$colnamesB.old = params$colnamesB
  params$colnamesB     = params$colnamesB.old[Bindicies]
  params$p.old         = params$p
  params$p             = length(Bindicies)
  params$BIndiciesKeep = Bindicies
  params$means         = params$means[Bindicies]
  params$sd            = params$sd[Bindicies]
  params$Bxty          = Bxty
  params <- add_to_log(params, "UpdateParamsLogistic.b3, UpdateDataLogistic.b3", read_time, read_size, 0, 0)
  return(params)
}


UpdateDataLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.a3\n\n")
  data$X = as.matrix(data$X[, params$AIndiciesKeep, drop = FALSE])
  return(data)
}


UpdateDataLogistic.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.b3\n\n")
  data$X = as.matrix(data$X[, params$BIndiciesKeep, drop = FALSE])
  return(data)
}


GetBetaALogistic.a3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.a3\n\n")
  converged       = NULL
  maxIterExceeded = NULL
  Abetas          = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasA.rdata"))
  read_size = sum(file.size(file.path(params$read_path[["T"]], c("converged.rdata",
                                                               "betasA.rdata"))))
  read_time = proc.time()[3] - read_time
  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  params$betas = Abetas
  params <- add_to_log(params, "GetBetaALogistic.a3", read_time, read_size, 0, 0)
  return(params)
}


GetBetaBLogistic.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.b3\n\n")
  converged       = NULL
  maxIterExceeded = NULL
  Bbetas          = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasB.rdata"))
  read_size = sum(file.size(file.path(params$read_path[["T"]], c("converged.rdata",
                                                               "betasB.rdata"))))
  read_time = proc.time()[3] - read_time
  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  params$betas = Bbetas
  params <- add_to_log(params, "GetBetaBLogistic.b3", read_time, read_size, 0, 0)
  return(params)
}


GetXAbetaALogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaALogistic.a3\n\n")
  XAbeta = data$X %*% params$betas

  write_time = proc.time()[3]
  save(XAbeta, file = file.path(params$write_path, "xabeta.rdata"))
  write_size = file.size(file.path(params$write_path, "xabeta.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetXAbetaALogistic.a3", 0, 0, write_time, write_size)
  return(params)
}


GetXBbetaBLogistic.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaBLogistic.b3\n\n")
  XBbeta = data$X %*% params$betas

  write_time = proc.time()[3]
  save(XBbeta, file = file.path(params$write_path, "xbbeta.rdata"))
  write_size = file.size(file.path(params$write_path, "xbbeta.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetXBbetaBLogistic.b3", 0, 0, write_time, write_size)
  return(params)
}


GetWeightsLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetWeightsLogistic.t3\n\n")
  XAbeta = NULL
  XBbeta = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["A"]], "xabeta.rdata"))  # Load XbetaB
  load(file.path(params$read_path[["B"]], "xbbeta.rdata"))  # Load XbetaB
  read_size = file.size(file.path(params$read_path[["A"]], "xabeta.rdata")) +
    file.size(file.path(params$read_path[["B"]], "xbbeta.rdata"))
  read_time = proc.time()[3] - read_time

  Xbeta = XAbeta + XBbeta
  pi_ = (1 + exp(-Xbeta))^(-1)
  params$pi_ = pi_

  write_time = proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi.rdata"))
  write_size = file.size(file.path(params$write_path, "pi.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "GetWeightsLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}


GetRVLogistic.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetRVLogistic.b3\n\n")
  pi_ = NULL
  write_time = 0
  write_size = 0
  read_time  = proc.time()[3]
  load(file.path(params$read_path[["T"]], "pi.rdata"))
  read_size  = file.size(file.path(params$read_path[["T"]], "pi.rdata"))
  read_time  = proc.time()[3] - read_time

  params$pi_ = pi_
  W = pi_ * (1 - params$pi_)
  XBTWXB = 0
  pbar = MakeProgressBar1(params$blocks$numBlocks, "R(I-Z*Z')W*XB", params$verbose)
  containerCt.RZ = 0
  containerCt.RV = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.RZ) {
      containerCt.RZ = containerCt.RZ + 1
      filename1 = paste0("crz_", containerCt.RZ, ".rdata")
      toRead = file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.RV) {
      containerCt.RV = containerCt.RV + 1
      filename2 = paste0("crv_", containerCt.RV, ".rdata")
      toWrite = file(file.path(params$write_path, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1

    Xblock  = data$X[strt:stp, , drop = FALSE]
    Wblock  = W[strt:stp]
    WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

    read_time = read_time - proc.time()[3]
    RZ = matrix(readBin(con = toRead, what = numeric(), n = n * n,
                        endian = "little"), nrow = n, ncol = n)
    read_time = read_time + proc.time()[3]

    RV = RZ %*% WXblock

    write_time = write_time - proc.time()[3]
    writeBin(as.vector(RV), con = toWrite, endian = "little")
    write_time = write_time + proc.time()[3]

    XBTWXB = XBTWXB + t(Xblock) %*% WXblock

    if ((i + 1) %in% params$container$filebreak.RZ || i == params$blocks$numBlocks) {
      close(toRead)
      read_size = read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.RV || i == params$blocks$numBlocks) {
      close(toWrite)
      write_size = write_size + file.size(file.path(params$write_path, filename2))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  write_time = write_time - proc.time()[3]
  save(XBTWXB, file = file.path(params$write_path, "xbtwxb.rdata"))
  write_size = write_size + sum(file.size(c(file.path(params$write_path, "xbtwxb.rdata"))))
  write_time = write_time + proc.time()[3]

  params <- add_to_log(params, "GetRVLogistic.b3", read_time, read_size, write_time, write_size)
  return(params)
}


ProcessVLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessVLogistic.t3\n\n")
  XBTWXB = NULL
  write_time = 0
  write_size = 0
  p2 = params$p2
  read_time = proc.time()[3]
  load(file.path(params$read_path[["B"]], "xbtwxb.rdata"))
  read_size = file.size(file.path(params$read_path[["B"]], "xbtwxb.rdata"))
  read_time = proc.time()[3] - read_time

  params$xbtwxb = XBTWXB

  numBlocks = params$blocks$numBlocks
  pbar = MakeProgressBar1(numBlocks, "(I-Z*Z')W*XB*R", params$verbose)

  containerCt.RV = 0
  containerCt.VR = 0
  for (i in 1:numBlocks) {
    if (i %in% params$container$filebreak.RV) {
      containerCt.RV = containerCt.RV + 1
      filename2 = paste0("crv_", containerCt.RV, ".rdata")
      toRead2 = file(file.path(params$read_path[["B"]], filename2), "rb")
    }
    if (i %in% params$container$filebreak.VR) {
      containerCt.VR = containerCt.VR + 1
      filename3 = paste0("cvr_", containerCt.VR, ".rdata")
      toWrite3 = file(file.path(params$write_path, filename3), "wb")
    }

    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1

    filename1 = paste0("r1_", i, ".rdata")
    filename4 = paste0("r3_", i, ".rdata")

    read_time = read_time - proc.time()[3]
    toRead1 = file(file.path(params$dplocalPath, filename1), "rb")
    R2 = matrix(readBin(con = toRead1, what = numeric(), n = n * n,
                        endian = "little"), nrow = n, ncol = n)
    read_size = read_size + file.size(file.path(params$dplocalPath, filename1))
    close(toRead1)
    RV = matrix(readBin(con = toRead2, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time = read_time + proc.time()[3]

    V = t(R2) %*% RV
    R3 = RandomOrthonomalMatrix(p2)
    VR = V %*% R3

    write_time = write_time - proc.time()[3]
    toWrite4 = file(file.path(params$dplocalPath, filename4), "wb")
    writeBin(as.vector(R3), con = toWrite4, endian = "little")
    close(toWrite4)
    write_size = write_size + file.size(file.path(params$dplocalPath, filename4))
    writeBin(as.vector(VR), con = toWrite3, endian = "little")
    write_time = write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.RV || i == numBlocks) {
      close(toRead2)
      read_size = read_size + file.size(file.path(params$dplocalPath, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.VR || i == numBlocks) {
      close(toWrite3)
      write_size = write_size + file.size(file.path(params$write_path, filename3))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "ProcessVLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}


GetXRLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXRLogistic.a3\n\n")
  pi_ = NULL
  p2 = params$p2
  write_time = 0
  write_size = 0
  read_time  = proc.time()[3]
  load(file.path(params$read_path[["T"]], "pi.rdata"))
  read_size  = file.size(file.path(params$read_path[["T"]], "pi.rdata"))
  read_time  = proc.time()[3] - read_time

  params$pi_ = pi_
  W = pi_ * (1 - params$pi_)
  XATWXA = 0

  pbar = MakeProgressBar1(params$blocks$numBlocks, "XA'(I-Z*Z')W*XB*R", params$verbose)

  containerCt.VR = 0
  containerCt.XR = 0
  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.RV) {
      containerCt.VR = containerCt.VR + 1
      filename1 = paste0("cvr_", containerCt.VR, ".rdata")
      toRead = file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.XR) {
      containerCt.XR = containerCt.XR + 1
      filename2 = paste0("cxr_", containerCt.XR, ".rdata")
      toWrite = file(file.path(params$write_path, filename2), "wb")
    }
    strt = params$blocks$starts[i]
    stp = params$blocks$stops[i]
    n = stp - strt + 1

    Xblock  = data$X[strt:stp, , drop = FALSE]
    Wblock  = W[strt:stp]
    WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

    read_time = read_time - proc.time()[3]
    VR = matrix(readBin(con = toRead, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time = read_time + proc.time()[3]

    XR = t(Xblock) %*% VR

    write_time = write_time - proc.time()[3]
    writeBin(as.vector(XR), con = toWrite, endian = "little")
    write_time = write_time + proc.time()[3]

    XATWXA = XATWXA + t(Xblock) %*% WXblock

    if ((i + 1) %in% params$container$filebreak.VR || i == params$blocks$numBlocks) {
      close(toRead)
      read_size = read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$numBlocks) {
      close(toWrite)
      write_size = write_size + file.size(file.path(params$write_path, filename2))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  write_time = write_time - proc.time()[3]
  save(XATWXA, file = file.path(params$write_path, "xatwxa.rdata"))
  write_size = write_size + sum(file.size(c(file.path(params$write_path, "xatwxa.rdata"))))
  write_time = write_time + proc.time()[3]

  params <- add_to_log(params, "GetXRLogistic.a3", read_time, read_size, write_time, write_size)
  return(params)
}


ProcessXtWXLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessXtWXLogistic.t3\n\n")
  XATWXA = NULL
  p1 = params$p1
  p2 = params$p2

  read_time = proc.time()[3]
  load(file.path(params$read_path[["A"]], "xatwxa.rdata"))
  read_size = file.size(file.path(params$read_path[["A"]], "xatwxa.rdata"))
  read_time = proc.time()[3] - read_time

  params$xatwxa = XATWXA

  pbar = MakeProgressBar1(params$blocks$numBlocks, "X'W*X", params$verbose)
  containerCt.XR = 0
  XATWXB = 0

  for (i in 1:params$blocks$numBlocks) {
    if (i %in% params$container$filebreak.XR) {
      containerCt.XR = containerCt.XR + 1
      filename1 = paste0("cxr_", containerCt.XR, ".rdata")
      toRead = file(file.path(params$read_path[["A"]], filename1), "rb")
    }

    filename2 = paste0("r3_", i, ".rdata")
    read_time = read_time - proc.time()[3]
    toRead1 = file(file.path(params$dplocalPath, filename2), "rb")
    R = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
                       endian = "little"), nrow = p2, ncol = p2)
    close(toRead1)
    XR = matrix(readBin(con = toRead, what = numeric(), n = p1 * p2,
                        endian = "little"), nrow = p1, ncol = p2)
    read_size = read_size + file.size(file.path(params$dplocalPath, filename2))
    read_time = read_time + proc.time()[3]

    XATWXB = XATWXB + XR %*% t(R)

    if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$numBlocks) {
      close(toRead)
      read_size = read_size + file.size(file.path(params$read_path[["A"]], filename1))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  xtwx = rbind(cbind(params$xatwxa, XATWXB), cbind(t(XATWXB), params$xbtwxb))
  params$xtwx = xtwx

  II = NULL
  tryCatch({
    II = solve(xtwx)
  },
  error = function(err) {
    II = NULL
  }
  )
  if (is.null(II)) {
    params$failed <- TRUE
    params$singular_matrix = TRUE
    params$error_message <-
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
    return(params)
  }
  params$II = II
  IIA = II[, 1:p1, drop = FALSE]
  IIB = II[, (p1 + 1):(p1 + p2), drop = FALSE]
  write_time = proc.time()[3]
  save(IIA, file = file.path(params$write_path, "IIA.rdata"))
  save(IIB, file = file.path(params$write_path, "IIB.rdata"))
  write_size = sum(file.size(file.path(params$write_path, c("IIA.rdata", "IIB.rdata"))))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "ProcessXtWXLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateBetaLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.a3\n\n")
  IIA = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "IIA.rdata"))
  read_size = file.size(file.path(params$read_path[["T"]], "IIA.rdata"))
  read_time = proc.time()[3] - read_time

  IA = params$Axty - t(data$X) %*% params$pi_
  AI = IIA %*% IA


  write_time = proc.time()[3]
  save(AI, file = file.path(params$write_path, "AI.rdata"))
  write_size = file.size(file.path(params$write_path, "AI.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "UpdateBetaLogistic.a3", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateBetaLogistic.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.b3\n\n")
  IIB = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "IIB.rdata"))
  read_size = file.size(file.path(params$read_path[["T"]], "IIB.rdata"))
  read_time = proc.time()[3] - read_time

  IB = params$Bxty - t(data$X) %*% params$pi_
  BI = IIB %*% IB

  write_time = proc.time()[3]
  save(BI, file = file.path(params$write_path, "BI.rdata"))
  write_size = file.size(file.path(params$write_path, "BI.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "UpdateBetaLogistic.b3", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateBetaLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.t3\n\n")
  AI = NULL
  BI = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["A"]], "AI.rdata"))
  load(file.path(params$read_path[["B"]], "BI.rdata"))
  read_size = file.size(file.path(params$read_path[["A"]], "AI.rdata")) +
    file.size(file.path(params$read_path[["B"]], "BI.rdata"))
  read_time = proc.time()[3] - read_time

  delta = AI + BI
  betas = params$betas + delta
  params$betas = betas
  converged = all(abs(delta) / (abs(betas) + .1) < params$cutoff)
  maxIterExceeded = (params$algIterationCounter >= params$maxIterations) && !converged

  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  Abetas = betas[1:params$p1]
  Bbetas = betas[(params$p1 + 1):(params$p1 + params$p2)]

  write_time = proc.time()[3]
  save(converged, maxIterExceeded,
       file = file.path(params$write_path, "converged.rdata"))
  save(Abetas, file = file.path(params$write_path, "betasA.rdata"))
  save(Bbetas, file = file.path(params$write_path, "betasB.rdata"))
  write_size = sum(file.size(file.path(params$write_path, c("betasA.rdata",
                                                          "betasB.rdata",
                                                          "converged.rdata"))))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "UpdateBetaLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}


GetFinalBetaLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.a3\n\n")
  betas = params$betas / params$sd
  offsetA = sum(betas[-1] * params$means[-1])
  AFinalFitted = t(params$sd * t(data$X) + params$means) %*% betas -
    t(params$sd[1] * t(data$X[, 1]) + params$means[1]) %*% betas[1]
  write_time = proc.time()[3]
  save(offsetA, AFinalFitted, file = file.path(params$write_path, "Afinalfitted.rdata"))
  write_size = file.size(file.path(params$write_path, "Afinalfitted.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "GetFinalBetaLogistic.a3", 0, 0, write_time, write_size)
  return(params)
}


GetFinalBetaLogistic.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.b3\n\n")
  betas = params$betas / params$sd
  offsetB = sum(betas * params$means)
  BFinalFitted = t(params$sd * t(data$X) + params$means) %*% betas

  write_time = proc.time()[3]
  save(offsetB, BFinalFitted, file = file.path(params$write_path, "Bfinalfitted.rdata"))
  write_size = file.size(file.path(params$write_path, "Bfinalfitted.rdata"))
  write_time = proc.time()[3] - write_time
  params <- add_to_log(params, "GetFinalBetaLogistic.b3", 0, 0, write_time, write_size)
  return(params)
}


GetFinalFittedLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetFinalFittedLogistic.t3\n\n")
  offsetA = NULL
  offsetB = NULL
  AFinalFitted = NULL
  BFinalFitted = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["A"]], "Afinalfitted.rdata"))
  load(file.path(params$read_path[["B"]], "Bfinalfitted.rdata"))
  read_size = file.size(file.path(params$read_path[["A"]], "Afinalfitted.rdata")) +
    file.size(file.path(params$read_path[["B"]], "Bfinalfitted.rdata"))
  read_time = proc.time()[3] - read_time

  betas = params$betas / c(params$sdA, params$sdB)
  betas[1] = betas[1] - offsetA - offsetB

  finalFitted = AFinalFitted + BFinalFitted + betas[1]

  params$betas = betas
  params$finalFitted = finalFitted

  write_time = proc.time()[3]
  save(finalFitted, file = file.path(params$write_path, "finalFitted.rdata"))
  write_size = file.size(file.path(params$write_path, "finalFitted.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "GetFinalBetaLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}

ComputeResultsLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.a3\n\n")
  finalFitted = NULL

  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "finalFitted.rdata"))
  read_size = file.size(file.path(params$read_path[["T"]], "finalFitted.rdata"))
  read_time = proc.time()[3] - read_time

  n = params$n
  ct      = sum(data$Y)
  params$FinalFitted = finalFitted
  resdev  = -2 * (sum(data$Y * finalFitted) - sum(log(1 + exp(finalFitted))))
  nulldev = -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))
  hoslem  = HoslemInternal(params, data)
  ROC     = RocInternal(params, data)

  write_time = proc.time()[3]
  save(resdev, nulldev, hoslem, ROC, file = file.path(params$write_path, "logisticstats.rdata"))
  write_size = file.size(file.path(params$write_path, "logisticstats.rdata"))
  write_time = proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeResultsLogistic.a3", read_time, read_size, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
ComputeResultsLogistic.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.t3\n\n")
  nulldev = NULL
  resdev  = NULL
  hoslem  = NULL
  ROC     = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["A"]], "logisticstats.rdata"))
  read_size = file.size(file.path(params$read_path[["A"]], "logisticstats.rdata"))
  read_time = proc.time()[3] - read_time

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
  Anames = params$colnamesA.old
  Bnames = params$colnamesB.old
  p1.old = params$p1.old
  p2.old = params$p2.old
  p.old  = params$p.old
  indicies = params$IndiciesKeep

  # If xtwx were singular, it would have been caught in GetII.A2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  cov1 = solve(params$xtwx)
  secoef = sqrt(diag(cov1)) / c(sdA, sdB)
  tmp = matrix(c(1, (-meansA / sdA)[-1], -meansB / sdB), ncol = 1)
  secoef[1] = sqrt(t(tmp) %*% cov1 %*% tmp)

  stats$party = c(rep("dp1", p1.old), rep("dp2", p2.old))
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
  stats$coefficients[indicies] = params$betas
  stats$secoef[indicies] = secoef
  tvals = params$betas / secoef
  pvals = 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[indicies] = tvals
  stats$pvals[indicies] = pvals
  stats$hoslem  = hoslem
  stats$ROC     = ROC
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

  params$stats      = stats

  params <- add_to_log(params, "ComputeResultsLogistic.t3", read_time, read_size, write_time, write_size)
  return(params)
}


GetResultsLogistic.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.a3\n\n")
  stats = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size = file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time = proc.time()[3] - read_time
  stats$Y           = data$Y # For Hoslem and ROC
  stats$FinalFitted = params$FinalFitted
  params$stats      = stats
  params <- add_to_log(params, "GetResultsLogistic.a3", read_time, read_size, 0, 0)
  return(params)
}


GetResultsLogistic.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.b3\n\n")
  stats = NULL
  read_time = proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size = file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time = proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsLogistic.b3", read_time, read_size, 0, 0)
  return(params)
}


############################### PARENT FUNCTIONS ###############################


PartyAProcess3Logistic = function(data,
                                  yname          = NULL,
                                  monitor_folder  = NULL,
                                  sleep_time      = 10,
                                  maxWaitingTime = 24 * 60 * 60,
                                  popmednet      = TRUE,
                                  trace          = FALSE,
                                  verbose        = TRUE) {

  params <- PrepareParams.3p("logistic", "A",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.3p(params)
  params <- InitializeStamps.3p(params)
  params <- InitializeTrackingTable.3p(params)
  Header(params)

  params   = PrepareFolderLinear.a3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLogistic.A23(params, data, yname)
  params <- add_to_log(params, "PrepareDataLogistic.A23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party A."
    MakeErrorMessage(params$write_path, message)
    files = c("error_message.rdata")
    params <- SendPauseQuit.3p(params, filesT = files, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- PrepareParamsLinear.a3(params, data)
  files = "pa.rdata"
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["T"]]))
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params$algIterationCounter = 1
  params <- PrepareBlocksLinear.a3(params)

  params <- GetZLinear.a3(params, data)
  files = SeqZW("cz_", length(params$container$filebreak.Z))
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- GetWRLinear.a3(params, data)
  files = c("xatxa.rdata", SeqZW("cpr_", length(params$container$filebreak.PR)))
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["T"]]))
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- UpdateParamsLogistic.a3(params)
  data = UpdateDataLogistic.a3(params, data)
  params <- GetBetaALogistic.a3(params)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetXAbetaALogistic.a3(params, data)
    files = c("xabeta.rdata")
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    params <- GetXRLogistic.a3(params, data)
    files = c("xatwxa.rdata", SeqZW("cxr_", length(params$container$filebreak.XR)))
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
      warning(ReadErrorMessage(params$read_path[["T"]]))
      params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
      return(params$stats)
    }

    params <- UpdateBetaLogistic.a3(params, data)
    files = c("ai.rdata")
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    params <- GetBetaALogistic.a3(params)
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- GetFinalBetaLogistic.a3(params, data)
  files = "Afinalfitted.rdata"
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  params <- ComputeResultsLogistic.a3(params, data)
  files = c("logisticstats.rdata")
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- GetResultsLogistic.a3(params, data)

  params <- SendPauseQuit.3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


PartyBProcess3Logistic = function(data,
                                  monitor_folder  = NULL,
                                  sleep_time      = 10,
                                  maxWaitingTime = 24 * 60 * 60,
                                  popmednet      = TRUE,
                                  trace          = FALSE,
                                  verbose        = TRUE) {
  params <- PrepareParams.3p("logistic", "B",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.3p(params)
  params <- InitializeStamps.3p(params)
  params <- InitializeTrackingTable.3p(params)

  Header(params)
  params   = PrepareFolderLinear.b3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  data = PrepareDataLogistic.B23(params, data)
  params <- add_to_log(params, "PrepareDataLogistic.B23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party B."
    MakeErrorMessage(params$write_path, message)
    files = c("error_message.rdata")
    params <- SendPauseQuit.3p(params, filesT = files, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- PrepareParamsLinear.b3(params, data)
  files = "pb.rdata"
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["T"]]))
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params$algIterationCounter = 1
  params <- PrepareBlocksLinear.b3(params)

  params <- GetRWLinear.b3(params, data)
  files = c("xbtxb.rdata", SeqZW("crw_", length(params$container$filebreak.RW)))
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["T"]]))
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- UpdateParamsLogistic.b3(params)
  data = UpdateDataLogistic.b3(params, data)
  params <- GetBetaBLogistic.b3(params)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetXBbetaBLogistic.b3(params, data)
    files = c("xbbeta.rdata")
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    params <- GetRVLogistic.b3(params, data)
    files = c("xbtwxb.rdata", SeqZW("crv_", length(params$container$filebreak.RV)))
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
      warning(ReadErrorMessage(params$read_path[["T"]]))
      params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
      return(params$stats)
    }

    params <- UpdateBetaLogistic.b3(params, data)
    files = c("bi.rdata")
    params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    params <- GetBetaBLogistic.b3(params)
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- GetFinalBetaLogistic.b3(params, data)
  files = "Bfinalfitted.rdata"
  params <- SendPauseContinue.3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  params <- GetResultsLogistic.b3(params)

  params <- SendPauseQuit.3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


PartyTProcess3Logistic = function(monitor_folder         = NULL,
                                  msreqid               = "v_default_0_000",
                                  blocksize             = 500,
                                  cutoff                = 1e-8,
                                  maxIterations         = 25,
                                  sleep_time             = 10,
                                  maxWaitingTime        = 24 * 60 * 60,
                                  popmednet             = TRUE,
                                  trace                 = FALSE,
                                  verbose               = TRUE) {
  params <- PrepareParams.3p("logistic", "T", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- InitializeLog.3p(params)
  params <- InitializeStamps.3p(params)
  params <- InitializeTrackingTable.3p(params)

  Header(params)
  params   = PrepareFolderLinear.t3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.3p(params, from = c("A", "B"), maxWaitingTime = maxWaitingTime)

  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata")) &&
      file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["A"]]))
    warning(ReadErrorMessage(params$read_path[["B"]]))
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["A"]]))
    file.copy(file.path(params$read_path[["A"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files = "error_message.rdata"
    params <- SendPauseContinue.3p(params, filesB = files, from = "B",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(ReadErrorMessage(params$read_path[["B"]]))
    file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files = "error_message.rdata"
    params <- SendPauseContinue.3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params   = PrepareParamsLinear.t3(params, cutoff, maxIterations)

  if (!params$failed) params <- PrepareBlocksLinear.t3(params, blocksize)

  if (params$failed) {
    warning(params$error_message)
    MakeErrorMessage(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- SendPauseContinue.3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files = "blocks.rdata"
  params <- SendPauseContinue.3p(params, filesA = files, from = "A",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params$algIterationCounter = 1
  params <- ProcessZLinear.t3(params)
  files = c("blocks.rdata", SeqZW("crz_", length(params$container$filebreak.RZ)))
  params <- SendPauseContinue.3p(params, filesB = files, from  = "B",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- ProcessWLinear.t3(params)
  files = c("p2.rdata", SeqZW("cwr_", length(params$container$filebreak.WR)))
  params <- SendPauseContinue.3p(params, filesA = files, from  = "A",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- GetProductsLinear.t3(params)

  params <- CheckColinearityLogistic.t3(params)

  if (params$failed) {
    warning(params$error_message)
    MakeErrorMessage(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- SendPauseContinue.3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params <- ComputeInitialBetasLogistic.t3(params)
  filesA = c("Aindicies.rdata", "betasA.rdata", "Axty.rdata", "converged.rdata")
  filesB = c("Bindicies.rdata", "betasB.rdata", "Bxty.rdata", "converged.rdata")
  params <- SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetWeightsLogistic.t3(params)
    files = "pi.rdata"
    params <- SendPauseContinue.3p(params, filesB = files, from  = "B",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    params <- ProcessVLogistic.t3(params)
    files = c("pi.rdata", SeqZW("cvr_", length(params$container$filebreak.RV)))
    params <- SendPauseContinue.3p(params, filesA = files, from  = "A",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    params <- ProcessXtWXLogistic.t3(params)

    if (params$failed) {
      warning(params$error_message)
      MakeErrorMessage(params$write_path, params$error_message)
      files = "error_message.rdata"
      params <- SendPauseContinue.3p(params, filesA = files, filesB = files,
                                    from = c("A", "B"),
                                    sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
      params <- SendPauseQuit.3p(params, sleep_time = sleep_time)
      SummarizeLog.3p(params)
      return(params$stats)
    }
    filesA = c("IIA.rdata")
    filesB = c("IIB.rdata")
    params <- SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    params <- UpdateBetaLogistic.t3(params)
    filesA = c("betasA.rdata", "converged.rdata")
    filesB = c("betasB.rdata", "converged.rdata")
    params <- SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- GetFinalFittedLogistic.t3(params)
  filesA = "finalfitted.rdata"
  params <- SendPauseContinue.3p(params, filesA = filesA, from  = "A",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- ComputeResultsLogistic.t3(params)
  files = "stats.rdata"
  params <- SendPauseContinue.3p(params, filesA = files, filesB = files, from  = c("A", "B"),
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- SendPauseQuit.3p(params, sleep_time = sleep_time)
  SummarizeLog.3p(params)
  return(params$stats)
}
