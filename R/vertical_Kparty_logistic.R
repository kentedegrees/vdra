################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ###################

GetProductsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLogistic.AC\n\n")
  read_time <- 0
  read_size = 0
  p = 0
  n = 0
  pi = c()

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
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "products.rdata"))
    load(file.path(params$readPathDP[id], "halfshare.rdata"))
    load(file.path(params$readPathDP[id], "colstats.rdata"))
    read_size = read_size + sum(file.size(file.path(params$readPathDP[id],
                                                  c("products.rdata",
                                                    "halfshare.rdata",
                                                    "colstats.rdata"))))
    read_time <- read_time + proc.time()[3]

    allproducts[[id]]  = products
    allhalfshare[[id]] = halfshare
    alltags[[id]]      = tags
    allcolmin          = c(allcolmin, colmin)
    allcolrange        = c(allcolrange, colrange)
    allcolsum          = c(allcolsum, colsum)
    allcolnames        = c(allcolnames, colnames)
    party              = c(party, rep(paste0("dp", id), length(colnames)))
    p = p + ncol(halfshare)
    pi = c(pi, ncol(halfshare))
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


  params$halfshare    = allhalfshare
  params$sts          = M[2:p, 2:p, drop = FALSE]
  params$sty          = M[2:p, 1, drop = FALSE]
  params$yty          = M[1, 1]
  params$meansy       = allcolsum[1] / n
  params$means        = allcolsum[-1] / n
  params$n            = n
  params$p            = p
  params$pi           = pi
  params$colmin       = allcolmin[-1]
  params$colrange     = allcolrange[-1]
  params$colsum       = allcolsum[-1]
  params$colnames     = allcolnames[-1]
  params$party        = party[-1]
  params$tags         = alltags

  params <- add_to_log(params, "GetProductsLogistic.AC", read_time, read_size, 0, 0)
  return(params)
}


CheckColinearityLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.AC\n\n")
  sts = params$sts
  sty = params$sty

  nrow = nrow(sts)
  indicies = c(1)
  for (i in 2:nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(sts[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }


  sts = sts[indicies, indicies, drop = FALSE]
  sty = sty[indicies, drop = FALSE]

  params$sts = sts
  params$sty = sty

  # Extract the indicies to keep for each party and check for errors.

  params$colmin       = params$colmin[indicies]
  params$colrange     = params$colrange[indicies]
  params$colsum       = params$colsum[indicies]
  params$fullindicies = indicies       # To be used when computing stats
  params$p            = params$p - 1   # Get rid of the response from the count

  indicies = indicies + 1  # take into account that pi still counts sty, which we removed earlier.

  params$indicies   = rep(list(list()), params$numDataPartners)
  tags              = rep(list(list()), params$numDataPartners)
  min = 1

  for (id in 1:params$numDataPartners) {
    max = min + params$pi[id] - 1
    idx = indicies[which(min <= indicies & indicies <= max)] - min + 1
    params$indicies[[id]] = idx
    if (id == 1) {
      idx = (idx - 1)[-1]
    }
    temp = params$tags[[id]]
    temp = temp[idx]
    tags[[id]] = temp
    min = max + 1
  }

  params$error_message = ""
  if ((length(unique(tags[[1]])) == 1) | (length(unique(tags[[1]])) >= 2 & !("numeric" %in% names(tags[[1]])))) {
    params$failed <- TRUE
    params$error_message = "Data Partner 1 must have no covariates or at least 2 covariates at least one of which is continuous.\n"
  }
  for (id in 2:params$numDataPartners) {
    if (length(unique(tags[[id]])) < 2) {
      params$failed <- TRUE
      params$error_message = paste0(params$error_message,
                                   paste("After removing colinear covariates, Data Partner", id, "has 1 or fewer covariates.\n"))
    } else if (!("numeric" %in% names(tags[[id]]))) {
      params$failed <- TRUE
      params$error_message = paste0(params$error_message,
                                   paste("After removing colinear covariates, Data Partner", id, "has no continuous covariates.\n"))
    }
  }

  indicies = params$indicies

  params$pReduct = c()
  for (id in 1:params$numDataPartners) {
    params$pReduct = c(params$pReduct, length(indicies[[id]]))
  }

  for (id in 1:params$numDataPartners) {
    params$halfshare[[id]] = params$halfshare[[id]][, indicies[[id]], drop = FALSE]
  }

  write_time <-proc.time()[3]
  save(indicies, file = file.path(params$write_path, "indicies.rdata"))
  write_size = file.size(file.path(params$write_path, "indicies.rdata"))
  write_time <-proc.time()[3] - write_time

  params <- add_to_log(params, "CheckColinearityLogistic.AC", 0, 0, write_time, write_size)

  return(params)
}


#' @importFrom stats runif
ComputeInitialBetasLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.AC\n\n")
  write_time <-0
  write_size = 0
  colsumS = (params$colsum - params$n * params$colmin) / params$colran
  beta = 4 * solve(params$sts) %*% (params$sty - 0.5 * colsumS)

  u = sum(runif(length(beta), min = 1, max = 5) * abs(beta))
  params$u = u
  start = 1
  for (id in 1:params$numDataPartners) {
    end = start + length(params$indicies[[id]]) - 1
    betas = beta[start:end]

    write_time <-write_time - proc.time()[3]
    save(u, betas, file = file.path(params$write_path, paste0("u_beta_", id, ".rdata")))
    write_size = write_size + file.size(file.path(params$write_path, paste0("u_beta_", id, ".rdata")))
    write_time <-write_time + proc.time()[3]
    start = end + 1
  }
  params <- add_to_log(params, "ComputeInitialBetasLogistic.AC", 0, 0, write_time, write_size)
  return(params)
}


UpdateParamsLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.DP\n\n")
  indicies = NULL
  u = NULL
  betas = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "indicies.rdata"))
  filename = paste0("u_beta_", params$data_partner_id, ".rdata")
  load(file.path(params$readPathAC, filename))
  read_size = file.size(file.path(params$readPathAC, "indicies.rdata")) +
    file.size(file.path(params$readPathAC, filename))
  read_time <- proc.time()[3] - read_time
  params$u = u
  params$betas = betas
  params$indicies = indicies
  params <- add_to_log(params, "UpdateParamsLogistic.DP", read_time, read_size, 0, 0)
  return(params)
}


UpdateDataLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.DP\n\n")
  if (params$data_partner_id == 1) {
    data$Y = data$X[, 1, drop = FALSE]
  }
  idx = params$indicies[[params$data_partner_id]]
  data$X = data$X[, idx, drop = FALSE]
  data$colmin = data$colmin[idx]
  data$colmax = data$colmax[idx]
  data$colsum = data$colsum[idx]
  data$colrange = data$colrange[idx]
  return(data)
}


#' @importFrom stats rnorm runif
ComputeSbetaLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSbetaLogistic.DP\n\n")
  set.seed(params$seed + params$algIterationCounter, kind = "Mersenne-Twister")
  V = matrix(rnorm(params$n, mean = runif(n = 1, min = -1, max = 1), sd = 10), ncol = 1)
  Vsum = 0
  for (id in 1:params$numDataPartners) {
    set.seed(params$seeds[id] + params$algIterationCounter, kind = "Mersenne-Twister")
    Vsum = Vsum + matrix(rnorm(params$n, mean = runif(n = 1, min = -1, max = 1), sd = 10), ncol = 1)
  }

  Sbeta = (data$X %*% params$betas + params$u) / (2 * params$u) + V - params$scaler / sum(params$scalers) * Vsum

  write_time <-proc.time()[3]
  save(Sbeta, file = file.path(params$write_path, "sbeta.rdata"))
  write_size = file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeSbetaLogistic.DP", 0, 0, write_time, write_size)
  return(params)
}


ComputeWeightsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeWeightsLogistic.AC\n\n")
  Sbeta = 0
  read_time <- 0
  read_size = 0
  sbeta = 0
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_time <- read_time + proc.time()[3]
    sbeta = sbeta + Sbeta
  }
  sbeta = 2 * params$u * sbeta - params$numDataPartners * params$u
  pi_ = 1 / (1 + exp(-sbeta))
  params$pi_ = pi_

  write_time <-proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi.rdata"))
  write_size = file.size(file.path(params$write_path, "pi.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "ComptueWeightsLogistic.AC", read_time, read_size, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm
ComputeStWSLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSLogistic.DP\n\n")
  pi_ = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "pi.rdata"))
  read_size = file.size(file.path(params$readPathAC, "pi.rdata"))
  read_time <- proc.time()[3] - read_time
  params$pi_ = pi_

  W = pi_ * (1 - pi_)
  C = rep(list(list()), params$numDataPartners)

  idx = params$indicies[[params$data_partner_id]]
  set.seed(params$seed, kind = "Mersenne-Twister")
  halfshare = matrix(rnorm(params$n * params$p, sd = 20),
                     nrow = params$n, ncol = params$p)[, idx, drop = FALSE]

  for (id in 1:params$numDataPartners) {
    if (id < params$data_partner_id) {
      set.seed(params$seeds[id], kind = "Mersenne-Twister")
      idx = params$indicies[[id]]
      halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
                           nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
      C[[id]] = params$scaler / (params$scaler + params$scalers[id]) *
        t(halfshareDP) %*% MultiplyDiagonalWTimesX(W, halfshare) +
        t(halfshareDP) %*% MultiplyDiagonalWTimesX(W, data$X - halfshare)
    } else if (id == params$data_partner_id) {
      C[[id]] = t(data$X) %*% MultiplyDiagonalWTimesX(W, data$X)
    } else {
      set.seed(params$seeds[id], kind = "Mersenne-Twister")
      idx = params$indicies[[id]]
      halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
                           nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
      C[[id]] = params$scaler / (params$scaler + params$scalers[id]) *
        t(halfshare) %*% MultiplyDiagonalWTimesX(W, halfshareDP) +
        t(data$X - halfshare) %*% MultiplyDiagonalWTimesX(W, halfshareDP)
    }
  }

  write_time <-proc.time()[3]
  save(C, file = file.path(params$write_path, "stwsshare.rdata"))
  write_size = file.size(file.path(params$write_path, "stwsshare.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeStWSLogistic.DP", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeStWSLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSLogistic.AC\n\n")
  read_time <- 0
  read_size = 0
  C        = NULL
  W = params$pi_ * (1 - params$pi_)
  StWS = matrix(0, sum(params$pReduct), sum(params$pReduct))

  for (id1 in 1:params$numDataPartners) {
    end = sum(params$pReduct[1:id1])
    start = end - params$pReduct[id1] + 1
    idx1 = start:end
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id1], "stwsshare.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id1], "stwsshare.rdata"))
    read_time <- read_time + proc.time()[3]
    for (id2 in 1:params$numDataPartners) {
      end = sum(params$pReduct[1:id2])
      start = end - params$pReduct[id2] + 1
      idx2 = start:end
      if (id1 < id2) {
        StWS[idx1, idx2] = StWS[idx1, idx2] + C[[id2]]
        StWS[idx2, idx1] = StWS[idx2, idx1] + t(C[[id2]])
      } else if (id1 == id2) {
        StWS[idx1, idx1] = C[[id1]]
      } else {
        StWS[idx2, idx1] = StWS[idx2, idx1] + C[[id2]]
        StWS[idx1, idx2] = StWS[idx1, idx2] + t(C[[id2]])
      }
    }
    if (id1 < params$numDataPartners) {
      for (id2 in (id1 + 1):params$numDataPartners) {
        end = sum(params$pReduct[1:id2])
        start = end - params$pReduct[id2] + 1
        idx2 = start:end
        temp = t(params$halfshare[[id1]]) %*% MultiplyDiagonalWTimesX(W, params$halfshare[[id2]])
        StWS[idx1, idx2] = StWS[idx1, idx2] + temp
        StWS[idx2, idx1] = StWS[idx2, idx1] + t(temp)
      }
    }
  }

  I = NULL
  tryCatch({
    I = solve(StWS)
  },
  error = function(err) {
    I = NULL
  }
  )
  if (is.null(I)) {
    params$failed <- TRUE
    params$singular_matrix = TRUE
    params$error_message =
      paste0("The matrix t(X)*W*X is not invertible.\n",
             "       This may be due to one of two possible problems.\n",
             "       1. Poor random initialization of the security matrices.\n",
             "       2. Near multicollinearity in the data\n",
             "SOLUTIONS: \n",
             "       1. Rerun the data analysis.\n",
             "       2. If the problem persists, check the variables for\n",
             "          duplicates for both parties and / or reduce the\n",
             "          number of variables used. Once this is done,\n",
             "          rerun the data analysis.")
    params <- add_to_log(params, "ComputeStWSLogistic.AC", read_time, read_size, 0, 0)
    return(params)
  }
  params$I = I
  halfshare = params$halfshare[[1]]
  for (id in 2:params$numDataPartners) {
    halfshare = cbind(halfshare, params$halfshare[[id]])
  }
  IDt = I %*% (params$sty - t(halfshare) %*% params$pi_)
  Itemp = I
  IDttemp = IDt

  write_time <-0
  write_size = 0
  start = 1
  stop  = params$pReduct[1]
  for (id in 1:params$numDataPartners) {
    I = Itemp[start:stop, , drop = FALSE]
    IDt = IDttemp[start:stop, , drop = FALSE]
    write_time <-write_time - proc.time()[3]
    save(I, IDt, file = file.path(params$write_path, paste0("ID", id, ".rdata")))
    write_size = write_size + file.size(file.path(params$write_path, paste0("ID", id, ".rdata")))
    write_time <-write_time + proc.time()[3]
    start = stop + 1
    stop = stop + params$pReduct[id + 1]
  }

  params <- add_to_log(params, "ComputeStWSLogistic.AC", read_time, read_size, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm runif
UpdateBetaLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.DP\n\n")
  I = IDt = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, paste0("ID", params$data_partner_id, ".rdata")))
  read_size = file.size(file.path(params$readPathAC, paste0("ID", params$data_partner_id, ".rdata")))
  read_time <- proc.time()[3] - read_time

  id = 1
  set.seed(params$seeds[id], kind = "Mersenne-Twister")
  idx = params$indicies[[id]]
  halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
                       nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
  for (id in 2:params$numDataPartners) {
    set.seed(params$seeds[id], kind = "Mersenne-Twister")
    idx = params$indicies[[id]]
    halfshareDP = cbind(halfshareDP,
                        matrix(rnorm(params$n * params$ps[id], sd = 20),
                               nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE])
  }

  D0 = t(halfshareDP) %*% params$pi_
  deltaBeta = IDt - I %*% D0
  params$betas = params$betas + deltaBeta
  maxdifference = max(abs(deltaBeta) / (abs(params$betas) + .1))
  utemp = sum(runif(length(deltaBeta), min = 1, max = 5) * abs(params$betas))

  write_time <-proc.time()[3]
  save(utemp, maxdifference, file = file.path(params$write_path, "u_converge.rdata"))
  write_size = file.size(file.path(params$write_path, "u_converge.rdata"))
  write_time <-proc.time()[3] - write_time

  params <- add_to_log(params, "UpdateBetaLogistic.DP", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeConvergeStatusLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeConvergeStatusLogistic.AC\n\n")
  read_time <- 0
  read_size = 0
  u = 0
  converged = TRUE
  utemp = NULL
  maxdifference = NULL
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "u_converge.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "u_converge.rdata"))
    read_time <- read_time + proc.time()[3]
    u = u + utemp
    converged = converged && (maxdifference < params$cutoff)
  }
  maxIterExceeded = params$algIterationCounter >= params$max_iterations
  params$maxIterExceeded = maxIterExceeded
  params$u = u
  params$converged = converged
  write_time <-proc.time()[3]
  save(u, converged, maxIterExceeded, file = file.path(params$write_path, "u_converge.rdata"))
  write_size = file.size(file.path(params$write_path, "u_converge.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeConvergeStatusLogistic.AC", read_time, read_size, write_time, write_size)
  return(params)
}


GetConvergeStatusLogistic.DP = function(params) {
  converged = NULL
  if (params$trace) cat(as.character(Sys.time()), "GetconvergeStatusLogistic.DP\n\n")
  u = converge = maxIterExceeded = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "u_converge.rdata"))
  read_size = file.size(file.path(params$readPathAC, "u_converge.rdata"))
  read_time <- proc.time()[3] - read_time
  params$u = u
  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  params <- add_to_log(params, "GetConvergeStatusLogistic.DP", read_time, read_size, 0, 0)
  return(params)
}

SendFinalBetasLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "SendFinalBetasLogistic.DP\n\n")
  betas = params$betas
  write_time <-proc.time()[3]
  save(betas, file = file.path(params$write_path, "finalbetas.rdata"))
  write_size = file.size(file.path(params$write_path, "finalbetas.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "SendFinalBetasLogistic.DP", 0, 0, write_time, write_size)
  return(params)
}

ComputeFinalSBetaLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeFinalSBetaLogistic.AC\n\n")
  Sbeta = 0
  read_time <- 0
  read_size = 0
  sbeta = 0
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_time <- read_time + proc.time()[3]
    sbeta = sbeta + Sbeta
  }
  sbeta = 2 * params$u * sbeta - params$numDataPartners * params$u

  write_time <-proc.time()[3]
  save(sbeta, file = file.path(params$write_path, "sbeta.rdata"))
  write_size = file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <-proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeFinalSBetaLogistic.AC", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeResultsLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.DP\n\n")
  sbeta = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "sbeta.rdata"))
  read_size = file.size(file.path(params$readPathAC, "sbeta.rdata"))
  read_time <- proc.time()[3] - read_time

  n       = params$n
  ct      = sum(data$Y)
  params$FinalFitted = sbeta
  resdev  = -2 * (sum(data$Y * sbeta) - sum(log(1 + exp(sbeta))))
  nulldev = -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

  hoslem  = HoslemInternal(params, data)
  ROC     = RocInternal(params, data)

  write_time <-proc.time()[3]
  save(resdev, nulldev, hoslem, ROC, file = file.path(params$write_path, "logisticstats.rdata"))
  write_size = file.size(file.path(params$write_path, "logisticstats.rdata"))
  write_time <-proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeResultsLogistic.DP", read_time, read_size, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
ComputeResultsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.AC\n\n")
  nulldev = NULL
  resdev  = NULL
  hoslem  = NULL
  ROC     = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "logisticstats.rdata"))
  read_size = file.size(file.path(params$readPathDP[1], "logisticstats.rdata"))
  read_time <- proc.time()[3] - read_time
  coefficients = c()
  p            = 0
  betas = NULL
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "finalbetas.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "finalbetas.rdata"))
    read_time <- read_time + proc.time()[3]
    coefficients = c(coefficients, betas)
    p            = p + length(params$indicies[[id]])
  }

  coefficients[2:p] = coefficients[2:p] / params$colran[2:p]
  coefficients[1] = coefficients[1] - sum(coefficients[2:p] * params$colmin[2:p])

  serror = rep(0, p)
  serror[2:p] = sqrt(diag(params$I)[2:p]) / params$colran[2:p]
  d1 = diag(c(1, params$colmin[-1] / params$colran[-1]))
  temp = d1 %*% params$I %*% d1
  serror[1] = sqrt(temp[1, 1] - 2 * sum(temp[1, 2:p]) + sum(temp[2:p, 2:p]))

  stats = params$stats
  stats$failed         = FALSE
  stats$converged      = params$converged


  # If xtwx were singular, it would have been caught in GetII.A2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  stats$party = params$party
  stats$coefficients = rep(NA, params$p)
  stats$secoef = rep(NA, params$p)
  stats$tvals  = rep(NA, params$p)
  stats$pvals  = rep(NA, params$p)
  stats$n  = params$n
  stats$nulldev = nulldev
  stats$resdev = resdev
  stats$aic = resdev + 2 * sum(params$pReduct)
  stats$bic = resdev + sum(params$pReduct) * log(params$n)
  stats$nulldev_df = params$n - 1
  stats$resdev_df = params$n - sum(params$pReduct)
  stats$coefficients[params$fullindicies] = coefficients
  stats$secoef[params$fullindicies] = serror
  tvals = coefficients / serror
  pvals = 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[params$fullindicies] = tvals
  stats$pvals[params$fullindicies] = pvals
  stats$hoslem  = hoslem
  stats$ROC     = ROC
  stats$iter    = params$algIterationCounter - 1
  names(stats$coefficients) = params$colnames
  names(stats$party) = params$colnames
  names(stats$secoef) = params$colnames
  names(stats$tvals) = params$colnames
  names(stats$pvals) = params$colnames

  write_time <-proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size = file.size(file.path(params$write_path, "stats.rdata"))
  write_time <-proc.time()[3] - write_time

  params$stats      = stats

  params <- add_to_log(params, "ComputeResultsLogistic.AC", read_time, read_size, write_time, write_size)
  return(params)
}


GetResultsLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.DP\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "stats.rdata"))
  read_size = file.size(file.path(params$readPathAC, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  if (params$data_partner_id == 1) {
    stats$Y           = data$Y # For Hoslem and ROC
    stats$FinalFitted = params$FinalFitted
  }
  params$stats      = stats
  params <- add_to_log(params, "GetResultsLogistic.DP", read_time, read_size, 0, 0)
  return(params)
}

############################## PARENT FUNCTIONS ###############################

DataPartnerKLogistic = function(data,
                                yname           = NULL,
                                numDataPartners = NULL,
                                data_partner_id   = NULL,
                                monitor_folder   = NULL,
                                sleep_time       = 10,
                                maxWaitingTime  = 24 * 60 * 60,
                                popmednet       = TRUE,
                                trace           = FALSE,
                                verbose         = TRUE) {

  params <- PrepareParams.kp("logistic", data_partner_id, numDataPartners, ac = FALSE,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- InitializeLog.kp(params)
  params <- InitializeStamps.kp(params)
  params <- InitializeTrackingTable.kp(params)
  Header(params)

  params   = PrepareFolder.ACDP(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  if (data_partner_id == 1) {
    data = PrepareDataLinLog.DP1(params, data, yname)
    params <- add_to_log(params, "PrepareParamsLinLog.DP1", 0, 0, 0, 0)
  } else {
    data = PrepareDataLinLog.DPk(params, data)
    params <- add_to_log(params, "PrepareParamsLinLog.DPk", 0, 0, 0, 0)
  }

  params <- add_to_log(params, "PrepareParamsLinear.DP", 0, 0, 0, 0)

  if (data$failed) {
    params$error_message = paste("Error processing data for data partner", params$data_partner_id, "\n")
    MakeErrorMessage(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
    params$error_message = ReadErrorMessage(params$readPathAC)
    warning(params$error_message)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- SendBasicInfo.DP(params, data)
  files = "n_analysis.rdata"
  params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  possibleError = ReceivedError.kp(params, from = "AC")
  if (possibleError$error) {
    params$error_message = possibleError$message
    warning(possibleError$message)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- PrepareParamsLinear.DP(params, data)
  files = "p_scaler_seed.rdata"
  params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  params <- PrepareSharesLinear.DP(params, data)
  files = c("products.rdata", "halfshare.rdata", "colstats.rdata")
  params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  possibleError = ReceivedError.kp(params, from = "AC")
  if (possibleError$error) {
    params$error_message = possibleError$message
    warning(possibleError$message)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- UpdateParamsLogistic.DP(params)

  data = UpdateDataLogistic.DP(params, data)
  params <- add_to_log(params, "UpdateDataLogistic.DP", 0, 0, 0, 0)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)

    params <- ComputeSbetaLogistic.DP(params, data)
    files = "Sbeta.rdata"
    params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    params <- ComputeStWSLogistic.DP(params, data)
    files = "stwsshare.rdata"
    params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

    possibleError = ReceivedError.kp(params, from = "AC")
    if (possibleError$error) {
      params$error_message = possibleError$message
      warning(possibleError$message)
      params <- SendPauseQuit.kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
      return(params$stats)
    }

    params <- UpdateBetaLogistic.DP(params)
    files = "u_converge.rdata"
    params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)


    params <- GetConvergeStatusLogistic.DP(params)

    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- ComputeSbetaLogistic.DP(params, data)

  params <- SendFinalBetasLogistic.DP(params)

  files = c("sbeta.rdata", "finalbetas.rdata")
  params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  if (data_partner_id == 1) {
    params <- ComputeResultsLogistic.DP(params, data)
    files = "logisticstats.rdata"
    params <- SendPauseContinue.kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
  }

  params <- GetResultsLogistic.DP(params, data)
  params <- SendPauseQuit.kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


AnalysisCenterKLogistic = function(numDataPartners = NULL,
                                   monitor_folder   = NULL,
                                   msreqid         = "v_default_0_000",
                                   cutoff          = 1E-8,
                                   max_iterations   = 25,
                                   sleep_time       = 10,
                                   maxWaitingTime  = 24 * 60 * 60,
                                   popmednet       = TRUE,
                                   trace           = FALSE,
                                   verbose         = TRUE) {

  filesList = rep(list(list()), numDataPartners)

  params <- PrepareParams.kp("logistic", 0, numDataPartners, msreqid, cutoff, max_iterations, ac = TRUE,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- InitializeLog.kp(params)
  params <- InitializeStamps.kp(params)
  params <- InitializeTrackingTable.kp(params)
  Header(params)

  params   = PrepareFolder.ACDP(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.kp(params, from = "DP", maxWaitingTime = maxWaitingTime)

  possibleError = ReceivedError.kp(params, from = "DP")
  if (possibleError$error) {
    params$error_message = possibleError$message
    warning(possibleError$message)
    MakeErrorMessage(params$write_path, possibleError$message)
    files = "error_message.rdata"
    params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- CheckAgreement.AC(params)

  if (params$failed) {
    MakeErrorMessage(params$write_path, params$error_message)
    files = "error_message.rdata"
    warning(params$error_message)
    params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  files = "empty.rdata"
  params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- GetProductsLogistic.AC(params)


  params <- CheckColinearityLogistic.AC(params)

  if (params$failed) {
    MakeErrorMessage(params$write_path, params$error_message)
    files = "error_message.rdata"
    warning(params$error_message)
    params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    params <- SendPauseQuit.kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- ComputeInitialBetasLogistic.AC(params)

  for (id in 1:numDataPartners) {
    filesList[[id]] = c(paste0("u_beta_", id, ".rdata"), "indicies.rdata")
  }

  params <- SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params$algIterationCounter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- ComputeWeightsLogistic.AC(params)
    files = "pi.rdata"
    params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    params <- ComputeStWSLogistic.AC(params)

    if (params$failed) {
      MakeErrorMessage(params$write_path, params$error_message)
      files = "error_message.rdata"
      warning(params$error_message)
      params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                    sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
      params <- SendPauseQuit.kp(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.kp(params)
      return(params$stats)
    }

    for (id in 1:numDataPartners) {
      filesList[[id]] = paste0("id", id, ".rdata")
    }
    params <- SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

    params <- ComputeConvergeStatusLogistic.AC(params)
    files = "u_converge.rdata"
    params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- ComputeFinalSBetaLogistic.AC(params)
  filesList = rep(list(list()), numDataPartners)
  filesList[[1]] = "sbeta.rdata"
  params <- SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)

  params <- ComputeResultsLogistic.AC(params)
  files = "stats.rdata"
  params <- SendPauseContinue.kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, maxWaitingTime = maxWaitingTime)
  params <- SendPauseQuit.kp(params, sleep_time = sleep_time)
  SummarizeLog.kp(params)
  return(params$stats)
}
