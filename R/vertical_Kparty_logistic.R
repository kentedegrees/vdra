################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

get_products_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "get_products_logistic_ac\n\n")
  read_time <- 0
  read_size <- 0
  p         <- 0
  n         <- 0
  pi        <- c()

  allproducts  <- rep(list(list()), params$num_data_partners)
  allhalfshare <- rep(list(list()), params$num_data_partners)
  alltags      <- rep(list(list()), params$num_data_partners)
  products  <- NULL
  halfshare <- NULL
  tags      <- NULL
  allcolmin <- allcolrange <- allcolsum <- allcolnames <- NULL
  colmin <- colrange <- colsum <- colnames <- NULL
  party <- NULL
  for (id in 1:params$num_data_partners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "products.rdata"))
    load(file.path(params$readPathDP[id], "halfshare.rdata"))
    load(file.path(params$readPathDP[id], "colstats.rdata"))
    read_size <- read_size + sum(file.size(file.path(params$readPathDP[id],
                                                     c("products.rdata",
                                                       "halfshare.rdata",
                                                       "colstats.rdata"))))
    read_time <- read_time + proc.time()[3]

    allproducts[[id]]  <- products
    allhalfshare[[id]] <- halfshare
    alltags[[id]]      <- tags
    allcolmin          <- c(allcolmin, colmin)
    allcolrange        <- c(allcolrange, colrange)
    allcolsum          <- c(allcolsum, colsum)
    allcolnames        <- c(allcolnames, colnames)
    party              <- c(party, rep(paste0("dp", id), length(colnames)))
    p <- p + ncol(halfshare)
    pi <- c(pi, ncol(halfshare))
    if (id == 1) n <- nrow(halfshare)
  }

  m <- matrix(0, p, p)
  colnames(m) <- allcolnames
  rownames(m) <- allcolnames
  offset1 <- 1
  params$pi <- rep(0, params$num_data_partners)
  for (id1 in 1:params$num_data_partners) {
    p1 <- ncol(allhalfshare[[id1]])
    params$pi[id1] <- p1
    offset2 <- offset1
    for (id2 in id1:params$num_data_partners) {
      p2 <- ncol(allhalfshare[[id2]])
      if (id1 == id2) {
        m[offset1:(offset1 + p1 - 1),
          offset2:(offset2 + p2 - 1)] <- allproducts[[id1]][[id2]]
      } else {
        temp <- allproducts[[id1]][[id2]] + allproducts[[id2]][[id1]] +
          t(allhalfshare[[id1]]) %*% allhalfshare[[id2]]
        m[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] <- temp
        m[offset2:(offset2 + p2 - 1), offset1:(offset1 + p1 - 1)] <- t(temp)
      }
      offset2 <- offset2 + p2
    }
    offset1 <- offset1 + p1
  }


  params$halfshare    <- allhalfshare
  params$sts          <- m[2:p, 2:p, drop = FALSE]
  params$sty          <- m[2:p, 1, drop = FALSE]
  params$yty          <- m[1, 1]
  params$means_y       <- allcolsum[1] / n
  params$means        <- allcolsum[-1] / n
  params$n            <- n
  params$p            <- p
  params$pi           <- pi
  params$colmin       <- allcolmin[-1]
  params$colrange     <- allcolrange[-1]
  params$colsum       <- allcolsum[-1]
  params$colnames     <- allcolnames[-1]
  params$party        <- party[-1]
  params$tags         <- alltags

  params <- add_to_log(params, "get_products_logistic_ac",
                       read_time, read_size, 0, 0)
  return(params)
}


check_colinearity_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "check_colinearity_logistic_ac\n\n")
  sts <- params$sts
  sty <- params$sty

  nrow <- nrow(sts)
  indicies <- c(1)
  for (i in 2:nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(sts[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }


  sts <- sts[indicies, indicies, drop = FALSE]
  sty <- sty[indicies, drop = FALSE]

  params$sts <- sts
  params$sty <- sty

  # Extract the indicies to keep for each party and check for errors.

  params$colmin       <- params$colmin[indicies]
  params$colrange     <- params$colrange[indicies]
  params$colsum       <- params$colsum[indicies]
  params$fullindicies <- indicies       # To be used when computing stats
  params$p            <- params$p - 1   # Get rid of the response from the count
  # take into account that pi still counts sty, which we removed earlier.
  indicies <- indicies + 1

  params$indicies   <- rep(list(list()), params$num_data_partners)
  tags              <- rep(list(list()), params$num_data_partners)
  min <- 1

  for (id in 1:params$num_data_partners) {
    max <- min + params$pi[id] - 1
    idx <- indicies[which(min <= indicies & indicies <= max)] - min + 1
    params$indicies[[id]] <- idx
    if (id == 1) {
      idx <- (idx - 1)[-1]
    }
    temp <- params$tags[[id]]
    temp <- temp[idx]
    tags[[id]] <- temp
    min <- max + 1
  }

  params$error_message <- ""
  if ((length(unique(tags[[1]])) == 1) ||
      (length(unique(tags[[1]])) >= 2 && !("numeric" %in% names(tags[[1]])))) {
    params$failed <- TRUE
    params$error_message <-
      paste("Data Partner 1 must have no covariates or at least 2",
            "covariates at least one of which is continuous.\n")
  }
  for (id in 2:params$num_data_partners) {
    if (length(unique(tags[[id]])) < 2) {
      params$failed <- TRUE
      params$error_message <-
        paste0(params$error_message,
               paste("After removing colinear covariates, Data Partner", id,
                     "has 1 or fewer covariates.\n"))
    } else if (!("numeric" %in% names(tags[[id]]))) {
      params$failed <- TRUE
      params$error_message <-
        paste0(params$error_message,
               paste("After removing colinear covariates, Data Partner", id,
                     "has no continuous covariates.\n"))
    }
  }

  indicies <- params$indicies

  params$p_reduct <- c()
  for (id in 1:params$num_data_partners) {
    params$p_reduct <- c(params$p_reduct, length(indicies[[id]]))
  }

  for (id in 1:params$num_data_partners) {
    params$halfshare[[id]] <-
      params$halfshare[[id]][, indicies[[id]], drop = FALSE]
  }

  write_time <- proc.time()[3]
  save(indicies, file = file.path(params$write_path, "indicies.rdata"))
  write_size <- file.size(file.path(params$write_path, "indicies.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "check_colinearity_logistic_ac",
                       0, 0, write_time, write_size)

  return(params)
}


#' @importFrom stats runif
comp_init_betas_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_init_betas_logistic_ac\n\n")
  write_time <- 0
  write_size <- 0
  colsum_s <- (params$colsum - params$n * params$colmin) / params$colran
  beta <- 4 * solve(params$sts) %*% (params$sty - 0.5 * colsum_s)

  u <- sum(runif(length(beta), min = 1, max = 5) * abs(beta))
  params$u <- u
  start <- 1
  for (id in 1:params$num_data_partners) {
    end <- start + length(params$indicies[[id]]) - 1
    betas <- beta[start:end]

    write_time <- write_time - proc.time()[3]
    save(u, betas, file = file.path(params$write_path,
                                    paste0("u_beta_", id, ".rdata")))
    write_size <- write_size +
      file.size(file.path(params$write_path, paste0("u_beta_", id, ".rdata")))
    write_time <- write_time + proc.time()[3]
    start <- end + 1
  }
  params <- add_to_log(params, "comp_init_betas_logistic_ac",
                       0, 0, write_time, write_size)
  return(params)
}


update_params_logistic_dp <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_params_logistic_dp\n\n")
  indicies <- NULL
  u <- NULL
  betas <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "indicies.rdata"))
  filename <- paste0("u_beta_", params$data_partner_id, ".rdata")
  load(file.path(params$readPathAC, filename))
  read_size <- file.size(file.path(params$readPathAC, "indicies.rdata")) +
    file.size(file.path(params$readPathAC, filename))
  read_time <- proc.time()[3] - read_time
  params$u <- u
  params$betas <- betas
  params$indicies <- indicies
  params <- add_to_log(params, "update_params_logistic_dp",
                       read_time, read_size, 0, 0)
  return(params)
}


update_data_logistic_dp <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_data_logistic_dp\n\n")
  if (params$data_partner_id == 1) {
    data$Y <- data$x[, 1, drop = FALSE]
  }
  idx <- params$indicies[[params$data_partner_id]]
  data$x <- data$x[, idx, drop = FALSE]
  data$colmin <- data$colmin[idx]
  data$colmax <- data$colmax[idx]
  data$colsum <- data$colsum[idx]
  data$colrange <- data$colrange[idx]
  return(data)
}


#' @importFrom stats rnorm runif
comp_s_beta_logistic_dp <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_s_beta_logistic_dp\n\n")
  set.seed(params$seed +
             params$alg_iteration_counter, kind = "Mersenne-Twister")
  v <- matrix(rnorm(params$n,
                    mean = runif(n = 1, min = -1, max = 1), sd = 10), ncol = 1)
  v_sum <- 0
  for (id in 1:params$num_data_partners) {
    set.seed(params$seeds[id] +
               params$alg_iteration_counter, kind = "Mersenne-Twister")
    v_sum <- v_sum + matrix(rnorm(params$n,
                                mean = runif(n = 1, min = -1, max = 1),
                                sd = 10), ncol = 1)
  }

  s_beta <- (data$x %*% params$betas + params$u) /
    (2 * params$u) + v - params$scaler / sum(params$scalers) * v_sum

  write_time <- proc.time()[3]
  save(s_beta, file = file.path(params$write_path, "sbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "comp_s_beta_logistic_dp",
                       0, 0, write_time, write_size)
  return(params)
}


comp_weights_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_weights_logistic_ac\n\n")
  s_beta <- 0
  read_time <- 0
  read_size <- 0
  sbeta <- 0
  for (id in 1:params$num_data_partners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_size <- read_size +
      file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_time <- read_time + proc.time()[3]
    sbeta <- sbeta + s_beta
  }
  sbeta <- 2 * params$u * sbeta - params$num_data_partners * params$u
  pi_ <- 1 / (1 + exp(-sbeta))
  params$pi_ <- pi_

  write_time <- proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi.rdata"))
  write_size <- file.size(file.path(params$write_path, "pi.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComptueWeightsLogistic.AC",
                       read_time, read_size, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm
comp_stws_logistic_dp <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_stws_logistic_dp\n\n")
  pi_ <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "pi.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "pi.rdata"))
  read_time <- proc.time()[3] - read_time
  params$pi_ <- pi_

  w <- pi_ * (1 - pi_)
  c_mat <- rep(list(list()), params$num_data_partners)

  idx <- params$indicies[[params$data_partner_id]]
  set.seed(params$seed, kind = "Mersenne-Twister")
  halfshare <- matrix(rnorm(params$n * params$p, sd = 20),
                      nrow = params$n, ncol = params$p)[, idx, drop = FALSE]

  for (id in 1:params$num_data_partners) {
    if (id < params$data_partner_id) {
      set.seed(params$seeds[id], kind = "Mersenne-Twister")
      idx <- params$indicies[[id]]
      halfshare_dp <- matrix(rnorm(params$n * params$ps[id], sd = 20),
                            nrow = params$n,
                            ncol = params$ps[id])[, idx, drop = FALSE]
      c_mat[[id]] <- params$scaler / (params$scaler + params$scalers[id]) *
        t(halfshare_dp) %*% MultiplyDiagonalWTimesX(w, halfshare) +
        t(halfshare_dp) %*% MultiplyDiagonalWTimesX(w, data$x - halfshare)
    } else if (id == params$data_partner_id) {
      c_mat[[id]] <- t(data$x) %*% MultiplyDiagonalWTimesX(w, data$x)
    } else {
      set.seed(params$seeds[id], kind = "Mersenne-Twister")
      idx <- params$indicies[[id]]
      halfshare_dp <- matrix(rnorm(params$n * params$ps[id], sd = 20),
                            nrow = params$n,
                            ncol = params$ps[id])[, idx, drop = FALSE]
      c_mat[[id]] <- params$scaler / (params$scaler + params$scalers[id]) *
        t(halfshare) %*% MultiplyDiagonalWTimesX(w, halfshare_dp) +
        t(data$x - halfshare) %*% MultiplyDiagonalWTimesX(w, halfshare_dp)
    }
  }

  write_time <- proc.time()[3]
  save(c_mat, file = file.path(params$write_path, "stwsshare.rdata"))
  write_size <- file.size(file.path(params$write_path, "stwsshare.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "comp_stws_logistic_dp",
                       read_time, read_size, write_time, write_size)
  return(params)
}


comp_stws_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_stws_logistic_ac\n\n")
  read_time <- 0
  read_size <- 0
  c_mat        <- NULL
  w <- params$pi_ * (1 - params$pi_)
  st_ws <- matrix(0, sum(params$p_reduct), sum(params$p_reduct))

  for (id1 in 1:params$num_data_partners) {
    end <- sum(params$p_reduct[1:id1])
    start <- end - params$p_reduct[id1] + 1
    idx1 <- start:end
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id1], "stwsshare.rdata"))
    read_size <- read_size +
      file.size(file.path(params$readPathDP[id1], "stwsshare.rdata"))
    read_time <- read_time + proc.time()[3]
    for (id2 in 1:params$num_data_partners) {
      end <- sum(params$p_reduct[1:id2])
      start <- end - params$p_reduct[id2] + 1
      idx2 <- start:end
      if (id1 < id2) {
        st_ws[idx1, idx2] <- st_ws[idx1, idx2] + c_mat[[id2]]
        st_ws[idx2, idx1] <- st_ws[idx2, idx1] + t(c_mat[[id2]])
      } else if (id1 == id2) {
        st_ws[idx1, idx1] <- c_mat[[id1]]
      } else {
        st_ws[idx2, idx1] <- st_ws[idx2, idx1] + c_mat[[id2]]
        st_ws[idx1, idx2] <- st_ws[idx1, idx2] + t(c_mat[[id2]])
      }
    }
    if (id1 < params$num_data_partners) {
      for (id2 in (id1 + 1):params$num_data_partners) {
        end <- sum(params$p_reduct[1:id2])
        start <- end - params$p_reduct[id2] + 1
        idx2 <- start:end
        temp <- t(params$halfshare[[id1]]) %*%
          MultiplyDiagonalWTimesX(w, params$halfshare[[id2]])
        st_ws[idx1, idx2] <- st_ws[idx1, idx2] + temp
        st_ws[idx2, idx1] <- st_ws[idx2, idx1] + t(temp)
      }
    }
  }

  i_mat <- NULL
  tryCatch({
    i_mat <- solve(st_ws)
  },
  error = function(err) {
    i_mat <- NULL
  }
  )
  if (is.null(i_mat)) {
    params$failed <- TRUE
    params$singular_matrix <- TRUE
    params$error_message <-
      paste0("The matrix t(x)*w*x is not invertible.\n",
             "       This may be due to one of two possible problems.\n",
             "       1. Poor random initialization of the security matrices.\n",
             "       2. Near multicollinearity in the data\n",
             "SOLUTIONS: \n",
             "       1. Rerun the data analysis.\n",
             "       2. If the problem persists, check the variables for\n",
             "          duplicates for both parties and / or reduce the\n",
             "          number of variables used. Once this is done,\n",
             "          rerun the data analysis.")
    params <- add_to_log(params, "comp_stws_logistic_ac",
                         read_time, read_size, 0, 0)
    return(params)
  }
  params$i_mat <- i_mat
  halfshare <- params$halfshare[[1]]
  for (id in 2:params$num_data_partners) {
    halfshare <- cbind(halfshare, params$halfshare[[id]])
  }
  IDt <- i_mat %*% (params$sty - t(halfshare) %*% params$pi_)
  Itemp <- i_mat
  IDttemp <- IDt

  write_time <- 0
  write_size <- 0
  start <- 1
  stop  <- params$p_reduct[1]
  for (id in 1:params$num_data_partners) {
    i_mat <- Itemp[start:stop, , drop = FALSE]
    IDt <- IDttemp[start:stop, , drop = FALSE]
    write_time <- write_time - proc.time()[3]
    save(i_mat, IDt,
         file = file.path(params$write_path, paste0("ID", id, ".rdata")))
    write_size <- write_size +
      file.size(file.path(params$write_path, paste0("ID", id, ".rdata")))
    write_time <- write_time + proc.time()[3]
    start <- stop + 1
    stop <- stop + params$p_reduct[id + 1]
  }

  params <- add_to_log(params, "comp_stws_logistic_ac",
                       read_time, read_size, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm runif
update_beta_logistic_dp <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_beta_logistic_dp\n\n")
  i_mat <- IDt <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC,
                 paste0("ID", params$data_partner_id, ".rdata")))
  read_size <-
    file.size(file.path(params$readPathAC,
                        paste0("ID", params$data_partner_id, ".rdata")))
  read_time <- proc.time()[3] - read_time

  id <- 1
  set.seed(params$seeds[id], kind = "Mersenne-Twister")
  idx <- params$indicies[[id]]
  halfshare_dp <- matrix(rnorm(params$n * params$ps[id], sd = 20),
                        nrow = params$n,
                        ncol = params$ps[id])[, idx, drop = FALSE]
  for (id in 2:params$num_data_partners) {
    set.seed(params$seeds[id], kind = "Mersenne-Twister")
    idx <- params$indicies[[id]]
    halfshare_dp <- cbind(halfshare_dp,
                         matrix(rnorm(params$n * params$ps[id], sd = 20),
                                nrow = params$n,
                                ncol = params$ps[id])[, idx, drop = FALSE])
  }

  D0 <- t(halfshare_dp) %*% params$pi_
  delta_beta <- IDt - i_mat %*% D0
  params$betas <- params$betas + delta_beta
  maxdifference <- max(abs(delta_beta) / (abs(params$betas) + .1))
  u_temp <- sum(runif(length(delta_beta),
                      min = 1, max = 5) * abs(params$betas))

  write_time <- proc.time()[3]
  save(u_temp, maxdifference,
       file = file.path(params$write_path, "u_converge.rdata"))
  write_size <- file.size(file.path(params$write_path, "u_converge.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "update_beta_logistic_dp",
                       read_time, read_size, write_time, write_size)
  return(params)
}


comp_conv_status_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_conv_status_logistic_ac\n\n")
  read_time <- 0
  read_size <- 0
  u <- 0
  converged <- TRUE
  u_temp <- NULL
  maxdifference <- NULL
  for (id in 1:params$num_data_partners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "u_converge.rdata"))
    read_size <- read_size +
      file.size(file.path(params$readPathDP[id], "u_converge.rdata"))
    read_time <- read_time + proc.time()[3]
    u <- u + u_temp
    converged <- converged && (maxdifference < params$cutoff)
  }
  maxIterExceeded <- params$alg_iteration_counter >= params$max_iterations
  params$maxIterExceeded <- maxIterExceeded
  params$u <- u
  params$converged <- converged
  write_time <- proc.time()[3]
  save(u, converged, maxIterExceeded,
       file = file.path(params$write_path, "u_converge.rdata"))
  write_size <- file.size(file.path(params$write_path, "u_converge.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "comp_conv_status_logistic_ac",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_conv_status_logistic_dp <- function(params) {
  converged <- NULL
  if (params$trace) cat(as.character(Sys.time()),
                        "get_conv_status_logistic_dp\n\n")
  u <- maxIterExceeded <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "u_converge.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "u_converge.rdata"))
  read_time <- proc.time()[3] - read_time
  params$u <- u
  params$converged <- converged
  params$maxIterExceeded <- maxIterExceeded
  params <- add_to_log(params, "get_conv_status_logistic_dp",
                       read_time, read_size, 0, 0)
  return(params)
}

SendFinalBetasLogistic.DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "SendFinalBetasLogistic.DP\n\n")
  betas <- params$betas
  write_time <- proc.time()[3]
  save(betas, file = file.path(params$write_path, "finalbetas.rdata"))
  write_size <- file.size(file.path(params$write_path, "finalbetas.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "SendFinalBetasLogistic.DP",
                       0, 0, write_time, write_size)
  return(params)
}

ComputeFinalSBetaLogistic.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "ComputeFinalSBetaLogistic.AC\n\n")
  s_beta <- 0
  read_time <- 0
  read_size <- 0
  sbeta <- 0
  for (id in 1:params$num_data_partners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_size <- read_size +
      file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_time <- read_time + proc.time()[3]
    sbeta <- sbeta + s_beta
  }
  sbeta <- 2 * params$u * sbeta - params$num_data_partners * params$u

  write_time <- proc.time()[3]
  save(sbeta, file = file.path(params$write_path, "sbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeFinalSBetaLogistic.AC",
                       read_time, read_size, write_time, write_size)
  return(params)
}


comp_results_logistic_dp <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_results_logistic_dp\n\n")
  sbeta <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "sbeta.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "sbeta.rdata"))
  read_time <- proc.time()[3] - read_time

  n       <- params$n
  ct      <- sum(data$Y)
  params$final_fitted <- sbeta
  resdev  <- -2 * (sum(data$Y * sbeta) - sum(log(1 + exp(sbeta))))
  nulldev <- -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

  hoslem  <- HoslemInternal(params, data)
  ROC     <- roc_internal(params, data)

  write_time <- proc.time()[3]
  save(resdev, nulldev, hoslem, ROC,
       file = file.path(params$write_path, "logisticstats.rdata"))
  write_size <- file.size(file.path(params$write_path, "logisticstats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "comp_results_logistic_dp",
                       read_time, read_size, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
comp_results_logistic_ac <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_results_logistic_ac\n\n")
  nulldev <- NULL
  resdev  <- NULL
  hoslem  <- NULL
  ROC     <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "logisticstats.rdata"))
  read_size <- file.size(file.path(params$readPathDP[1], "logisticstats.rdata"))
  read_time <- proc.time()[3] - read_time
  coefficients <- c()
  p            <- 0
  betas <- NULL
  for (id in 1:params$num_data_partners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "finalbetas.rdata"))
    read_size <- read_size +
      file.size(file.path(params$readPathDP[id], "finalbetas.rdata"))
    read_time <- read_time + proc.time()[3]
    coefficients <- c(coefficients, betas)
    p            <- p + length(params$indicies[[id]])
  }

  coefficients[2:p] <- coefficients[2:p] / params$colran[2:p]
  coefficients[1] <- coefficients[1] -
    sum(coefficients[2:p] * params$colmin[2:p])

  serror <- rep(0, p)
  serror[2:p] <- sqrt(diag(params$i_mat)[2:p]) / params$colran[2:p]
  d1 <- diag(c(1, params$colmin[-1] / params$colran[-1]))
  temp <- d1 %*% params$i_mat %*% d1
  serror[1] <- sqrt(temp[1, 1] - 2 * sum(temp[1, 2:p]) + sum(temp[2:p, 2:p]))

  stats <- params$stats
  stats$failed         <- FALSE
  stats$converged      <- params$converged


  # If xtwx were singular, it would have been caught in GetII.A2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  stats$party <- params$party
  stats$coefficients <- rep(NA, params$p)
  stats$secoef <- rep(NA, params$p)
  stats$tvals  <- rep(NA, params$p)
  stats$pvals  <- rep(NA, params$p)
  stats$n <- params$n
  stats$nulldev <- nulldev
  stats$resdev <- resdev
  stats$aic <- resdev + 2 * sum(params$p_reduct)
  stats$bic <- resdev + sum(params$p_reduct) * log(params$n)
  stats$nulldev_df <- params$n - 1
  stats$resdev_df <- params$n - sum(params$p_reduct)
  stats$coefficients[params$fullindicies] <- coefficients
  stats$secoef[params$fullindicies] <- serror
  tvals <- coefficients / serror
  pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[params$fullindicies] <- tvals
  stats$pvals[params$fullindicies] <- pvals
  stats$hoslem  <- hoslem
  stats$ROC     <- ROC
  stats$iter    <- params$alg_iteration_counter - 1
  names(stats$coefficients) <- params$colnames
  names(stats$party) <- params$colnames
  names(stats$secoef) <- params$colnames
  names(stats$tvals) <- params$colnames
  names(stats$pvals) <- params$colnames

  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params$stats      <- stats

  params <- add_to_log(params, "comp_results_logistic_ac",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_results_logistic_dp <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_logistic_dp\n\n")
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "stats.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  if (params$data_partner_id == 1) {
    stats$Y           <- data$Y # For Hoslem and ROC
    stats$final_fitted <- params$final_fitted
  }
  params$stats      <- stats
  params <- add_to_log(params, "get_results_logistic_dp",
                       read_time, read_size, 0, 0)
  return(params)
}

############################## PARENT FUNCTIONS ###############################

data_partner_k_logistic <- function(data,
                                 y_name           = NULL,
                                 num_data_partners = NULL,
                                 data_partner_id   = NULL,
                                 monitor_folder   = NULL,
                                 sleep_time       = 10,
                                 max_waiting_time  = 24 * 60 * 60,
                                 popmednet       = TRUE,
                                 trace           = FALSE,
                                 verbose         = TRUE) {

  params <- prepare_params_kp("logistic",
                              data_partner_id,
                              num_data_partners,
                              ac = FALSE,
                              popmednet = popmednet,
                              trace = trace,
                              verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- initialize_log_kp(params)
  params <- initialize_time_stamps_kp(params)
  params <- initialize_tracking_table_kp(params)
  header(params)

  params <- prepare_folder_acdp(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  if (data_partner_id == 1) {
    data <- prepare_data_linlog_dp1(params, data, y_name)
    params <- add_to_log(params, "PrepareParamsLinLog.DP1", 0, 0, 0, 0)
  } else {
    data <- prepare_data_linlog_dpk(params, data)
    params <- add_to_log(params, "PrepareParamsLinLog.DPk", 0, 0, 0, 0)
  }

  params <- add_to_log(params, "prepare_params_linear_dp", 0, 0, 0, 0)

  if (data$failed) {
    params$error_message <- paste("Error processing data for data partner",
                                  params$data_partner_id, "\n")
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_kp(params,
                                     filesAC = files,
                                     from = "AC",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)
    params$error_message <- read_error_message(params$readPathAC)
    warning(params$error_message)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- send_basic_info_dp(params, data)
  files <- "n_analysis.rdata"
  params <- send_pause_continue_kp(params,
                                   filesAC = files,
                                   from = "AC",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  possible_error <- ReceivedError.kp(params, from = "AC")
  if (possible_error$error) {
    params$error_message <- possible_error$message
    warning(possible_error$message)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_dp(params, data)
  files <- "p_scaler_seed.rdata"
  params <- send_pause_continue_kp(params,
                                   filesDP = files,
                                   from = "DP",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  params <- prepare_shares_linear_dp(params, data)
  files <- c("products.rdata", "halfshare.rdata", "colstats.rdata")
  params <- send_pause_continue_kp(params,
                                   filesAC = files,
                                   from = "AC",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  possible_error <- ReceivedError.kp(params, from = "AC")
  if (possible_error$error) {
    params$error_message <- possible_error$message
    warning(possible_error$message)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- update_params_logistic_dp(params)

  data <- update_data_logistic_dp(params, data)
  params <- add_to_log(params, "update_data_logistic_dp", 0, 0, 0, 0)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)

    params <- comp_s_beta_logistic_dp(params, data)
    files <- "s_beta.rdata"
    params <- send_pause_continue_kp(params,
                                     filesAC = files,
                                     from = "AC",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    params <- comp_stws_logistic_dp(params, data)
    files <- "stwsshare.rdata"
    params <- send_pause_continue_kp(params,
                                     filesAC = files,
                                     from = "AC",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    possible_error <- ReceivedError.kp(params, from = "AC")
    if (possible_error$error) {
      params$error_message <- possible_error$message
      warning(possible_error$message)
      params <- send_pause_quit_kp(params,
                                   sleep_time = sleep_time,
                                   wait_for_turn = TRUE)
      return(params$stats)
    }

    params <- update_beta_logistic_dp(params)
    files <- "u_converge.rdata"
    params <- send_pause_continue_kp(params,
                                     filesAC = files,
                                     from = "AC",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)


    params <- get_conv_status_logistic_dp(params)

    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }

  params <- comp_s_beta_logistic_dp(params, data)

  params <- SendFinalBetasLogistic.DP(params)

  files <- c("sbeta.rdata", "finalbetas.rdata")
  params <- send_pause_continue_kp(params,
                                   filesAC = files,
                                   from = "AC",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  if (data_partner_id == 1) {
    params <- comp_results_logistic_dp(params, data)
    files <- "logisticstats.rdata"
    params <- send_pause_continue_kp(params,
                                     filesAC = files,
                                     from = "AC",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
  }

  params <- get_results_logistic_dp(params, data)
  params <- send_pause_quit_kp(params,
                               sleep_time = sleep_time,
                               wait_for_turn = TRUE)
  return(params$stats)
}


analysis_center_k_logistic <- function(num_data_partners = NULL,
                                    monitor_folder   = NULL,
                                    msreqid         = "v_default_0_000",
                                    cutoff          = 1E-8,
                                    max_iterations   = 25,
                                    sleep_time       = 10,
                                    max_waiting_time  = 24 * 60 * 60,
                                    popmednet       = TRUE,
                                    trace           = FALSE,
                                    verbose         = TRUE) {

  files_list <- rep(list(list()), num_data_partners)

  params <- prepare_params_kp("logistic", 0,
                              num_data_partners,
                              msreqid,
                              cutoff,
                              max_iterations,
                              ac = TRUE,
                              popmednet = popmednet,
                              trace = trace,
                              verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- initialize_log_kp(params)
  params <- initialize_time_stamps_kp(params)
  params <- initialize_tracking_table_kp(params)
  header(params)

  params <- prepare_folder_acdp(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.kp(params, from = "DP",
                             max_waiting_time = max_waiting_time)

  possible_error <- ReceivedError.kp(params, from = "DP")
  if (possible_error$error) {
    params$error_message <- possible_error$message
    warning(possible_error$message)
    make_error_message(params$write_path, possible_error$message)
    files <- "error_message.rdata"
    params <- send_pause_continue_kp(params,
                                     filesDP = files,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- check_agreement_ac(params)

  if (params$failed) {
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    warning(params$error_message)
    params <- send_pause_continue_kp(params,
                                     filesDP = files,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  files <- "empty.rdata"
  params <- send_pause_continue_kp(params,
                                   filesDP = files,
                                   from = "DP",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_products_logistic_ac(params)


  params <- check_colinearity_logistic_ac(params)

  if (params$failed) {
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    warning(params$error_message)
    params <- send_pause_continue_kp(params,
                                     filesDP = files,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- comp_init_betas_logistic_ac(params)

  for (id in 1:num_data_partners) {
    files_list[[id]] <- c(paste0("u_beta_", id, ".rdata"), "indicies.rdata")
  }

  params <- send_pause_continue_kp(params,
                                   filesDP = files_list,
                                   from = "DP",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- comp_weights_logistic_ac(params)
    files <- "pi.rdata"
    params <- send_pause_continue_kp(params,
                                     filesDP = files,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    params <- comp_stws_logistic_ac(params)

    if (params$failed) {
      make_error_message(params$write_path, params$error_message)
      files <- "error_message.rdata"
      warning(params$error_message)
      params <- send_pause_continue_kp(params,
                                       filesDP = files,
                                       from = "DP",
                                       sleep_time = sleep_time,
                                       max_waiting_time = max_waiting_time)
      params <- send_pause_quit_kp(params,
                                   sleep_time = sleep_time,
                                   job_failed = TRUE)
      SummarizeLog.kp(params)
      return(params$stats)
    }

    for (id in 1:num_data_partners) {
      files_list[[id]] <- paste0("id", id, ".rdata")
    }
    params <- send_pause_continue_kp(params,
                                     filesDP = files_list,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    params <- comp_conv_status_logistic_ac(params)
    files <- "u_converge.rdata"
    params <- send_pause_continue_kp(params,
                                     filesDP = files,
                                     from = "DP",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }

  params <- ComputeFinalSBetaLogistic.AC(params)
  files_list <- rep(list(list()), num_data_partners)
  files_list[[1]] <- "sbeta.rdata"
  params <- send_pause_continue_kp(params,
                                   filesDP = files_list,
                                   from = "DP1",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- comp_results_logistic_ac(params)
  files <- "stats.rdata"
  params <- send_pause_continue_kp(params,
                                   filesDP = files,
                                   from = "DP",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)
  params <- send_pause_quit_kp(params,
                               sleep_time = sleep_time)
  SummarizeLog.kp(params)
  return(params$stats)
}
