################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

check_colinearity_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "check_colinearity_logistic_t3\n\n")
  xtx <- params$xtx
  xty <- params$xty

  nrow <- nrow(xtx)
  indicies <- c(1)
  for (i in 2:nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }

  xtx <- xtx[indicies, indicies, drop = FALSE]
  xty <- xty[indicies, drop = FALSE]

  a_names               <- params$colnamesA
  b_names               <- params$colnamesB
  a_index               <- which(indicies <= length(a_names))
  params$IndiciesKeep  <- indicies
  params$a_indicies_keep <- indicies[a_index]
  params$b_indicies_keep <- indicies[-a_index] - length(a_names)
  a_names_keep           <- a_names[params$a_indicies_keep]
  b_names_keep           <- b_names[params$b_indicies_keep]
  params$colnames_a_old <- params$colnamesA
  params$colnamesB_old <- params$colnamesB
  params$colnamesA     <- a_names_keep
  params$colnamesB     <- b_names_keep
  params$p1_old        <- params$p1
  params$p2_old        <- params$p2
  params$p1            <- length(a_names_keep)
  params$p2            <- length(b_names_keep)
  params$p_old         <- params$p1_old + params$p2_old
  params$p             <- params$p1 + params$p2
  params$means_a        <- params$means_a[params$a_indicies_keep]
  params$means_b        <- params$means_b[params$b_indicies_keep]
  params$sda           <- params$sda[params$a_indicies_keep]
  params$sdb           <- params$sdb[params$b_indicies_keep]
  params$xtx           <- xtx
  params$xty           <- xty

  a_indicies <- params$a_indicies_keep
  b_indicies <- params$b_indicies_keep

  write_time <- proc.time()[3]
  save(a_indicies, file = file.path(params$write_path, "Aindicies.rdata"))
  save(b_indicies, file = file.path(params$write_path, "Bindicies.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("Aindicies.rdata",
                                                             "Bindicies.rdata"))))
  write_time <- proc.time()[3] - write_time

  b_tags <- params$b_tags[params$b_indicies_keep]
  a_tags <- params$a_tags[params$a_indicies_keep][-1]

  if ((length(unique(a_tags)) == 1) ||
      (length(unique(a_tags)) >= 2 && !("numeric" %in% names(a_tags)))) {
    params$failed <- TRUE
    params$error_message <-
      paste("A must have no covariates or at least 2",
            "covariates at least one of which is continuous.")
  } else if (length(unique(b_tags)) < 2) {
    params$failed <- TRUE
    params$error_message <-
      "After removing colinear covariates, Party B has 1 or fewer covariates."
  } else if (!("numeric" %in% names(b_tags))) {
    params$failed <- TRUE
    params$error_message <-
      paste("After removing colinear covariates,",
            "Party B has no continuous covariates.")
  }

  params <- add_to_log(params, "check_colinearity_logistic_t3",
                       0, 0, write_time, write_size)
  return(params)
}


compute_initial_betas_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_initial_betas_logistic_t3\n\n")
  # de-standardize xty
  p1     <- params$p1
  p2     <- params$p2
  xty    <- params$xty
  xtx    <- params$xtx

  betas <- 4 * solve(xtx) %*% xty

  a_betas   <- betas[1:p1]
  b_betas   <- betas[(p1 + 1):(p1 + p2)]
  a_xty     <- xty[1:p1]
  b_xty     <- xty[(p1 + 1):(p1 + p2)]

  params$a_xty      <- a_xty
  params$b_xty      <- b_xty
  params$betas     <- betas
  params$betas_a    <- a_betas
  params$betas_a_old  <- matrix(0, p1, 1)
  params$betas_b    <- b_betas

  params$alg_iteration_counter      <- 1
  params$delta_beta <- Inf
  params$converged <- FALSE
  converged <- FALSE
  max_iter_exceeded <- FALSE

  write_time <- proc.time()[3]
  save(a_betas, file = file.path(params$write_path, "betasA.rdata"))
  save(p2, a_xty,   file = file.path(params$write_path, "a_xty.rdata"))
  save(b_betas, file = file.path(params$write_path, "betasB.rdata"))
  save(b_xty,   file = file.path(params$write_path, "b_xty.rdata"))
  save(converged, max_iter_exceeded,
       file = file.path(params$write_path, "converged.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("betasA.rdata",
                                                             "betasB.rdata",
                                                             "a_xty.rdata",
                                                             "b_xty.rdata",
                                                             "converged.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_initial_betas_logistic_t3",
                       0, 0, write_time, write_size)

  return(params)
}


update_params_logistic_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_params_logistic_a3\n\n")
  a_indicies <- NULL
  a_xty      <- NULL
  p2        <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Aindicies.rdata"))
  load(file.path(params$read_path[["T"]], "a_xty.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]],
                                       c("Aindicies.rdata", "a_xty.rdata"))))

  read_time <- proc.time()[3] - read_time
  params$colnames_a_old <- params$colnamesA
  params$colnamesA     <- params$colnames_a_old[a_indicies]
  params$p_old         <- params$p
  params$p             <- length(a_indicies)
  params$p2            <- p2
  params$a_indicies_keep <- a_indicies
  params$means         <- params$means[a_indicies]
  params$sd            <- params$sd[a_indicies]
  params$a_xty          <- a_xty
  params <- add_to_log(params,
                       "update_params_logistic_a3, update_data_logistic_a3",
                       read_time, read_size, 0, 0)
  return(params)
}


update_params_logistic_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_params_logistic_b3\n\n")
  b_indicies <- NULL
  b_xty      <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Bindicies.rdata"))
  load(file.path(params$read_path[["T"]], "b_xty.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]],
                                       c("Bindicies.rdata", "b_xty.rdata"))))

  read_time <- proc.time()[3] - read_time
  params$colnamesB_old <- params$colnamesB
  params$colnamesB     <- params$colnamesB_old[b_indicies]
  params$p_old         <- params$p
  params$p             <- length(b_indicies)
  params$b_indicies_keep <- b_indicies
  params$means         <- params$means[b_indicies]
  params$sd            <- params$sd[b_indicies]
  params$b_xty          <- b_xty
  params <- add_to_log(params,
                       "update_params_logistic_b3, update_data_logistic_b3",
                       read_time, read_size, 0, 0)
  return(params)
}


update_data_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_logistic_a3\n\n")
  data$x <- as.matrix(data$x[, params$a_indicies_keep, drop = FALSE])
  return(data)
}


update_data_logistic_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_logistic_b3\n\n")
  data$x <- as.matrix(data$x[, params$b_indicies_keep, drop = FALSE])
  return(data)
}


get_beta_a_logistic_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.a3\n\n")
  converged       <- NULL
  max_iter_exceeded <- NULL
  a_betas          <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasA.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]],
                                       c("converged.rdata", "betasA.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$converged <- converged
  params$max_iter_exceeded <- max_iter_exceeded
  params$betas <- a_betas
  params <- add_to_log(params, "get_beta_a_logistic_a3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_beta_b_logistic_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.b3\n\n")
  converged       <- NULL
  max_iter_exceeded <- NULL
  b_betas          <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasB.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]],
                                       c("converged.rdata", "betasB.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$converged <- converged
  params$max_iter_exceeded <- max_iter_exceeded
  params$betas <- b_betas
  params <- add_to_log(params, "get_beta_b_logistic_b3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_xa_beta_a_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaALogistic.a3\n\n")
  xa_beta <- data$x %*% params$betas

  write_time <- proc.time()[3]
  save(xa_beta, file = file.path(params$write_path, "xabeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "xabeta.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_xa_beta_a_logistic_a3",
                       0, 0, write_time, write_size)
  return(params)
}


get_xb_beta_b_logistic_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaBLogistic.b3\n\n")
  xb_beta <- data$x %*% params$betas

  write_time <- proc.time()[3]
  save(xb_beta, file = file.path(params$write_path, "xbbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbbeta.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_xb_beta_b_logistic_b3",
                       0, 0, write_time, write_size)
  return(params)
}


get_weights_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_weights_logistic_t3\n\n")
  xa_beta <- NULL
  xb_beta <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "xabeta.rdata"))
  load(file.path(params$read_path[["B"]], "xbbeta.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "xabeta.rdata")) +
    file.size(file.path(params$read_path[["B"]], "xbbeta.rdata"))
  read_time <- proc.time()[3] - read_time

  x_beta <- xa_beta + xb_beta
  pi_ <- (1 + exp(-x_beta))^(-1)
  params$pi_ <- pi_

  write_time <- proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi.rdata"))
  write_size <- file.size(file.path(params$write_path, "pi.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "get_weights_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_rv_logistic_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_rv_logistic_b3\n\n")
  pi_ <- NULL
  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "pi.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "pi.rdata"))
  read_time <- proc.time()[3] - read_time

  params$pi_ <- pi_
  w <- pi_ * (1 - params$pi_)
  xb_t_w_xb <- 0
  pbar <- make_progress_bar_1(params$blocks$num_blocks,
                              "R(I-z*z')w*XB", params$verbose)
  container_ct_rz <- 0
  container_ct_rv <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak_rz) {
      container_ct_rz <- container_ct_rz + 1
      filename1 <- paste0("crz_", container_ct_rz, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak_rv) {
      container_ct_rv <- container_ct_rv + 1
      filename2 <- paste0("crv_", container_ct_rv, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    x_block  <- data$x[strt:stp, , drop = FALSE]
    w_block  <- w[strt:stp]
    wx_block <- MultiplyDiagonalWTimesX(w_block, x_block)

    read_time <- read_time - proc.time()[3]
    rz <- matrix(readBin(con = to_read, what = numeric(), n = n * n,
                         endian = "little"), nrow = n, ncol = n)
    read_time <- read_time + proc.time()[3]

    rv <- rz %*% wx_block

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(rv), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    xb_t_w_xb <- xb_t_w_xb + t(x_block) %*% wx_block

    if ((i + 1) %in% params$container$filebreak_rz ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak_rv ||
        i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  write_time <- write_time - proc.time()[3]
  save(xb_t_w_xb, file = file.path(params$write_path, "xbtwxb.rdata"))
  write_size <- write_size +
    sum(file.size(c(file.path(params$write_path, "xbtwxb.rdata"))))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "get_rv_logistic_b3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


process_v_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "process_v_logistic_t3\n\n")
  xb_t_w_xb <- NULL
  write_time <- 0
  write_size <- 0
  p2 <- params$p2
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["B"]], "xbtwxb.rdata"))
  read_size <- file.size(file.path(params$read_path[["B"]], "xbtwxb.rdata"))
  read_time <- proc.time()[3] - read_time

  params$xbtwxb <- xb_t_w_xb

  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "(I-z*z')w*XB*R", params$verbose)

  container_ct_rv <- 0
  container_ct_vr <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak_rv) {
      container_ct_rv <- container_ct_rv + 1
      filename2 <- paste0("crv_", container_ct_rv, ".rdata")
      to_read_2 <- file(file.path(params$read_path[["B"]], filename2), "rb")
    }
    if (i %in% params$container$filebreak_vr) {
      container_ct_vr <- container_ct_vr + 1
      filename3 <- paste0("cvr_", container_ct_vr, ".rdata")
      to_write3 <- file(file.path(params$write_path, filename3), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    filename1 <- paste0("r1_", i, ".rdata")
    filename4 <- paste0("r3_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    to_read_1 <- file(file.path(params$dplocalPath, filename1), "rb")
    r2  <- matrix(readBin(con = to_read_1, what = numeric(), n = n * n,
                          endian = "little"), nrow = n, ncol = n)
    read_size <- read_size + file.size(file.path(params$dplocalPath, filename1))
    close(to_read_1)
    rv  <- matrix(readBin(con = to_read_2, what = numeric(), n = n * p2,
                          endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    v <- t(r2) %*% rv
    r_3 <- random_orthonormal_matrix(p2)
    vr <- v %*% r_3

    write_time <- write_time - proc.time()[3]
    to_write4 <- file(file.path(params$dplocalPath, filename4), "wb")
    writeBin(as.vector(r_3), con = to_write4, endian = "little")
    close(to_write4)
    write_size <- write_size +
      file.size(file.path(params$dplocalPath, filename4))
    writeBin(as.vector(vr), con = to_write3, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak_rv || i == num_blocks) {
      close(to_read_2)
      read_size <- read_size +
        file.size(file.path(params$dplocalPath, filename1))
    }
    if ((i + 1) %in% params$container$filebreak_vr || i == num_blocks) {
      close(to_write3)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename3))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "process_v_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_xr_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_xr_logistic_a3\n\n")
  pi_ <- NULL
  p2 <- params$p2
  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "pi.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "pi.rdata"))
  read_time <- proc.time()[3] - read_time

  params$pi_ <- pi_
  w <- pi_ * (1 - params$pi_)
  xa_t_w_xa <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks,
                              "XA'(I-z*z')w*XB*R", params$verbose)

  container_ct_vr <- 0
  container_ct_xr <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak_rv) {
      container_ct_vr <- container_ct_vr + 1
      filename1 <- paste0("cvr_", container_ct_vr, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak_xr) {
      container_ct_xr <- container_ct_xr + 1
      filename2 <- paste0("cxr_", container_ct_xr, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    x_block  <- data$x[strt:stp, , drop = FALSE]
    w_block  <- w[strt:stp]
    wx_block <- MultiplyDiagonalWTimesX(w_block, x_block)

    read_time <- read_time - proc.time()[3]
    vr  <- matrix(readBin(con = to_read, what = numeric(), n = n * p2,
                          endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    xr <- t(x_block) %*% vr

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(xr), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    xa_t_w_xa <- xa_t_w_xa + t(x_block) %*% wx_block

    if ((i + 1) %in% params$container$filebreak_vr ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak_xr ||
        i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  write_time <- write_time - proc.time()[3]
  save(xa_t_w_xa, file = file.path(params$write_path, "xatwxa.rdata"))
  write_size <- write_size +
    sum(file.size(c(file.path(params$write_path, "xatwxa.rdata"))))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "get_xr_logistic_a3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


process_xt_wx_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "process_xt_wx_logistic_t3\n\n")
  xa_t_w_xa <- NULL
  p1 <- params$p1
  p2 <- params$p2

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "xatwxa.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "xatwxa.rdata"))
  read_time <- proc.time()[3] - read_time

  params$xatwxa <- xa_t_w_xa

  pbar <- make_progress_bar_1(params$blocks$num_blocks, "x'w*x", params$verbose)
  container_ct_xr <- 0
  xa_t_w_xb <- 0

  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak_xr) {
      container_ct_xr <- container_ct_xr + 1
      filename1 <- paste0("cxr_", container_ct_xr, ".rdata")
      to_read <- file(file.path(params$read_path[["A"]], filename1), "rb")
    }

    filename2 <- paste0("r3_", i, ".rdata")
    read_time <- read_time - proc.time()[3]
    to_read_1 <- file(file.path(params$dplocalPath, filename2), "rb")
    r  <- matrix(readBin(con = to_read_1, what = numeric(), n = p2 * p2,
                         endian = "little"), nrow = p2, ncol = p2)
    close(to_read_1)
    xr  <- matrix(readBin(con = to_read, what = numeric(), n = p1 * p2,
                          endian = "little"), nrow = p1, ncol = p2)
    read_size <- read_size + file.size(file.path(params$dplocalPath, filename2))
    read_time <- read_time + proc.time()[3]

    xa_t_w_xb <- xa_t_w_xb + xr %*% t(r)

    if ((i + 1) %in% params$container$filebreak_xr ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["A"]], filename1))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  xtwx <- rbind(cbind(params$xatwxa, xa_t_w_xb),
                cbind(t(xa_t_w_xb), params$xbtwxb))
  params$xtwx <- xtwx

  ii <- NULL
  tryCatch({
    ii <- solve(xtwx)
  },
  error = function(err) {
    ii <- NULL
  }
  )
  if (is.null(ii)) {
    params$failed <- TRUE
    params$singular_matrix <- TRUE
    params$error_message <-
      paste0("The matrix t(x)*w*x is not invertible.\n",
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
  params$ii <- ii
  ii_a <- ii[, 1:p1, drop = FALSE]
  ii_b <- ii[, (p1 + 1):(p1 + p2), drop = FALSE]
  write_time <- proc.time()[3]
  save(ii_a, file = file.path(params$write_path, "IIA.rdata"))
  save(ii_b, file = file.path(params$write_path, "IIB.rdata"))
  write_size <- sum(file.size(file.path(params$write_path,
                                        c("IIA.rdata", "IIB.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "process_xt_wx_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


update_beta_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_beta_logistic_a3\n\n")
  ii_a <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "IIA.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "IIA.rdata"))
  read_time <- proc.time()[3] - read_time

  i_a <- params$a_xty - t(data$x) %*% params$pi_
  ai <- ii_a %*% i_a


  write_time <- proc.time()[3]
  save(ai, file = file.path(params$write_path, "AI.rdata"))
  write_size <- file.size(file.path(params$write_path, "AI.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "update_beta_logistic_a3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


update_beta_logistic_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_beta_logistic_b3\n\n")
  ii_b <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "IIB.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "IIB.rdata"))
  read_time <- proc.time()[3] - read_time

  i_b <- params$b_xty - t(data$x) %*% params$pi_
  bi <- ii_b %*% i_b

  write_time <- proc.time()[3]
  save(bi, file = file.path(params$write_path, "BI.rdata"))
  write_size <- file.size(file.path(params$write_path, "BI.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "update_beta_logistic_b3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


update_beta_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "update_beta_logistic_t3\n\n")
  ai <- NULL
  bi <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "AI.rdata"))
  load(file.path(params$read_path[["B"]], "BI.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "AI.rdata")) +
    file.size(file.path(params$read_path[["B"]], "BI.rdata"))
  read_time <- proc.time()[3] - read_time

  delta <- ai + bi
  betas <- params$betas + delta
  params$betas <- betas
  converged <- all(abs(delta) / (abs(betas) + .1) < params$cutoff)
  max_iter_exceeded <- (params$alg_iteration_counter >=
                          params$max_iterations) && !converged

  params$converged <- converged
  params$max_iter_exceeded <- max_iter_exceeded
  a_betas <- betas[1:params$p1]
  b_betas <- betas[(params$p1 + 1):(params$p1 + params$p2)]

  write_time <- proc.time()[3]
  save(converged, max_iter_exceeded,
       file = file.path(params$write_path, "converged.rdata"))
  save(a_betas, file = file.path(params$write_path, "betasA.rdata"))
  save(b_betas, file = file.path(params$write_path, "betasB.rdata"))
  write_size <- sum(file.size(file.path(params$write_path,
                                        c("betasA.rdata",
                                          "betasB.rdata",
                                          "converged.rdata"))))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "update_beta_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_final_beta_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.a3\n\n")
  betas <- params$betas / params$sd
  offset_a <- sum(betas[-1] * params$means[-1])
  a_final_fitted <- t(params$sd * t(data$x) + params$means) %*% betas -
    t(params$sd[1] * t(data$x[, 1]) + params$means[1]) %*% betas[1]
  write_time <- proc.time()[3]
  save(offset_a, a_final_fitted,
       file = file.path(params$write_path, "Afinalfitted.rdata"))
  write_size <- file.size(file.path(params$write_path, "Afinalfitted.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "get_final_beta_logistic_a3",
                       0, 0, write_time, write_size)
  return(params)
}


get_final_beta_logistic_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.b3\n\n")
  betas <- params$betas / params$sd
  offset_b <- sum(betas * params$means)
  b_final_fitted <- t(params$sd * t(data$x) + params$means) %*% betas

  write_time <- proc.time()[3]
  save(offset_b, b_final_fitted,
       file = file.path(params$write_path, "Bfinalfitted.rdata"))
  write_size <- file.size(file.path(params$write_path, "Bfinalfitted.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "get_final_beta_logistic_b3",
                       0, 0, write_time, write_size)
  return(params)
}


get_final_fitted_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "get_final_fitted_logistic_t3\n\n")
  offset_a <- NULL
  offset_b <- NULL
  a_final_fitted <- NULL
  b_final_fitted <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "Afinalfitted.rdata"))
  load(file.path(params$read_path[["B"]], "Bfinalfitted.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]],
                                   "Afinalfitted.rdata")) +
    file.size(file.path(params$read_path[["B"]], "Bfinalfitted.rdata"))
  read_time <- proc.time()[3] - read_time

  betas <- params$betas / c(params$sda, params$sdb)
  betas[1] <- betas[1] - offset_a - offset_b

  final_fitted <- a_final_fitted + b_final_fitted + betas[1]

  params$betas <- betas
  params$final_fitted <- final_fitted

  write_time <- proc.time()[3]
  save(final_fitted, file = file.path(params$write_path, "finalFitted.rdata"))
  write_size <- file.size(file.path(params$write_path, "finalFitted.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_final_beta_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}

compute_results_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_results_logistic_a3\n\n")
  final_fitted <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "finalFitted.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]],
                                   "finalFitted.rdata"))
  read_time <- proc.time()[3] - read_time

  n <- params$n
  ct      <- sum(data$Y)
  params$final_fitted <- final_fitted
  resdev  <- -2 * (sum(data$Y * final_fitted) - sum(log(1 + exp(final_fitted))))
  nulldev <- -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))
  hoslem  <- HoslemInternal(params, data)
  roc     <- roc_internal(params, data)

  write_time <- proc.time()[3]
  save(resdev, nulldev, hoslem, roc,
       file = file.path(params$write_path, "logisticstats.rdata"))
  write_size <- file.size(file.path(params$write_path, "logisticstats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_results_logistic_a3",
                       read_time, read_size, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
compute_results_logistic_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_results_logistic_t3\n\n")
  nulldev <- NULL
  resdev  <- NULL
  hoslem  <- NULL
  roc     <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "logisticstats.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]],
                                   "logisticstats.rdata"))
  read_time <- proc.time()[3] - read_time

  stats <- params$stats
  stats$failed         <- FALSE
  stats$converged      <- params$converged

  n      <- params$n
  p1     <- params$p1
  p2     <- params$p2
  sda    <- params$sda
  sdb    <- params$sdb
  means_a <- params$means_a
  means_b <- params$means_b
  a_names <- params$colnames_a_old
  b_names <- params$colnamesB_old
  p1_old <- params$p1_old
  p2_old <- params$p2_old
  p_old  <- params$p_old
  indicies <- params$IndiciesKeep

  # If xtwx were singular, it would have been caught in GetII.A2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  cov1 <- solve(params$xtwx)
  secoef <- sqrt(diag(cov1)) / c(sda, sdb)
  tmp  <- matrix(c(1, (-means_a / sda)[-1], -means_b / sdb), ncol = 1)
  secoef[1] <- sqrt(t(tmp) %*% cov1 %*% tmp)

  stats$party <- c(rep("dp1", p1_old), rep("dp2", p2_old))
  stats$coefficients <- rep(NA, p_old)
  stats$secoef <- rep(NA, p_old)
  stats$tvals  <- rep(NA, p_old)
  stats$pvals  <- rep(NA, p_old)
  stats$n  <- n
  stats$nulldev <- nulldev
  stats$resdev <- resdev
  stats$aic <- resdev + 2 * (p1 + p2)
  stats$bic <- resdev + (p1 + p2) * log(n)
  stats$nulldev_df <- n - 1
  stats$resdev_df <- n - (p1 + p2)
  stats$coefficients[indicies] <- params$betas
  stats$secoef[indicies] <- secoef
  tvals <- params$betas / secoef
  pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[indicies] <- tvals
  stats$pvals[indicies] <- pvals
  stats$hoslem  <- hoslem
  stats$roc     <- roc
  stats$iter    <- params$alg_iteration_counter - 1
  names_old <- c(a_names, b_names)
  names(stats$coefficients) <- names_old
  names(stats$party) <- names_old
  names(stats$secoef) <- names_old
  names(stats$tvals) <- names_old
  names(stats$pvals) <- names_old

  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params$stats      <- stats

  params <- add_to_log(params, "compute_results_logistic_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_results_logistic_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_logistic_a3\n\n")
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  stats$Y           <- data$Y # For Hoslem and ROC
  stats$final_fitted <- params$final_fitted
  params$stats      <- stats
  params <- add_to_log(params, "get_results_logistic_a3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_results_logistic_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_logistic_b3\n\n")
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats <- stats
  params <- add_to_log(params, "get_results_logistic_b3",
                       read_time, read_size, 0, 0)
  return(params)
}


############################### PARENT FUNCTIONS ###############################


party_a_process_3_logistic <- function(data,
                                       y_name          = NULL,
                                       monitor_folder  = NULL,
                                       sleep_time      = 10,
                                       max_waiting_time = 24 * 60 * 60,
                                       popmednet      = TRUE,
                                       trace          = FALSE,
                                       verbose        = TRUE) {

  params <- prepare_params_3p("logistic", "A",
                              popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)
  header(params)

  params   <- prepare_folder_linear_a3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_logistic_a23(params, data, y_name)
  params <- add_to_log(params, "prepare_data_logistic_a23", 0, 0, 0, 0)

  if (data$failed) {
    message <- "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params, files_t = files,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_a3(params, data)
  files <- "pa.rdata"
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params$alg_iteration_counter <- 1
  params <- prepare_blocks_linear_a3(params)

  params <- get_z_linear_a3(params, data)
  files <- seq_zw("cz_", length(params$container$file_break_z))
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_wr_linear_a3(params, data)
  files <- c("xatxa.rdata",
             seq_zw("cpr_", length(params$container$filebreak_pr)))
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- update_params_logistic_a3(params)
  data <- update_data_logistic_a3(params, data)
  params <- get_beta_a_logistic_a3(params)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- get_xa_beta_a_logistic_a3(params, data)
    files <- c("xabeta.rdata")
    params <- send_pause_continue_3p(params,
                                     files_t = files,
                                     from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    params <- get_xr_logistic_a3(params, data)
    files <- c("xatwxa.rdata",
               seq_zw("cxr_", length(params$container$filebreak_xr)))
    params <- send_pause_continue_3p(params,
                                     files_t = files,
                                     from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    if (file.exists(file.path(params$read_path[["T"]],
                              "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["T"]]))
      params <- send_pause_quit_3p(params,
                                   sleep_time = sleep_time,
                                   job_failed = TRUE,
                                   wait_for_turn = TRUE)
      return(params$stats)
    }

    params <- update_beta_logistic_a3(params, data)
    files <- c("ai.rdata")
    params <- send_pause_continue_3p(params,
                                     files_t = files,
                                     from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    params <- get_beta_a_logistic_a3(params)
    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }

  params <- get_final_beta_logistic_a3(params, data)
  files <- "Afinalfitted.rdata"
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  params <- compute_results_logistic_a3(params, data)
  files <- c("logisticstats.rdata")
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_results_logistic_a3(params, data)

  params <- send_pause_quit_3p(params,
                               sleep_time = sleep_time,
                               wait_for_turn = TRUE)
  return(params$stats)
}


party_b_process_3_logistic <- function(data,
                                       monitor_folder  = NULL,
                                       sleep_time      = 10,
                                       max_waiting_time = 24 * 60 * 60,
                                       popmednet      = TRUE,
                                       trace          = FALSE,
                                       verbose        = TRUE) {
  params <- prepare_params_3p("logistic", "B",
                              popmednet = popmednet,
                              trace = trace,
                              verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)

  header(params)
  params   <- prepare_folder_linear_b3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  data <- prepare_data_logistic_b23(params, data)
  params <- add_to_log(params, "prepare_data_logistic_b23", 0, 0, 0, 0)

  if (data$failed) {
    message <- "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params,
                                 files_t = files,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_b3(params, data)
  files <- "pb.rdata"
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params$alg_iteration_counter <- 1
  params <- prepare_blocks_linear_b3(params)

  params <- get_rw_linear_b3(params, data)
  files <- c("xbtxb.rdata",
             seq_zw("crw_", length(params$container$filebreak_RW)))
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE,
                                 wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- update_params_logistic_b3(params)
  data <- update_data_logistic_b3(params, data)
  params <- get_beta_b_logistic_b3(params)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- get_xb_beta_b_logistic_b3(params, data)
    files <- c("xbbeta.rdata")
    params <- send_pause_continue_3p(params,
                                     files_t = files
                                     , from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    params <- get_rv_logistic_b3(params, data)
    files <- c("xbtwxb.rdata",
               seq_zw("crv_", length(params$container$filebreak_rv)))
    params <- send_pause_continue_3p(params,
                                     files_t = files,
                                     from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    if (file.exists(file.path(params$read_path[["T"]],

                              "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["T"]]))
      params <- send_pause_quit_3p(params,
                                   sleep_time = sleep_time,
                                   job_failed = TRUE,
                                   wait_for_turn = TRUE)
      return(params$stats)
    }

    params <- update_beta_logistic_b3(params, data)
    files <- c("bi.rdata")
    params <- send_pause_continue_3p(params,
                                     files_t = files,
                                     from = "T",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time,
                                     wait_for_turn = TRUE)

    params <- get_beta_b_logistic_b3(params)
    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }

  params <- get_final_beta_logistic_b3(params, data)
  files <- "Bfinalfitted.rdata"
  params <- send_pause_continue_3p(params,
                                   files_t = files,
                                   from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  params <- get_results_logistic_b3(params)

  params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                               wait_for_turn = TRUE)
  return(params$stats)
}


party_t_process_3_logistic <- function(monitor_folder   = NULL,
                                       msreqid          = "v_default_0_000",
                                       blocksize        = 500,
                                       cutoff           = 1e-8,
                                       max_iterations   = 25,
                                       sleep_time       = 10,
                                       max_waiting_time = 24 * 60 * 60,
                                       popmednet        = TRUE,
                                       trace            = FALSE,
                                       verbose          = TRUE) {
  params <- prepare_params_3p("logistic", "T",
                              msreqid = msreqid,
                              popmednet = popmednet,
                              trace = trace,
                              verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)

  header(params)
  params   <- prepare_folder_linear_t3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.3p(params, from = c("A", "B"),
                             max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata")) &&
      file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["A"]]))
    warning(read_error_message(params$read_path[["B"]]))
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["A"]]))
    file.copy(file.path(params$read_path[["A"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params,
                                     files_b = files,
                                     from = "B",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["B"]]))
    file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params,
                                     files_a = files,
                                     from = "A",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params   <- prepare_params_linear_t3(params, cutoff, max_iterations)

  if (!params$failed) params <- prepare_blocks_linear_t3(params, blocksize)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params,
                                     files_a = files,
                                     files_b = files,
                                     from = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files <- "blocks.rdata"
  params <- send_pause_continue_3p(params,
                                   files_a = files,
                                   from = "A",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params$alg_iteration_counter <- 1
  params <- process_z_linear_t3(params)
  files <- c("blocks.rdata",
             seq_zw("crz_", length(params$container$filebreak_rz)))
  params <- send_pause_continue_3p(params,
                                   files_b = files,
                                   from  = "B",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- process_w_linear_t3(params)
  files <- c("p2.rdata", seq_zw("cwr_", length(params$container$filebreak_wr)))
  params <- send_pause_continue_3p(params,
                                   files_a = files,
                                   from  = "A",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_products_linear_t3(params)

  params <- check_colinearity_logistic_t3(params)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params,
                                     files_a = files,
                                     files_b = files,
                                     from = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params <- compute_initial_betas_logistic_t3(params)
  files_a <-
    c("Aindicies.rdata", "betasA.rdata", "a_xty.rdata", "converged.rdata")
  files_b <-
    c("Bindicies.rdata", "betasB.rdata", "b_xty.rdata", "converged.rdata")
  params <- send_pause_continue_3p(params,
                                   files_a = files_a,
                                   files_b = files_b,
                                   from  = c("A", "B"),
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- get_weights_logistic_t3(params)
    files <- "pi.rdata"
    params <- send_pause_continue_3p(params,
                                     files_b = files,
                                     from  = "B",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    params <- process_v_logistic_t3(params)
    files <- c("pi.rdata",

               seq_zw("cvr_", length(params$container$filebreak_rv)))
    params <- send_pause_continue_3p(params,
                                     files_a = files,
                                     from  = "A",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    params <- process_xt_wx_logistic_t3(params)

    if (params$failed) {
      warning(params$error_message)
      make_error_message(params$write_path, params$error_message)
      files <- "error_message.rdata"
      params <- send_pause_continue_3p(params,
                                       files_a = files,
                                       files_b = files,
                                       from = c("A", "B"),
                                       sleep_time = sleep_time,
                                       max_waiting_time = max_waiting_time)
      params <- send_pause_quit_3p(params, sleep_time = sleep_time)
      SummarizeLog.3p(params)
      return(params$stats)
    }
    files_a <- c("IIA.rdata")
    files_b <- c("IIB.rdata")
    params <- send_pause_continue_3p(params,
                                     files_a = files_a,
                                     files_b = files_b,
                                     from  = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)

    params <- update_beta_logistic_t3(params)
    files_a <- c("betasA.rdata", "converged.rdata")
    files_b <- c("betasB.rdata", "converged.rdata")
    params <- send_pause_continue_3p(params,
                                     files_a = files_a,
                                     files_b = files_b,
                                     from  = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }

  params <- get_final_fitted_logistic_t3(params)
  files_a <- "finalfitted.rdata"
  params <- send_pause_continue_3p(params,
                                   files_a = files_a,
                                   from  = "A",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- compute_results_logistic_t3(params)
  files <- "stats.rdata"
  params <- send_pause_continue_3p(params,
                                   files_a = files,
                                   files_b = files,
                                   from  = c("A", "B"),
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- send_pause_quit_3p(params, sleep_time = sleep_time)
  SummarizeLog.3p(params)
  return(params$stats)
}
