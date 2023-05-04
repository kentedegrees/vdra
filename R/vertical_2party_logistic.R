################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

prepare_folder_logistic_a2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_folder_logistic_a2\n\n")
  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "inputfiles")
  params$read_path       <- file.path(monitor_folder, "msoc1")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  if (!is.character(monitor_folder)) {
    warning("monitor_folder directory is not valid. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }
  params$error_message <- NULL
  if (!create_io_location(monitor_folder, "dplocal")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dp_local_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "rprograms")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$r_programs_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "macros")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macros_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "inputfiles")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc1")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  params <- add_to_log(params,
                       "prepare_data_logistic_a23, prepare_folder_logistic_a2",
                       0, 0, 0, 0)
  return(params)
}


prepare_folder_logistic_b2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_folder_logistic_b2\n\n")

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "msoc")
  params$read_path       <- file.path(monitor_folder, "inputfiles")

  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  if (!is.character(monitor_folder)) {
    warning("monitor_folder directory is not valid. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }
  params$error_message <- NULL
  if (!create_io_location(monitor_folder, "dplocal")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dp_local_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "rprograms")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$r_programs_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "macros")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$macros_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "inputfiles")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "ould not create directory",
                                paste0(params$read_path, "."),
                                "Check the path and restart the program.")
  }

  Sys.sleep(1)
  delete_trigger("files_done.ok", params$read_path)

  params <- add_to_log(params,
                       "PrepareDataLogisitc.b23, prepare_folder_logistic_b2",
                       0, 0, 0, 0)

  return(params)
}


#' @importFrom stats model.matrix
prepare_data_logistic_a23 <- function(params, data, y_name = NULL) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_data_logistic_a23\n\n")

  workdata <- list()
  workdata$failed <- FALSE

  workdata$failed <- check_data_format(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data <- data.frame(data) # convert to a clean data.frame

  response_index <- check_response(params, data, y_name)

  if (is.null(response_index)) {
    workdata$failed <- TRUE
    return(workdata)
  }
  covariate_index <- setdiff(seq_len(ncol(data)), response_index)
  workdata$tags <- create_model_matrix_tags(data[, covariate_index,
                                                 drop = FALSE])
  workdata$tags <- c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] <- "numeric"

  x <- model.matrix(~ ., data[, c(response_index, covariate_index),
                              drop = FALSE])
  rownames(x) <- NULL
  covariate_index <- setdiff(seq_len(ncol(x)), 2)

  means <- apply(x, 2, mean)
  sd    <- apply(x, 2, sd)
  sd    <- sapply(sd, function(x) {
    ifelse(x > 0, x, 1)
  })
  workdata$Y      <- x[, 2, drop = FALSE]
  workdata$x      <- x[, covariate_index, drop = FALSE]
  workdata$means  <- means[covariate_index]
  workdata$sd     <- sd[covariate_index]
  workdata$yty    <- t(workdata$Y) %*% workdata$Y

  if (ncol(workdata$x) >= 2) {
    for (i in 2:ncol(workdata$x)) {
      workdata$x[, i] <- (workdata$x[, i] - workdata$means[i]) / workdata$sd[i]
    }
  }

  return(workdata)
}

#' @importFrom stats model.matrix sd
prepare_data_logistic_b23 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_data_logistic_b23\n\n")

  workdata <- list()
  workdata$failed <- FALSE

  workdata$failed <- check_data_format(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data <- data.frame(data) # convert to a clean data.frame

  workdata$tags <- create_model_matrix_tags(data)
  if (ncol(data) < 2 || !("numeric" %in% names(workdata$tags))) {
    warning("The data partner that does not have the response must have ",
            "at least 2 covariates at least one of which must be numeric.")
    workdata$failed <- TRUE
    return(workdata)
  }
  workdata$x <- model.matrix(~ ., data)
  rownames(workdata$x) <- NULL
  workdata$x <- workdata$x[, -1, drop = FALSE]
  workdata$means <- apply(workdata$x, 2, mean)
  workdata$sd    <- apply(workdata$x, 2, sd)
  workdata$sd    <- sapply(workdata$sd, function(x) {
    ifelse(x > 0, x, 1)
  })

  for (i in seq_len(ncol(workdata$x))) {
    workdata$x[, i] <- (workdata$x[, i] - workdata$means[i]) / workdata$sd[i]
  }


  return(workdata)
}

prepare_params_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_params_logistic_b2\n\n")
  params$failed         <- FALSE
  params$converged      <- FALSE
  params$halted         <- FALSE

  params$n             <- nrow(data$x)
  params$num_events     <- 0
  params$p1            <- 0
  params$p2            <- ncol(data$x)
  params$p             <- params$p1 + params$p2
  params$p1_old        <- 0
  params$p2_old        <- params$p2
  params$a_col_names     <- c("")
  params$b_col_names     <- colnames(data$x)
  params$y_name         <- ""
  params$a_col_names_old <- c("")
  params$b_col_names_old <- c("")
  params$cutoff        <- 1
  params$max_iterations <- 1

  params$means_a        <- 0
  params$sda           <- 0
  params$means_b        <- data$means
  params$sdb           <- data$sd
  params$yty           <- 0

  pb          <- list()
  pb$p2       <- params$p2
  pb$n        <- params$n
  pb$means    <- data$means
  pb$sd       <- data$sd
  pb$analysis <- params$analysis
  pb$b_col_names <- params$b_col_names
  pb$tags      <- data$tags

  write_time <- proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_logistic_b2",
                       0, 0, write_time, write_size)
  return(params)
}


prepare_params_logistic_a2 <- function(params, data,
                                       cutoff = 0.01, max_iterations = 25) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_params_logistic_a2\n\n")

  params$converged       <- FALSE
  params$halted          <- FALSE
  params$pmn_step_counter  <- 1
  pb                     <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb, b_col_names
  read_size <- sum(file.size(file.path(params$read_path, "pb.rdata")))
  read_time <- proc.time()[3] - read_time

  if (params$analysis != pb$analysis) {
    params$error_message <-
      paste("Party A is running", params$analysis,
            "regression and Party B is running", pb$analysis, "regression.")
    warning(params$error_message)
    params$failed <- TRUE
    return(params)
  }

  params$n <- nrow(data$x)
  if (pb$n != params$n) {
    params$error_message <-
      paste("Party A has", params$n,
            "observations and Party B has", pb$n, "observations.")
    warning(params$error_message)
    params$failed <- TRUE
  }

  params$p1 <- ncol(data$x)
  params$p2 <- pb$p2
  params$p  <- params$p1 + params$p2
  params$p1_old <- params$p1
  params$p2_old <- params$p2

  params$a_col_names <- colnames(data$x)
  params$b_col_names <- pb$b_col_names
  params$y_name     <- colnames(data$Y)
  params$a_col_names_old <- c("")
  params$b_col_names_old <- c("")
  params$a_tags         <- data$tags
  params$b_tags         <- pb$tags

  if (cutoff <= 0) cutoff <- 0.01
  if (cutoff >= 1) cutoff <- 0.05
  params$cutoff           <- cutoff

  if (max_iterations < 1) max_iterations <- 1
  params$max_iterations <- max_iterations

  params$means_a <- data$means
  params$sda    <- data$sd
  params$means_b <- pb$means
  params$sdb    <- pb$sd
  params$yty    <- data$yty

  pa               <- list()
  pa$p1            <- params$p1
  pa$means         <- data$means
  pa$sd            <- data$sd
  pa$yty           <- data$yty
  pa$y_name         <- data$y_name
  pa$cutoff        <- params$cutoff
  pa$max_iterations <- params$max_iterations
  pa$a_col_names     <- params$a_col_names
  pa$tags          <- data$tags

  write_time <- proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_logistic_a2",
                       read_time, read_size,
                    write_time, write_size)

  return(params)
}


prepare_blocks_logistic_a2 <- function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_blocks_logistic_a2\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n <- params$n
  p1 <- params$p1
  p2 <- params$p2

  minimum_block_size <- get_block_size(p1, p2)
  if (n < minimum_block_size) {
    max_a_covariates <- trunc(sqrt(p2 * n) - p2 - 1)

    params$error_message <-
      paste("The minimum secure blocksize of", minimum_block_size,
            "is larger than the number of observations", paste0(n, ".\n"),
            "Your options are:\n",
            "Increase the number of observations to at least",
            paste0(minimum_block_size, ".\n"),
            "Decrease the number of A covariates to",
            max_a_covariates, "or less.")

    b <- n - 2 * p1 - 2
    discrim <- b^2 - 4 * (p1 + 1)^2
    if (discrim >= 0) {
      min_b_covariates <- trunc(1 + (b - sqrt(discrim)) / 2)
      max_b_covariates <- trunc((b + sqrt(discrim)) / 2)
      params$error_message <-
        paste0(params$error_message,
               "\nSet the number of B covariates to be between ",
               min_b_covariates, "and",
               paste0(max_b_covariates, "."))
    }
    warning(params$error_message)
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_blocks_logistic_a2", 0, 0, 0, 0)
    return(params)
  }

  if (is.null(blocksize)) {
    blocksize <- minimum_block_size
  }
  if (blocksize < minimum_block_size) {
    message(paste("Block size of", blocksize,
                  "is too small. Proceeding with minimum blocksize of",
                  paste0(minimum_block_size, ".")))
    blocksize <- minimum_block_size
  } else if (n < blocksize) {
    message(paste("Block size of", blocksize,
                  "is larger than size of data.  Proceeding with blocksize of",
                  paste0(n, ".")))
  }

  params$blocks    <- create_blocks(p1, p2, n, blocksize)
  params$container <- create_containers(p1, p2, params$blocks)
  write_time <- proc.time()[3]
  save(blocksize, file = file.path(params$write_path, "blocksize.rdata"))
  write_size <- file.size(file.path(params$write_path, "blocksize.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_blocks_logistic_a2",
                       0, 0, write_time, write_size)
  return(params)
}


get_z_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_z_logistic_a2\n\n")
  write_time <- 0
  write_size <- 0

  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "z", params$verbose)
  container_ct_z <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename <- paste0("cz_", container_ct_z, ".rdata")
      to_write <- file(file.path(params$write_path, filename), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1
    g <- params$blocks$g[i]
    z <- FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$x[strt:stp, ]), g)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(z), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "get_z_logistic_a2",
                       0, 0, write_time, write_size)
  return(params)
}


finalize_params_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "finalize_params_logistic_b2\n\n")
  pa <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # Load pa, a_col_names
  read_size <- sum(file.size(file.path(params$read_path, "pa.rdata")))
  read_time <- proc.time()[3] - read_time
  params$p1            <- pa$p1
  params$p             <- params$p1 + params$p2
  params$means_a        <- pa$means
  params$sda           <- pa$sd
  params$yty           <- pa$yty
  params$y_name         <- pa$y_name
  params$cutoff        <- pa$cutoff
  params$max_iterations <- pa$max_iterations
  params$a_col_names     <- pa$a_col_names
  params$a_tags         <- pa$tags
  params$b_tags         <- data$tags

  params <- add_to_log(params, "finalize_params_logistic_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


prepare_blocks_logistic_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_blocks_logistic_b2\n\n")
  blocksize <- NULL
  # For now, assuming that p1 > 0 and p2 > 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "blocksize.rdata")) # load blocksize
  read_size <- file.size(file.path(params$read_path, "blocksize.rdata"))
  read_time <- proc.time()[3] - read_time
  params$blocks    <- create_blocks(params$p1, params$p2, params$n, blocksize)
  params$container <- create_containers(params$p1, params$p2, params$blocks)
  params <- add_to_log(params, "prepare_blocks_logistic_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


get_w_logistics_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_w_logistics_b2\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks,
                              "(I-z*z')XB", params$verbose)

  xb_t_xb <- t(data$x) %*% data$x

  container_ct_z <- 0
  container_ct_w <- 0

  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    if (i %in% params$container$filebreak_w) {
      container_ct_w <- container_ct_w + 1
      filename2 <- paste0("cw_", container_ct_w, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1
    g1 <- params$blocks$g[i]

    read_time <- read_time - proc.time()[3]
    z <- matrix(readBin(con = to_read, what = numeric(), n = n2 * g1,
                       endian = "little"), nrow = n2, ncol = g1)
    read_time <- read_time + proc.time()[3]

    w <- data$x[strt:stp, ] - z %*% (t(z) %*% data$x[strt:stp, ])

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(w), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z ||
        i == params$blocks$num_blocks) {
      close(to_read)
    }
    if ((i + 1) %in% params$container$filebreak_w ||
        i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  write_time <- write_time - proc.time()[3]
  save(xb_t_xb, file = file.path(params$write_path, "xbtxb.rdata"))
  write_size <- write_size +
    file.size(file.path(params$write_path, "xbtxb.rdata"))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "get_w_logistics_b2",
                       read_time, read_size, write_time, write_size)

  return(params)
}


check_colinearity_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "check_colinearity_logistic_a2\n\n")
  p2 <- params$p2
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  xb_t_xb     <- NULL

  read_time <- read_time - proc.time()[3]
  load(file.path(params$read_path, "xbtxb.rdata")) # load xb_t_xb
  read_size <- file.size(file.path(params$read_path, "xbtxb.rdata"))
  read_time <- read_time + proc.time()[3]
  xa_t_xa <- t(data$x) %*% data$x
  xa_t_xb <- 0
  xa_t_y  <- t(data$x) %*% data$Y
  y_t_xb  <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks, "X'X", params$verbose)

  container_ct_w <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak_w) {
      container_ct_w <- container_ct_w + 1
      filename <- paste0("cw_", container_ct_w, ".rdata")
      to_read <- file(file.path(params$read_path, filename), "rb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1

    read_time <- read_time - proc.time()[3]
    w  <- matrix(readBin(con = to_read, what = numeric(), n = n2 * p2,
                       endian = "little"), nrow = n2, ncol = p2)
    read_time <- read_time + proc.time()[3]

    xa_t_xb <- xa_t_xb + t(data$x[strt:stp, ]) %*% w
    y_t_xb  <- y_t_xb  + t(data$Y[strt:stp, ]) %*% w

    if ((i + 1) %in% params$container$filebreak_w ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  xtx <- rbind(cbind(xa_t_xa, xa_t_xb), cbind(t(xa_t_xb), xb_t_xb))
  x_t_y <- rbind(xa_t_y, t(y_t_xb))

  nrow <- nrow(xtx)
  indicies <- c(1)
  for (i in 2:nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }

  xtx <- xtx[indicies, indicies]
  x_t_y  <- matrix(x_t_y[indicies], ncol = 1)

  a_names   <- params$a_col_names
  b_names   <- params$b_col_names
  a_index   <- which(indicies <= length(a_names))
  params$IndiciesKeep  <- indicies
  params$a_indicies_keep <- indicies[a_index]
  params$b_indicies_keep <- indicies[-a_index] - length(a_names)

  a_names_keep <- a_names[params$a_indicies_keep]
  b_names_keep <- b_names[params$b_indicies_keep]
  params$a_col_names_old <- params$a_col_names
  params$b_col_names_old <- params$b_col_names
  params$a_col_names     <- a_names_keep
  params$b_col_names     <- b_names_keep
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
  params$xty           <- x_t_y

  a_indicies <- params$a_indicies_keep
  b_indicies <- params$b_indicies_keep

  write_time <- write_time - proc.time()[3]
  save(a_indicies, b_indicies,
       file = file.path(params$write_path, "indicies.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "indicies.rdata")))
  write_time <- write_time + proc.time()[3]

  tags <- params$b_tags[params$b_indicies_keep]

  if (length(unique(tags)) < 2) {
    params$failed <- TRUE
    params$error_message <- paste("After removing colinear covariates,",
                                  "Party B has 1 or fewer covariates.\n\n")
  } else if (!("numeric" %in% names(tags))) {
    params$failed <- TRUE
    params$error_message <- paste("After removing colinear covariates,",
                                  "Party B has no continuous covariates.\n\n")
  }

  params <- add_to_log(params, "check_colinearity_logistic_a2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


update_data_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_logistic_a2\n\n")
  data$x <- as.matrix(data$x[, params$a_indicies_keep, drop = FALSE])
  return(data)
}


compute_initial_betas_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_initial_betas_logistic_a2\n\n")
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

  write_time <- proc.time()[3]
  save(b_betas, b_xty, file = file.path(params$write_path, "Bbetas_xty.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "Bbetas_xty.rdata")))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_initial_betas_logistic_a2",
                       0, 0, write_time, write_size)

  return(params)
}


update_params_logistic_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "update_params_logistic_b2\n\n")
  a_indicies <- NULL
  b_indicies <- NULL
  b_betas    <- NULL
  b_xty      <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "indicies.rdata"))
  load(file.path(params$read_path, "Bbetas_xty.rdata"))
  read_size <- sum(file.size(file.path(params$read_path,
                                       c("indicies.rdata",
                                         "Bbetas_xty.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$a_col_names_old <- params$a_col_names
  params$b_col_names_old <- params$b_col_names
  params$a_col_names     <- params$a_col_names_old[a_indicies]
  params$b_col_names     <- params$b_col_names_old[b_indicies]
  params$p1_old <- params$p1
  params$p2_old <- params$p2
  params$p1     <- length(a_indicies)
  params$p2     <- length(b_indicies)
  params$p_old  <- params$p
  params$p      <- params$p1 + params$p2
  params$b_indicies_keep <- b_indicies
  params$a_indicies_keep <- a_indicies
  params$betas_b    <- b_betas
  params$betas_b_old  <- matrix(0, params$p2, 1)
  params$means_b <- params$means_b[b_indicies]
  params$sdb    <- params$sdb[b_indicies]
  params$b_xty   <- b_xty
  params <- add_to_log(params, "update_params_logistic_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


update_data_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_logistic_b2\n\n")
  data$x <- as.matrix(data$x[, params$b_indicies_keep, drop = FALSE])
  return(data)
}


get_x_beta_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_x_beta_logistic_b2\n\n")
  x_beta_b <- data$x %*% params$betas_b

  write_time <- proc.time()[3]
  save(x_beta_b, file = file.path(params$write_path, "xbetab.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbetab.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_x_beta_logistic_b2",
                       0, 0, write_time, write_size)
  return(params)
}


get_weights_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_weights_logistic_a2\n\n")
  n      <- params$n
  p1     <- params$p1
  x_beta_b <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "xbetab.rdata"))  # Load x_beta_b
  read_size <- file.size(file.path(params$read_path, "xbetab.rdata"))
  read_time <- proc.time()[3] - read_time

  x_beta_a <- data$x %*% params$betas_a
  x_beta <- x_beta_a + x_beta_b
  pi_ <- (1 + exp(-x_beta))^(-1)
  params$pi_ <- pi_

  write_time <- proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi_.rdata"))
  write_size <- file.size(file.path(params$write_path, "pi_.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "get_weights_logistic_a2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_v_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_v_logistic_b2\n\n")
  pi_       <- NULL
  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pi_.rdata"))
  read_size <- file.size(file.path(params$read_path, "pi_.rdata"))
  read_time <- proc.time()[3] - read_time

  params$pi_ <- pi_
  w <- pi_ * (1 - params$pi_)

  xb_t_w_xb <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks,
                              "(I - z*z')w*XB", params$verbose)

  container_ct_z <- 0
  container_ct_v <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
    }
    if (i %in% params$container$filebreak_v) {
      container_ct_v <- container_ct_v + 1
      filename2 <- paste0("cv_", container_ct_v, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1
    g <- params$blocks$g[i]

    x_block  <- data$x[strt:stp, ]
    w_block  <- w[strt:stp]
    wx_block <- MultiplyDiagonalWTimesX(w_block, x_block)

    read_time <- read_time - proc.time()[3]
    z <- matrix(readBin(con = to_read, what = numeric(), n = n * g,
                       endian = "little"), nrow = n, ncol = g)
    read_time <- read_time + proc.time()[3]

    v <- wx_block - z %*% (t(z) %*% wx_block)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(v), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    xb_t_w_xb <- xb_t_w_xb + t(x_block) %*% wx_block
    if ((i + 1) %in% params$container$file_break_z ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    if ((i + 1) %in% params$container$filebreak_v ||
        i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }
  # sums of each column in WX_B
  sums_w_x_b <- apply(MultiplyDiagonalWTimesX(w, data$x), 2, sum)
  # This information needs to be shared in order to get the intercept term

  write_time <- write_time - proc.time()[3]
  save(sums_w_x_b, xb_t_w_xb,
       file = file.path(params$write_path, "sumswx_xbtwxb.rdata"))
  write_size <- write_size +
    sum(file.size(c(file.path(params$write_path, "sumswx_xbtwxb.rdata"))))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "get_v_logistic_b2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_ii_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_ii_logistic_a2\n\n")
  p1 <- params$p1
  p2 <- params$p2
  sums_w_x_b <- NULL
  xb_t_w_xb  <- NULL

  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "sumswx_xbtwxb.rdata"))
  read_size <- sum(file.size(file.path(params$read_path,
                                       "sumswx_xbtwxb.rdata")))
  read_time <- proc.time()[3] - read_time

  params$sums_w_x_b <- sums_w_x_b

  i_a <- params$a_xty - t(data$x) %*% params$pi_
  w <- params$pi_ * (1 - params$pi_)
  sums_w_xa <- apply(MultiplyDiagonalWTimesX(w, data$x), 2, sum)[-1]
  params$sums_w_xa <- sums_w_xa

  xa_t_w_xa <- t(data$x) %*% MultiplyDiagonalWTimesX(w, data$x)

  pbar <- make_progress_bar_1(params$blocks$num_blocks, "X'W*X", params$verbose)

  xa_t_w_xb <- 0
  container_ct_v <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak_v) {
      container_ct_v <- container_ct_v + 1
      filename1 <- paste0("cv_", container_ct_v, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    read_time <- read_time - proc.time()[3]
    v  <- matrix(readBin(con = to_read, what = numeric(),
                       n = n * p2, endian = "little"), n, p2)
    read_time <- read_time + proc.time()[3]
    xa_t_w_xb <- xa_t_w_xb + t(data$x[strt:stp, ]) %*% v

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
    if ((i + 1) %in% params$container$filebreak_v ||
        i == params$blocks$num_blocks) {
      close(to_read)
    }
  }

  x_t_w_x <- rbind(cbind(xa_t_w_xa, xa_t_w_xb), cbind(t(xa_t_w_xb), xb_t_w_xb))

  params$xtwx <- x_t_w_x

  ii <- NULL
  tryCatch({
    ii <- solve(params$xtwx)
  }, # dims are 1 + p1 + p2
  error = function(err) {
    ii <- NULL
  }
  )
  if (is.null(ii)) {
    params$singular_matrix <- TRUE
    params$failed <- TRUE
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
    warning(params$error_message)
  } else {
    params$ii <- ii
    params$i_a <- i_a

    a21i1 <- ii[(p1 + 1):(p1 + p2), 1:p1] %*% matrix(i_a, p1, 1)
    a11i1 <- ii[1:p1, 1:p1] %*% matrix(i_a, p1, 1)
    params$a11i1 <- a11i1

    write_time <- proc.time()[3]
    save(a21i1, x_t_w_x, file = file.path(params$write_path,
                                           "a21i1_xtwx.rdata"))
    write_size <- sum(file.size(file.path(params$write_path,
                                          "a21i1_xtwx.rdata")))
    write_time <- proc.time()[3] - write_time
  }
  params <- add_to_log(params, "get_ii_logistic_a2",
                       read_time, read_size, write_time, write_size)

  return(params)
}


get_coef_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_coef_logistic_b2\n\n")
  p1 <- params$p1
  p2 <- params$p2
  x_t_w_x  <- NULL
  a21i1 <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "a21i1_xtwx.rdata"))
  read_size <- sum(file.size(file.path(params$read_path, "a21i1_xtwx.rdata")))
  read_time <- proc.time()[3] - read_time

  i_b <- params$b_xty - t(data$x) %*% params$pi_

  ii <- solve(x_t_w_x)

  params$ii <- ii
  params$i_b <- i_b

  a22i2 <- ii[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2), drop = FALSE] %*% i_b
  a12i2 <- ii[1:p1, (p1 + 1):(p1 + p2), drop = FALSE] %*% i_b
  params$a22i2 <- a22i2

  params$betas_b_old <- params$betas_b
  params$betas_b <- params$betas_b + a21i1 + a22i2

  delta_beta_b <- max(abs(params$betas_b -
                            params$betas_b_old) / (abs(params$betas_b) + 0.1))

  write_time <- proc.time()[3]
  save(a12i2, delta_beta_b, file = file.path(params$write_path,
                                             "a12_deltabetaB.rdata"))
  write_size <- sum(file.size(file.path(params$write_path,
                                        "a12_deltabetaB.rdata")))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_coef_logistic_b2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_coef_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_coef_logistic_a2\n\n")
  a12i2      <- NULL
  delta_beta_b <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "a12_deltabetaB.rdata"))
  read_size <- sum(file.size(file.path(params$read_path,
                                       "a12_deltabetaB.rdata")))
  read_time <- proc.time()[3] - read_time

  params$betas_a_old <- params$beta_a
  params$betas_a <- params$beta_a + params$a11i1 + a12i2

  deltabeta <- max(abs(params$betas_a -
                         params$betas_a_old) / (abs(params$betas_a) + 0.1),
                   delta_beta_b)

  if (deltabeta < params$cutoff)  {
    params$converged <- TRUE
  } else if (params$alg_iteration_counter >= params$max_iterations) {
    params$max_iter_exceeded <- TRUE
    warning(paste("Failed to converged in",
                  params$max_iterations, "iterations."))
  }

  write_time <- proc.time()[3]
  save(deltabeta, file = file.path(params$write_path, "deltabeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "deltabeta.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "get_coef_logistic_a2",
                       read_time, read_size, write_time, write_size)


  return(params)
}


get_converged_status_logistic_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "GetconvergedStatusLogistic.b2\n\n")
  deltabeta <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "deltabeta.rdata"))
  read_size <- file.size(file.path(params$read_path, "deltabeta.rdata"))
  read_time <- proc.time()[3] - read_time

  if (deltabeta < params$cutoff)  {
    params$converged <- TRUE
  } else if (params$alg_iteration_counter >= params$max_iterations) {
    params$max_iter_exceeded <- TRUE
    warning(paste("Failed to converged in",
                  params$max_iterations, "iterations."))
  }

  params <- add_to_log(params, "get_converged_status_logistic_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


get_final_coef_logistic_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "get_final_coef_logistic_b2\n\n")
  betas_b <- params$betas_b / params$sdb
  offset_b <- sum(betas_b * params$means_b)
  b_final_fitted <- t(params$sdb * t(data$x) + params$means_b) %*% betas_b
  write_time <- proc.time()[3]
  save(betas_b, b_final_fitted, offset_b,
       file = file.path(params$write_path, "b_final.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "b_final.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "get_final_coef_logistic_b2",
                       0, 0, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
compute_results_logistic_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_results_logistic_a2\n\n")
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
  a_names <- params$a_col_names_old
  b_names <- params$b_col_names_old
  p1_old <- params$p1_old
  p2_old <- params$p2_old
  p_old  <- params$p_old
  indicies <- params$IndiciesKeep


  betas_b <- NULL
  offset_b <- NULL
  b_final_fitted <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "b_final.rdata"))
  read_size <- sum(file.size(file.path(params$read_path, "b_final.rdata")))
  read_time <- proc.time()[3] - read_time
  betas_a <- params$betas_a / sda
  offset_a <- sum(betas_a[-1] * params$means_a[-1])
  betas_a[1] <- betas_a[1] - offset_a - offset_b
  betas <- c(betas_a, betas_b)

  a_final_fitted <- t(sda * t(data$x) + means_a) %*% betas_a -
    t(sda[1] * t(data$x[, 1]) + means_a[1]) %*% betas_a[1]
  final_fitted <- a_final_fitted + b_final_fitted + betas[1]
  params$final_fitted <- final_fitted

  n <- params$n
  ct      <- sum(data$Y)
  resdev  <- -2 * (sum(data$Y * final_fitted) - sum(log(1 + exp(final_fitted))))
  nulldev <- -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

  # If xtwx were singular, it would have been caught in GetII.a2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  cov1 <- solve(params$xtwx)
  secoef <- sqrt(diag(cov1)) / c(sda, sdb)
  tmp <- matrix(c(1, (-means_a / sda)[-1], -means_b / sdb), ncol = 1)
  secoef[1] <- sqrt(t(tmp) %*% cov1 %*% tmp)


  stats$party <- c(rep("dp0", p1_old), rep("dp1", p2_old))
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
  stats$coefficients[indicies] <- betas
  stats$secoef[indicies] <- secoef
  tvals <- betas / secoef
  pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)
  stats$tvals[indicies] <- tvals
  stats$pvals[indicies] <- pvals

  stats$nulldev <- nulldev
  stats$resdev  <- resdev
  stats$hoslem  <- HoslemInternal(params, data)
  stats$ROC     <- roc_internal(params, data)
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

  stats$Y           <- data$Y # For Hoslem and ROC
  stats$final_fitted <- final_fitted
  params$stats      <- stats

  params <- add_to_log(params, "compute_results_logistic_b2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_results_logistic_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_logistic_b2\n\n")
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size <- file.size(file.path(params$read_path, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats <- stats
  params <- add_to_log(params, "get_results_logistic_b2",
                       read_time, read_size, 0, 0)
  return(params)
}



############################### PARENT FUNCTIONS ###############################


party_a_process_2_logistic <- function(data,
                                  y_name                 = NULL,
                                  monitor_folder         = NULL,
                                  msreqid               = "v_default_00_000",
                                  blocksize             = 500,
                                  cutoff                = 1e-8,
                                  max_iterations         = 25,
                                  sleep_time             = 10,
                                  max_waiting_time        = 24 * 60 * 60,
                                  popmednet             = TRUE,
                                  trace                 = FALSE,
                                  verbose               = TRUE) {
  params <- prepare_params_2p("logistic", "A", msreqid = msreqid,
                            popmednet = popmednet,
                            trace = trace, verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params   <- prepare_folder_logistic_a2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_logistic_a23(params, data, y_name)

  params <- PauseContinue.2p(params,  max_waiting_time)
  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed <- TRUE
    warning(read_error_message(params$read_path))
    params$pmn_step_counter <- 1
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
    params$completed <- TRUE
    message <- "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params$pmn_step_counter <- 1
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- prepare_params_logistic_a2(params, data, cutoff, max_iterations)

  if (params$failed) {   # Check for failed from prepare_params_logistic_a2()
    params$completed <- TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- prepare_blocks_logistic_a2(params, blocksize)

  if (params$failed) { # Check for failed from prepare_blocks_logistic_a2()
    params$completed <- TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- get_z_logistic_a2(params, data)

  files <- c("pa.rdata", "blocksize.rdata",
            seq_zw("cz_", length(params$container$file_break_z)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- check_colinearity_logistic_a2(params, data)

  if (params$failed) { # Check for check_colinearity_logistic_a2() failed
    params$completed <- TRUE
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }
  data <- update_data_logistic_a2(params, data)
  params <- add_to_log(params, "update_data_logistic_a2", 0, 0, 0, 0)
  params <- compute_initial_betas_logistic_a2(params, data)

  files <- c("indicies.rdata", "Bbetas_xty.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- get_weights_logistic_a2(params, data)
    files <- c("pi_.rdata")
    params <- send_pause_continue_2p(params, files,
                                     sleep_time, max_waiting_time)
    params <- get_ii_logistic_a2(params, data)

    if (params$failed) { # Check for failed from ComputeInverseLogistic.a2()
      params$completed <- TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    files <- c("a21i1_xtwx.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    params <- get_coef_logistic_a2(params, data)
    files <- "deltabeta.rdata"
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }
  params$completed <- TRUE

  params <- compute_results_logistic_a2(params, data)

  files <- c("stats.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)
  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  SummarizeLog.2p(params)
  return(invisible(params$stats))
}

party_b_process_2_logistic <- function(data,
                                  monitor_folder      = "v_default_00_000",
                                  sleep_time          = 10,
                                  max_waiting_time    = 24 * 60 * 60,
                                  popmednet           = TRUE,
                                  trace               = FALSE,
                                  verbose             = TRUE) {
  params <- prepare_params_2p("logistic", "B",
                            popmednet = popmednet, trace = trace,
                            verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params   <- prepare_folder_logistic_b2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_logistic_b23(params, data)

  if (data$failed) { # Check for Error from prepare_data_logistic_b2()
    params$completed <- TRUE
    message <- "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_2p(params, files, sleep_time = sleep_time,
                                 job_failed = TRUE)
    return(params$stats)
  }

  params   <- prepare_params_logistic_b2(params, data)

  files <- c("pb.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed <- TRUE
    warning(read_error_message(params$read_path))
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    return(params$stats)
  }

  params <- finalize_params_logistic_b2(params, data)
  params <- prepare_blocks_logistic_b2(params)
  params <- get_w_logistics_b2(params, data)

  files <- c("xbtxb.rdata", seq_zw("cw_", length(params$container$filebreak_w)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed <- TRUE
    warning(read_error_message(params$read_path))
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    return(params$stats)
  }

  params <- update_params_logistic_b2(params)
  data <- update_data_logistic_b2(params, data)
  params <- add_to_log(params, "update_data_logistic_b2", 0, 0, 0, 0)

  params$alg_iteration_counter <- 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- get_x_beta_logistic_b2(params, data)

    files <- c("xbetab.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    params <- get_v_logistic_b2(params, data)
    files <- c("sumswx_xbtwxb.rdata",
              seq_zw("cv_", length(params$container$filebreak_v)))
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$completed <- TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }

    params <- get_coef_logistic_b2(params, data)
    files <- c("a12_deltabetaB.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    params <- get_converged_status_logistic_b2(params)

    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }
  params$completed <- TRUE

  params <- get_final_coef_logistic_b2(params, data)
  files <- c("b_final.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- get_results_logistic_b2(params)
  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  return(invisible(params$stats))
}
