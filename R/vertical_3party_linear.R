################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

prepare_folder_linear_a3 <- function(params, monitor_folder = NULL) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_folder_linear_a3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  if (!is.character(monitor_folder)) {
    warning("monitor_folder directory is not valid. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "msoc")
  params$read_path       <- c(file.path(monitor_folder, "inputfiles"),
                              NA,
                              file.path(monitor_folder, "msoc2"))
  names(params$read_path) <- c("T", "A", "B")

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
                                  paste0(params$read_path[["T"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc2")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$read_path[["B"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$write_path, "."),
                                  "Check the path and restart the program.")
  }

  Sys.sleep(1)
  delete_trigger("files_done.ok", params$read_path[1])
  delete_trigger("files_done.ok", params$read_path[3])

  params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)

  return(params)
}


prepare_folder_linear_b3 <- function(params, monitor_folder = NULL) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_folder_linear_b3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_b3", 0, 0, 0, 0)
    return(params)
  }
  if (!is.character(monitor_folder)) {
    warning("monitor_folder directory is not valid. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_b3", 0, 0, 0, 0)
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "msoc")
  params$read_path       <- c(file.path(monitor_folder, "inputfiles"),
                              file.path(monitor_folder, "msoc1"),
                              NA)
  names(params$read_path) <- c("T", "A", "B")

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
                                  paste0(params$read_path[["T"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc1")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$read_path[["A"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$write_path, "."),
                                  "Check the path and restart the program.")
  }

  Sys.sleep(1)
  delete_trigger("files_done.ok", params$read_path[1])
  delete_trigger("files_done.ok", params$read_path[2])

  params <- add_to_log(params, "prepare_folder_linear_b3", 0, 0, 0, 0)

  return(params)
}


prepare_folder_linear_t3 <- function(params, monitor_folder = NULL) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_folder_linear_t3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  if (!is.character(monitor_folder)) {
    warning("monitor_folder directory is not valid. ",
            "Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  while (!dir.exists(monitor_folder)) {
    Sys.sleep(1)
  }

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "inputfiles")
  params$read_path       <- c(NA,
                              file.path(monitor_folder, "msoc1"),
                              file.path(monitor_folder, "msoc2"))
  names(params$read_path) <- c("T", "A", "B")

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
  if (!create_io_location(monitor_folder, "msoc1")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$read_path[["A"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc2")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$read_path[["B"]], "."),
                                  "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$write_path, "."),
                                  "Check the path and restart the program.")
  }

  params <- add_to_log(params, "prepare_folder_linear_t3", 0, 0, 0, 0)

  return(params)
}


prepare_params_linear_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_params_linear_a3\n\n")
  params$n     <- nrow(data$x)
  params$p     <- ncol(data$x)
  params$means <- data$means
  params$sd    <- data$sd

  pa          <- list()
  pa$analysis <- params$analysis
  pa$n        <- params$n
  pa$p        <- params$p
  pa$means    <- data$means
  pa$sd       <- data$sd
  pa$yty      <- data$yty
  pa$means_y   <- data$means_y
  pa$sdy       <- data$sdy
  pa$y_name    <- data$y_name
  pa$colnames  <- colnames(data$x)
  pa$tags      <- data$tags

  write_time <- proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size <- file.size(file.path(params$write_path, "pa.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_linear_a3",
                       0, 0, write_time, write_size)
  return(params)
}


prepare_params_linear_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_params_linear_b3\n\n")
  params$n     <- nrow(data$x)
  params$p     <- ncol(data$x)
  params$means <- data$means
  params$sd    <- data$sd

  pb           <- list()
  pb$analysis  <- params$analysis
  pb$n         <- params$n
  pb$p         <- params$p
  pb$means     <- data$means
  pb$sd        <- data$sd
  pb$colnames  <- colnames(data$x)
  pb$tags      <- data$tags

  write_time <- proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size <- file.size(file.path(params$write_path, "pb.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_linear_b3",
                       0, 0, write_time, write_size)
  return(params)
}


prepare_params_linear_t3 <- function(params, cutoff = 1e-8,
                                     max_iterations = 25) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_params_linear_t3\n\n")
  pa <- NULL
  pb <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "pa.rdata"))
  load(file.path(params$read_path[["B"]], "pb.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "pa.rdata")) +
    file.size(file.path(params$read_path[["B"]], "pb.rdata"))
  read_time <- proc.time()[3] - read_time
  if (length(table(c(pa$analysis, pb$analysis, params$analysis))) > 1) {
    params$failed <- TRUE
    params$error_message <- paste("Party A specified", pa$analysis,
                                  "regression, ",
                                  "Party B specified", pb$analysis,
                                  "regression, ",
                                  "and Party T specified", params$analysis,
                                  "regression. ")
  }
  if (pa$n != pb$n) {
    params$failed <- TRUE
    params$error_message <- paste0(params$error_message,
                                   paste("Party A has", pa$n,
                                         "observtions and Party B has", pb$n,
                                         "observations."))
  }
  params$analysis      <- pa$analysis
  params$n             <- pa$n
  params$p1            <- pa$p
  params$p2            <- pb$p
  params$p1_old        <- params$p1
  params$p2_old        <- params$p2
  params$p             <- pa$p + pb$p
  params$means_a       <- pa$means
  params$sda           <- pa$sd
  params$means_b       <- pb$means
  params$sdb           <- pb$sd
  params$means_y        <- pa$means_y
  params$sdy           <- pa$sdy
  params$yty           <- pa$yty
  params$colnamesA     <- pa$colnames
  params$colnamesB     <- pb$colnames
  params$a_tags         <- pa$tags
  params$b_tags         <- pb$tags
  params$y_name         <- pa$y_name
  params$cutoff        <- cutoff
  params$max_iterations <- max_iterations

  params <- add_to_log(params, "prepare_params_linear_t3",
                       read_time, read_size, 0, 0)
  return(params)
}


prepare_blocks_linear_t3 <- function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_blocks_linear_t3\n\n")
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
    params <- add_to_log(params, "prepare_blocks_linear_t3", 0, 0, 0, 0)
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
  blocks     <- params$blocks
  containers <- params$container
  write_time <- proc.time()[3]
  save(blocks, containers, file = file.path(params$write_path, "blocks.rdata"))
  write_size <- file.size(file.path(params$write_path, "blocks.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_blocks_linear_t3",
                       0, 0, write_time, write_size)
  return(params)
}


prepare_blocks_linear_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_blocks_linear_a3\n\n")
  blocks     <- NULL
  containers <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_time <- proc.time()[3] - read_time

  params$blocks <- blocks
  params$containers <- containers
  params <- add_to_log(params, "prepare_blocks_linear_a3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_z_linear_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_z_linear_a3\n\n")
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
  params <- add_to_log(params, "get_z_linear_a3", 0, 0, write_time, write_size)
  return(params)
}


process_z_linear_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "process_z_linear_t3\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "R(I-z*z')", params$verbose)
  container_ct_z <- 0
  container_ct_rz <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      to_read <- file(file.path(params$read_path[["A"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak_rz) {
      container_ct_rz <- container_ct_rz + 1
      filename2 <- paste0("crz_", container_ct_rz, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    filename3 <- paste0("r1_", i, ".rdata")

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1
    g <- params$blocks$g[i]

    read_time <- read_time - proc.time()[3]
    z <- matrix(readBin(con = to_read, what = numeric(), n = n * g,
                        endian = "little"), nrow = n, ncol = g)
    read_time <- read_time + proc.time()[3]
    r <- random_orthonormal_matrix(n)
    rz <- r - (r %*% z) %*% t(z)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(rz), con = to_write, endian = "little")
    to_write2 <- file(file.path(params$dp_local_path, filename3), "wb")
    writeBin(as.vector(r), con = to_write2, endian = "little")
    close(to_write2)
    write_size <- write_size +
      file.size(file.path(params$dp_local_path, filename3))
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["A"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak_rz || i == num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "process_z_linear_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


prepare_blocks_linear_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "prepare_blocks_linear_b3\n\n")
  blocks     <- NULL
  containers <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_time <- proc.time()[3] - read_time

  params$blocks <- blocks
  params$containers <- containers
  params <- add_to_log(params, "prepare_blocks_linear_b3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_rw_linear_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_rw_linear_b3\n\n")

  read_time <- 0
  read_size <- 0

  num_blocks <- params$blocks$num_blocks

  xb_t_xb <- t(data$x) %*% data$x
  write_time <- proc.time()[3]
  save(xb_t_xb, file = file.path(params$write_path, "xbtxb.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbtxb.rdata"))
  write_time <- proc.time()[3] - write_time

  pbar <- make_progress_bar_1(num_blocks, "R(I-z*z')XB", params$verbose)
  container_ct_rz <- 0
  container_ct_rw <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak_rz) {
      container_ct_rz <- container_ct_rz + 1
      filename1 <- paste0("crz_", container_ct_rz, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak_rw) {
      container_ct_rw <- container_ct_rw + 1
      filename2 <- paste0("crw_", container_ct_rw, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp  <- params$blocks$stops[i]
    n    <- stp - strt + 1
    g <- params$blocks$g[i]

    xb <- data$x[strt:stp, , drop = FALSE]
    read_time <- read_time - proc.time()[3]
    rz <- matrix(readBin(con = to_read, what = numeric(), n = n * n,
                         endian = "little"), nrow = n, ncol = n)
    read_time <- read_time + proc.time()[3]

    rw <- rz %*% xb

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(rw), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak_rz || i == num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak_rw || i == num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "get_rw_linear_b3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


process_w_linear_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "process_w_linear_t3\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  p2 <- params$p2

  write_time <- proc.time()[3]
  save(p2, file = file.path(params$write_path, "p2.rdata"))
  write_size <- file.size(file.path(params$write_path, "p2.rdata"))
  write_time <- proc.time()[3] - write_time

  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "(I-z*z')XB*R", params$verbose)

  container_ct_rw <- 0
  container_ct_wr <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak_rw) {
      container_ct_rw <- container_ct_rw + 1
      filename2 <- paste0("crw_", container_ct_rw, ".rdata")
      to_read_2 <- file(file.path(params$read_path[["B"]], filename2), "rb")
    }
    if (i %in% params$container$filebreak_wr) {
      container_ct_wr <- container_ct_wr + 1
      filename3 <- paste0("cwr_", container_ct_wr, ".rdata")
      to_write3 <- file(file.path(params$write_path, filename3), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    filename1 <- paste0("r1_", i, ".rdata")
    filename4 <- paste0("r2_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    to_read_1 <- file(file.path(params$dp_local_path, filename1), "rb")
    r1  <- matrix(readBin(con = to_read_1, what = numeric(), n = n * n,
                          endian = "little"), nrow = n, ncol = n)
    read_size <- read_size +
      file.size(file.path(params$dp_local_path, filename1))
    close(to_read_1)
    rw  <- matrix(readBin(con = to_read_2, what = numeric(), n = n * p2,
                          endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    w <- t(r1) %*% rw
    r2 <- random_orthonormal_matrix(p2)
    wr2 <- w %*% r2

    write_time <- write_time - proc.time()[3]
    to_write4 <- file(file.path(params$dp_local_path, filename4), "wb")
    writeBin(as.vector(r2), con = to_write4, endian = "little")
    close(to_write4)
    write_size <- write_size +
      file.size(file.path(params$dp_local_path, filename4))
    writeBin(as.vector(wr2), con = to_write3, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak_rw || i == num_blocks) {
      close(to_read_2)
      read_size <- read_size +
        file.size(file.path(params$read_path[["B"]], filename2))
    }
    if ((i + 1) %in% params$container$filebreak_wr || i == num_blocks) {
      close(to_write3)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename3))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "process_w_linear_t3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_wr_linear_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_wr_linear_a3\n\n")
  xa_t_xa <- t(data$x) %*% data$x
  xa_t_y  <- t(data$x) %*% data$Y
  write_time <- proc.time()[3]
  save(xa_t_xa, xa_t_y, file = file.path(params$write_path, "xatxa.rdata"))
  write_size <- file.size(file.path(params$write_path, "xatxa.rdata"))
  write_time <- proc.time()[3] - write_time

  p2 <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "p2.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "p2.rdata"))
  read_time <- proc.time()[3] - read_time

  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "XA'(I-z*z')XB*R", params$verbose)

  container_ct_wr <- 0
  container_ct_pr <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak_wr) {
      container_ct_wr <- container_ct_wr + 1
      filename1 <- paste0("cwr_", container_ct_wr, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak_pr) {
      container_ct_pr <- container_ct_pr + 1
      filename2 <- paste0("cpr_", container_ct_pr, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n <- stp - strt + 1

    read_time <- read_time - proc.time()[3]
    wr  <- matrix(readBin(con = to_read, what = numeric(), n = n * p2,
                          endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    y_xa <- cbind(data$Y[strt:stp, ], data$x[strt:stp, ])
    pr <- t(y_xa) %*% wr
    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(pr), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak_wr || i == num_blocks) {
      close(to_read)
      read_size <- read_size +
        file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak_pr || i == num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "get_wr_linear_a3",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_products_linear_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_products_linear_t3\n\n")
  n <- params$n
  p1 <- params$p1
  p2 <- params$p2
  xa_t_xa <- 0
  xb_t_xb <- 0
  xa_t_y  <- 0
  y_xa_t_xb <- 0

  num_blocks <- params$blocks$num_blocks
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["B"]], "xbtxb.rdata"))
  load(file.path(params$read_path[["A"]], "xatxa.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["B"]], "xbtxb.rdata")),
                   file.size(file.path(params$read_path[["A"]], "xatxa.rdata")))
  read_time <- proc.time()[3] - read_time

  pbar <- make_progress_bar_1(num_blocks, "X'X", params$verbose)

  container_ct_pr <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak_pr) {
      container_ct_pr <- container_ct_pr + 1
      filename1 <- paste0("cpr_", container_ct_pr, ".rdata")
      to_read <- file(file.path(params$read_path[["A"]], filename1), "rb")
      read_size <- read_size +
        file.size(file.path(params$read_path[["A"]], filename1))
    }

    filename1 <- paste0("r2_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    to_read_1 <- file(file.path(params$dp_local_path, filename1), "rb")
    r2  <- matrix(readBin(con = to_read_1, what = numeric(), n = p2 * p2,
                          endian = "little"), p2, p2)
    read_size <- read_size +
      file.size(file.path(params$dp_local_path, filename1))
    close(to_read_1)
    pr  <- matrix(readBin(con = to_read, what = numeric(), n = (p1 + 1) * p2,
                          endian = "little"), p1 + 1, p2)
    read_time <- read_time + proc.time()[3]

    y_xa_t_xb <- y_xa_t_xb + pr %*% t(r2)

    if ((i + 1) %in% params$container$filebreak_pr || i == num_blocks) {
      close(to_read)
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  y_t_xb <- y_xa_t_xb[1, , drop = FALSE]
  xa_t_xb <- y_xa_t_xb[-1, , drop = FALSE]
  xtx <- rbind(cbind(xa_t_xa, xa_t_xb), cbind(t(xa_t_xb), xb_t_xb))
  x_t_y <- rbind(xa_t_y, t(y_t_xb))

  x_t_x_lasso <- xtx / (n - 1)
  x_t_y_lasso <- params$sdy * x_t_y / sqrt(n - 1)

  params$xtx <- xtx
  params$xty <- x_t_y
  params$xtxLasso <- x_t_x_lasso
  params$xtyLasso <- x_t_y_lasso

  params$converged <- TRUE

  params <- add_to_log(params, "get_products_linear_t3",
                       read_time, read_size, 0, 0)
  return(params)
}


#' @importFrom  stats pf pt
comp_results_linear_t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()),
                        "comp_results_linear_t3\n\n")
  stats    <- params$stats
  stats$converged <- params$converged
  stats$failed    <- FALSE
  a_names   <- params$colnamesA
  b_names   <- params$colnamesB
  n        <- params$n
  yty      <- params$yty
  xty      <- params$xty
  xtx      <- params$xtx
  sdy      <- params$sdy
  sda      <- params$sda
  sdb      <- params$sdb
  means_y   <- params$means_y
  means_a   <- params$means_a
  means_b   <- params$means_b

  # First we de-standardize.
  xtx <- diag(c(sda, sdb)) %*% xtx %*% diag(c(sda, sdb))
  offset  <- matrix(c(means_a, means_b), ncol = 1) %*%
    matrix(c(means_a, means_b), nrow = 1) * n
  offset[1, 1] <- 0
  xtx <- xtx + offset

  xty <- diag(c(sda, sdb)) %*% xty * sdy
  offset <- n * means_y * matrix(c(means_a, means_b), ncol = 1)
  xty <- xty + offset

  # Now, we check for colinearity
  nrow <- nrow(xtx)
  indicies <- c(1)
  for (i in 2:nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }


  a_index        <- which(indicies <= length(a_names))
  a_indicies_keep <- indicies[a_index]
  b_indicies_keep <- indicies[-a_index] - length(a_names)
  names_old     <- c(a_names, b_names)
  p             <- length(indicies)
  xtx_old       <- xtx
  xty_old       <- xty
  xtx           <- xtx[indicies, indicies, drop = FALSE]
  xty            <- matrix(xty[indicies, ], ncol = 1)

  invxtx <- solve(xtx)
  betas  <- drop(invxtx %*% xty)

  num_covariates <- p - 1

  #   # If true sse is approximately 0, random variations could cause this
  #   # calculation to be less than 0
  #   # If calculated sse is less than 0, we set it equal to 0.
  sse     <- max(drop(yty - 2 * t(xty) %*% betas +
                        (t(betas) %*% xtx) %*% betas), 0)
  rstderr <- drop(sqrt(sse / (n - num_covariates - 1)))
  sst     <- drop(yty - means_y^2 * n)
  ssr     <- sst - sse
  df1     <- num_covariates
  df2     <- n - num_covariates - 1
  if (sse == 0) {
    Fstat <- Inf
  } else {
    Fstat <- (ssr / df1) / (sse / df2)
  }
  Fpval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  if (sse == 0) {
    r_sq <- 1
  } else {
    r_sq <- drop(1 - sse / sst)
  }
  adj_r_sq <- drop(1 - (n - 1) / (n - num_covariates - 1) * (1 - r_sq))
  if (rstderr == 0) {
    tvals <- rep(Inf, num_covariates + 1)
  } else {
    tvals   <- betas / (rstderr * sqrt(diag(invxtx)))
  }

  secoef  <- tvals^-1 * betas
  pvals   <- 2 * pt(abs(tvals), n - num_covariates - 1, lower.tail = FALSE)
  stats$party                  <- c(rep("dp1", params$p1_old),
                                    rep("dp2", params$p2_old))
  stats$responseParty          <- "dp1"
  stats$coefficients           <- rep(NA, params$p)
  stats$tvals                  <- rep(NA, params$p)
  stats$secoef                 <- rep(NA, params$p)
  stats$pvals                  <- rep(NA, params$p)

  stats$sse                    <- sse
  stats$coefficients[indicies] <- betas
  stats$tvals[indicies]        <- tvals
  stats$secoef[indicies]       <- secoef
  stats$pvals[indicies]        <- pvals
  stats$rstderr                <- rstderr
  stats$rsquare                <- r_sq
  stats$adjrsquare             <- adj_r_sq
  stats$Fstat                  <- Fstat
  stats$Fpval <- Fpval
  stats$df1                    <- df1
  stats$df2                    <- df2
  stats$n                      <- params$n
  stats$xtx                    <- xtx_old
  stats$xty                    <- xty_old
  stats$yty                    <- yty
  stats$means_y                 <- means_y
  stats$means                  <- c(means_a, means_b)

  names(stats$party)           <- names_old
  names(stats$coefficients)    <- names_old
  names(stats$secoef)          <- names_old
  names(stats$tvals)           <- names_old
  names(stats$pvals)           <- names_old

  colnames(stats$xtx)          <- names_old
  rownames(stats$xtx)          <- names_old
  colnames(stats$xty)          <- colnames(params$xty)
  rownames(stats$xty)          <- names_old

  params$stats <- stats

  tags <- params$b_tags[b_indicies_keep]

  if (length(unique(tags)) < 2) {
    params$failed <- TRUE
    params$error_message <-
      paste("After removing colinear covariates,",
            "Party B has 1 or fewer covariates.")
  } else if (!("numeric" %in% names(tags))) {
    params$failed <- TRUE
    params$error_message <-
      paste("After removing colinear covariates,",
            "Party B has no continuous covariates.")
  }

  stats$failed <- params$failed

  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "comp_results_linear_t3",
                       0, 0, write_time, write_size)
  return(params)
}


get_results_linear_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_linear_a3\n\n")
  params$converged <- TRUE
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats <- stats

  params <- add_to_log(params, "get_results_linear_a3",
                       read_time, read_size, 0, 0)
  return(params)
}


get_results_linear_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_linear_b3\n\n")
  params$converged <- TRUE
  stats <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats <- stats

  params <- add_to_log(params, "get_results_linear_b3",
                       read_time, read_size, 0, 0)
  return(params)
}


############################## PARENT FUNCTIONS ###############################


party_a_process_3_linear <- function(data,
                                     y_name          = NULL,
                                     monitor_folder  = NULL,
                                     sleep_time      = 10,
                                     max_waiting_time = 24 * 60 * 60,
                                     popmednet      = TRUE,
                                     trace          = FALSE,
                                     verbose        = TRUE) {

  params <- prepare_params_3p("linear", "A",
                              popmednet = popmednet,
                              trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)
  header(params)

  params   <- prepare_folder_linear_a3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_linear_a23(params, data, y_name)
  params <- add_to_log(params, "prepare_data_linear_a23", 0, 0, 0, 0)

  if (data$failed) {
    message <- "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params, files_t = files,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_a3(params, data)
  files <- "pa.rdata"
  params <- send_pause_continue_3p(params, files_t = files, from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_blocks_linear_a3(params)
  params <- get_z_linear_a3(params, data)
  files <- seq_zw("cz_", length(params$container$file_break_z))
  params <- send_pause_continue_3p(params, files_t = files, from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_wr_linear_a3(params, data)
  files <- c("xatxa.rdata",
             seq_zw("cpr_", length(params$container$filebreak_pr)))
  params <- send_pause_continue_3p(params, files_t = files, from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- get_results_linear_a3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                               wait_for_turn = TRUE)
  return(params$stats)
}


party_b_process_3_linear <- function(data,
                                     monitor_folder  = NULL,
                                     sleep_time      = 10,
                                     max_waiting_time = 24 * 60 * 60,
                                     popmednet      = TRUE,
                                     trace          = FALSE,
                                     verbose        = TRUE) {
  params <- prepare_params_3p("linear", "B",
                              popmednet = popmednet, trace = trace,
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

  data <- prepare_data_linear_b23(params, data)
  params <- add_to_log(params, "prepare_data_linear_b23", 0, 0, 0, 0)

  if (data$failed) {
    message <- "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params, files_t = files,
                                 sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_b3(params, data)
  files <- "pb.rdata"
  params <- send_pause_continue_3p(params, files_t = files, from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time,
                                   wait_for_turn = TRUE)


  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- prepare_blocks_linear_b3(params)
  params <- get_rw_linear_b3(params, data)
  files <- c("xbtxb.rdata",
             seq_zw("crw_", length(params$container$filebreak_rw)))
  params <- send_pause_continue_3p(params, files_t = files, from = "T",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE, wait_for_turn = TRUE)
    return(params$stats)
  }

  params <- get_results_linear_b3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                               wait_for_turn = TRUE)
  return(params$stats)
}


party_t_process_3_linear <- function(monitor_folder = NULL,
                                     msreqid             = "v_default_0_000",
                                     blocksize           = 500,
                                     sleep_time          = 10,
                                     max_waiting_time    = 24 * 60 * 60,
                                     popmednet           = TRUE,
                                     trace               = FALSE,
                                     verbose             = TRUE) {
  params <- prepare_params_3p("linear", "T", msreqid = msreqid,
                              popmednet = popmednet, trace = trace,
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
    warning(paste(read_error_message(params$read_path[["A"]]), "\n",
                  read_error_message(params$read_path[["B"]])))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["A"]]))
    file.copy(file.path(params$read_path[["A"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, files_b = files, from = "B",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["B"]]))
    file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, files_a = files, from = "A",
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params   <- prepare_params_linear_t3(params)
  if (!params$failed) params <- prepare_blocks_linear_t3(params, blocksize)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, files_a = files, files_b = files,
                                     from = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files <- "blocks.rdata"
  params <- send_pause_continue_3p(params, files_a = files, from = "A",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- process_z_linear_t3(params)
  files <- c("blocks.rdata",
             seq_zw("crz_", length(params$container$filebreak_rz)))
  params <- send_pause_continue_3p(params, files_b = files, from  = "B",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- process_w_linear_t3(params)
  files <- c("p2.rdata", seq_zw("cwr_", length(params$container$filebreak_wr)))
  params <- send_pause_continue_3p(params, files_a = files, from  = "A",
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- get_products_linear_t3(params)
  params <- comp_results_linear_t3(params)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, files_a = files, files_b = files,
                                     from = c("A", "B"),
                                     sleep_time = sleep_time,
                                     max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files <- "stats.rdata"
  params <- send_pause_continue_3p(params, files_a = files,
                                   files_b = files, from  = c("A", "B"),
                                   sleep_time = sleep_time,
                                   max_waiting_time = max_waiting_time)

  params <- send_pause_quit_3p(params, sleep_time = sleep_time)
  SummarizeLog.3p(params)
  return(params$stats)
}
