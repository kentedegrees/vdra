################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

prepare_folder_linear_a3 <- function(params, monitor_folder = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_folder_linear_a3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
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
  if (params$trace) cat(as.character(Sys.time()), "prepare_folder_linear_b3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_b3", 0, 0, 0, 0)
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
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
  if (params$trace) cat(as.character(Sys.time()), "prepare_folder_linear_t3\n\n")
  if (is.null(monitor_folder)) {
    warning("monitor_folder must be specified.  Please use the same monitor_folder as the DataMart Client.")
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_folder_linear_a3", 0, 0, 0, 0)
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning("monitor_folder directory is not valid.  Please use the same monitor_folder as the DataMart Client.")
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
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_linear_a3\n\n")
  params$n     = nrow(data$X)
  params$p     = ncol(data$X)
  params$means = data$means
  params$sd    = data$sd

  pa          = list()
  pa$analysis = params$analysis
  pa$n        = params$n
  pa$p        = params$p
  pa$means    = data$means
  pa$sd       = data$sd
  pa$yty      = data$yty
  pa$meansy   = data$meansy
  pa$sdy      = data$sdy
  pa$y_name    = data$y_name
  pa$colnames = colnames(data$X)
  pa$tags     = data$tags

  write_time <- proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size <- file.size(file.path(params$write_path, "pa.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_linear_a3", 0, 0, write_time, write_size)
  return(params)
}


prepare_params_linear_b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_linear_b3\n\n")
  params$n     = nrow(data$X)
  params$p     = ncol(data$X)
  params$means = data$means
  params$sd    = data$sd

  pb           = list()
  pb$analysis  = params$analysis
  pb$n         = params$n
  pb$p         = params$p
  pb$means     = data$means
  pb$sd        = data$sd
  pb$colnames  = colnames(data$X)
  pb$tags      = data$tags

  write_time <- proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size <- file.size(file.path(params$write_path, "pb.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_linear_b3", 0, 0, write_time, write_size)
  return(params)
}


prepare_params_linear_t3 <- function(params, cutoff = 1e-8, max_iterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_linear_t3\n\n")
  pa = NULL
  pb = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "pa.rdata"))
  load(file.path(params$read_path[["B"]], "pb.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "pa.rdata")) +
    file.size(file.path(params$read_path[["B"]], "pb.rdata"))
  read_time <- proc.time()[3] - read_time
  if (length(table(c(pa$analysis, pb$analysis, params$analysis))) > 1) {
    params$failed <- TRUE
    params$error_message <- paste("Party A specified", pa$analysis, "regression, ",
                                "Party B specified", pb$analysis, "regression, ",
                                "and Party T specified", params$analysis, "regression. ")
  }
  if (pa$n != pb$n) {
    params$failed <- TRUE
    params$error_message <- paste0(params$error_message,
                                 paste("Party A has", pa$n,
                                       "observtions and Party B has", pb$n,
                                       "observations."))
  }
  params$analysis      = pa$analysis
  params$n             = pa$n
  params$p1            = pa$p
  params$p2            = pb$p
  params$p1.old        = params$p1
  params$p2.old        = params$p2
  params$p             = pa$p + pb$p
  params$meansA        = pa$means
  params$sdA           = pa$sd
  params$meansB        = pb$means
  params$sdB           = pb$sd
  params$meansy        = pa$meansy
  params$sdy           = pa$sdy
  params$yty           = pa$yty
  params$colnamesA     = pa$colnames
  params$colnamesB     = pb$colnames
  params$Atags         = pa$tags
  params$Btags         = pb$tags
  params$y_name         = pa$y_name
  params$cutoff        = cutoff
  params$max_iterations = max_iterations

  params <- add_to_log(params, "prepare_params_linear_t3", read_time, read_size, 0, 0)
  return(params)
}


prepare_blocks_linear_t3 <- function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_blocks_linear_t3\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  = params$n
  p1 = params$p1
  p2 = params$p2


  minimum_block_size = get_block_size(p1, p2)
  if (n < minimum_block_size) {
    max_a_covariates = trunc(sqrt(p2 * n) - p2 - 1)

    params$error_message <-
      paste("The minimum secure blocksize of", minimum_block_size,
            "is larger than the number of observations", paste0(n, ".\n"),
            "Your options are:\n",
            "Increase the number of observations to at least",
            paste0(minimum_block_size, ".\n"),
            "Decrease the number of A covariates to", max_a_covariates, "or less.")

    b = n - 2 * p1 - 2
    discrim = b^2 - 4 * (p1 + 1)^2
    if (discrim >= 0) {
      min_b_covariates = trunc(1 + (b - sqrt(discrim)) / 2)
      max_b_covariates = trunc((b + sqrt(discrim)) / 2)
      params$error_message <-
        paste0(params$error_message,
               "\nSet the number of B covariates to be between ", min_b_covariates, "and",
               paste0(max_b_covariates, "."))
    }
    warning(params$error_message)
    params$failed <- TRUE
    params <- add_to_log(params, "prepare_blocks_linear_t3", 0, 0, 0, 0)
    return(params)
  }

  if (is.null(blocksize)) {
    blocksize = minimum_block_size
  }
  if (blocksize < minimum_block_size) {
    message(paste("Block size of", blocksize,
                  "is too small. Proceeding with minimum blocksize of",
                  paste0(minimum_block_size, ".")))
    blocksize = minimum_block_size
  } else if (n < blocksize) {
    message(paste("Block size of", blocksize,
                  "is larger than size of data.  Proceeding with blocksize of",
                  paste0(n, ".")))
  }

  params$blocks    = CreateBlocks(p1, p2, n, blocksize)
  params$container = CreateContainers(p1, p2, params$blocks)
  blocks     = params$blocks
  containers = params$container
  write_time <- proc.time()[3]
  save(blocks, containers, file = file.path(params$write_path, "blocks.rdata"))
  write_size <- file.size(file.path(params$write_path, "blocks.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_blocks_linear_t3", 0, 0, write_time, write_size)
  return(params)
}


prepare_blocks_linear_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_blocks_linear_a3\n\n")
  blocks     = NULL
  containers = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_time <- proc.time()[3] - read_time

  params$blocks = blocks
  params$containers = containers
  params <- add_to_log(params, "prepare_blocks_linear_a3", read_time, read_size, 0, 0)
  return(params)
}


get_z_linear_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_z_linear_a3\n\n")
  write_time <- 0
  write_size <- 0

  num_blocks = params$blocks$num_blocks
  pbar = MakeProgressBar1(num_blocks, "Z", params$verbose)
  container_ct_z <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename = paste0("cz_", container_ct_z, ".rdata")
      to_write <- file(file.path(params$write_path, filename), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]
    Z = FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$X[strt:stp, ]), g)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(Z), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename))
    }
    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "get_z_linear_a3", 0, 0, write_time, write_size)
  return(params)
}


ProcessZLinear.t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessZLinear.t3\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  num_blocks = params$blocks$num_blocks
  pbar = MakeProgressBar1(num_blocks, "R(I-Z*Z')", params$verbose)
  container_ct_z <- 0
  container_ct_RZ = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      toRead = file(file.path(params$read_path[["A"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.RZ) {
      container_ct_RZ = container_ct_RZ + 1
      filename2 <- paste0("crz_", container_ct_RZ, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    filename3 = paste0("r1_", i, ".rdata")

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]

    read_time <- read_time - proc.time()[3]
    Z = matrix(readBin(con = toRead, what = numeric(), n = n * g,
                       endian = "little"), nrow = n, ncol = g)
    read_time <- read_time + proc.time()[3]
    R = RandomOrthonomalMatrix(n)
    RZ = R - (R %*% Z) %*% t(Z)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(RZ), con = to_write, endian = "little")
    to_write2 = file(file.path(params$dp_local_path, filename3), "wb")
    writeBin(as.vector(R), con = to_write2, endian = "little")
    close(to_write2)
    write_size <- write_size + file.size(file.path(params$dp_local_path, filename3))
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(toRead)
      read_size <- read_size + file.size(file.path(params$read_path[["A"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.RZ || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "ProcessZLinear.t3", read_time, read_size, write_time, write_size)
  return(params)
}


prepare_blocks_linear_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_blocks_linear_b3\n\n")
  blocks     = NULL
  containers = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "blocks.rdata"))
  read_time <- proc.time()[3] - read_time

  params$blocks = blocks
  params$containers = containers
  params <- add_to_log(params, "prepare_blocks_linear_b3", read_time, read_size, 0, 0)
  return(params)
}


GetRWLinear.b3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetRWLinear.b3\n\n")

  read_time <- 0
  read_size <- 0

  num_blocks = params$blocks$num_blocks

  xb_t_xb <- t(data$X) %*% data$X
  write_time <- proc.time()[3]
  save(xb_t_xb, file = file.path(params$write_path, "xbtxb.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbtxb.rdata"))
  write_time <- proc.time()[3] - write_time

  pbar = MakeProgressBar1(num_blocks, "R(I-Z*Z')XB", params$verbose)
  container_ct_RZ = 0
  container_ct_RW = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.RZ) {
      container_ct_RZ = container_ct_RZ + 1
      filename1 <- paste0("crz_", container_ct_RZ, ".rdata")
      toRead = file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.RW) {
      container_ct_RW = container_ct_RW + 1
      filename2 <- paste0("crw_", container_ct_RW, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp  = params$blocks$stops[i]
    n    = stp - strt + 1
    g = params$blocks$g[i]

    XB = data$X[strt:stp, , drop = FALSE]
    read_time <- read_time - proc.time()[3]
    RZ = matrix(readBin(con = toRead, what = numeric(), n = n * n,
                        endian = "little"), nrow = n, ncol = n)
    read_time <- read_time + proc.time()[3]

    RW = RZ %*% XB

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(RW), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.RZ || i == num_blocks) {
      close(toRead)
      read_size <- read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.RW || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "GetRWLinear.b3", read_time, read_size, write_time, write_size)
  return(params)
}


ProcessWLinear.t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessWLinear.t3\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  p2 = params$p2

  write_time <- proc.time()[3]
  save(p2, file = file.path(params$write_path, "p2.rdata"))
  write_size <- file.size(file.path(params$write_path, "p2.rdata"))
  write_time <- proc.time()[3] - write_time

  num_blocks = params$blocks$num_blocks
  pbar = MakeProgressBar1(num_blocks, "(I-Z*Z')XB*R", params$verbose)

  container_ct_RW = 0
  container_ct_wr = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.RW) {
      container_ct_RW = container_ct_RW + 1
      filename2 <- paste0("crw_", container_ct_RW, ".rdata")
      toRead2 = file(file.path(params$read_path[["B"]], filename2), "rb")
    }
    if (i %in% params$container$filebreak.WR) {
      container_ct_wr = container_ct_wr + 1
      filename3 = paste0("cwr_", container_ct_wr, ".rdata")
      to_write3 = file(file.path(params$write_path, filename3), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1

    filename1 <- paste0("r1_", i, ".rdata")
    filename4 = paste0("r2_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    toRead1 = file(file.path(params$dp_local_path, filename1), "rb")
    R1 = matrix(readBin(con = toRead1, what = numeric(), n = n * n,
                        endian = "little"), nrow = n, ncol = n)
    read_size <- read_size + file.size(file.path(params$dp_local_path, filename1))
    close(toRead1)
    RW = matrix(readBin(con = toRead2, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    W = t(R1) %*% RW
    R2 = RandomOrthonomalMatrix(p2)
    WR2 = W %*% R2

    write_time <- write_time - proc.time()[3]
    to_write4 = file(file.path(params$dp_local_path, filename4), "wb")
    writeBin(as.vector(R2), con = to_write4, endian = "little")
    close(to_write4)
    write_size <- write_size + file.size(file.path(params$dp_local_path, filename4))
    writeBin(as.vector(WR2), con = to_write3, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.RW || i == num_blocks) {
      close(toRead2)
      read_size <- read_size + file.size(file.path(params$read_path[["B"]], filename2))
    }
    if ((i + 1) %in% params$container$filebreak.WR || i == num_blocks) {
      close(to_write3)
      write_size <- write_size + file.size(file.path(params$write_path, filename3))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "ProcessWLinear.t3", read_time, read_size, write_time, write_size)
  return(params)
}


get_wr_linear_a3 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_wr_linear_a3\n\n")
  XATXA = t(data$X) %*% data$X
  XATY  = t(data$X) %*% data$Y
  write_time <- proc.time()[3]
  save(XATXA, XATY, file = file.path(params$write_path, "xatxa.rdata"))
  write_size <- file.size(file.path(params$write_path, "xatxa.rdata"))
  write_time <- proc.time()[3] - write_time

  p2 = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "p2.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "p2.rdata"))
  read_time <- proc.time()[3] - read_time

  num_blocks = params$blocks$num_blocks
  pbar = MakeProgressBar1(num_blocks, "XA'(I-Z*Z')XB*R", params$verbose)

  container_ct_wr = 0
  container_ct_PR = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.WR) {
      container_ct_wr = container_ct_wr + 1
      filename1 <- paste0("cwr_", container_ct_wr, ".rdata")
      toRead = file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.PR) {
      container_ct_PR = container_ct_PR + 1
      filename2 <- paste0("cpr_", container_ct_PR, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1

    read_time <- read_time - proc.time()[3]
    WR = matrix(readBin(con = toRead, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    YXA = cbind(data$Y[strt:stp, ], data$X[strt:stp, ])
    PR = t(YXA) %*% WR
    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(PR), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.WR || i == num_blocks) {
      close(toRead)
      read_size <- read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.PR || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "get_wr_linear_a3", read_time, read_size, write_time, write_size)
  return(params)
}


GetProductsLinear.t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLinear.t3\n\n")
  n  = params$n
  p1 = params$p1
  p2 = params$p2
  XATXA = 0
  xb_t_xb <- 0
  XATY  = 0
  YXATXB = 0

  num_blocks = params$blocks$num_blocks
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["B"]], "xbtxb.rdata"))
  load(file.path(params$read_path[["A"]], "xatxa.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["B"]], "xbtxb.rdata")),
                 file.size(file.path(params$read_path[["A"]], "xatxa.rdata")))
  read_time <- proc.time()[3] - read_time

  pbar = MakeProgressBar1(num_blocks, "X'X", params$verbose)

  container_ct_PR = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.PR) {
      container_ct_PR = container_ct_PR + 1
      filename1 <- paste0("cpr_", container_ct_PR, ".rdata")
      toRead = file(file.path(params$read_path[["A"]], filename1), "rb")
      read_size <- read_size + file.size(file.path(params$read_path[["A"]], filename1))
    }

    filename1 <- paste0("r2_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    toRead1 = file(file.path(params$dp_local_path, filename1), "rb")
    R2 = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
                        endian = "little"), p2, p2)
    read_size <- read_size + file.size(file.path(params$dp_local_path, filename1))
    close(toRead1)
    PR = matrix(readBin(con = toRead, what = numeric(), n = (p1 + 1) * p2,
                        endian = "little"), p1 + 1, p2)
    read_time <- read_time + proc.time()[3]

    YXATXB = YXATXB + PR %*% t(R2)

    if ((i + 1) %in% params$container$filebreak.PR || i == num_blocks) {
      close(toRead)
    }

    pbar = MakeProgressBar2(i, pbar, params$verbose)
  }

  YTXB = YXATXB[1, , drop = FALSE]
  XATXB = YXATXB[-1, , drop = FALSE]
  xtx = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), xb_t_xb))
  XTY = rbind(XATY, t(YTXB))

  XTXLasso = xtx / (n - 1)
  XTYLasso = params$sdy * XTY / sqrt(n - 1)

  params$xtx = xtx
  params$xty = XTY
  params$xtxLasso = XTXLasso
  params$xtyLasso = XTYLasso

  params$converged = TRUE

  params <- add_to_log(params, "GetProductsLinear.t3", read_time, read_size, 0, 0)
  return(params)
}


#' @importFrom  stats pf pt
ComputeResultsLinear.t3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLinear.t3\n\n")
  stats    = params$stats
  stats$converged = params$converged
  stats$failed    = FALSE
  a_names   = params$colnamesA
  b_names   = params$colnamesB
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


  a_index        = which(indicies <= length(a_names))
  a_indicies_keep = indicies[a_index]
  b_indicies_keep = indicies[-a_index] - length(a_names)
  names_old     = c(a_names, b_names)
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
  stats$party                  = c(rep("dp1", params$p1.old),
                                   rep("dp2", params$p2.old))
  stats$responseParty          = "dp1"
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

  names(stats$party)           = names_old
  names(stats$coefficients)    = names_old
  names(stats$secoef)          = names_old
  names(stats$tvals)           = names_old
  names(stats$pvals)           = names_old

  colnames(stats$xtx)          = names_old
  rownames(stats$xtx)          = names_old
  colnames(stats$xty)          = colnames(params$xty)
  rownames(stats$xty)          = names_old

  params$stats = stats

  tags = params$Btags[b_indicies_keep]

  if (length(unique(tags)) < 2) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has 1 or fewer covariates."
  } else if (!("numeric" %in% names(tags))) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has no continuous covariates."
  }

  stats$failed = params$failed

  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeResultsLinear.t3", 0, 0, write_time, write_size)
  return(params)
}


get_results_linear_a3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_linear_a3\n\n")
  params$converged = TRUE
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats

  params <- add_to_log(params, "get_results_linear_a3", read_time, read_size, 0, 0)
  return(params)
}


get_results_linear_b3 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_linear_b3\n\n")
  params$converged = TRUE
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats

  params <- add_to_log(params, "get_results_linear_b3", read_time, read_size, 0, 0)
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
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)
  header(params)

  params   = prepare_folder_linear_a3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = prepare_data_linear_a23(params, data, y_name)
  params <- add_to_log(params, "prepare_data_linear_a23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_a3(params, data)
  files <- "pa.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_blocks_linear_a3(params)
  params <- get_z_linear_a3(params, data)
  files <- seq_zw("cz_", length(params$container$file_break_z))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- get_wr_linear_a3(params, data)
  files <- c("xatxa.rdata", seq_zw("cpr_", length(params$container$filebreak.PR)))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- get_results_linear_a3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


PartyBProcess3Linear <- function(data,
                                monitor_folder  = NULL,
                                sleep_time      = 10,
                                max_waiting_time = 24 * 60 * 60,
                                popmednet      = TRUE,
                                trace          = FALSE,
                                verbose        = TRUE) {
  params <- prepare_params_3p("linear", "B",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)

  header(params)
  params   = prepare_folder_linear_b3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  data = prepare_data_linear_B23(params, data)
  params <- add_to_log(params, "prepare_data_linear_B23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_b3(params, data)
  files <- "pb.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)


  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_blocks_linear_b3(params)
  params <- GetRWLinear.b3(params, data)
  files <- c("xbtxb.rdata", seq_zw("crw_", length(params$container$filebreak.RW)))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- get_results_linear_b3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


PartyTProcess3Linear <- function(monitor_folder         = NULL,
                                msreqid               = "v_default_0_000",
                                blocksize             = 500,
                                sleep_time             = 10,
                                max_waiting_time        = 24 * 60 * 60,
                                popmednet             = TRUE,
                                trace                 = FALSE,
                                verbose               = TRUE) {
  params <- prepare_params_3p("linear", "T", msreqid = msreqid,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)

  header(params)
  params   = prepare_folder_linear_t3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.3p(params, from = c("A", "B"), max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata")) &&
      file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(paste(read_error_message(params$read_path[["A"]]), "\n",
                  read_error_message(params$read_path[["B"]])))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["A"]]))
    file.copy(file.path(params$read_path[["A"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, filesB = files, from = "B",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["B"]]))
    file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params   = prepare_params_linear_t3(params)
  if (!params$failed) params <- prepare_blocks_linear_t3(params, blocksize)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files <- "blocks.rdata"
  params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- ProcessZLinear.t3(params)
  files <- c("blocks.rdata", seq_zw("crz_", length(params$container$filebreak.RZ)))
  params <- send_pause_continue_3p(params, filesB = files, from  = "B",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- ProcessWLinear.t3(params)
  files <- c("p2.rdata", seq_zw("cwr_", length(params$container$filebreak.WR)))
  params <- send_pause_continue_3p(params, filesA = files, from  = "A",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- GetProductsLinear.t3(params)
  params <- ComputeResultsLinear.t3(params)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files <- "stats.rdata"
  params <- send_pause_continue_3p(params, filesA = files, filesB = files, from  = c("A", "B"),
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- send_pause_quit_3p(params, sleep_time = sleep_time)
  SummarizeLog.3p(params)
  return(params$stats)
}
