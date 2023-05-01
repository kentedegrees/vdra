################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

PrepareFolderLogistic.a2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.a2\n\n")
  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "inputfiles")
  params$read_path       <- file.path(monitor_folder, "msoc1")

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

  params <- add_to_log(params, "PrepareDataLogistic.a23, PrepareFolderLogistic.a2", 0, 0, 0, 0)
  return(params)
}


PrepareFolderLogistic.b2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLogistic.b2\n\n")

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "msoc")
  params$read_path       <- file.path(monitor_folder, "inputfiles")

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

  params <- add_to_log(params, "PrepareDataLogisitc.b23, PrepareFolderLogistic.b2", 0, 0, 0, 0)

  return(params)
}


#' @importFrom stats model.matrix
PrepareDataLogistic.a23 <- function(params, data, y_name = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.a23\n\n")

  workdata = list()
  workdata$failed = FALSE

  workdata$failed = CheckDataFormat(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data = data.frame(data) # convert to a clean data.frame

  response_index = CheckResponse(params, data, y_name)

  if (is.null(response_index)) {
    workdata$failed = TRUE
    return(workdata)
  }
  covariate_index = setdiff(1:ncol(data), response_index)
  workdata$tags = CreateModelMatrixTags(data[, covariate_index, drop = FALSE])
  workdata$tags = c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] = "numeric"

  x = model.matrix(~ ., data[, c(response_index, covariate_index), drop = FALSE])
  rownames(x) = NULL
  covariate_index = setdiff(1:ncol(x), 2)

  means = apply(x, 2, mean)
  sd    = apply(x, 2, sd)
  sd    = sapply(sd, function(x) {
    ifelse(x > 0, x, 1)
  })
  workdata$Y      = x[, 2, drop = FALSE]
  workdata$x      = x[, covariate_index, drop = FALSE]
  workdata$means  = means[covariate_index]
  workdata$sd     = sd[covariate_index]
  workdata$yty    = t(workdata$Y) %*% workdata$Y

  if (ncol(workdata$x) >= 2) {
    for (i in 2:ncol(workdata$x)) {
      workdata$x[, i] = (workdata$x[, i] - workdata$means[i]) / workdata$sd[i]
    }
  }

  return(workdata)
}

#' @importFrom stats model.matrix sd
PrepareDataLogistic.b23 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.b23\n\n")

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
  workdata$x = model.matrix(~ ., data)
  rownames(workdata$x) = NULL
  workdata$x = workdata$x[, -1, drop = FALSE]
  workdata$means = apply(workdata$x, 2, mean)
  workdata$sd    = apply(workdata$x, 2, sd)
  workdata$sd    = sapply(workdata$sd, function(x) {
    ifelse(x > 0, x, 1)
  })

  for (i in 1:ncol(workdata$x)) {
    workdata$x[, i] = (workdata$x[, i] - workdata$means[i]) / workdata$sd[i]
  }


  return(workdata)
}

PrepareParamsLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLogistic.b2\n\n")
  params$failed         = FALSE
  params$converged      = FALSE
  params$halted         = FALSE

  params$n             = nrow(data$x)
  params$num_events     = 0
  params$p1            = 0
  params$p2            = ncol(data$x)
  params$p             = params$p1 + params$p2
  params$p1_old        = 0
  params$p2_old        = params$p2
  params$a_col_names     = c("")
  params$b_col_names     = colnames(data$x)
  params$y_name         = ""
  params$a_col_names_old = c("")
  params$b_col_names_old = c("")
  params$cutoff        = 1
  params$max_iterations = 1

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
  pb$b_col_names = params$b_col_names
  pb$tags      = data$tags

  write_time <- proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLogistic.b2", 0, 0, write_time, write_size)
  return(params)
}


PrepareParamsLogistic.a2 <- function(params, data, cutoff = 0.01, max_iterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLogistic.a2\n\n")

  params$converged       = FALSE
  params$halted          = FALSE
  params$pmnStepCounter  = 1
  pb                     = NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb, b_col_names
  read_size <- sum(file.size(file.path(params$read_path, "pb.rdata")))
  read_time <- proc.time()[3] - read_time

  if (params$analysis != pb$analysis) {
    params$error_message <-
      paste("Party A is running", params$analysis, "regression and Party B is running", pb$analysis, "regression.")
    warning(params$error_message)
    params$failed <- TRUE
    return(params)
  }

  params$n  = nrow(data$x)
  if (pb$n != params$n) {
    params$error_message <-
      paste("Party A has", params$n, "observations and Party B has", pb$n, "observations.")
    warning(params$error_message)
    params$failed <- TRUE
  }

  params$p1 = ncol(data$x)
  params$p2 = pb$p2
  params$p  = params$p1 + params$p2
  params$p1_old = params$p1
  params$p2_old = params$p2

  params$a_col_names = colnames(data$x)
  params$b_col_names = pb$b_col_names
  params$y_name     = colnames(data$Y)
  params$a_col_names_old = c("")
  params$b_col_names_old = c("")
  params$Atags         = data$tags
  params$Btags         = pb$tags

  if (cutoff <= 0) cutoff = 0.01
  if (cutoff >= 1) cutoff = 0.05
  params$cutoff           = cutoff

  if (max_iterations < 1) max_iterations = 1
  params$max_iterations = max_iterations

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
  pa$y_name         = data$y_name
  pa$cutoff        = params$cutoff
  pa$max_iterations = params$max_iterations
  pa$a_col_names     = params$a_col_names
  pa$tags          = data$tags

  write_time <- proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareParamsLogistic.a2", read_time, read_size,
                    write_time, write_size)

  return(params)
}


PrepareBlocksLogistic.a2 <- function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLogistic.a2\n\n")
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
    params <- add_to_log(params, "PrepareBlocksLogistic.a2", 0, 0, 0, 0)
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
  write_time <- proc.time()[3]
  save(blocksize, file = file.path(params$write_path, "blocksize.rdata"))
  write_size <- file.size(file.path(params$write_path, "blocksize.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "PrepareBlocksLogistic.a2", 0, 0, write_time, write_size)
  return(params)
}


GetZLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLogistic.a2\n\n")
  write_time <- 0
  write_size <- 0

  num_blocks = params$blocks$num_blocks
  pbar <- MakeProgressBar1(num_blocks, "z", params$verbose)
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
    z <- FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$x[strt:stp, ]), g)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(z), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "GetZLogistic.a2", 0, 0, write_time, write_size)
  return(params)
}


FinalizeParamsLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "FinalizeParamsLogistic.b2\n\n")
  pa = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # Load pa, a_col_names
  read_size <- sum(file.size(file.path(params$read_path, "pa.rdata")))
  read_time <- proc.time()[3] - read_time
  params$p1            = pa$p1
  params$p             = params$p1 + params$p2
  params$meansA        = pa$means
  params$sdA           = pa$sd
  params$yty           = pa$yty
  params$y_name         = pa$y_name
  params$cutoff        = pa$cutoff
  params$max_iterations = pa$max_iterations
  params$a_col_names     = pa$a_col_names
  params$Atags         = pa$tags
  params$Btags         = data$tags

  params <- add_to_log(params, "FinalizeParamsLogistic.b2", read_time, read_size, 0, 0)
  return(params)
}


PrepareBlocksLogistic.b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLogistic.b2\n\n")
  blocksize = NULL
  # For now, assuming that p1 > 0 and p2 > 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "blocksize.rdata")) # load blocksize
  read_size <- file.size(file.path(params$read_path, "blocksize.rdata"))
  read_time <- proc.time()[3] - read_time
  params$blocks    = CreateBlocks(params$p1, params$p2, params$n, blocksize)
  params$container = CreateContainers(params$p1, params$p2, params$blocks)
  params <- add_to_log(params, "PrepareBlocksLogistic.b2", read_time, read_size, 0, 0)
  return(params)
}


GetWLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWLogistic.b2\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "(I-z*z')XB", params$verbose)

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
    if (i %in% params$container$filebreak.w) {
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

    w = data$x[strt:stp, ] - z %*% (t(z) %*% data$x[strt:stp, ])

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(w), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z || i == params$blocks$num_blocks) {
      close(to_read)
    }
    if ((i + 1) %in% params$container$filebreak.w || i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }

  write_time <- write_time - proc.time()[3]
  save(xb_t_xb, file = file.path(params$write_path, "xbtxb.rdata"))
  write_size <- write_size + file.size(file.path(params$write_path, "xbtxb.rdata"))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "GetWLogistic.b2", read_time, read_size, write_time, write_size)

  return(params)
}


CheckColinearityLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.a2\n\n")
  p2 = params$p2
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  xb_t_xb     = NULL

  read_time <- read_time - proc.time()[3]
  load(file.path(params$read_path, "xbtxb.rdata")) # load xb_t_xb
  read_size <- file.size(file.path(params$read_path, "xbtxb.rdata"))
  read_time <- read_time + proc.time()[3]
  XATXA = t(data$x) %*% data$x
  XATXB = 0
  XATY  = t(data$x) %*% data$Y
  YTXB  = 0

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "X'X", params$verbose)

  container_ct_w <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.w) {
      container_ct_w <- container_ct_w + 1
      filename = paste0("cw_", container_ct_w, ".rdata")
      to_read <- file(file.path(params$read_path, filename), "rb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1

    read_time <- read_time - proc.time()[3]
    w = matrix(readBin(con = to_read, what = numeric(), n = n2 * p2,
                       endian = "little"), nrow = n2, ncol = p2)
    read_time <- read_time + proc.time()[3]

    XATXB = XATXB + t(data$x[strt:stp, ]) %*% w
    YTXB  = YTXB  + t(data$Y[strt:stp, ]) %*% w

    if ((i + 1) %in% params$container$filebreak.w || i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }

  xtx = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), xb_t_xb))
  XTY = rbind(XATY, t(YTXB))

  nrow = nrow(xtx)
  indicies = c(1)
  for (i in 2:nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  xtx = xtx[indicies, indicies]
  XTY = matrix(XTY[indicies], ncol = 1)

  a_names   = params$a_col_names
  b_names   = params$b_col_names
  a_index   = which(indicies <= length(a_names))
  params$IndiciesKeep  = indicies
  params$a_indicies_keep = indicies[a_index]
  params$b_indicies_keep = indicies[-a_index] - length(a_names)

  AnamesKeep = a_names[params$a_indicies_keep]
  BnamesKeep = b_names[params$b_indicies_keep]
  params$a_col_names_old = params$a_col_names
  params$b_col_names_old = params$b_col_names
  params$a_col_names     = AnamesKeep
  params$b_col_names     = BnamesKeep
  params$p1_old        = params$p1
  params$p2_old        = params$p2
  params$p1            = length(AnamesKeep)
  params$p2            = length(BnamesKeep)
  params$p_old         = params$p1_old + params$p2_old
  params$p             = params$p1 + params$p2
  params$meansA        = params$meansA[params$a_indicies_keep]
  params$meansB        = params$meansB[params$b_indicies_keep]
  params$sdA           = params$sdA[params$a_indicies_keep]
  params$sdB           = params$sdB[params$b_indicies_keep]
  params$xtx           = xtx
  params$xty           = XTY

  a_indicies = params$a_indicies_keep
  b_indicies = params$b_indicies_keep

  write_time <- write_time - proc.time()[3]
  save(a_indicies, b_indicies, file = file.path(params$write_path, "indicies.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "indicies.rdata")))
  write_time <- write_time + proc.time()[3]

  tags = params$Btags[params$b_indicies_keep]

  if (length(unique(tags)) < 2) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has 1 or fewer covariates.\n\n"
  } else if (!("numeric" %in% names(tags))) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has no continuous covariates.\n\n"
  }

  params <- add_to_log(params, "CheckColinearityLogistic.a2", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateDataLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.a2\n\n")
  data$x = as.matrix(data$x[, params$a_indicies_keep, drop = FALSE])
  return(data)
}


ComputeInitialBetasLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.a2\n\n")
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

  params$alg_iteration_counter      = 1
  params$deltabeta = Inf

  write_time <- proc.time()[3]
  save(Bbetas, Bxty, file = file.path(params$write_path, "Bbetas_xty.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "Bbetas_xty.rdata")))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeInitialBetasLogistic.a2", 0, 0, write_time, write_size)

  return(params)
}


UpdateParamsLogistic.b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.b2\n\n")
  a_indicies = NULL
  b_indicies = NULL
  Bbetas    = NULL
  Bxty      = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "indicies.rdata")) # load a_indicies, b_indicies
  load(file.path(params$read_path, "Bbetas_xty.rdata"))     # Load Bbetas
  read_size <- sum(file.size(file.path(params$read_path, c("indicies.rdata",
                                                        "Bbetas_xty.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$a_col_names_old = params$a_col_names
  params$b_col_names_old = params$b_col_names
  params$a_col_names     = params$a_col_names_old[a_indicies]
  params$b_col_names     = params$b_col_names_old[b_indicies]
  params$p1_old = params$p1
  params$p2_old = params$p2
  params$p1     = length(a_indicies)
  params$p2     = length(b_indicies)
  params$p_old  = params$p
  params$p      = params$p1 + params$p2
  params$b_indicies_keep = b_indicies
  params$a_indicies_keep = a_indicies
  params$betasB    = Bbetas
  params$betasBold = matrix(0, params$p2, 1)
  params$meansB = params$meansB[b_indicies]
  params$sdB    = params$sdB[b_indicies]
  params$Bxty   = Bxty
  params <- add_to_log(params, "UpdateParamsLogistic.b2", read_time, read_size, 0, 0)
  return(params)
}


UpdateDataLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.b2\n\n")
  data$x = as.matrix(data$x[, params$b_indicies_keep, drop = FALSE])
  return(data)
}


GetXBetaLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaLogistic.b2\n\n")
  XbetaB = data$x %*% params$betasB

  write_time <- proc.time()[3]
  save(XbetaB, file = file.path(params$write_path, "xbetab.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbetab.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "GetXBetaLogistic.b2", 0, 0, write_time, write_size)
  return(params)
}


GetWeightsLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWeightsLogistic.a2\n\n")
  n      = params$n
  p1     = params$p1
  XbetaB = NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "xbetab.rdata"))  # Load XbetaB
  read_size <- file.size(file.path(params$read_path, "xbetab.rdata"))
  read_time <- proc.time()[3] - read_time

  XbetaA = data$x %*% params$betasA
  Xbeta = XbetaA + XbetaB
  pi_ = (1 + exp(-Xbeta))^(-1)
  params$pi_ = pi_

  write_time <- proc.time()[3]
  save(pi_, file = file.path(params$write_path, "pi_.rdata"))
  write_size <- file.size(file.path(params$write_path, "pi_.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetWeightsLogistic.a2", read_time, read_size, write_time, write_size)
  return(params)
}


GetVLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetVLogistic.b2\n\n")
  pi_       = NULL
  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pi_.rdata"))
  read_size <- file.size(file.path(params$read_path, "pi_.rdata"))
  read_time <- proc.time()[3] - read_time

  params$pi_ = pi_
  w = pi_ * (1 - params$pi_)

  XBTWXB = 0

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "(I - z*z')w*XB", params$verbose)

  container_ct_z <- 0
  container_ct_V = 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
    }
    if (i %in% params$container$filebreak.V) {
      container_ct_V = container_ct_V + 1
      filename2 <- paste0("cv_", container_ct_V, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1
    g = params$blocks$g[i]

    Xblock  = data$x[strt:stp, ]
    Wblock  = w[strt:stp]
    WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

    read_time <- read_time - proc.time()[3]
    z <- matrix(readBin(con = to_read, what = numeric(), n = n * g,
                       endian = "little"), nrow = n, ncol = g)
    read_time <- read_time + proc.time()[3]

    V = WXblock - z %*% (t(z) %*% WXblock)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(V), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    XBTWXB = XBTWXB + t(Xblock) %*% WXblock
    if ((i + 1) %in% params$container$file_break_z || i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }
  # sums of each column in WX_B
  sumsWXB = apply(MultiplyDiagonalWTimesX(w, data$x), 2, sum)
  # This information needs to be shared in order to get the intercept term

  write_time <- write_time - proc.time()[3]
  save(sumsWXB, XBTWXB, file = file.path(params$write_path, "sumswx_xbtwxb.rdata"))
  write_size <- write_size + sum(file.size(c(file.path(params$write_path, "sumswx_xbtwxb.rdata"))))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "GetVLogistic.b2", read_time, read_size, write_time, write_size)
  return(params)
}


GetIILogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetIILogistic.a2\n\n")
  p1 = params$p1
  p2 = params$p2
  sumsWXB = NULL
  XBTWXB  = NULL

  write_time <- 0
  write_size <- 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "sumswx_xbtwxb.rdata")) # load sumsWXB, XBTWXB
  read_size <- sum(file.size(file.path(params$read_path, "sumswx_xbtwxb.rdata")))
  read_time <- proc.time()[3] - read_time

  params$sumsWXB = sumsWXB

  IA = params$Axty - t(data$x) %*% params$pi_
  w = params$pi_ * (1 - params$pi_)
  sumsWXA = apply(MultiplyDiagonalWTimesX(w, data$x), 2, sum)[-1]
  params$sumsWXA = sumsWXA

  XATWXA = t(data$x) %*% MultiplyDiagonalWTimesX(w, data$x)

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "X'W*X", params$verbose)

  XATWXB = 0
  container_ct_V = 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.V) {
      container_ct_V = container_ct_V + 1
      filename1 <- paste0("cv_", container_ct_V, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1

    read_time <- read_time - proc.time()[3]
    V = matrix(readBin(con = to_read, what = numeric(),
                       n = n * p2, endian = "little"), n, p2)
    read_time <- read_time + proc.time()[3]
    XATWXB = XATWXB + t(data$x[strt:stp, ]) %*% V

    pbar <- MakeProgressBar2(i, pbar, params$verbose)
    if ((i + 1) %in% params$container$filebreak.V || i == params$blocks$num_blocks) {
      close(to_read)
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
    params$singular_matrix = TRUE
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
    params$II = II
    params$IA = IA

    a21i1 = II[(p1 + 1):(p1 + p2), 1:p1] %*% matrix(IA, p1, 1)
    a11i1 = II[1:p1, 1:p1] %*% matrix(IA, p1, 1)
    params$a11i1 = a11i1

    write_time <- proc.time()[3]
    save(a21i1, XTWX, file = file.path(params$write_path, "a21i1_xtwx.rdata"))
    write_size <- sum(file.size(file.path(params$write_path, "a21i1_xtwx.rdata")))
    write_time <- proc.time()[3] - write_time
  }
  params <- add_to_log(params, "GetIILogistic.a2", read_time, read_size, write_time, write_size)

  return(params)
}


GetCoefLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.b2\n\n")
  p1 = params$p1
  p2 = params$p2
  XTWX  = NULL
  a21i1 = NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "a21i1_xtwx.rdata"))   # load a21i1, XTWX
  read_size <- sum(file.size(file.path(params$read_path, "a21i1_xtwx.rdata")))
  read_time <- proc.time()[3] - read_time

  IB = params$Bxty - t(data$x) %*% params$pi_

  II = solve(XTWX)

  params$II = II
  params$IB = IB

  a22i2 = II[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  a12i2 = II[1:p1, (p1 + 1):(p1 + p2), drop = FALSE] %*% IB
  params$a22i2 = a22i2

  params$betasBold = params$betasB
  params$betasB = params$betasB + a21i1 + a22i2

  deltabetaB = max(abs(params$betasB - params$betasBold) / (abs(params$betasB) + 0.1))

  write_time <- proc.time()[3]
  save(a12i2, deltabetaB, file = file.path(params$write_path, "a12_deltabetaB.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "a12_deltabetaB.rdata")))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "GetCoefLogistic.b2", read_time, read_size, write_time, write_size)
  return(params)
}


GetCoefLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetCoefLogistic.a2\n\n")
  a12i2      = NULL
  deltabetaB = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "a12_deltabetaB.rdata")) # load  a12i2, deltabetab
  read_size <- sum(file.size(file.path(params$read_path, "a12_deltabetaB.rdata")))
  read_time <- proc.time()[3] - read_time

  params$betasAold = params$betasA
  params$betasA = params$betasA + params$a11i1 + a12i2

  deltabeta = max(abs(params$betasA - params$betasAold) / (abs(params$betasA) + 0.1), deltabetaB)

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$alg_iteration_counter >= params$max_iterations) {
    params$maxIterExceeded = TRUE
    warning(paste("Failed to converged in", params$max_iterations, "iterations."))
  }

  write_time <- proc.time()[3]
  save(deltabeta, file = file.path(params$write_path, "deltabeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "deltabeta.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "GetCoefLogistic.a2", read_time, read_size, write_time, write_size)


  return(params)
}


GetConvergedStatusLogistic.b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetconvergedStatusLogistic.b2\n\n")
  deltabeta = NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "deltabeta.rdata")) # load deltabeta.rdata
  read_size <- file.size(file.path(params$read_path, "deltabeta.rdata"))
  read_time <- proc.time()[3] - read_time

  if (deltabeta < params$cutoff)  {
    params$converged = TRUE
  } else if (params$alg_iteration_counter >= params$max_iterations) {
    params$maxIterExceeded = TRUE
    warning(paste("Failed to converged in", params$max_iterations, "iterations."))
  }

  params <- add_to_log(params, "GetConvergedStatusLogistic.b2", read_time, read_size, 0, 0)
  return(params)
}


GetFinalCoefLogistic.b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetFinalCoefLogistic.b2\n\n")
  betasB = params$betasB / params$sdB
  offsetB = sum(betasB * params$meansB)
  BFinalFitted = t(params$sdB * t(data$x) + params$meansB) %*% betasB
  write_time <- proc.time()[3]
  save(betasB, BFinalFitted, offsetB, file = file.path(params$write_path, "b_final.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "b_final.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetFinalCoefLogistic.b2", 0, 0, write_time, write_size)
  return(params)
}

#' @importFrom stats pnorm
ComputeResultsLogistic.a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.a2\n\n")
  stats <- params$stats
  stats$failed         <- FALSE
  stats$converged      <- params$converged

  n      <- params$n
  p1     <- params$p1
  p2     <- params$p2
  sdA    <- params$sdA
  sdB    <- params$sdB
  meansA <- params$meansA
  meansB <- params$meansB
  a_names <- params$a_col_names_old
  b_names <- params$b_col_names_old
  p1_old <- params$p1_old
  p2_old <- params$p2_old
  p_old  <- params$p_old
  indicies <- params$IndiciesKeep


  betasB <- NULL
  offsetB <- NULL
  BFinalFitted <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "b_final.rdata"))  # betasB, offsetB, BFinalFitted
  read_size <- sum(file.size(file.path(params$read_path, "b_final.rdata")))
  read_time <- proc.time()[3] - read_time
  betasA <- params$betasA / sdA
  offsetA <- sum(betasA[-1] * params$meansA[-1])
  betasA[1] <- betasA[1] - offsetA - offsetB
  betas <- c(betasA, betasB)

  AFinalFitted <- t(sdA * t(data$x) + meansA) %*% betasA -
    t(sdA[1] * t(data$x[, 1]) + meansA[1]) %*% betasA[1]
  FinalFitted <- AFinalFitted + BFinalFitted + betas[1]
  params$FinalFitted <- FinalFitted

  n <- params$n
  ct      <- sum(data$Y)
  resdev  <- -2 * (sum(data$Y * FinalFitted) - sum(log(1 + exp(FinalFitted))))
  nulldev <- -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

  # If xtwx were singular, it would have been caught in GetII.a2(), so we may
  # assume that xtwx is NOT singular and so we do not have to do a check.
  cov1 <- solve(params$xtwx)
  secoef <- sqrt(diag(cov1)) / c(sdA, sdB)
  tmp <- matrix(c(1, (-meansA / sdA)[-1], -meansB / sdB), ncol = 1)
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
  stats$ROC     <- RocInternal(params, data)
  stats$iter    <- params$alg_iteration_counter - 1

  names_old <- c(a_names, b_names)
  names(stats$coefficients) <- names_old
  names(stats$party) <- names_old
  names(stats$secoef) <- names_old
  names(stats$tvals) <- names_old
  names(stats$pvals) <- names_old

  write_time <- proc.time()[3]
  save(stats, file <- file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  stats$Y           <- data$Y # For Hoslem and ROC
  stats$FinalFitted <- FinalFitted
  params$stats      <- stats

  params <- add_to_log(params, "ComputeResultsLogistic.b2", read_time, read_size, write_time, write_size)
  return(params)
}


GetResultsLogistic.b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.b2\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size <- file.size(file.path(params$read_path, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsLogistic.b2", read_time, read_size, 0, 0)
  return(params)
}



############################### PARENT FUNCTIONS ###############################


PartyAProcess2Logistic <- function(data,
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
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params   = PrepareFolderLogistic.a2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLogistic.a23(params, data, y_name)

  params <- PauseContinue.2p(params,  max_waiting_time)
  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(read_error_message(params$read_path))
    params$pmnStepCounter = 1
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
    params$completed = TRUE
    message = "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params$pmnStepCounter = 1
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- PrepareParamsLogistic.a2(params, data, cutoff, max_iterations)

  if (params$failed) {   # Check for failed from PrepareParamsLogistic.a2()
    params$completed = TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- PrepareBlocksLogistic.a2(params, blocksize)

  if (params$failed) { # Check for failed from PrepareBlocksLogistic.a2()
    params$completed = TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- GetZLogistic.a2(params, data)

  files <- c("pa.rdata", "blocksize.rdata",
            seq_zw("cz_", length(params$container$file_break_z)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- CheckColinearityLogistic.a2(params, data)

  if (params$failed) { # Check for CheckColinearityLogistic.a2() failed
    params$completed = TRUE
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }
  data = UpdateDataLogistic.a2(params, data)
  params <- add_to_log(params, "UpdateDataLogistic.a2", 0, 0, 0, 0)
  params <- ComputeInitialBetasLogistic.a2(params, data)

  files <- c("indicies.rdata", "Bbetas_xty.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetWeightsLogistic.a2(params, data)
    files <- c("pi_.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)
    params <- GetIILogistic.a2(params, data)

    if (params$failed) { # Check for failed from ComputeInverseLogistic.a2()
      params$completed = TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
      params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    files <- c("a21i1_xtwx.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

    params <- GetCoefLogistic.a2(params, data)
    files <- "deltabeta.rdata"
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

    EndingIteration(params)
    params$alg_iteration_counter = params$alg_iteration_counter + 1
  }
  params$completed = TRUE

  params <- ComputeResultsLogistic.a2(params, data)

  files <- c("stats.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)
  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  SummarizeLog.2p(params)
  return(invisible(params$stats))
}

PartyBProcess2Logistic <- function(data,
                                  monitor_folder       = "v_default_00_000",
                                  sleep_time           = 10,
                                  max_waiting_time      = 24 * 60 * 60,
                                  popmednet           = TRUE,
                                  trace               = FALSE,
                                  verbose             = TRUE) {
  params <- prepare_params_2p("logistic", "B",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params   = PrepareFolderLogistic.b2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = PrepareDataLogistic.b23(params, data)

  if (data$failed) { # Check for Error from PrepareDataLogistic.b2()
    params$completed = TRUE
    message = "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_2p(params, files, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params   = PrepareParamsLogistic.b2(params, data)

  files <- c("pb.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(read_error_message(params$read_path))
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params <- FinalizeParamsLogistic.b2(params, data)
  params <- PrepareBlocksLogistic.b2(params)
  params <- GetWLogistic.b2(params, data)

  files <- c("xbtxb.rdata", seq_zw("cw_", length(params$container$filebreak.w)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$completed = TRUE
    warning(read_error_message(params$read_path))
    params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
    return(params$stats)
  }

  params <- UpdateParamsLogistic.b2(params)
  data = UpdateDataLogistic.b2(params, data)
  params <- add_to_log(params, "UpdateDataLogistic.b2", 0, 0, 0, 0)

  params$alg_iteration_counter = 1
  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    params <- GetXBetaLogistic.b2(params, data)

    files <- c("xbetab.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

    params <- GetVLogistic.b2(params, data)
    files <- c("sumswx_xbtwxb.rdata",
              seq_zw("cv_", length(params$container$filebreak.V)))
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$completed = TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time, job_failed = TRUE)
      return(params$stats)
    }

    params <- GetCoefLogistic.b2(params, data)
    files <- c("a12_deltabetaB.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

    params <- GetConvergedStatusLogistic.b2(params)

    EndingIteration(params)
    params$alg_iteration_counter = params$alg_iteration_counter + 1
  }
  params$completed = TRUE

  params <- GetFinalCoefLogistic.b2(params, data)
  files <- c("b_final.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- GetResultsLogistic.b2(params)
  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  return(invisible(params$stats))
}
