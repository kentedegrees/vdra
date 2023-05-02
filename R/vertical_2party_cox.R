#################### DISTRIBUTED COX REGRESSION FUNCTIONS ####################

prepare_folder_cox_a2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_folder_cox_a2\n\n")
  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "inputfiles")
  params$read_path       <- file.path(monitor_folder, "msoc1")

  if (is.null(monitor_folder)) {
    warning(paste("monitor_folder must be specified.",
                  "Please use the same monitor_folder as the DataMart Client."))
    params$failed <- TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning(paste("monitor_folder directory is not valid.",
                  "Please use the same monitor_folder as the DataMart Client."))
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

  params <- add_to_log(params, "prepare_data_cox_a23, prepare_folder_cox_a2",
                       0, 0, 0, 0)
  return(params)
}


prepare_folder_cox_b2 <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_folder_cox_b2\n\n")

  params$dp_local_path   <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  params$write_path      <- file.path(monitor_folder, "msoc")
  params$read_path       <- file.path(monitor_folder, "inputfiles")

  if (is.null(monitor_folder)) {
    warning(paste("monitor_folder must be specified.",
                  "Please use the same monitor_folder as the DataMart Client."))
    params$failed <- TRUE
    return(params)
  }
  if (class(monitor_folder) != "character") {
    warning(paste("monitor_folder directory is not valid.",
                  "Please use the same monitor_folder as the DataMart Client."))
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
                                  "Could not create directory",
                                  paste0(params$read_path, "."),
                                  "Check the path and restart the program.")
  }

  Sys.sleep(1)
  delete_trigger("files_done.ok", params$read_path)

  params <- add_to_log(params, "prepare_data_cox_b23, prepare_folder_cox_b2",
                       0, 0, 0, 0)

  return(params)
}


extract_strata <- function(params, data, stratas, mask) {
  strata <- list()
  strata$failed <- FALSE
  if (!is.null(stratas)) {
    if (!("character" %in% class(stratas))) {
      warning("Strata is not a valid variable name(s).")
      strata$failed <- TRUE
      return(strata)
    }
    if (length(stratas) > 0) {
      idx <- stratas %in% colnames(data)
      if (!is.null(params$party_name) && params$party_name == "A")  {
        strata$strata_from_a <- stratas[idx]
        strata$strata_from_b <- stratas[!idx]
      } else {
        strata$strata_from_b <- stratas[idx]
        strata$strata_from_a <- stratas[!idx]
      }
      if (!is.null(params$data_partner_id) && params$data_partner_id == "1") {
        strata$strata_from_me     <- stratas[idx]
        strata$strata_from_others <- stratas[!idx]
      } else {
        strata$strata_from_me     <- stratas[idx]
        strata$strata_from_others <- stratas[!idx]
      }
      strata$strata_index <- which(colnames(data) %in% stratas)
      if (length(strata$strata_index) > 0) {
        strata$x <- data[, strata$strata_index, drop = FALSE]
        strata$legend <- list()
        for (i in seq_len(ncol(strata$x))) {
          levels <- unique(strata$x[, i])
          strata$legend[[colnames(strata$x)[i]]] <- levels
          if (mask) {
            levels <- sample(levels)
            strata$legend[[colnames(strata$x)[i]]] <- rep("NA", length(levels))
          }
          strata$x[, i] <- sapply(strata$x[, i],
                                  function(x) {
                                    which(levels %in% x)
                                  })
        }
      }
    }
  }
  return(strata)
}


#' @importFrom stats model.matrix
prepare_data_cox_23 <- function(params, data, y_name, strata, mask) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_data_cox_23\n\n")

  workdata <- list()
  workdata$failed <- check_data_format(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data <- data.frame(data) # convert to a clean data.frame

  workdata$strata <- extract_strata(params, data, strata, mask)
  if (workdata$strata$failed) {
    workdata$failed <- TRUE
    return(workdata)
  }
  strata_index <- workdata$strata$strata_index
  response_index <- numeric()

  if (params$party_name == "A") {
    response_index <- check_response(params, data, y_name)

    if (is.null(response_index)) {
      workdata$failed <- TRUE
      return(workdata)
    }

    workdata$survival        <- list()
    workdata$survival$rank   <- data[, response_index[1]]
    workdata$survival$status <- data[, response_index[2]]
    if (length(intersect(strata_index, response_index)) > 0) {
      warning("Response and strata share a variable.")
      workdata$failed <- TRUE
      return(workdata)
    }
  }

  covariate_index <- setdiff(seq_len(ncol(data)), union(strata_index, response_index))

  if (length(covariate_index) == 0) {
    if (params$party_name == "A") {
      workdata$x  <- matrix(0, nrow = nrow(data), ncol = 0)
    } else {
      warning(paste("After removing strata, data is empty.",
                    "Party B must supply at least one non-strata covariate."))
      workdata$failed <- TRUE
      return(workdata)
    }
  } else {
    workdata$tags <- create_model_matrix_tags(data[, covariate_index,
                                                   drop = FALSE])
    if (params$party_name == "B" && (ncol(data) < 2 ||
                                     !("numeric" %in% names(workdata$tags)))) {
      warning(paste("The data partner that does not have the response",
                    "must have at least 2 covariates at least one of",
                    "which must be numeric."))
      workdata$failed <- TRUE
      return(workdata)
    }
    workdata$x <- scale(model.matrix(~ ., data[, covariate_index,
                                               drop = FALSE]),
                        center = TRUE, scale = FALSE)
    workdata$x <- workdata$x[, -1, drop = FALSE]
  }

  return(workdata)
}


prepare_params_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_b2\n\n")
  params$failed          <- FALSE
  params$converged       <- FALSE
  params$halted          <- FALSE
  params$singular_matrix <- FALSE

  params$n <- nrow(data$x)
  params$num_events <- 0
  params$p1 <- 0
  params$p2 <- ncol(data$x)
  params$p  <- params$p1 + params$p2
  params$p1_old <- params$p1
  params$p2_old <- params$p2
  params$a_col_names <- c("")
  params$b_col_names <- colnames(data$x)
  params$a_col_names_old <- c("")
  params$b_col_names_old <- c("")
  params$cutoff        <- 1
  params$max_iterations <- 1

  params$survival_installed <- requireNamespace("survival", quietly = TRUE)
  if (params$survival_installed && !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  pb           <- list()
  pb$p2        <- params$p2
  pb$n         <- params$n
  pb$analysis  <- params$analysis
  pb$b_col_names <- params$b_col_names
  pb$strata_b   <- list()
  pb$strata_b$strata_from_a <- data$strata$strata_from_a
  pb$strata_b$strata_from_b <- data$strata$strata_from_b
  pb$tags      <- data$tags

  write_time <- proc.time()[3]
  save(pb, file <- file.path(params$write_path, "pb.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_cox_b2", 0, 0,
                       write_time, write_size)
  return(params)
}

check_strata_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "check_strata_cox_a2\n\n")
  pb <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb
  read_size <- sum(file.size(file.path(params$read_path, "pb.rdata")))
  read_time <- proc.time()[3] - read_time

  if (length(pb$strata_b$strata_from_a) == length(data$strata$strata_from_a) &&
      length(pb$strata_b$strata_from_b) == length(data$strata$strata_from_b) &&
      ifelse(length(pb$strata_b$strata_from_a) == 0, TRUE,
             order(pb$strata_b$strata_from_a) ==
             order(data$strata$strata_from_a)) &&
      ifelse(length(pb$strata_b$strata_from_b) == 0, TRUE,
             order(pb$strata_b$strata_from_b) ==
             order(data$strata$strata_from_b))) {
    params$get_strata_from_b <- length(data$strata$strata_from_b) > 0
  } else {
    params$get_strata_from_b <- FALSE
    a_cap_b <- intersect(data$strata$strata_from_a, pb$strata_b$strata_from_b)
    b_cap_a <- intersect(data$strata$strata_from_b, pb$strata_b$strata_from_a)
    if (length(a_cap_b) > 0) {
      params$error_message <-
        paste("Party A and Party B have",
              length(a_cap_b),
              "variable(s) with the same name which are used in the strata.",
              "These variable(s) are <", paste0(a_cap_b, collapse = ", "), ">.",
              "Make sure the variables from each party have distinct names_")
    } else if (length(b_cap_a) > 0) {
      params$error_message <-
        paste("Party A and Party B have specified",
              length(b_cap_a),
              "variable(s) for the strata which are not found in the data.",
              "These variable(s) are <", paste0(b_cap_a, collapse = ", "), ">.",
              "Check the spelling of the variables names and / or",
              "remove them from the strata.")
    } else {
      params$error_message <-
        paste("Party A and Party B have specified different strata.",
              "Verify that both parties specify the same strata.")
    }
    warning(params$error_message)
    params$failed <- TRUE
  }
  params <- add_to_log(params, "CheckStrata.a2", read_time, read_size, 0, 0)
  return(params)
}


send_strata_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "send_strata_cox_b2\n\n")
  strata_temp <- data$strata
  write_time <- proc.time()[3]
  save(strata_temp, file = file.path(params$write_path, "strata.rdata"))
  write_size <- file.size(file.path(params$write_path, "strata.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "send_strata_cox_b2", 0, 0,
                       write_time, write_size)
  return(params)
}


prepare_strata_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_strata_cox_a2\n\n")
  # I will need to update both data$survival and params$logs, etc in this
  # function.
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  strata_temp <- NULL

  if (params$get_strata_from_b) {
    read_time <- proc.time()[3]
    load(file = file.path(params$read_path, "strata.rdata"))
    read_size <- file.size(file.path(params$read_path, "strata.rdata"))
    read_time <- proc.time()[3] - read_time
  }
  if (length(data$strata$strata_from_a) == 0 &&
      length(data$strata$strata_from_b) == 0) {
    strata_temp$x <- data.frame(const__ = rep(1, params$n))
    strata_temp$legend <- FALSE
  } else if (length(data$strata$strata_from_a) > 0 &&
             length(data$strata$strata_from_b) == 0) {
    strata_temp$x <- data$strata$x
    strata_temp$legend <- data$strata$legend
  } else if (length(data$strata$strata_from_a) > 0 &&
             length(data$strata$strata_from_b) > 0) {
    strata_temp$x <- cbind(data$strata$x, strata_temp$x)
    strata_temp$legend <- c(data$strata$legend, strata_temp$legend)
  }

  sorted <- do.call("order", cbind(strata_temp$x,
                                   data$survival$rank,
                                   data$survival$status))
  strata_temp$x <- strata_temp$x[sorted, , drop = FALSE]
  data$survival$rank   <- data$survival$rank[sorted]
  data$survival$status <- data$survival$status[sorted]
  data$x <- data$x[sorted, , drop = FALSE]
  data$survival$sorted <- sorted
  ranks <- which(apply(abs(apply(strata_temp$x, 2, diff)), 1, sum) > 0)
  ranks <- c(ranks, nrow(strata_temp$x))
  names(ranks) <- NULL
  strata <- rep(list(list()), length(ranks))
  if (length(ranks) == 1 && colnames(strata_temp$x)[1] == "const__") {
    strata[[1]]$start <- 1
    strata[[1]]$end   <- as.integer(nrow(data$x))
    strata[[1]]$label <- ""
  } else {
    start <- 1
    for (i in  seq_along(ranks)) {
      strata[[i]]$start <- start
      strata[[i]]$end   <- as.integer(ranks[i])
      label <- ""
      for (j in seq_len(ncol(strata_temp$x))) {
        temp  <- colnames(strata_temp$x)[j]
        label <- paste0(label, temp, "=",
                        strata_temp$legend[[temp]][strata_temp$x[start, j]])
        if (j < ncol(strata_temp$x)) {
          label <- paste0(label, ", ")
        }
      }
      strata[[i]]$label <- label
      start <- as.numeric(ranks[i]) + 1
    }
  }

  empty_strata <- c()
  data$x <- cbind(matrix(0, nrow = nrow(data$x), ncol = length(strata)), data$x)
  for (i in seq_along(strata)) {
    idx <- strata[[i]]$start:strata[[i]]$end
    data$x[idx, i] <- 1
    temp <- table(data$survival$rank[idx])
    m <- length(temp)
    # number of unique observed times, including where no one fails
    # Count the number of 0's and 1's for each observed time
    temp0 <- table(data$survival$rank[idx], data$survival$status[idx])
    # Check if there are all 1's or all 0's .  If so, add them into the table.
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 <- cbind(temp0, 0)
        colnames(temp0) <- c("0", "1")
      } else {
        temp0 <- cbind(0, temp0)
        colnames(temp0) <- c("0", "1")
      }
    }
    # number of distinct failure times
    strata[[i]]$J <- as.integer(sum(temp0[, 2] > 0))
    if (strata[[i]]$J == 0) {
      empty_strata <- c(empty_strata, i)
    }
    # The number of failures at each rank which has a failure
    strata[[i]]$nfails <- as.numeric(temp0[which(temp0[, 2] > 0), 2])
    # The first index of the ranks for which the number of failures is > 0.
    strata[[i]]$start0 <- c(1, (cumsum(temp)[1:(m - 1)] +
                                  1))[which(temp0[, 2] > 0)]
    # The first index of a failure for each rank which has a failure
    strata[[i]]$start1 <- strata[[i]]$start0 + temp0[which(temp0[, 2] > 0), 1]
    # The last index of a failure for each rank which has a failure
    strata[[i]]$stop1  <- as.numeric(strata[[i]]$start1 +
                                       strata[[i]]$nfails - 1)
    strata[[i]]$start0 <- as.numeric(strata[[i]]$start0 +
                                       strata[[i]]$start - 1)
    strata[[i]]$start1 <- as.numeric(strata[[i]]$start1 +
                                       strata[[i]]$start - 1)
    strata[[i]]$stop1  <- as.numeric(strata[[i]]$stop1  +
                                       strata[[i]]$start - 1)
  }

  data$survival$strata <- strata

  write_time <- proc.time()[3]
  survival <- data$survival
  save(survival, file = file.path(params$write_path, "survival.rdata"))
  write_size <- file.size(file.path(params$write_path, "survival.rdata"))
  write_time <- proc.time()[3] - write_time
  data$read_time <- read_time
  data$read_size <- read_size
  data$write_time <- write_time
  data$write_size <- write_size

  return(data)
}


prepare_params_cox_a2 <- function(params,
                                  data,
                                  cutoff = 0.01,
                                  max_iterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_a2\n\n")
  params$converged       <- FALSE
  params$halted          <- FALSE
  params$singular_matrix <- FALSE
  params$pmn_step_counter  <- 1
  pb <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pb.rdata")) # load pb
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
  params$num_events <- sum(data$survival$status)
  if (pb$n != params$n) {
    params$error_message <-
      paste("Party A has", params$n, "observations and Party B has",
            pb$n, "observations.")
    warning(params$error_message)
    params$failed <- TRUE
    return(params)
  }

  params$p1 <- ncol(data$x)
  params$p2 <- pb$p2
  params$p  <- params$p1 + params$p2
  params$p1_old <- params$p1
  params$p2_old <- params$p2

  params$a_col_names <- colnames(data$x)
  params$b_col_names <- pb$b_col_names
  params$a_col_names_old <- c("")
  params$b_col_names_old <- c("")

  params$a_tags         <- data$tags
  params$b_tags         <- pb$tags

  if (cutoff <= 0) cutoff <- 0.01
  if (cutoff >= 1) cutoff <- 0.05
  params$cutoff        <- cutoff

  if (max_iterations < 1) max_iterations <- 1
  params$max_iterations <- max_iterations

  params$survival_installed <- requireNamespace("survival", quietly = TRUE)
  if (params$survival_installed && !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  pa <- list()
  pa$p1 <- params$p1
  pa$cutoff <- params$cutoff
  pa$max_iterations <- params$max_iterations
  pa$a_col_names <- params$a_col_names
  pa$tags <- data$tags

  write_time <- proc.time()[3]
  save(pa, file <- file.path(params$write_path, "pa.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_cox_a2", read_time, read_size,
                       write_time, write_size)
  return(params)
}


prepare_blocks_cox_a2 <- function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_blocks_cox_a2\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  <- params$n
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
    params <- add_to_log(params, "prepare_blocks_cox_a2", 0, 0, 0, 0)
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
  params <- add_to_log(params, "prepare_blocks_cox_a2", 0, 0,
                       write_time, write_size)
  return(params)
}


get_z_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_z_cox_a2\n\n")
  write_time <- 0
  write_size <- 0

  num_blocks <- params$blocks$num_blocks
  pbar <- make_progress_bar_1(num_blocks, "z Files", params$verbose)
  container_ct_z <- 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename <- paste0("cz_", container_ct_z, ".rdata")
      to_write <- file(file.path(params$write_path, filename), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    g <- params$blocks$g[i]
    z <- FindOrthogonalVectors(data$x[strt:stp, ], g)
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
  params <- add_to_log(params, "get_z_cox_a2", 0, 0, write_time, write_size)
  return(params)
}


sort_data_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "sort_data_cox_b2\n\n")
  survival <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "survival.rdata"))
  read_size <- file.size(file.path(params$read_path, "survival.rdata"))
  read_time <- proc.time()[3] - read_time
  data$x <- data$x[survival$sorted, , drop = FALSE]
  data$survival <- survival
  data$read_time <- read_time
  data$read_size <- read_size
  return(data)
}


finalize_params_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "finalize_params_cox_b2\n\n")
  pa <- NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # read pa
  read_size <- sum(file.size(file.path(params$read_path, "pa.rdata")))
  read_time             <- proc.time()[3] - read_time
  params$num_events     <- sum(data$survival$status)
  params$p1            <- pa$p1
  params$cutoff        <- pa$cutoff
  params$max_iterations <- pa$max_iterations
  params$p             <- params$p1 + params$p2
  params$a_col_names     <- pa$a_col_names
  params <- add_to_log(params, "finalize_params_cox_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


prepare_blocks_cox_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_blocks_cox_b2\n\n")
  blocksize <- NULL
  # For now, assuming that p1 > 0 and p2 > 0
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "blocksize.rdata")) # load blocksize
  read_size <- file.size(file.path(params$read_path, "blocksize.rdata"))
  read_time <- proc.time()[3] - read_time
  params$blocks    <- create_blocks(params$p1, params$p2, params$n, blocksize)
  params$container <- create_containers(params$p1, params$p2, params$blocks)
  params <- add_to_log(params, "prepare_blocks_cox_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


get_w_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_w_cox_b2\n\n")
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks,
                              "(I-z*z')x", params$verbose)

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

    w <- data$x[strt:stp, ] - z %*% (t(z) %*% data$x[strt:stp, ])

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(w), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z ||
        i == params$blocks$num_blocks) {
      close(to_read)
    }
    if ((i + 1) %in% params$container$filebreak.w ||
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

  params <- add_to_log(params, "get_w_cox_b2",
                       read_time, read_size, write_time, write_size)

  return(params)
}


check_colinearity_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "check_colinearity_cox_a2\n\n")
  p2 <- params$p2
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  xb_t_xb <- NULL

  read_time <- read_time - proc.time()[3]
  load(file.path(params$read_path, "xbtxb.rdata")) # load xb_t_xb
  read_size <- file.size(file.path(params$read_path, "xbtxb.rdata"))
  read_time <- read_time + proc.time()[3]
  xa_t_xa <- t(data$x) %*% data$x
  xa_t_xb <- 0

  pbar <- make_progress_bar_1(params$blocks$num_blocks, "X'X", params$verbose)

  container_ct_w <- 0
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.w) {
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

    if ((i + 1) %in% params$container$filebreak.w ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename))
    }
    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }

  xtx <- rbind(cbind(xa_t_xa, xa_t_xb), cbind(t(xa_t_xb), xb_t_xb))

  nrow <- nrow(xtx)
  indicies <- seq_along(data$survival$strata)
  for (i in (1 + length(data$survival$strata)):nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }

  num_strata <- length(data$survival$strata)
  a_names <- params$a_col_names
  b_names <- params$b_col_names
  # Get rid of the strata indicators
  indicies <- indicies[-(seq_along(data$survival$strata))] - num_strata
  a_index <- which(indicies < length(a_names))
  b_index <- which(indicies > length(a_names))
  params$a_indicies_keep <- indicies[a_index]
  params$b_indicies_keep <- indicies[b_index] - length(a_names)
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

  a_indicies <- params$a_indicies_keep
  b_indicies <- params$b_indicies_keep


  write_time <- write_time - proc.time()[3]
  save(a_indicies, b_indicies,
       file = file.path(params$write_path, "indicies.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "indicies.rdata")))
  write_time <- write_time + proc.time()[3]
  tags <- params$b_tags[b_indicies]

  if (length(unique(tags)) == 0) {
    params$failed <- TRUE
    params$error_message <-
      "After removing colinear covariates, Party B has no covariates."
  }
  params <- add_to_log(params, "check_colinearity_cox_a2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


update_data_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_cox_a2\n\n")
  num_strata <- length(data$survival$strata)
  data$x <- as.matrix(data$x[, params$a_indicies_keep + num_strata,
                             drop = FALSE])
  return(data)
}


prepare_loop_cox_a2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_loop_cox_a2\n\n")
  params$betas_a     <- matrix(0, params$p1, 1)
  params$betas_a_old  <- matrix(0, params$p1, 1)
  params$alg_iteration_counter      <- 1
  params$delta_beta <- Inf
  params <- add_to_log(params, "prepare_loop_cox_a2", 0, 0, 0, 0)
  return(params)
}


compute_log_likelihood_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_log_likelihood_cox_a2\n\n")
  n <- params$n
  p1 <- params$p1
  read_time <- 0
  read_size <- 0
  xb_betas_b <- NULL

  if (params$alg_iteration_counter == 1) {
    params$x_betas_old  <- matrix(0, n, 1)
    x_betas  <- matrix(0, n, 1)
    params$loglikelihood_old <- -Inf
  } else {
    read_time <- -proc.time()[3]
    # load xb_betas_b
    load(file.path(params$read_path, "XB_betasB.rdata"))
    read_size <- file.size(file.path(params$read_path, "XB_betasB.rdata"))
    read_time <- read_time + proc.time()[3]
    params$x_betas_old <- params$x_betas
    x_betas <- params$xa_betas_a + xb_betas_b
    params$loglikelihood_old <- params$loglikelihood
  }

  num_events <- sum(data$survival$status)

  step_size <- 1
  w <- exp(x_betas)

  while (max(w) == Inf) {
    x_betas <- (x_betas + params$x_betas_old) * 0.5
    step_size <- step_size * 0.5
    w <- exp(x_betas)
  }

  compute_log_likelihood <- TRUE

  while (compute_log_likelihood) {
    step_counter <- 0
    pbar <- make_progress_bar_1(num_events, "Loglikelihood", params$verbose)
    loglikelihood <- 0
    for (i in seq_along(data$survival$strata)) {
      if (data$survival$strata[[i]]$J > 0) {
        for (j in 1:data$survival$strata[[i]]$J) {
          nj <- data$survival$strata[[i]]$nfails[j]
          y_start <- data$survival$strata[[i]]$start0[j]
          y_end   <- data$survival$strata[[i]]$end
          y_index <- y_start:y_end
          z_start <- data$survival$strata[[i]]$start1[j]
          z_end   <- data$survival$strata[[i]]$stop1[j]
          z_index <- z_start:z_end
          a_j1 <- sum(w[y_index])
          a_j2 <- sum(w[z_index]) / nj
          loglikelihood <- loglikelihood + sum(log(w[z_index]))
          for (r in 0:(nj - 1)) {
            a_jr <- a_j1 - r * a_j2
            loglikelihood <- loglikelihood - log(a_jr)
          }
          step_counter <- step_counter + nj
          pbar <- make_progress_bar_2(step_counter, pbar, params$verbose)
        }
      }
    }
    if (loglikelihood > params$loglikelihood_old || step_size < 0.5^6) {
      compute_log_likelihood <- FALSE
    } else {
      if (params$verbose) cat("Step Halving\n\n")
      x_betas <- (x_betas + params$x_betas_old) * 0.5
      step_size <- step_size * 0.5
      w <- exp(x_betas)
    }
  }

  deltal <- as.numeric(data$survival$status)
  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  deltal[1] <- deltal[1]
  w_xa  <- matrix(0, n, p1)
  step_counter <- 0

  .Call("compute_cox", data$survival$strata, data$x, w, deltal, w_xa,
        as.integer(n), as.integer(p1), as.integer(num_events),
        as.integer(params$verbose))

  params$loglikelihood <- loglikelihood
  params$deltal <- deltal
  params$t_xa_w_xa <- t(data$x) %*% w_xa
  params$t_xa_delta_l <- t(data$x) %*% deltal

  if (params$alg_iteration_counter == 1) {
    params$nullloglikelihood <- loglikelihood
    params$null_score <- params$t_xa_delta_l
  }
  params$x_betas <- x_betas
  params$step_size <- step_size

  write_time <- proc.time()[3]
  save(x_betas, step_size, file = file.path(params$write_path,
                                            "Xbetas_ss.rdata"))
  write_size <- sum(file.size(file.path(params$write_path,
                                        "Xbetas_ss.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "compute_log_likelihood_cox_a2",
                       read_time, read_size,
                       write_time, write_size)
  return(params)
}


update_params_cox_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "update_params_cox_b2\n\n")
  a_indicies <- NULL
  b_indicies <- NULL

  read_time <- proc.time()[3]
  # load a_indicies, b_indicies
  load(file.path(params$read_path, "indicies.rdata"))
  read_size <- sum(file.size(file.path(params$read_path, "indicies.rdata")))
  read_time <- proc.time()[3] - read_time
  params$a_col_names_old <- params$a_col_names
  params$b_col_names_old <- params$b_col_names
  params$a_col_names     <- params$a_col_names_old[a_indicies]
  params$b_col_names     <- params$b_col_names_old[b_indicies]
  params$p1_old <- params$p1
  params$p2_old <- params$p2
  params$p1     <- length(a_indicies)
  params$p2     <- length(b_indicies)
  params$p_old  <- params$p - 1
  params$p      <- params$p1 + params$p2
  params$b_indicies_keep <- b_indicies
  params$a_indicies_keep <- a_indicies
  params$betas_b     <- matrix(0, params$p2, 1)
  params$betas_b_old  <- matrix(0, params$p2, 1)
  params <- add_to_log(params, "update_params_cox_b2",
                       read_time, read_size, 0, 0)
  return(params)
}

update_data_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_cox_b2\n\n")
  data$x <- as.matrix(data$x[, params$b_indicies_keep, drop = FALSE])
  return(data)
}


compute_log_likelihood_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_log_likelihood_cox_b2\n\n")
  n <- params$n
  p2  <- params$p2
  step_size <- NULL
  x_betas  <- NULL

  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0

  read_time <- read_time - proc.time()[3]
  # load x_betas, step_size
  load(file.path(params$read_path, "Xbetas_ss.rdata"))
  read_size <- sum(file.size(file.path(params$read_path, "Xbetas_ss.rdata")))
  read_time <- read_time + proc.time()[3]

  params$step_size <- step_size

  w <- exp(x_betas)
  deltal <- as.numeric(data$survival$status)
  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  deltal[1] <- deltal[1]
  num_events <- sum(data$survival$status)
  w_xb  <- matrix(0, n, p2)

  .Call("compute_cox", data$survival$strata, data$x, w, deltal, w_xb,
        as.integer(n), as.integer(p2), as.integer(num_events),
        as.integer(params$verbose))

  pbar <- make_progress_bar_1(params$blocks$num_blocks, "(I-z*z')WX",
                              params$verbose)
  container_ct_z <- 0
  container_ct_cox <- 0

  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$file_break_z) {
      container_ct_z <- container_ct_z + 1
      filename1 <- paste0("cz_", container_ct_z, ".rdata")
      to_read <- file(file.path(params$read_path, filename1), "rb")
    }
    if (i %in% params$container$filebreak.Cox) {
      container_ct_cox <- container_ct_cox + 1
      filename2 <- paste0("cCox_", container_ct_cox, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1
    g <- params$blocks$g[i]

    read_time <- read_time - proc.time()[3]
    z <- matrix(readBin(con = to_read, what = numeric(), n = n2 * g,
                        endian = "little"), nrow = n2, ncol = g)
    read_time <- read_time + proc.time()[3]

    iz_tz_w_xbtemp <- w_xb[strt:stp, , drop = FALSE] -
      z %*% (t(z) %*% w_xb[strt:stp, , drop = FALSE])

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(iz_tz_w_xbtemp), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$file_break_z ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename1))
    }
    if ((i + 1) %in% params$container$filebreak.Cox ||
        i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size +
        file.size(file.path(params$write_path, filename2))
    }

    pbar <- make_progress_bar_2(i, pbar, params$verbose)
  }
  params$deltal <- deltal
  params$txb_w_xb   <- t(data$x) %*% w_xb
  params$t_xb_delta_l <- t(data$x) %*% deltal
  if (params$alg_iteration_counter == 1) {
    params$null_score <- params$t_xb_delta_l
  }

  txb_w_xb <- params$txb_w_xb
  write_time <- write_time - proc.time()[3]
  save(txb_w_xb, file = file.path(params$write_path,
                                  "txb_w_xb.rdata"))
  write_size <- write_size + file.size(file.path(params$write_path,
                                                 "txb_w_xb.rdata"))
  write_time <- write_time + proc.time()[3]
  params <- add_to_log(params, "compute_log_likelihood_cox_b2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


compute_inverse_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_inverse_cox_a2\n\n")
  p1 <- params$p1
  p2 <- params$p2
  read_time <- 0
  read_size <- 0
  write_time <- 0
  write_size <- 0
  txb_w_xb <- 0
  txa_w_xb <- 0

  read_time <- read_time - proc.time()[3]
  # load txb_w_xb
  load(file.path(params$read_path, "txb_w_xb.rdata"))
  read_size <- read_size + file.size(file.path(params$read_path,
                                               "txb_w_xb.rdata"))
  read_time <- read_time + proc.time()[3]

  container_ct_cox <- 0

  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.Cox) {
      container_ct_cox <- container_ct_cox + 1
      filename <- paste0("cCox_", container_ct_cox, ".rdata")
      to_read <- file(file.path(params$read_path, filename), "rb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1

    read_time <- read_time - proc.time()[3]
    iz_tz_w_xb  <- matrix(readBin(con = to_read, what = numeric(), n = n2 * p2,
                                  endian = "little"), nrow = n2, ncol = p2)
    read_time <- read_time + proc.time()[3]

    txa_w_xb <- txa_w_xb + t(data$x[strt:stp, ]) %*% iz_tz_w_xb
    if ((i + 1) %in% params$container$filebreak.Cox ||
        i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path, filename))
    }
  }

  m <- rbind(cbind(params$t_xa_w_xa, txa_w_xb),
             cbind(t(txa_w_xb), txb_w_xb))

  if (params$alg_iteration_counter == 1) {
    params$nullHessian <- m
  }
  params$XtWX <- m


  inv <- NULL
  tryCatch({
    inv <- solve(m)
  },
  error = function(err) {
    inv <- NULL
  })
  m <- inv
  if (is.null(m)) {
    params$failed <- TRUE
    params$error_message <-
      paste("The matrix XWX is singular.",
            "This is probably due to divergence of the coefficients.")
    warning(params$error_message)

    betas <- rep(NA, length(params$a_col_names_old))
    betas[params$a_indicies_keep] <- params$beta_a
    betas <- data.frame(betas)
    rownames(betas) <- params$a_col_names_old
    params <- add_to_log(params, "compute_inverse_cox_a2", read_time, read_size, write_time, write_size)
    return(params)
  }
  m3 <- m[(p1 + 1):(p1 + p2), 1:p1]
  params$m <- m
  params$m3_txa_delta_l <- m3 %*% params$t_xa_delta_l
  m3_txa_delta_l <- params$m3_txa_delta_l
  write_time <- write_time - proc.time()[3]
  save(m, m3_txa_delta_l, file = file.path(params$write_path, "M.rdata"))
  write_size <- write_size + sum(file.size(file.path(params$write_path, "M.rdata")))
  write_time <- write_time + proc.time()[3]
  params <- add_to_log(params, "compute_inverse_cox_a2", read_time, read_size, write_time, write_size)

  return(params)
}


compute_beta_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_beta_cox_b2\n\n")
  p1 <- params$p1
  p2 <- params$p2
  m  <- 0
  m3_txa_delta_l <- 0

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "M.rdata")) # load m, m3_txa_delta_l
  read_size <- sum(file.size(file.path(params$read_path, "M.rdata")))
  read_time <- proc.time()[3] - read_time

  if (params$step_size < 1) {  # Is this in the wrong spot?
    params$betas_b <- params$betas_b_old + (params$betas_b - params$betas_b_old) * params$step_size
  }

  params$betas_b_old <- params$betas_b
  m4 <- m[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
  m2 <- m[1:p1, (p1 + 1):(p1 + p2)]

  params$m2_txb_deta_l <- m2 %*% params$t_xb_delta_l
  params$betas_b <- params$betas_b + m3_txa_delta_l + m4 %*% params$t_xb_delta_l

  params$xb_betas_b <- data$x %*% params$betas_b

  m2_txb_deta_l <- params$m2_txb_deta_l
  xb_betas_b <- params$xb_betas_b
  write_time <- proc.time()[3]
  save(m2_txb_deta_l, file = file.path(params$write_path, "M2_tXB_deltal.rdata"))
  save(xb_betas_b, file = file.path(params$write_path, "XB_betasB.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("M2_tXB_deltal.rdata",
                                                             "XB_betasB.rdata"))))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "compute_beta_cox_b2", read_time, read_size, write_time, write_size)

  return(params)
}


compute_beta_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_beta_cox_a2\n\n")
  p1 <- params$p1
  m2_txb_deta_l <- 0

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "M2_tXB_deltal.rdata")) # load m2_txb_deta_l
  read_size <- sum(file.size(file.path(params$read_path, "M2_tXB_deltal.rdata")))
  read_time <- proc.time()[3] - read_time

  if (params$step_size < 1) { # Is this in the wrong spot?
    params$betas_a <- params$betas_a_old + (params$betas_a - params$betas_a_old) * params$step_size
  }

  params$betas_a_old <- params$betas_a
  params$betas_a <- params$betas_a + params$m[1:p1, 1:p1] %*% params$t_xa_delta_l + m2_txb_deta_l

  converged <- abs(params$loglikelihood - params$loglikelihood_old) /
    (abs(params$loglikelihood) + 0.1) < params$cutoff
  params$converged <- converged

  params$xa_betas_a <- data$x %*% params$betas_a

  if (params$alg_iteration_counter >= params$max_iterations) {
    params$halted <- TRUE
    warning(paste("Failed to converged in", params$max_iterations, "iterations."))
  }

  write_time <- proc.time()[3]
  save(converged, file = file.path(params$write_path, "converged.rdata"))
  write_size <- file.size(file.path(params$write_path, "converged.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_beta_cox_a2", read_time, read_size, write_time, write_size)
  return(params)
}


get_converged_status_cox_b2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_converged_status_cox_b2\n\n")
  converged <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "converged.rdata")) # load deltabeta
  read_size <- file.size(file.path(params$read_path, "converged.rdata"))
  read_time <- proc.time()[3] - read_time
  params$converged <- converged
  if (params$alg_iteration_counter > params$max_iterations) {
    params$halted <- TRUE
    warning(paste("Failed to converged in", params$max_iterations, "iterations."))
  }
  write_time <- 0
  write_size <- 0
  if (params$converged || params$halted) {
    betas_b <- params$betas_b
    null_score_b <- params$null_score
    write_time <- proc.time()[3]
    save(betas_b, null_score_b, file = file.path(params$write_path, "B_betas_ns.rdata"))
    write_size <- sum(file.size(file.path(params$write_path, "B_betas_ns.rdata")))
    write_time <- proc.time()[3] - write_time
  }
  params <- add_to_log(params, "get_converged_status_cox_b2", read_time, read_size, write_time, write_size)
  return(params)
}


survfit_cox_a2 <- function(params, survival, pred) {
  if (params$trace) cat(as.character(Sys.time()), "survfit_cox_a2\n\n")
  surv <- rep(1, length(survival$rank))
  for (i in seq_along(survival$strata)) {
    if (survival$strata[[i]]$J > 0) {
      start   <- survival$strata[[i]]$start
      end     <- survival$strata[[i]]$end
      risk    <- exp(pred[start:end])
      dtime   <- survival$rank[start:end]
      status  <- survival$status[start:end]
      death   <- status == 1                                  # times where death happened
      time    <- sort(unique(dtime))                          # (A) get unique event times
      rcumsum <- function(x) rev(cumsum(rev(x)))
      nevent  <- as.vector(rowsum(as.numeric(death), dtime))  # (A) Count the number of deaths at each event time
      ndeath  <- rowsum(status, dtime)                        # (A) number of deaths at each unique event time
      nrisk   <- rcumsum(rowsum(risk, dtime))                 # (A) rowsum = sum of risk at each time, then reverse cum sum, sorted by time
      erisk   <- rowsum(risk * death, dtime)                  # (A) risk score sums of death at each unique event time
      n       <- length(nevent)
      sum1    <- double(n)   # a vector of 0's, length number of unique event times
      for (i in 1:n) {
        d <- ndeath[i]
        if (d == 1) {
          sum1[i] <- 1 / nrisk[i]
        } else if (d > 1) {
          for (j in 0:(d - 1)) {
            sum1[i] <- sum1[i] + 1 / (d * nrisk[i] - erisk[i] * j)
          }
        }
      }
      temp <- exp(-cumsum(nevent * sum1))
      for (i in start:end) {
        surv[i] <- temp[which(time == survival$rank[i])]
      }
    }
  }
  return(surv)
}

#' @importFrom  stats pchisq pnorm qnorm
compute_results_cox_a2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_results_cox_a2\n\n")
  stats <- params$stats
  stats$converged <- params$converged
  stats$failed    <- FALSE
  read_time <- 0
  read_size <- 0
  betas_b     <- NULL
  null_score_b <- NULL
  xb_betas_b  <- NULL

  read_time <- proc.time()[3]
  load(file.path(params$read_path, "B_betas_ns.rdata")) # load betas_b, nullscoreB
  load(file.path(params$read_path, "XB_betasB.rdata")) # load xb_betas_b
  read_size <- sum(file.size(file.path(params$read_path, c("B_betas_ns.rdata",
                                                           "XB_betasB.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$betas_b <- betas_b
  params$null_score <- rbind(params$null_score, null_score_b)

  names_new          <- c(params$a_col_names, params$b_col_names)
  names_old          <- c(params$a_col_names_old, params$b_col_names_old)
  idx_a               <- params$a_indicies_keep
  idx_b               <- params$b_indicies_keep
  idx                <- c(idx_a, idx_b + length(params$a_col_names_old))
  stats$party        <- c(rep("dp0", length(params$a_col_names_old)),
                          rep("dp1", length(params$b_col_names_old)))
  stats$coefficients <- rep(NA, length(stats$party))
  tempcoefs          <- c(params$betas_a, params$betas_b)
  stats$coefficients[idx] <- tempcoefs
  stats$expcoef      <- exp(stats$coefficients)  # hazard ratios
  stats$expncoef     <- exp(-stats$coefficients)
  tempvar            <- solve(params$XtWX)
  stats$var          <- matrix(0, length(names_old), length(names_old))
  stats$var[idx, idx] <- tempvar
  stats$secoef       <- rep(NA, length(names_old))
  stats$secoef[idx]  <- sqrt(diag(tempvar))  # standard error

  stats$zvals        <- stats$coefficients / stats$secoef  # z values
  stats$pvals        <- 2 * pnorm(abs(stats$zvals), lower.tail = FALSE)   # pvals
  stats$stars        <- matrix(sapply(stats$pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  }))
  stats$lower95      <- exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      <- exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  stats$loglik       <- c(params$nullLoglikelihood, params$loglikelihood)
  stats$n            <- params$n
  stats$nevent       <- params$num_events
  stats$iter         <- params$alg_iteration_counter - 1
  stats$df           <- params$p
  stats$score        <- t(params$null_score) %*%
    solve(params$nullHessian) %*%
    params$null_score
  stats$score        <- c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       <- "efron"
  stats$lrt          <- 2 * (stats$loglik[2] - stats$loglik[1])
  stats$lrt          <- c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      <- c(1 - exp(-stats$lrt[1] / stats$n),
                          1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    <- t(tempcoefs) %*% params$XtWX %*% tempcoefs
  stats$wald.test    <- c(stats$wald.test,
                          1 - pchisq(stats$wald.test, stats$df))
  pred <- -params$xa_betas_a - xb_betas_b
  if (params$survival_installed) {
    surv <- survival::Surv(data$survival$rank, data$survival$status)
    strat <- rep(0, length(surv))
    for (i in seq_along(data$survival$strata)) {
      strat[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] <- i
    }
    results <- survival::concordance(surv ~ pred + strata(strat))
    if (is.matrix(results$stats)) {
      # more than one strata
      stats$concordance <- c(apply(results$count, 2, sum)[1:4], results$concordance, sqrt(results$var))
    } else {
      # only one strata, so a numeric vector
      stats$concordance <- c(results$count[1:4], results$concordance, sqrt(results$var))
    }
  } else {
    stats$concordance <- c(NA, NA, NA, NA, NA, NA)
  }

  stats$survival <- data.frame(
    rank   = data$survival$rank,
    status = data$survival$status,
    sorted = data$survival$sorted,
    surv   = survfit_cox_a2(params, data$survival, pred)
  )

  stats$strata <- as.data.frame(matrix(0, length(data$survival$strata), 3))
  stats$strata$label <- ""
  colnames(stats$strata) <- c("start", "end", "events", "label")
  for (i in seq_along(data$survival$strata)) {
    stats$strata$start[i]  <- data$survival$strata[[i]]$start
    stats$strata$end[i]    <- data$survival$strata[[i]]$end
    stats$strata$events[i] <- sum(data$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  <- data$survival$strata[[i]]$label
  }

  names(stats$party)           <- names_old
  names(stats$coefficients)    <- names_old
  names(stats$expcoef)         <- names_old
  names(stats$expncoef)        <- names_old
  rownames(stats$var)          <- names_old
  colnames(stats$var)          <- names_old
  names(stats$secoef)          <- names_old
  names(stats$zvals)           <- names_old
  names(stats$pvals)           <- names_old
  names(stats$stars)           <- names_old
  names(stats$lower95)         <- names_old
  names(stats$upper95)         <- names_old
  names(stats$loglik)          <- c("loglikelihood", "null loglikelihood")
  names(stats$score)           <- c("score", "p-value")
  names(stats$lrt)             <- c("likelihood ratio", "p-value")
  names(stats$rsquare)         <- c("r-square", "max possible")
  names(stats$wald.test)       <- c("wald", "p-value")
  names(stats$concordance)     <- c("concordant", "discordant",
                                    "tied.risk", "tied.time",
                                    "concordance", "stderr")

  params$stats <- stats
  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_results_cox_a2",
                       read_time, read_size, write_time, write_size)
  return(params)
}


get_results_cox_b2 <- function(params) {
  stats <- NULL
  if (params$trace) cat(as.character(Sys.time()), "get_results_cox_b2\n\n")
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size <- file.size(file.path(params$read_path, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats <- stats
  params <- add_to_log(params, "get_results_cox_b2", read_time, read_size, 0, 0)
  return(params)
}


####################### REGRESSION BY B ONLY FUNCTIONS #######################

finalize_params_2_cox_b2 <- function(params, data) {
  pa <- NULL
  if (params$trace) cat(as.character(Sys.time()),
                        "finalize_params_2_cox_b2\n\n")
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "pa.rdata")) # read pa
  read_size <- file.size(file.path(params$read_path, "pa.rdata"))
  read_time <- proc.time()[3] - read_time
  params$num_events <- sum(data$survival$status)
  params$p1 <- pa$p1
  params$cutoff <- pa$cutoff
  params$max_iterations <- pa$max_iterations
  params$p <- params$p1 + params$p2
  params <- add_to_log(params, "finalize_params_2_cox_b2",
                       read_time, read_size, 0, 0)
  return(params)
}


check_colinearity_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "check_colinearity_cox_b2\n\n")
  # Add in the strata here
  num_strata <- length(data$survival$strata)
  x <- cbind(matrix(0, nrow = nrow(data$x), ncol = num_strata), data$x)
  for (i in 1:num_strata) {
    x[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] <- 1
  }

  xtx <- t(x) %*% x
  nrow <- nrow(xtx)
  indicies <- 1:num_strata
  for (i in (num_strata + 1):nrow) {
    temp_indicies <- c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies <- c(indicies, i)
    }
  }
  # Get rid of the Strata
  indicies <- indicies[-(1:num_strata)] - num_strata

  params$a_indicies_keep <- c()
  params$a_col_names_old <- c()
  params$a_col_names     <- c()
  params$p1_old          <- params$p1
  params$p1              <- 0

  b_names                <- params$b_col_names
  b_names_keep           <- b_names[indicies]
  params$b_indicies_keep <- indicies
  params$b_col_names_old <- params$b_col_names
  params$b_col_names     <- b_names_keep
  params$p2_old          <- params$p2
  params$p2              <- length(b_names_keep)
  params$p               <- params$p1 + params$p2

  if (params$p2 == 0) {
    params$failed <- TRUE
    params$error_message <-
      "Party A has no covariates and all of Party B's covariates are linear."
    warning(params$error_message)
  }
  params <- add_to_log(params, "check_colinearity_cox_b2", 0, 0, 0, 0)

  return(params)
}

#' @importFrom stats as.formula
compute_cox_from_survival_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()),
                        "compute_cox_from_survival_b2\n\n")
  # We have loaded survival previously

  strata <- rep(0, nrow(data$x))
  for (i in seq_along(data$survival$strata)) {
    strata[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] <- i
  }

  colnames(data$x) <- paste0("V", seq_len(ncol(data$x)))
  f <- paste(c("Surv(rank, status) ~ strata(strata)",
               paste0("V", seq_len(ncol(data$x)))), collapse = " + ")

  error <- tryCatch(
    {
      fit <- survival::coxph(as.formula(f),
                             data <- data.frame(rank = data$survival$rank,
                                                status = data$survival$status,
                                                strata = strata,
                                                data$x),
                             iter.max = params$max_iterations)
    },
    error = function(e) {
      return(TRUE)
    },
    warning = function(e) {
      return(FALSE)
    }
  )

  if ((class(error) == "logical" && error)) {
    params$converged <- FALSE
    params$failed    <- TRUE
    params$error_message <- "Coxph in the survival package failed to converge."
    warning(params$error_message)
  } else {
    params$converged <- TRUE
    if (class(error) == "logical") {
      fit <- suppressWarnings(survival::coxph(
        as.formula(f),
        data <- data.frame(rank = data$survival$rank,
                           status = data$survival$status,
                           strata = strata,
                           data$x),
        iter.max = params$max_iterations))
      params$converged <- FALSE
    }
    fit$linear.predictors <- NULL
    fit$residuals <- NULL
    fit$y <- NULL
    fit$formula <- NULL
    fit$call <- NULL
    fit$assign <- NULL
    fit$terms <- NULL
    fit$means <- NULL
    params$fit <- fit
  }
  params <- add_to_log(params, "compute_cox_from_survival_b2", 0, 0, 0, 0)
  return(params)
}


compute_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_cox_b2\n\n")
  n           <- params$n
  p2          <- params$p2
  params$alg_iteration_counter <- 1
  x_betas_old  <- matrix(0, n, 1)
  x_betas      <- matrix(0, n, 1)
  betas_b       <- matrix(0, p2, 1)
  betas_b_old   <- betas_b
  loglikelihood_old <- -Inf


  while (params$alg_iteration_counter <= params$max_iterations &&
         !params$converged) {
    BeginningIteration(params)
    loglikelihood <- 0
    step_size <- 1
    w <- exp(x_betas)
    while (max(w) == Inf) {
      if (params$verbose) cat("Step Halving\n\n")
      x_betas <- (x_betas + x_betas_old) * 0.5
      step_size <- step_size * 0.5
      w <- exp(x_betas)
    }
    compute_log_likelihood <- TRUE
    while (compute_log_likelihood) {
      num_events <- sum(data$survival$status)
      step_counter <- 0
      pbar <- make_progress_bar_1(num_events, "Loglikelihood", params$verbose)
      loglikelihood <- 0
      for (i in seq_along(data$survival$strata)) {
        if (data$survival$strata[[i]]$J > 0) {
          for (j in 1:data$survival$strata[[i]]$J) {
            nj <- data$survival$strata[[i]]$nfails[j]
            y_start <- data$survival$strata[[i]]$start0[j]
            y_end   <- data$survival$strata[[i]]$end
            z_start <- data$survival$strata[[i]]$start1[j]
            z_end   <- data$survival$strata[[i]]$stop1[j]
            y_index <- y_start:y_end
            z_index <- z_start:z_end
            a_j1 <- sum(w[y_index])
            a_j2 <- sum(w[z_index]) / nj
            loglikelihood <- loglikelihood + sum(log(w[z_index]))
            for (r in 0:(nj - 1)) {
              a_jr <- a_j1 - r * a_j2
              loglikelihood <- loglikelihood - log(a_jr)
            }
            step_counter <- step_counter + nj
            pbar <- make_progress_bar_2(step_counter, pbar, params$verbose)
          }
        }
      }
      if (loglikelihood > loglikelihood_old || step_size < 0.5^6) {
        compute_log_likelihood <- FALSE
      } else {
        if (params$verbose) cat("Step Halving\n\n")
        x_betas <- (x_betas + x_betas_old) * 0.5
        step_size <- step_size * 0.5
        w <- exp(x_betas)
      }
    }
    num_events <- sum(data$survival$status)
    deltal <- as.numeric(data$survival$status)
    # This is to force R to make a copy since we are exploiting
    # a pass by reference with the C call.
    deltal[1] <- deltal[1]
    w_xb  <- matrix(0, n, p2)

    .Call("compute_cox", data$survival$strata, data$x, w, deltal, w_xb,
          as.integer(n), as.integer(p2), as.integer(num_events),
          as.integer(params$verbose))

    m <- t(data$x) %*% w_xb

    params$XtWX <- m
    if (params$alg_iteration_counter == 1) {
      params$nullHessian <- m
    }

    inv <- NULL
    tryCatch({
      inv <- solve(m)
    },
    error = function(err) {
      inv <- NULL
    })
    m <- inv
    if (is.null(m)) {
      params$failed <- TRUE
      params$error_message <-
        paste("The matrix t(X)WX is singular.",
              "This is probably due to divergence of the coefficients.")
      warning(params$error_message)

      betas <- rep(NA, length(params$b_col_names_old))
      betas[params$b_indicies_keep] <- betas_b
      betas <- data.frame(betas)
      rownames(betas) <- params$b_col_names_old
      params <- add_to_log(params, "compute_cox_b2", 0, 0, 0, 0)
      return(params)
    }
    delta_beta  <- m %*% t(data$x) %*% deltal
    betas_b     <- betas_b_old + (betas_b - betas_b_old) * step_size
    betas_b_old <- betas_b
    betas_b     <- betas_b + delta_beta
    x_betas     <- data$x %*% betas_b

    converged <- abs(loglikelihood - loglikelihood_old) /
      (abs(loglikelihood) + 0.1) < params$cutoff
    params$converged <- converged

    if (params$alg_iteration_counter == 1) {
      params$null_score         <- t(data$x) %*% deltal
      params$nullloglikelihood <- loglikelihood
    }
    loglikelihood_old <- loglikelihood
    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }
  params$loglikelihood <- loglikelihood
  params$betas_b <- betas_b
  params$x_betas <- x_betas
  if (!params$converged) {
    params$failed <- TRUE
    params$error_message <-
      "Cox failed to converge in the specified number of iterations."
    warning(params$error_message)
  }
  params <- add_to_log(params, "compute_cox_b2", 0, 0, 0, 0)
  return(params)
}

#' @importFrom  stats pchisq pnorm qnorm
compute_results_cox_b2 <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_results_cox_b2\n\n")
  stats <- params$stats
  stats$converged <- params$converged
  stats$party_name <- params$party_name
  stats$failed    <- FALSE

  fit_exists <- !is.null(params$fit)
  names_old          <- c(params$a_col_names_old, params$b_col_names_old)
  idx_a               <- params$a_indicies_keep
  idx_b               <- params$b_indicies_keep
  idx                <- c(idx_a, idx_b + length(params$a_col_names_old))
  stats$party        <- c(rep("dp0", length(params$a_col_names_old)),
                          rep("dp1", length(params$b_col_names_old)))
  stats$coefficients <- rep(NA, length(stats$party))
  if (fit_exists) {
    stats$coefficients[idx] <- params$fit$coefficients
    tempvar          <- params$fit$var
  } else {
    stats$coefficients[idx] <- params$betas_b
    tempvar            <- solve(params$XtWX)
  }
  stats$expcoef      <- exp(stats$coefficients)  # hazard ratios
  stats$expncoef     <- exp(-stats$coefficients)
  stats$var          <- matrix(0, length(names_old), length(names_old))
  stats$var[idx, idx] <- tempvar
  stats$secoef       <- rep(NA, length(names_old))
  stats$secoef[idx]  <- sqrt(diag(tempvar))  # standard error

  stats$zvals        <- stats$coefficients / stats$secoef  # z values
  stats$pvals        <- 2 * pnorm(abs(stats$zvals), lower.tail = FALSE) # pvals
  stats$stars        <- matrix(sapply(stats$pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  }))
  stats$lower95      <- exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      <- exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  if (fit_exists) {
    stats$loglik     <- params$fit$loglik
    stats$n          <- params$fit$n
    stats$nevent     <- params$fit$nevent
    stats$iter       <- params$fit$iter
    stats$score      <- params$fit$score
    stats$wald.test  <- params$fit$wald.test
    stats$concordance <- params$fit$concordance[c(2, 1, 3, 4, 6, 7)]
  } else {
    stats$loglik       <- c(params$nullLoglikelihood, params$loglikelihood)
    stats$n            <- params$n
    stats$nevent       <- params$num_events
    stats$iter         <- params$alg_iteration_counter - 1
    stats$score        <- t(params$null_score) %*% solve(params$nullHessian) %*%
      params$null_score
    stats$wald.test    <- t(params$betas_b) %*% params$XtWX %*% params$betas_b
    stats$concordance <- c(NA, NA, NA, NA, NA, NA)
  }
  stats$df           <- params$p
  stats$score        <- c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       <- "efron"
  stats$lrt          <- 2 * (stats$loglik[2] - stats$loglik[1])
  stats$lrt          <- c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      <- c(1 - exp(-stats$lrt[1] / stats$n),
                          1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    <- c(stats$wald.test,
                          1 - pchisq(stats$wald.test, stats$df))

  pred <- data$x %*% stats$coefficients[idx]
  stats$survival <- data.frame(
    rank   = data$survival$rank,
    status = data$survival$status,
    sorted = data$survival$sorted,
    surv   = survfit_cox_a2(params, data$survival, pred)
  )
  stats$strata <- as.data.frame(matrix(0, length(data$survival$strata), 3))
  stats$strata$label <- ""
  colnames(stats$strata) <- c("start", "end", "events", "label")
  for (i in seq_along(data$survival$strata)) {
    stats$strata$start[i]  <- data$survival$strata[[i]]$start
    stats$strata$end[i]    <- data$survival$strata[[i]]$end
    start                  <- stats$strata$start[i]
    end                    <- stats$strata$end[i]
    stats$strata$events[i] <- sum(data$survival$status[start:end])
    stats$strata$label[i]  <- data$survival$strata[[i]]$label
  }

  names(stats$party)           <- names_old
  names(stats$coefficients)    <- names_old
  names(stats$expcoef)         <- names_old
  names(stats$expncoef)        <- names_old
  rownames(stats$var)          <- names_old
  colnames(stats$var)          <- names_old
  names(stats$secoef)          <- names_old
  names(stats$zvals)           <- names_old
  names(stats$pvals)           <- names_old
  names(stats$stars)           <- names_old
  names(stats$lower95)         <- names_old
  names(stats$upper95)         <- names_old
  names(stats$loglik)          <- c("loglikelihood", "null loglikelihood")
  names(stats$score)           <- c("score", "p-value")
  names(stats$lrt)             <- c("likelihood ratio", "p-value")
  names(stats$rsquare)         <- c("r-square", "max possible")
  names(stats$wald.test)       <- c("wald", "p-value")
  names(stats$concordance)     <- c("concordant", "discordant",
                                    "tied.risk", "tied.time",
                                    "concordance", "stderr")

  params$stats <- stats
  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "compute_results_cox_b2", 0, 0,
                       write_time, write_size)
  return(params)
}


get_results_cox_a2 <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_cox_a2\n\n")
  read_time <- proc.time()[3]
  load(file.path(params$read_path, "stats.rdata"))
  read_size <- file.size(file.path(params$read_path, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  stats$party_name <- params$party_name
  params$stats <- stats
  params <- add_to_log(params, "get_results_cox_a2", read_time, read_size, 0, 0)
  return(params)
}


############################## PARENT FUNCTIONS ##############################


party_a_process_2_cox <- function(data,
                                  y_name          = NULL,
                                  strata         = NULL,
                                  mask           = TRUE,
                                  monitor_folder  = NULL,
                                  msreqid        = "v_default_00_000",
                                  blocksize      = 500,
                                  cutoff         = 1e-8,
                                  max_iterations  = 25,
                                  sleep_time     = 10,
                                  max_waiting_time = 24 * 60 * 60,
                                  popmednet      = TRUE,
                                  trace          = FALSE,
                                  verbose        = TRUE) {
  params <- prepare_params_2p("cox", "A", msreqid = msreqid,
                              popmednet = popmednet, trace = trace,
                              verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params   <- prepare_folder_cox_a2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_cox_23(params, data, y_name, strata, mask)

  params <- PauseContinue.2p(params, max_waiting_time)
  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$complete <- TRUE
    warning(read_error_message(params$read_path))
    params$pmn_step_counter <- 1
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (data$failed) {
    params$complete <- TRUE
    message <- "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params$pmn_step_counter <- 1
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- check_strata_cox_a2(params, data)
  if (params$get_strata_from_b) {
    files <- c()
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)
  }

  params <- prepare_params_cox_a2(params, data, cutoff, max_iterations)

  if (params$failed) {   # Check for failed from prepare_params_cox_a2()
    params$complete <- TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }


  data <- prepare_strata_cox_a2(params, data)
  params <- add_to_log(params, "prepare_strata_cox_a2", data$read_time,
                       data$read_size, data$write_time, data$write_size)

  # Check for $p1 == 0 => no covariates, only strata
  if (params$p1 == 0) {
    MakeTransferMessage(params$write_path)
    files <- c("transferControl.rdata", "pa.rdata", "survival.rdata")
    params <- send_pause_continue_2p(params, files,
                                     sleep_time, max_waiting_time)
    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$complete <- TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    params <- get_results_cox_a2(params)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- prepare_blocks_cox_a2(params, blocksize)

  if (params$failed) { # Check for failed from prepare_blocks_cox_a2()
    params$complete <- TRUE
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  params <- get_z_cox_a2(params, data)

  # This works even in the case that we send no blocks over.  Just set
  # file_break_z to c()
  files <- c("pa.rdata", "blocksize.rdata", "survival.rdata",
             seq_zw("cz_", length(params$container$file_break_z)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- check_colinearity_cox_a2(params, data)

  if (params$failed) { # Check for check_colinearity_cox_a2() failed
    params$complete <- TRUE
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files <- c("error_message.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  if (params$p1 == 0) { # No covariates left.  All colinear with Strata
    MakeTransferMessage(params$write_path)
    files <- c("transferControl.rdata", "indicies.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$complete <- TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }
    params <- get_results_cox_a2(params)
    params <- send_pause_quit_2p(params, sleep_time = sleep_time)
    SummarizeLog.2p(params)
    return(params$stats)
  }

  data <- update_data_cox_a2(params, data)
  params <- add_to_log(params, "update_data_cox_a2", 0, 0, 0, 0)
  params <- prepare_loop_cox_a2(params)

  while (!params$converged && !params$halted) {
    BeginningIteration(params)

    params <- compute_log_likelihood_cox_a2(params, data)

    files <- c("Xbetas_ss.rdata")
    if (params$alg_iteration_counter == 1) {
      files <- c("indicies.rdata", files)
    } else {
      files <- c("converged.rdata", files)
    }

    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    params <- compute_inverse_cox_a2(params, data)
    # Check for failed from compute_inverse_cox_a2()
    if (params$failed) {
      params$complete <- TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      SummarizeLog.2p(params)
      return(params$stats)
    }

    files <- c("M.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    params <- compute_beta_cox_a2(params, data)

    EndingIteration(params)
    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }
  params$complete <- TRUE
  files <- c("converged.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- compute_results_cox_a2(params, data)

  files <- c("stats.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time = sleep_time)
  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  SummarizeLog.2p(params)
  return(params$stats)
}


party_b_process_2_cox <- function(data,
                                  strata              = NULL,
                                  mask                = TRUE,
                                  monitor_folder      = NULL,
                                  sleep_time          = 10,
                                  max_waiting_time    = 24 * 60 * 60,
                                  popmednet           = TRUE,
                                  trace               = FALSE,
                                  verbose             = TRUE) {
  params <- prepare_params_2p("cox", "B", popmednet = popmednet, trace = trace,
                              verbose = verbose)
  params <- initialize_log_2p(params)
  params <- initialize_time_stamps_2p(params)
  params <- initialize_tracking_table_2p(params)
  header(params)
  params <- prepare_folder_cox_b2(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data <- prepare_data_cox_23(params, data, NULL, strata, mask)

  if (data$failed) { # Check for Error from prepare_data_cox_b2()
    params$complete <- TRUE
    message <- "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files <- c("error_message.rdata")
    params <- send_pause_quit_2p(params, files, sleep_time = sleep_time,
                                 job_failed = TRUE)
    return(params$stats)
  }

  params <- prepare_params_cox_b2(params, data)
  files <- c("pb.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
    params$complete <- TRUE
    warning(read_error_message(params$read_path))
    params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                 job_failed = TRUE)
    return(params$stats)
  }

  if (!is.null(data$strata$strata_from_b)) {
    params <- send_strata_cox_b2(params, data)
    files <- c("strata.rdata")
    params <- send_pause_continue_2p(params, files, sleep_time,
                                     max_waiting_time)

    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$complete <- TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }
  }

  data <- sort_data_cox_b2(params, data)
  params <- add_to_log(params, "sort_data_cox_b2", data$read_time,
                       data$read_size, 0, 0)

  if (file.exists(file.path(params$read_path, "transferControl.rdata"))) {
    params <- finalize_params_2_cox_b2(params, data)
    params <- check_colinearity_cox_b2(params, data)

    if (params$failed) {  # Happens if pb_new == 0
      params$complete <- TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_quit_2p(params, files, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }
    data <- update_data_cox_b2(params, data)
    params <- add_to_log(params, "update_data_cox_b2", 0, 0, 0, 0)
    if (params$survival_installed) {
      params <- compute_cox_from_survival_b2(params, data)
    } else {
      params <- compute_cox_b2(params, data)
    }
    # We could get a job_failed here from coefficient explosion
    if (params$failed) {
      params$complete <- TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_quit_2p(params, files, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }
    params <- compute_results_cox_b2(params, data)
    stats <- params$stats
    save(stats, file = file.path(params$write_path, "stats.rdata"))
    files <- c("stats.rdata")
    params <- send_pause_quit_2p(params, files, sleep_time = sleep_time)
    return(params$stats)
  }

  params <- finalize_params_cox_b2(params, data)
  params <- prepare_blocks_cox_b2(params)
  params <- get_w_cox_b2(params, data)

  files <- c("xbtxb.rdata", seq_zw("cw_", length(params$container$filebreak.w)))
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  if (file.exists(file.path(params$read_path, "transferControl.rdata"))) {
    params <- update_params_cox_b2(params)
    data <- update_data_cox_b2(params, data)
    params <- add_to_log(params, "update_data_cox_b2", 0, 0, 0, 0)
    if (params$survival_installed) {
      params <- compute_cox_from_survival_b2(params, data)
    } else {
      params <- compute_cox_b2(params, data)
    }

    # We could get a job_failed here from coefficient explosion
    if (params$failed) {
      params$complete <- TRUE
      make_error_message(params$write_path, params$error_message)
      files <- c("error_message.rdata")
      params <- send_pause_quit_2p(params, files, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }
    params <- compute_results_cox_b2(params, data)
    stats <- params$stats
    save(stats, file = file.path(params$write_path, "stats.rdata"))
    files <- c("stats.rdata")
    params <- send_pause_quit_2p(params, files, sleep_time = sleep_time)
    return(params$stats)
  }

  params$alg_iteration_counter <- 1
  repeat {
    if (params$alg_iteration_counter == 1) {
      if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
        params$complete <- TRUE
        warning(read_error_message(params$read_path))
        params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                     job_failed = TRUE)
        return(params$stats)
      }
      params <- update_params_cox_b2(params)
      data <- update_data_cox_b2(params, data)
      params <- add_to_log(params, "update_data_cox_b2", 0, 0, 0, 0)
    } else {
      params <- get_converged_status_cox_b2(params)
      if (params$converged || params$halted) {
        break
      }
    }

    BeginningIteration(params)

    params <- compute_log_likelihood_cox_b2(params, data)

    files <- c("txb_w_xb.rdata",
               seq_zw("cCox_", length(params$container$filebreak.Cox)))
    params <- send_pause_continue_2p(params, files,
                                     sleep_time, max_waiting_time)

    if (file.exists(file.path(params$read_path, "error_message.rdata"))) {
      params$complete <- TRUE
      warning(read_error_message(params$read_path))
      params <- send_pause_quit_2p(params, sleep_time = sleep_time,
                                   job_failed = TRUE)
      return(params$stats)
    }

    params <- compute_beta_cox_b2(params, data)

    files <- c("XB_betasB.rdata", "M2_tXB_deltal.rdata")
    params <- send_pause_continue_2p(params, files,
                                     sleep_time, max_waiting_time)

    EndingIteration(params)

    params$alg_iteration_counter <- params$alg_iteration_counter + 1
  }
  params$complete <- TRUE

  files <- c("B_betas_ns.rdata")
  params <- send_pause_continue_2p(params, files, sleep_time, max_waiting_time)

  params <- get_results_cox_b2(params)

  params <- send_pause_quit_2p(params, sleep_time = sleep_time)
  return(params$stats)
}
