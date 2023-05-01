################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

PrepareFolder.ACDP <- function(params, monitor_folder) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolder.ACDP\n\n")
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

  params$dplocalPath     <- file.path(monitor_folder, "dplocal")
  params$r_programs_path <- file.path(monitor_folder, "rprograms")
  params$macros_path     <- file.path(monitor_folder, "macros")
  if (params$data_partner_id == 0) {
    params$write_path   <- file.path(monitor_folder, "inputfiles")
  } else {
    params$write_path   <- file.path(monitor_folder, "msoc")
  }
  params$readPathAC    <- file.path(monitor_folder, "inputfiles")
  params$readPathDP    <- file.path(monitor_folder, paste0("msoc", 1:params$numDataPartners))

  if (!create_io_location(monitor_folder, "dplocal")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$dplocalPath, "."),
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
                                paste0(params$readPathAC, "."),
                                "Check the path and restart the program.")
  }
  if (!create_io_location(monitor_folder, "msoc")) {
    params$failed <- TRUE
    params$error_message <- paste(params$error_message,
                                "Could not create directory",
                                paste0(params$write_path, "."),
                                "Check the path and restart the program.")
  }
  for (id in 1:params$numDataPartners) {
    if (!create_io_location(monitor_folder, paste0("msoc", id))) {
      params$failed <- TRUE
      params$error_message <- paste(params$error_message,
                                  "Could not create directory",
                                  paste0(params$readPathDP[id], "."),
                                  "Check the path and restart the program.")
    }
  }

  if (params$data_partner_id != 0) {
    Sys.sleep(1)
    delete_trigger("files_done.ok", params$readPathAC)
    for (id in 1:params$numDataPartners) {
      delete_trigger("files_done.ok", params$readpathDP[id])
    }
  }

  empty = NULL
  write_time <- proc.time()[3]
  save(empty, file = file.path(params$write_path, "empty.rdata"))
  write_size = file.size(file.path(params$write_path, "empty.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "PrepareFolder.ACDP", 0, 0, write_time, write_size)
  return(params)
}


#' @importFrom stats model.matrix
PrepareDataLinLog.DP1 <- function(params, data, y_name = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinLog.DP1\n\n")

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
  covariate_index <- setdiff(1:ncol(data), response_index)
  workdata$tags <- create_model_matrix_tags(data[, covariate_index, drop = FALSE])
  workdata$tags = c("(Intercept)", workdata$tags)
  names(workdata$tags)[1] = "numeric"
  x = model.matrix(~ ., data[, c(response_index, covariate_index), drop = FALSE])
  rownames(x) = NULL
  covariate_index <- setdiff(1:ncol(x), 2)
  workdata$x = x[, c(2, covariate_index), drop = FALSE]

  workdata$n        = nrow(workdata$x)
  workdata$colmin   = apply(workdata$x, 2, min)
  workdata$colmax   = apply(workdata$x, 2, max)
  workdata$colsum   = apply(workdata$x, 2, sum)
  workdata$colrange <- workdata$colmax - workdata$colmin
  for (i in 1:ncol(workdata$x)) {
    if (workdata$colmin[i] == workdata$colmax[i]) {
      workdata$colmin[i] = 0
      workdata$colrange[i] = 1
    }
    workdata$x[, i] = (workdata$x[, i] - workdata$colmin[i]) / workdata$colrange[i]
  }

  return(workdata)
}

#' @importFrom stats model.matrix
PrepareDataLinLog.DPk <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinLog.DPk\n\n")

  workdata <- list()
  workdata$failed <- FALSE

  workdata$failed <- check_data_format(params, data)

  if (workdata$failed) {
    return(workdata)
  }

  data <- data.frame(data) # convert to a clean data.frame

  workdata$tags <- create_model_matrix_tags(data)
  workdata$x = model.matrix(~ ., data)
  rownames(workdata$x) = NULL
  workdata$x <- workdata$x[, -1, drop = FALSE]

  workdata$n        = nrow(workdata$x)
  workdata$colmin   = apply(workdata$x, 2, min)
  workdata$colmax   = apply(workdata$x, 2, max)
  workdata$colsum   = apply(workdata$x, 2, sum)
  workdata$colrange <- workdata$colmax - workdata$colmin
  for (i in 1:ncol(workdata$x)) {
    if (workdata$colmin[i] == workdata$colmax[i]) {
      workdata$colmin[i] = 0
      workdata$colrange[i] = 1
    }
    workdata$x[, i] = (workdata$x[, i] - workdata$colmin[i]) / workdata$colrange[i]
  }

  return(workdata)
}

SendBasicInfo.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendbasicInfo.DP\n\n")
  n = data$n
  params$n <- n
  analysis = params$analysis
  data_partner_id = params$data_partner_id
  write_time <- proc.time()[3]
  save(analysis, n, data_partner_id, file = file.path(params$write_path, "n_analysis.rdata"))
  write_size = file.size(file.path(params$write_path, "n_analysis.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "SendBasicInfo.DP", 0, 0, write_time, write_size)
  return(params)
}

CheckAgreement.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckAgreement.AC\n\n")
  read_time <- 0
  read_size = 0
  analysisAll = rep("", params$numDataPartners)
  nAll        = rep(0, params$numDataPartners)
  ndata_partner_id = rep(0, params$numDataPartners)
  message1    = NULL
  message2    = NULL
  n           = NULL
  analysis    = NULL
  data_partner_id = NULL
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "n_analysis.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "n_analysis.rdata"))
    read_time <- read_time + proc.time()[3]
    analysisAll[id] = analysis
    nAll[id]        = n
    ndata_partner_id[id] = data_partner_id
  }

  if (any(params$analysis != analysisAll)) {
    params$failed <- TRUE
    message1 = "Different regressions have been specified.\n"
    message1 = paste(message1, "Analysis center specified", params$analysis, "regression.\n")
    for (id in 1:params$numDataPartners) {
      message1 = paste(message1, "Data partner", id, "specified", analysisAll[id], "regression.\n")
    }
  }

  if (min(nAll) < max(nAll)) {
    params$failed <- TRUE
    message2 = "Data partners provided different numbers of observations.\n"
    for (id in 1:params$numDataPartners) {
      message2 = paste(message2, "Data partner", id, "has", nAll[id], "observations.\n")
    }
  }

  message3error = FALSE
  message3 = ""
  for (i in 1:params$numDataPartners) {
    if (i != ndata_partner_id[i]) {
      message3error = TRUE
      params$failed <- TRUE
      message3 = paste0(message3, "Data Partner ", i, " reports its ID as ", ndata_partner_id[i], "\n")
    }
  }
  if (message3error) {
    message3 = paste0("Check PopMedNet DataMart setup.\n", message3)
  }

  if (params$failed) {
    params$error_message <- paste0(message1, message2, message3)
  }

  params <- add_to_log(params, "CheckAgreement.AC", read_time, read_size, 0, 0)
  return(params)
}


#' @importFrom stats runif
prepare_params_linear_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_linear_DP\n\n")
  params$n          = nrow(data$x)
  params$p          = ncol(data$x)
  temp = as.numeric(Sys.time())
  set.seed((temp - trunc(temp)) * .Machine$integer.max)
  params$seed       = floor(runif(1) * .Machine$integer.max)
  params$scaler     = 1 + runif(1)

  p = params$p
  seed = params$seed
  scaler = params$scaler

  write_time <- proc.time()[3]
  save(p, scaler, seed, file = file.path(params$write_path, "p_scaler_seed.rdata"))
  write_size = file.size(file.path(params$write_path, "p_scaler_seed.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "prepare_params_linear_DP", 0, 0, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm
PrepareSharesLinear.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareSharesLinear.DP\n\n")
  read_time <- 0
  read_size = 0
  p = seed <- scaler = NULL

  set.seed(params$seed, kind = "Mersenne-Twister")
  halfshare = matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)

  products = rep(list(list()), params$numDataPartners)

  params$ps = c()
  params$scalers = c()
  params$seeds = c()

  for (id in 1:params$numDataPartners) {
    if (id == params$data_partner_id) {
      products[[id]] = t(data$x) %*% data$x
      params$ps      = c(params$ps, params$p)
      params$scalers = c(params$scalers, params$scaler)
      params$seeds   = c(params$seeds, params$seed)
      next
    }
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
    read_size = read_size + file.size(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
    read_time <- read_time + proc.time()[3]
    params$ps      = c(params$ps, p)
    params$scalers = c(params$scalers, scaler)
    params$seeds   = c(params$seeds, seed)

    set.seed(seed, kind = "Mersenne-Twister")
    halfShare2 = matrix(rnorm(params$n * p, sd = 20), nrow = params$n, ncol = p)

    if (id < params$data_partner_id) {
      products[[id]] = t(halfShare2) %*% (data$x - scaler / (scaler + params$scaler) * halfshare)
    }

    if (id > params$data_partner_id) {
      products[[id]] = t(data$x - scaler / (scaler + params$scaler) * halfshare) %*% halfShare2
    }
  }

  halfshare = data$x - halfshare
  colmin    = data$colmin
  colrange  = data$colrange
  colsum    = data$colsum
  colnames  = colnames(data$x)
  tags      = data$tags

  write_time <- proc.time()[3]
  save(products, file = file.path(params$write_path, "products.rdata"))
  save(halfshare, file = file.path(params$write_path, "halfshare.rdata"))
  save(colmin, colrange, colsum, colnames, tags, file = file.path(params$write_path, "colstats.rdata"))
  write_size = sum(file.size(file.path(params$write_path, c("products.rdata",
                                                          "halfshare.rdata",
                                                          "colstats.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "PrepareSharesLinear.DP", read_time, read_size, write_time, write_size)
  return(params)
}


get_products_linear_AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_products_linear_AC\n\n")
  read_time <- 0
  read_size = 0
  p = 0
  n = 0

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
    if (id == 1) n = nrow(halfshare)
  }

  m = matrix(0, p, p)
  colnames(m) = allcolnames
  rownames(m) = allcolnames
  offset1 = 1
  params$pi = rep(0, params$numDataPartners)
  for (id1 in 1:params$numDataPartners) {
    p1 = ncol(allhalfshare[[id1]])
    params$pi[id1] = p1
    offset2 = offset1
    for (id2 in id1:params$numDataPartners) {
      p2 = ncol(allhalfshare[[id2]])
      if (id1 == id2) {
        m[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = allproducts[[id1]][[id2]]
      } else {
        temp = allproducts[[id1]][[id2]] + allproducts[[id2]][[id1]] +
          t(allhalfshare[[id1]]) %*% allhalfshare[[id2]]
        m[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = temp
        m[offset2:(offset2 + p2 - 1), offset1:(offset1 + p1 - 1)] = t(temp)
      }
      offset2 = offset2 + p2
    }
    offset1 = offset1 + p1
  }

  m = diag(allcolrange) %*% m %*% diag(allcolrange) +
    outer(allcolmin, allcolsum) + outer(allcolsum, allcolmin) -
    n * outer(allcolmin, allcolmin)

  params$xtx          = m[2:p, 2:p, drop = FALSE]
  params$xty          = m[2:p, 1, drop = FALSE]
  params$yty          = m[1, 1]
  params$means_y       = allcolsum[1] / n
  params$means        = allcolsum[-1] / n
  params$n            = n
  params$p            = p
  params$colnames     = allcolnames[-1]
  params$party        = party[-1]
  params$converged    = TRUE
  params$tags         = alltags

  params <- add_to_log(params, "get_products_linear_AC", read_time, read_size, 0, 0)
  return(params)
}


#' @importFrom  stats pf pt
compute_results_linear_AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "compute_results_linear_AC\n\n")
  stats           = params$stats
  stats$converged <- params$converged
  n        = params$n
  yty      = params$yty
  xty      = params$xty
  xtx      = params$xtx
  means_y   = params$means_y

  # First we de-standardize.

  nrow = nrow(xtx)
  indicies = c(1)
  for (i in 2:nrow) {
    temp_indicies = c(indicies, i)
    if (rcond(xtx[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  tags = params$tags
  min = 1
  for (id in 1:params$numDataPartners) {
    max = min + params$pi[id] - 1
    if (id == 1) {
      max = max - 1
    }
    idx <- indicies[which(min <= indicies & indicies <= max)] - min + 1
    temp = tags[[id]]
    temp = temp[idx]
    tags[[id]] = temp
    min = max + 1
  }

  params$error_message <- ""
  numeric_found = FALSE
  for (id in 2:params$numDataPartners) {
    if (length(unique(tags[[id]])) == 0) {
      params$failed <- TRUE
      params$error_message <- paste0(params$error_message,
                                   paste("After removing colinear covariates, Data Partner",
                                         id, "has no covariates."))
    } else {
      numeric_found = numeric_found | "numeric" %in% names(tags[[id]])
    }
  }
  if (!numeric_found) {
    params$failed <- TRUE
    params$error_message <- paste0(params$error_message,
                                 paste("After removing colinear covariates, no Data Partner > DP1 has a numeric covariate."))
  }

  stats$failed    = params$failed

  p             = length(indicies)
  p1            = ncol(xtx)
  xtx_old       = xtx
  xty_old       = xty
  xtx           = xtx[indicies, indicies, drop = FALSE]
  xty           = xty[indicies, , drop = FALSE]

  invxtx = solve(xtx)
  betas  = drop(invxtx %*% xty)

  num_covariates <- p - 1

  #   # If true sse is approximately 0, random variations could cause this
  #   # calculation to be less than 0
  #   # If calculated sse is less than 0, we set it equal to 0.
  sse     = max(drop(yty - 2 * t(xty) %*% betas + (t(betas) %*% xtx) %*% betas), 0)
  rstderr = drop(sqrt(sse / (n - num_covariates - 1)))
  sst     = drop(yty - means_y^2 * n)
  ssr     = sst - sse
  df1     = num_covariates
  df2     = n - num_covariates - 1
  if (sse == 0) {
    f_stat <- Inf
  } else {
    f_stat <- (ssr / df1) / (sse / df2)
  }
  f_pval <- pf(f_stat, df1, df2, lower.tail = FALSE)
  if (sse == 0) {
    r_sq <- 1
  } else {
    r_sq <- drop(1 - sse / sst)
  }
  adj_r_sq <- drop(1 - (n - 1) / (n - num_covariates - 1) * (1 - r_sq))
  if (rstderr == 0) {
    tvals = rep(Inf, num_covariates + 1)
  } else {
    tvals   = betas / (rstderr * sqrt(diag(invxtx)))
  }
  secoef  = tvals^-1 * betas
  pvals   = 2 * pt(abs(tvals), n - num_covariates - 1, lower.tail = FALSE)
  stats$party                  = params$party
  stats$responseParty          = "dp1"
  stats$coefficients           = rep(NA, p1)
  stats$tvals                  = rep(NA, p1)
  stats$secoef                 = rep(NA, p1)
  stats$pvals                  = rep(NA, p1)

  stats$sse                    = sse
  stats$coefficients[indicies] = betas
  stats$tvals[indicies]        = tvals
  stats$secoef[indicies]       = secoef
  stats$pvals[indicies]        = pvals
  stats$rstderr                = rstderr
  stats$rsquare                = r_sq
  stats$adjrsquare             = adj_r_sq
  stats$f_stat                  = f_stat
  stats$f_pval <- f_pval
  stats$df1                    = df1
  stats$df2                    = df2
  stats$n                      = params$n
  stats$xtx                    = xtx_old
  stats$xty                    = xty_old
  stats$yty                    = yty
  stats$means_y                 = means_y
  stats$means                  = params$means

  names(stats$party)           = params$colnames
  names(stats$coefficients)    = params$colnames
  names(stats$secoef)          = params$colnames
  names(stats$tvals)           = params$colnames
  names(stats$pvals)           = params$colnames
  colnames(stats$xtx)          = params$colnames
  rownames(stats$xtx)          = params$colnames
  colnames(stats$xty)          = colnames(params$xty)
  rownames(stats$xty)          = params$colnames

  class(stats) = "vdralinear"

  params$stats = stats
  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size = file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "compute_results_linear_AC", 0, 0, write_time, write_size)
  return(params)
}


get_results_linear_DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_linear_DP\n\n")
  params$converged = TRUE
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "stats.rdata"))
  read_size = file.size(file.path(params$readPathAC, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats

  params <- add_to_log(params, "get_results_linear_DP", read_time, read_size, 0, 0)
  return(params)
}


############################## PARENT FUNCTIONS ###############################


DataPartnerKLinear <- function(data,
                              y_name           = NULL,
                              numDataPartners = NULL,
                              data_partner_id   = NULL,
                              monitor_folder   = NULL,
                              sleep_time       = 10,
                              max_waiting_time  = 24 * 60 * 60,
                              popmednet      = TRUE,
                              trace          = FALSE,
                              verbose        = TRUE) {

  params <- prepare_params_kp("linear", data_partner_id, numDataPartners, ac = FALSE,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- initialize_log_kp(params)
  params <- initialize_time_stamps_kp(params)
  params <- initialize_tracking_table_kp(params)
  header(params)

  params   = PrepareFolder.ACDP(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  if (data_partner_id == 1) {
    data = PrepareDataLinLog.DP1(params, data, y_name)
    params <- add_to_log(params, "PrepareDataLinLog.DP1", 0, 0, 0, 0)
  } else {
    data = PrepareDataLinLog.DPk(params, data)
    params <- add_to_log(params, "PrepareDataLinLog.DP2", 0, 0, 0, 0)
  }


  if (data$failed) {
    params$error_message <- paste("Error processing data for data partner", params$data_partner_id, "\n")
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)
    params$error_message <- read_error_message(params$readPathAC)
    warning(params$error_message)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- SendBasicInfo.DP(params, data)
  files <- "n_analysis.rdata"
  params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  possibleError = ReceivedError.kp(params, from = "AC")
  if (possibleError$error) {
    params$error_message <- possibleError$message
    warning(possibleError$message)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_linear_DP(params, data)
  files <- "p_scaler_seed.rdata"
  params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  params <- PrepareSharesLinear.DP(params, data)
  files <- c("products.rdata", "halfshare.rdata", "colstats.rdata")
  params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  possibleError = ReceivedError.kp(params, from = "AC")
  if (possibleError$error) {
    params$error_message <- possibleError$message
    warning(possibleError$message)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  } else {
    params <- get_results_linear_DP(params)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }
}


AnalysisCenterKLinear <- function(numDataPartners = NULL,
                                 monitor_folder   = NULL,
                                 msreqid         = "v_default_0_000",
                                 sleep_time       = 10,
                                 max_waiting_time  = 24 * 60 * 60,
                                 popmednet       = TRUE,
                                 trace           = FALSE,
                                 verbose         = TRUE) {
  params <- prepare_params_kp("linear", 0, numDataPartners, msreqid, ac = TRUE,
                            popmednet = popmednet, trace = trace, verbose = verbose)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  params <- initialize_log_kp(params)
  params <- initialize_time_stamps_kp(params)
  params <- initialize_tracking_table_kp(params)
  header(params)

  params   = PrepareFolder.ACDP(params, monitor_folder)

  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.kp(params, from = "DP", max_waiting_time = max_waiting_time)

  possibleError = ReceivedError.kp(params, from = "DP")
  if (possibleError$error) {
    params$error_message <- possibleError$message
    warning(possibleError$message)
    make_error_message(params$write_path, possibleError$message)
    files <- "error_message.rdata"
    params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- CheckAgreement.AC(params)

  if (params$failed) {
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    warning(params$error_message)
    params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  files <- "empty.rdata"
  params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- get_products_linear_AC(params)
  params <- compute_results_linear_AC(params)

  if (params$failed) {
    make_error_message(params$write_path, params$error_message)
    files <- "error_message.rdata"
    warning(params$error_message)
    params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  } else {
    files <- "stats.rdata"
    params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time)
    SummarizeLog.kp(params)
    return(params$stats)
  }
}
