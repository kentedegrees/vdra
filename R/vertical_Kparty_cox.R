################### DISTRIBUTED COX REGRESSION FUNCTIONS ##################

#' @importFrom stats model.matrix
prepare_data_cox_DP <- function(params, data, y_name, strata, mask) {
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
  if (params$data_partner_id == 1) {
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
  covariate_index <- setdiff(1:ncol(data), union(strata_index, response_index))
  workdata$n = nrow(data)
  if (length(covariate_index) == 0) {
    if (params$data_partner_id == 1) {
      workdata$x  <- matrix(0, nrow = nrow(data), ncol = 0)
    } else {
      warning("After removing strata, data is empty.  Party B must supply at least one non-strata covariate.")
      workdata$failed <- TRUE
      return(workdata)
    }
  } else {
    workdata$tags <- create_model_matrix_tags(data[, covariate_index, drop = FALSE])
    workdata$x <- model.matrix(~ ., data[, covariate_index, drop = FALSE])
    workdata$x <- workdata$x[, -1, drop = FALSE]
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
  }

  return(workdata)
}


SendStrataNamesCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataNamesCox.DP\n\n")
  strataNames <- list()
  strataNames$strataFromMe <- data$strata$strataFromMe
  strataNames$strataFromOthers <- data$strata$strataFromOthers
  write_time <- proc.time()[3]
  save(strataNames, file = file.path(params$write_path, "strata_names.rdata"))
  write_size <- file.size(file.path(params$write_path, "strata_names.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "SendStrataNamesCox.DP", 0, 0, write_time, write_size)
  return(params)
}


check_strata_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "check_strata_cox_DP\n\n")
  read_time <- 0
  read_size <- 0
  strataNames          = NULL
  strataClaimed        = rep(list(list()), params$numDataPartners)
  strataUnclaimed      = rep(list(list()), params$numDataPartners)
  if (length(data$strata$strataFromMe) > 0) {
    strataClaimed[[1]]   <- data$strata$strataFromMe
  }
  if (length(data$strata$strataFromOthers) > 0) {
    strataUnclaimed[[1]] <- data$strata$strataFromOthers
  }
  for (i in 2:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[i], "strata_names.rdata"))
    read_size <- read_size + file.size(file.path(params$readPathDP[i], "strata_names.rdata"))
    read_time <- read_time + proc.time()[3]
    if (length(strataNames$strataFromMe) > 0) {
      strataClaimed[[i]] = strataNames$strataFromMe
    }
    if (length(strataNames$strataFromOthers) > 0) {
      strataUnclaimed[[i]] = strataNames$strataFromOthers
    }
  }

  params$strataClaimed = strataClaimed

  claimed = c()
  unclaimed = c()
  specified = rep(list(list()), params$numDataPartners)
  for (i in 1:params$numDataPartners) {
    temp = c(as.character(strataClaimed[[i]]), as.character(strataUnclaimed[[i]]))
    if (length(temp) > 0) {
      specified[[i]] = sort(temp)
    }
    if (length(strataClaimed[[i]]) > 0) {
      claimed = c(claimed, strataClaimed[[i]])
    }
    if (length(strataUnclaimed[[i]]) > 0) {
      unclaimed = union(unclaimed, strataUnclaimed[[i]])
    }
  }

  if (length(claimed) > 0) {
    claimed   = sort(claimed)
  }
  if (length(unclaimed) > 0) {
    unclaimed = sort(unclaimed)
  }

  passed = TRUE
  for (i in 2:params$numDataPartners) {
    passed = passed && length(specified[[1]]) == length(specified[[i]]) && all(specified[[1]] == specified[[i]])
  }

  if (!passed) {
    params$failed <- TRUE
    params$error_message <- "Data partners specified different strata:\n"
    for (i in 1:params$numDataPartners) {
      temp = NULL
      if (length(specified[[i]] > 0)) {
        temp = paste0(specified[[i]], collapse = ", ")
      }
      params$error_message <- paste0(params$error_message,
                                   paste("Data Partner", i, "specified strata:", temp, "\n"))
    }
    params <- add_to_log(params, "check_strata_cox_DP", read_time, read_size, 0, 0)
    return(params)
  }
  passed = TRUE
  if (length(claimed) > 0) {
    tab = table(claimed)
    passed = all(tab == 1)
  }

  if (!passed) {
    params$failed <- TRUE
    params$error_message <- "The following strata are claimed by two or more data partners: "
    params$error_message <- paste0(params$error_message, paste0(names(tab[which(tab > 1)]), collapse = ", "), "\n")
    params$error_message <- paste0(params$error_message, "Make Sure that strata covariate names are unique to each data partner.")
    params <- add_to_log(params, "check_strata_cox_DP", read_time, read_size, 0, 0)
    return(params)
  }

  passed = TRUE
  claimed1 = sort(unique(claimed))
  passed = length(claimed1) == length(unclaimed) && all(claimed1 == unclaimed)

  if (!passed) {
    params$failed <- TRUE
    params$error_message <- "No data partner has the following specified strata: "
    params$error_message <- paste0(params$error_message, paste0(unclaimed[which(!(unclaimed %in% claimed1))], collapse = ", "), ".")
    params <- add_to_log(params, "check_strata_cox_DP", read_time, read_size, 0, 0)
    return(params)
  }

  params <- add_to_log(params, "check_strata_cox_DP", read_time, read_size, 0, 0)
  return(params)
}


send_strata_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "send_strata_cox_DP\n\n")
  strata <- list()
  strata$x <- data$strata$x
  strata$legend <- data$strata$legend
  write_time <- proc.time()[3]
  save(strata, file = file.path(params$write_path, "strata.rdata"))
  write_size <- file.size(file.path(params$write_path, "strata.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "send_strata_cox_DP", 0, 0, write_time, write_size)
  return(params)
}


prepare_strata_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_strata_cox_DP\n\n")
  strata = NULL
  survival = data$survival
  strataTemp   <- list()
  totalStrata = 0
  read_size <- 0
  read_time <- 0
  for (id in 1:params$numDataPartners) {
    totalStrata = totalStrata + length(params$strataClaimed[[id]])
  }

  if (totalStrata == 0) {
    strataTemp$x <- data.frame(const__ = rep(1, params$n))
    strataTemp$legend = FALSE
  } else {
    first = TRUE
    strataTemp$legend <- list()
    for (id in 1:params$numDataPartners) {
      if (length(params$strataClaimed[[id]]) > 0) {
        if (id == 1) {
          strataTemp$x <- data$strata$x
          strataTemp$legend <- data$strata$legend
          first = FALSE
        } else {
          read_time <- read_time - proc.time()[3]
          load(file.path(params$readPathDP[[id]], "strata.rdata"))
          read_size <- read_size + file.size(file.path(params$readPathDP[[id]], "strata.rdata"))
          read_time <- read_time + proc.time()[3]
          if (first) {
            strataTemp$x = strata$x
            strataTemp$legend = strata$legend
            first = FALSE
          } else {
            strataTemp$x <- cbind(strataTemp$x, strata$x)
            strataTemp$legend = c(strataTemp$legend, strata$legend)
          }
        }
      }
    }
  }

  sorted <- do.call("order", cbind(strataTemp$x, survival$rank, survival$status))
  unsort = order(sorted)
  strataTemp$x = strataTemp$x[sorted, , drop = FALSE]
  survival$rank   = survival$rank[sorted]
  survival$status = survival$status[sorted]
  survival$sortedIdx = sorted
  survival$unsortIdx = unsort
  ranks <- which(apply(abs(apply(strataTemp$x, 2, diff)), 1, sum) > 0)
  ranks = c(ranks, nrow(strataTemp$x))
  names(ranks) = NULL
  strata = rep(list(list()), length(ranks))
  if (length(ranks) == 1 && colnames(strataTemp$x)[1] == "const__") {
    strata[[1]]$start = 1
    strata[[1]]$end   = as.integer(length(survival$rank))
    strata[[1]]$label <- ""
  } else {
    start = 1
    for (i in  1:length(ranks)) {
      strata[[i]]$start = start
      strata[[i]]$end   = as.integer(ranks[i])
      label <- ""
      for (j in 1:ncol(strataTemp$x)) {
        temp = colnames(strataTemp$x)[j]
        label <- paste0(label, temp, "=", strataTemp$legend[[temp]][strataTemp$x[start, j]])
        if (j < ncol(strataTemp$x)) {
          label <- paste0(label, ", ")
        }
      }
      strata[[i]]$label <- label
      start = as.numeric(ranks[i]) + 1
    }
  }
  for (i in 1:length(strata)) {
    idx <- strata[[i]]$start:strata[[i]]$end
    temp  = table(survival$rank[idx])
    m = length(temp)   # number of unique observed times, including where no one fails
    # Count the number of 0's and 1's for each observed time
    temp0 <- table(survival$rank[idx], survival$status[idx])
    # Check if there are all 1's or all 0's .  If so, add them into the table.
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 <- cbind(temp0, 0)
        colnames(temp0) = c("0", "1")
      } else {
        temp0 <- cbind(0, temp0)
        colnames(temp0) = c("0", "1")
      }
    }
    strata[[i]]$J = as.integer(sum(temp0[, 2] > 0))    # number of distinct failure times
    # The number of failures at each rank which has a failure
    strata[[i]]$nfails = as.numeric(temp0[which(temp0[, 2] > 0), 2])
    # The first index of the ranks for which the number of failures is > 0.
    strata[[i]]$start0 = c(1, (cumsum(temp)[1:(m - 1)] + 1))[which(temp0[, 2] > 0)]
    # The first index of a failure for each rank which has a failure
    strata[[i]]$start1 = strata[[i]]$start0 + temp0[which(temp0[, 2] > 0), 1]
    # The last index of a failure for each rank which has a failure
    strata[[i]]$stop1  = as.numeric(strata[[i]]$start1 + strata[[i]]$nfails  - 1)
    strata[[i]]$start0 = as.numeric(strata[[i]]$start0 + strata[[i]]$start - 1)
    strata[[i]]$start1 = as.numeric(strata[[i]]$start1 + strata[[i]]$start - 1)
    strata[[i]]$stop1  = as.numeric(strata[[i]]$stop1  + strata[[i]]$start - 1)
  }

  survival$strata = strata
  survival$num_events <- sum(survival$status)
  params$survival = survival

  if (totalStrata == 0) {
    params$strata = as.matrix(strataTemp$x)
  } else {
    ptemp = 1
    for (i in 1:ncol(strataTemp$x)) {
      ptemp = ptemp + max(strataTemp$x[, i] - min(strataTemp$x[, i]))
    }
    params$strata  <- matrix(0, nrow = params$n, ncol = ptemp)
    colnames(params$strata) = paste0("strata:", 1:ptemp)
    params$strata[, 1] = 1
    idx <- 2
    for (i in 1:ncol(strataTemp$x)) {
      min1 = min(strataTemp$x[, i])
      max1 = max(strataTemp$x[, i])
      if (min1 < max1) {
        for (j in (min1 + 1):max1) {
          params$strata[which(strataTemp$x[, i] == j), idx] = 1
          idx <- idx + 1
        }
      }
    }
  }

  pStrata = ncol(params$strata)
  params$pStrata = pStrata

  write_time <- proc.time()[3]
  save(pStrata, survival, file = file.path(params$write_path, "survival.rdata"))
  write_size <- file.size(file.path(params$write_path, "survival.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "prepare_strata_cox_DP", read_time, read_size, write_time, write_size)

  return(params)
}


AddStrataToDataCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "AddStrataToDataCox.DP\n\n")
  data$x <- cbind(params$strata, data$x)

  colmin1 = apply(params$strata, 2, min)
  colmax1 = apply(params$strata, 2, max)
  colrange1 = colmax1 - colmin1
  colsum1   = apply(params$strata, 2, sum)
  data$colmax = c(colmax1, data$colmax)
  data$colmin = c(colmin1, data$colmin)
  data$colmin[1] = 0
  data$colrange = c(colrange1, data$colrange)
  data$colsum = c(colsum1, data$colsum)

  return(data)
}

#' @importFrom stats runif
prepare_params_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_DP\n\n")
  params$strata     = NULL
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
  write_size <- file.size(file.path(params$write_path, "p_scaler_seed.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "prepare_params_cox_DP", 0, 0, write_time, write_size)
  return(params)
}

#' @importFrom stats rnorm
PrepareSharesCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareSharesCox.DP\n\n")
  read_time <- 0
  read_size <- 0
  n <- params$n
  p = params$p
  scaler = NULL
  seed   = NULL
  set.seed(params$seed, kind = "Mersenne-Twister")
  halfshare.L   <- matrix(rnorm(n * p, sd = 20), nrow = n, ncol = p)

  halfshare.R  = data$x - halfshare.L

  products = rep(list(list()), params$numDataPartners)

  params$ps = c()
  params$scalers = c()
  params$seeds = c()

  for (id in 1:params$numDataPartners) {
    if (id == params$data_partner_id) {
      products[[id]] <- t(data$x) %*% data$x
      params$ps      = c(params$ps, params$p)
      params$scalers = c(params$scalers, params$scaler)
      params$seeds   = c(params$seeds, params$seed)
      next
    }
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
    read_size <- read_size + file.size(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
    read_time <- read_time + proc.time()[3]
    params$ps      = c(params$ps, p)
    params$scalers = c(params$scalers, scaler)
    params$seeds   = c(params$seeds, seed)

    set.seed(seed, kind = "Mersenne-Twister")
    halfShare2  <- matrix(rnorm(params$n * p, sd = 20), nrow = params$n, ncol = p)

    if (id < params$data_partner_id) {
      products[[id]] = t(halfShare2) %*% (data$x - scaler / (scaler + params$scaler) * halfshare.L)
    }

    if (id > params$data_partner_id) {
      products[[id]] = t(data$x - scaler / (scaler + params$scaler) * halfshare.L) %*% halfShare2
    }
  }

  colmin    = data$colmin
  colrange  = data$colrange
  colsum    = data$colsum
  colnames  = colnames(data$x)
  tags      = data$tags

  write_time <- proc.time()[3]
  save(products, file = file.path(params$write_path, "products.rdata"))
  save(halfshare.R, file = file.path(params$write_path, "halfshare.rdata"))
  save(colmin, colrange, colsum, colnames, tags, file = file.path(params$write_path, "colstats.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("products.rdata",
                                                          "halfshare.rdata",
                                                          "colstats.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "PrepareSharesCox.DP", read_time, read_size, write_time, write_size)

  return(params)
}


GetStrata.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetStrata.AC\n\n")
  survival = NULL
  pStrata  = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "survival.rdata"))
  read_size <- file.size(file.path(params$readPathDP[1], "survival.rdata"))
  read_time <- proc.time()[3] - read_time
  params$survival = survival
  params$pStrata = pStrata
  params <- add_to_log(params, "GetStrata.AC", read_time, read_size, 0, 0)
  return(params)
}


get_products_cox_AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_products_cox_AC\n\n")
  read_time <- 0
  read_size <- 0
  p = 0
  n = 0

  allproducts  = rep(list(list()), params$numDataPartners)
  allhalfshare = rep(list(list()), params$numDataPartners)
  alltags      = rep(list(list()), params$numDataPartners)
  products    = NULL
  halfshare.R = NULL
  tags        = NULL
  allcolmin = allcolrange = allcolsum = allcolnames = NULL
  colmin = colrange = colsum = colnames = NULL
  party = NULL
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "products.rdata"))
    load(file.path(params$readPathDP[id], "halfshare.rdata"))
    load(file.path(params$readPathDP[id], "colstats.rdata"))
    read_size <- read_size + sum(file.size(file.path(params$readPathDP[id],
                                                  c("products.rdata",
                                                    "halfshare.rdata",
                                                    "colstats.rdata"))))
    read_time <- read_time + proc.time()[3]

    allproducts[[id]]  = products
    allhalfshare[[id]] = halfshare.R
    alltags[[id]]      = tags
    allcolmin          = c(allcolmin, colmin)
    allcolrange        = c(allcolrange, colrange)
    allcolsum          = c(allcolsum, colsum)
    allcolnames        = c(allcolnames, colnames)
    party              = c(party, rep(paste0("dp", id), length(colnames)))
    p = p + ncol(halfshare.R)
    if (id == 1) n = nrow(halfshare.R)
  }

  m  <- matrix(0, p, p)
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

  params$sts          = m
  params$n            = n
  params$p            = p
  params$colmin       = allcolmin
  params$colrange     = allcolrange
  params$colsum       = allcolsum
  params$colnames     = allcolnames
  params$party        = party
  params$tags         = alltags

  params$halfshare    = allhalfshare

  params <- add_to_log(params, "get_products_cox_AC", read_time, read_size, 0, 0)
  return(params)
}

check_colinearity_cox_AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "check_colinearity_cox_AC\n\n")
  sts = params$sts

  nrow = nrow(sts)
  indicies = c(1)
  for (i in 2:nrow) {
    temp_indicies = c(indicies, i)
    if (rcond(sts[temp_indicies, temp_indicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  if (max(indicies) > params$pStrata) {
    indicies = indicies[which(indicies > params$pStrata)]
  } else {
    indicies = c()
  }


  sts = sts[indicies, indicies, drop = FALSE]

  params$sts = sts

  # Extract the indicies to keep for each party and check for errors.

  params$colmin       = params$colmin[indicies]
  params$colrange     = params$colrange[indicies]
  params$colsum       = params$colsum[indicies]
  params$fullindicies = indicies       # To be used when computing stats
  params$p            = params$p - 1   # Get rid of the response from the count

  params$indicies = rep(list(c()), params$numDataPartners)
  params$idx      = rep(list(c()), params$numDataPartners)
  tags            = rep(list(c()), params$numDataPartners)

  start = 1
  min = 1
  for (id in 1:params$numDataPartners) {
    max = min + params$pi[id] - 1
    idx <- which(min <= indicies & indicies <= max)
    if (length(idx) > 0) {
      idx_1 = indicies[idx] - min + 1
      end   = start + length(idx) - 1

      params$indicies[[id]] = idx_1
      params$idx[[id]] = start:end

      start = end + 1
      if (id >= 2) {
        temp = params$tags[[id]]
        temp = temp[idx_1]
        tags[[id]] = temp
      }

    }
    min = max + 1
  }

  params$error_message <- ""
  numeric_found = FALSE
  for (id in 2:params$numDataPartners) {
    if (length(unique(tags[[id]])) == 0) {
      params$failed <- TRUE
      params$error_message <- paste0(params$error_message,
                                   paste("After removing colinear covariates, Data Partner", id, "has no covariates."))
    }
  }

  if (params$failed) {
    params <- add_to_log(params, "check_colinearity_logistic_AC", 0, 0, 0, 0)
  }
  indicies = params$indicies
  idx      = params$idx

  params$pReduct = c()
  for (id in 1:params$numDataPartners) {
    params$pReduct = c(params$pReduct, length(indicies[[id]]))
  }

  for (id in 1:params$numDataPartners) {
    params$halfshare[[id]] = params$halfshare[[id]][, indicies[[id]], drop = FALSE]
  }
  params$halfshare <- do.call(cbind, params$halfshare)

  params$halfsharecolsum = apply(params$halfshare, 2, sum)

  cutoff = params$cutoff
  pReduct = params$pReduct

  write_time <- proc.time()[3]
  save(cutoff, idx, indicies, pReduct, file = file.path(params$write_path, "indicies.rdata"))
  write_size <- file.size(file.path(params$write_path, "indicies.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "check_colinearity_logistic_AC", 0, 0, write_time, write_size)
  return(params)
}

ComputeUCox.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeUCox.AC\n\n")
  read_time <- 0
  read_size <- 0
  if (params$alg_iteration_counter == 1) {
    u = 1
  } else {
    uTemp = 0
    for (id in 1:params$numDataPartners) {
      read_time <- read_time - proc.time()[3]
      load(file.path(params$readPathDP[id], "u.rdata"))
      read_size <- read_size + file.size(file.path(params$readPathDP[id], "u.rdata"))
      read_time <- read_time + proc.time()[3]
      uTemp = uTemp + u
    }
    u = uTemp
  }
  params$u = u
  write_time <- proc.time()[3]
  save(u, file = file.path(params$write_path, "u.rdata"))
  write_size <- file.size(file.path(params$write_path, "u.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeUCox.AC", read_time, read_size, write_time, write_size)
  return(params)
}

update_params_cox_DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "update_params_cox_DP\n\n")
  indicies = NULL
  idx      = NULL
  u        = NULL
  cutoff   = NULL
  pReduct  = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "indicies.rdata"))
  load(file.path(params$readPathAC, "u.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "indicies.rdata")) +
    file.size(file.path(params$readPathAC, "u.rdata"))
  read_time <- proc.time()[3] - read_time
  betas  <- matrix(0, nrow = length(indicies[[params$data_partner_id]]), ncol = 1)
  params$u             = u
  params$idx           = idx
  params$betas         = betas
  params$delta_beta_old = betas
  params$indicies      = indicies
  params$pReduct       = pReduct
  if (params$data_partner_id == 1) {
    params$loglikelihood <- -Inf
    params$sBeta_old = rep(0, params$n)
    params$cutoff    = cutoff
  }
  params <- add_to_log(params, "update_params_cox_DP", read_time, read_size, 0, 0)
  return(params)
}

#' @importFrom stats rnorm
update_data_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "update_data_cox_DP\n\n")
  idx <- params$indicies[[params$data_partner_id]]
  data$x <- data$x[, idx, drop = FALSE]
  x = data$x
  data$colmin = data$colmin[idx]
  data$colmax = data$colmax[idx]
  data$colsum = data$colsum[idx]
  data$colrange = data$colrange[idx]

  if (params$data_partner_id == 1) {
    halfshare = rep(list(list()), params$numDataPartners)
    for (id in 1:params$numDataPartners) {
      set.seed(params$seeds[id], kind = "Mersenne-Twister")
      halfshare[[id]]  <- matrix(rnorm(params$n * params$ps[id], sd = 20), nrow = params$n, ncol = params$ps[id])
      halfshare[[id]] = halfshare[[id]][, params$indicies[[id]], drop = FALSE]
    }
    data$halfshare <- do.call(cbind, halfshare)
  }
  return(data)
}

#' @importFrom stats rnorm runif
ComputeSBetaCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSBetaCox.DP\n\n")
  u = NULL
  n <- params$n
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "u.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "u.rdata"))
  read_time <- proc.time()[3] - read_time

  sBetaPart = (data$x %*% params$betas + u) / (2 * u)
  v = 0
  for (id in 1:params$numDataPartners) {
    set.seed(params$seeds[id] + params$alg_iteration_counter)
    v = rnorm(n, mean = runif(n = 1, min = -1, max = 1), sd = 20)
    v = v + v
    if (id == params$data_partner_id) {
      sBetaPart = sBetaPart + v
    }
  }

  sBetaPart = sBetaPart - params$scalers[params$data_partner_id] / sum(params$scalers) * v

  write_time <- proc.time()[3]
  save(sBetaPart, file = file.path(params$write_path, "sbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeSBetaCox.DP", read_time, read_size, write_time, write_size)
  return(params)
}

GetSBetaCox.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetSBetaCox.AC\n\n")
  sBeta = 0
  sBetaPart = NULL
  read_time <- 0
  read_size <- 0
  for (id in 1:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_size <- read_size + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
    read_time <- read_time + proc.time()[3]
    sBeta = sBeta + sBetaPart
  }
  sBeta = 2 * params$u * sBeta - params$u * params$numDataPartners
  params$sBeta = sBeta
  write_time <- proc.time()[3]
  save(sBeta, file = file.path(params$write_path, "sbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetBetaCox.AC", read_time, read_size, write_time, write_size)
  return(params)
}

compute_log_likelihood_cox_DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "compute_log_likelihood_cox_DP\n\n")

  sBeta = 0
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "sbeta.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "sbeta.rdata"))
  read_time <- proc.time()[3] - read_time
  sBeta = sBeta[params$survival$sortedIdx]
  loglikelihood_old = params$loglikelihood
  stephalving = TRUE
  scale = 1
  while (stephalving) {
    w = exp(sBeta)
    loglikelihood <- 0
    step_counter <- 0
    pbar = MakeProgressBar1(params$survival$num_events, "loglikelihood", params$verbose)
    for (i in 1:length(params$survival$strata)) {
      if (params$survival$strata[[i]]$J > 0) {
        for (j in 1:params$survival$strata[[i]]$J) {
          nj = params$survival$strata[[i]]$nfails[j]
          y_index = params$survival$strata[[i]]$start0[j]:params$survival$strata[[i]]$end
          zIndex = params$survival$strata[[i]]$start1[j]:params$survival$strata[[i]]$stop1[j]
          a_j1 = sum(w[y_index])
          a_j2 = sum(w[zIndex]) / nj
          loglikelihood <- loglikelihood + sum(log(w[zIndex]))
          for (r in 0:(nj - 1)) {
            a_jr = a_j1 - r * a_j2
            loglikelihood <- loglikelihood - log(a_jr)
          }
          step_counter <- step_counter + nj
          pbar = MakeProgressBar2(step_counter, pbar, params$verbose)
        }
      }
    }
    if (loglikelihood < loglikelihood_old && scale > 0.5^6) {
      sBeta = 0.5 * (sBeta + params$sBeta)
      scale <- scale / 2
      if (params$verbose) cat("Step Halving\n\n")
    } else {
      stephalving = FALSE
    }
  }

  if (params$alg_iteration_counter == 1) {
    params$nullloglikelihood <- loglikelihood
  }

  converged = abs(loglikelihood - loglikelihood_old) / (abs(loglikelihood) + 0.1) < params$cutoff
  params$converged = converged
  params$scale     <- scale
  params$sBeta     = sBeta
  params$loglikelihood <- loglikelihood

  write_time <- proc.time()[3]
  save(sBeta, file = file.path(params$write_path, "sbeta.rdata"))
  save(scale, converged, file = file.path(params$write_path, "converged.rdata"))
  write_size <- file.size(file.path(params$write_path, "converged.rdata")) +
    file.size(file.path(params$write_path, "sbeta.rdata"))
  write_time <- proc.time()[3] - proc.time()[3]

  params <- add_to_log(params, "compute_log_likelihood_cox_DP", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeSDelLCox.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSDeltaCox.AC\n\n")
  sBeta = 0
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "sbeta.rdata"))
  read_size <- file.size(file.path(params$readPathDP[1], "sbeta.rdata"))
  read_time <- proc.time()[3] - read_time
  halfshare = params$halfshare[params$survival$sortedIdx, , drop = FALSE]
  p = ncol(halfshare)
  n <- params$n
  w = exp(sBeta)

  deltal = as.numeric(params$survival$status)
  deltal[1] <- deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  w_s_r  <- matrix(0, n, p)
  num_events <- params$survival$num_events

  .Call("ComputeCox", params$survival$strata, halfshare, w, deltal, w_s_r,
        as.integer(n), as.integer(p), as.integer(num_events),
        as.integer(params$verbose))

  w_s_r = w_s_r[params$survival$unsortIdx, , drop = FALSE]
  tS.deltal.R = t(halfshare) %*% deltal
  w_s_r_1 = w_s_r[, params$idx[[1]], drop = FALSE]
  params$w_s_r = w_s_r
  params$tS.deltal.R = tS.deltal.R

  write_time <- proc.time()[3]
  save(w_s_r_1, file = file.path(params$write_path, "wsr1.rdata"))
  write_size <- file.size(file.path(params$write_path, "wsr1.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeSDelLCox.AC", read_time, read_size, write_time, write_size)
  return(params)
}


#' @importFrom stats rnorm
ComputeSDelLCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSDelLCox.DP\n\n")
  halfshare = data$halfshare[params$survival$sortedIdx, , drop = FALSE]
  p = ncol(halfshare)
  n <- params$n
  w = exp(params$sBeta)

  deltal = as.numeric(params$survival$status)
  deltal[1] <- deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  w_s_l  <- matrix(0, n, p)
  num_events <- params$survival$num_events

  .Call("ComputeCox", params$survival$strata, halfshare, w, deltal, w_s_l,
        as.integer(n), as.integer(p), as.integer(num_events),
        as.integer(params$verbose))

  w_s_l = w_s_l[params$survival$unsortIdx, , drop = FALSE]
  tS.deltal.L = t(halfshare) %*% deltal
  params$w_s_l = w_s_l
  params$tS.deltal.L = tS.deltal.L

  colmin.w_s_l   = apply(w_s_l, 2, min)
  colmax.w_s_l   = apply(w_s_l, 2, max)
  colrange.w_s_l = colmax.w_s_l - colmin.w_s_l

  for (i in 1:ncol(w_s_l)) {
    if (colmin.w_s_l[i] == colmax.w_s_l[i]) {
      colmin.w_s_l[i] = 0
    }
    if (colrange.w_s_l[i] == 0) {
      colrange.w_s_l[i] = 1
    }
  }

  scaled.w_s_l = w_s_l

  for (i in 1:ncol(scaled.w_s_l)) {
    scaled.w_s_l[, i] = (scaled.w_s_l[, i] - colmin.w_s_l[i]) / colrange.w_s_l[i]
  }

  temp = as.numeric(Sys.time())
  set.seed((temp - trunc(temp)) * .Machine$integer.max)
  scaled.w_s_l.L  <- matrix(rnorm(n * p, sd = 20), nrow = n, ncol = p)
  scaled.w_s_l.R <- scaled.w_s_l - scaled.w_s_l.L

  params$colmin.w_s_l = colmin.w_s_l
  params$colrange.w_s_l = colrange.w_s_l
  params$scaled.w_s_l.L <- scaled.w_s_l.L

  write_time <- proc.time()[3]
  save(colmin.w_s_l, colrange.w_s_l, scaled.w_s_l.R,
       file = file.path(params$write_path, "scaledwslr.rdata"))
  save(tS.deltal.L, file = file.path(params$write_path, "tsdeltal.rdata"))
  save(scaled.w_s_l.L, file = file.path(params$write_path, "scaledwsll.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "scaledwslr.rdata"),
                            file.path(params$write_path, "tsdeltal.rdata"),
                            file.path(params$write_path, "scaledwsll.rdata")))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeSDelLCox.DP", 0, 0, write_time, write_size)

  return(params)
}


#' @importFrom stats rnorm
ComputeProductsCox.DP <- function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeProductsCox.DP\n\n")
  if (params$data_partner_id == 1) {
    w_s_r_1 = NULL
    read_time <- proc.time()[3]
    load(file.path(params$readPathAC, "wsr1.rdata"))
    read_size <- file.size(file.path(params$readPathAC, "wsr1.rdata"))
    read_time <- proc.time()[3] - read_time
    E = rep(list(list()), params$numDataPartners)
    E1 = rep(list(matrix(0, 0, 0)), params$numDataPartners)
    p = length(params$idx[[1]])
    if (p == 0) {
      for (id2 in 1:params$numDataPartners) {
        E1[[id2]]  <- matrix(0, nrow = 0, ncol = length(params$idx[[id2]]))
      }
    } else {
      E1[[1]] <- t(data$x) %*% (params$w_s_l[, params$idx[[1]], drop = FALSE] + w_s_r_1)
      D = diag(params$colrange.w_s_l[params$idx[[1]]], ncol = p, nrow = p)
      for (id2 in 2:params$numDataPartners) {
        set.seed(params$seeds[id2], kind = "Mersenne-Twister")
        halfshare.RL  <- matrix(rnorm(params$n * params$ps[id2], sd = 20),
                              nrow = params$n, ncol = params$ps[id2])  # needed to get randomization to right spot
        halfshare.RL  <- matrix(rnorm(params$n * params$ps[id2], sd = 20),
                              nrow = params$n, ncol = params$ps[id2])[, params$indicies[[id2]], drop = FALSE]
        E1[[id2]] = solve(D) %*% t(data$x) %*% params$w_s_l[, params$idx[[id2]], drop = FALSE] +
          params$scaler / (params$scaler + params$scalers[id2]) *
          t(params$scaled.w_s_l.L[, params$idx[[1]], drop = FALSE]) %*% halfshare.RL
      }
    }
    E[[1]] = E1
    for (id1 in 2:params$numDataPartners) {
      E1 = rep(list(matrix(0, 0, 0)), params$numDataPartners)
      p = length(params$idx[[id1]])
      D = diag(params$colrange.w_s_l[params$idx[[id1]]], ncol = p, nrow = p)
      set.seed(params$seeds[id1], kind = "Mersenne-Twister")
      halfshare.L  <- matrix(rnorm(params$n * params$ps[id1], sd = 20),
                           nrow = params$n, ncol = params$ps[id1])[, params$indicies[[id1]], drop = FALSE]
      for (id2 in 2:params$numDataPartners) {
        set.seed(params$seeds[id2], kind = "Mersenne-Twister")
        halfshare.RL  <- matrix(rnorm(params$n * params$ps[id2], sd = 20),
                              nrow = params$n, ncol = params$ps[id2])
        halfshare.RL  <- matrix(rnorm(params$n * params$ps[id2], sd = 20),
                              nrow = params$n, ncol = params$ps[id2])[, params$indicies[[id2]], drop = FALSE]
        E1[[id2]] = 0.5 * solve(D) %*% t(halfshare.L) %*% params$w_s_l[, params$idx[[id2]], drop = FALSE] +
          params$scaler / (params$scaler + params$scalers[id2]) *
          t(params$scaled.w_s_l.L[, params$idx[[id1]], drop = FALSE]) %*% halfshare.RL
      }
      E[[id1]] = E1
    }
    write_time <- proc.time()[3]
    save(E, file = file.path(params$write_path, "products.rdata"))
    write_size <- file.size(file.path(params$write_path, "products.rdata"))
    write_time <- proc.time()[3] - write_time
  } else {
    scaled.w_s_l.L = NULL
    read_time <- proc.time()[3]
    load(file.path(params$readPathDP[1], "scaledwsll.rdata"))
    read_size <- file.size(file.path(params$readPathDP[1], "scaledwsll.rdata"))
    read_time <- proc.time()[3] - read_time

    F1 = rep(list(list()), params$numDataPartners)
    set.seed(params$seed, kind = "Mersenne-Twister")
    halfshare.L  <- matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)[, params$indicies[[params$data_partner_id]], drop = FALSE]
    halfshare.R = data$x - halfshare.L
    halfshare.R.L  <- matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)[, params$indicies[[params$data_partner_id]], drop = FALSE]
    halfshare.R.R = halfshare.R - halfshare.R.L  # This is halfshare.R.L and halfshare.R.R
    for (id in 1:params$numDataPartners) {
      F1[[id]] = params$scaler / (params$scalers[1] + params$scaler) * t(scaled.w_s_l.L[, params$idx[[id]], drop = FALSE]) %*% halfshare.R.L +
        t(scaled.w_s_l.L[, params$idx[[id]], drop = FALSE]) %*% halfshare.R.R
    }
    write_time <- proc.time()[3]
    save(F1, file = file.path(params$write_path, "products.rdata"))
    write_size <- file.size(file.path(params$write_path, "products.rdata"))
    write_time <- proc.time()[3] - write_time
  }
  params <- add_to_log(params, "ComputeProductsCox.DP", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateConvergeStatus.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateConvergeStatus.AC\n\n")
  converged = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "converged.rdata"))
  read_size <- file.size(file.path(params$readPathDP[1], "converged.rdata"))
  read_time <- proc.time()[3] - read_time
  params$converged = converged
  params <- add_to_log(params, "UpdateConvergeStatus.AC", read_time, read_size, 0, 0)
}


ComputeStWSCox.AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSCox.AC\n\n")
  E = NULL
  F1 = NULL
  scaled.w_s_l.R = NULL
  colmin.w_s_l = NULL
  colrange.w_s_l = NULL
  F1 = rep(list(list()), params$numDataPartners)
  read_time <- proc.time()[3]
  load(file.path(params$readPathDP[1], "scaledwslr.rdata"))
  load(file.path(params$readPathDP[1], "products.rdata"))
  read_size <- file.size(file.path(params$readPathDP[1], "scaledwslr.rdata")) +
    file.size(file.path(params$readPathDP[1], "products.rdata"))
  read_time <- proc.time()[3] - read_time
  for (id in 2:params$numDataPartners) {
    read_time <- read_time - proc.time()[3]
    load(file.path(params$readPathDP[id], "products.rdata"))
    read_size <- read_size + file.size(file.path(params$readPathDP[id], "products.rdata"))
    read_time <- read_time + proc.time()[3]
    F1[[id]] = F1
  }
  p = 0
  for (id in 1:params$numDataPartners) {
    p = p + length(params$idx[[id]])
  }
  m  <- matrix(0, p, p)
  startrow = 1
  for (id1 in 1:params$numDataPartners) {
    p1 = length(params$idx[[id1]])
    if (p1 == 0) next
    endrow = startrow + p1 - 1
    startcol = startrow
    D = diag(x = colrange.w_s_l[params$idx[[id1]]], nrow = p1, ncol = p1)
    for (id2  in id1:params$numDataPartners) {
      p2 = length(params$idx[[id2]])
      endcol = startcol + p2 - 1
      if (p2 == 0) next
      if (id1 == 1 && id2 == 1) {
        m[startrow:endrow, startcol:endcol] = E[[1]][[1]]
      } else if (id1 == id2) {
        idx <- params$idx[[id1]]
        G = D %*% (E[[id1]][[id1]] + F1[[id1]][[id1]] + t(scaled.w_s_l.R[, idx, drop = FALSE]) %*%
                     params$halfshare[, idx, drop = FALSE]) +
          outer(colmin.w_s_l[idx], params$halfsharecolsum[idx])
        m[startrow:endrow, startcol:endcol] = G + t(G) +
          t(params$halfshare[, idx, drop = FALSE]) %*% params$w_s_r[, idx, drop = FALSE]
      } else {
        idx1 = params$idx[[id1]]
        idx2 = params$idx[[id2]]
        if (id1 == 1) {
          temp = D %*% (E[[id1]][[id2]] + F1[[id2]][[id1]] + t(scaled.w_s_l.R[, idx1, drop = FALSE]) %*%
                          params$halfshare[, idx2, drop = FALSE]) +
            t(params$halfshare[, idx1, drop = FALSE]) %*% params$w_s_r[, idx2, drop = FALSE] +
            outer(colmin.w_s_l[idx1], params$halfsharecolsum[idx2])
          m[startrow:endrow, startcol:endcol] = temp
          m[startcol:endcol, startrow:endrow] = t(temp)
        } else {
          p1 = length(params$idx[[id2]])
          D2 = diag(x = colrange.w_s_l[params$idx[[id2]]], nrow = p2, ncol = p2)
          G23 = D %*% (E[[id1]][[id2]] + F1[[id2]][[id1]] + t(scaled.w_s_l.R[, idx1, drop = FALSE]) %*%
                         params$halfshare[, idx2, drop = FALSE]) +
            outer(colmin.w_s_l[idx1], params$halfsharecolsum[idx2])
          G32 = D2 %*% (E[[id2]][[id1]] + F1[[id1]][[id2]] + t(scaled.w_s_l.R[, idx2, drop = FALSE]) %*%
                          params$halfshare[, idx1, drop = FALSE]) +
            outer(colmin.w_s_l[idx2], params$halfsharecolsum[idx1])
          temp = G23 + t(G32) + t(params$halfshare[, idx1, drop = FALSE]) %*% params$w_s_r[, idx2, drop = FALSE]
          m[startrow:endrow, startcol:endcol] = temp
          m[startcol:endcol, startrow:endrow] = t(temp)
        }
      }
      startcol = endcol + 1
    }
    startrow = endrow + 1
  }

  I = NULL
  tryCatch({
    I = solve(m)
  },
  error = function(err) {
    I = NULL
  }
  )
  if (is.null(I)) {
    params$failed <- TRUE
    params$singular_matrix = TRUE
    params$error_message <-
      paste0("The matrix t(S)*w*S is not invertible.\n",
             "       This may be due to one of two possible problems.\n",
             "       1. Poor random initialization of the security halfshares.\n",
             "       2. Near multicollinearity in the data\n",
             "SOLUTIONS: \n",
             "       1. Rerun the data analysis.\n",
             "       2. If the problem persists, check the variables for\n",
             "          duplicates for both parties and / or reduce the\n",
             "          number of variables used. Once this is done,\n",
             "          rerun the data analysis.")
    params <- add_to_log(params, "computeStWSCox.AC", read_time, read_size, 0, 0)
    return(params)
  }

  params$I = I

  IDt = I %*% params$tS.deltal.R

  if (params$alg_iteration_counter == 1) {
    params$score = t(params$tS.deltal.R) %*% I %*% params$tS.deltal.R
  }
  params$max_iter_exceeded = params$alg_iteration_counter > params$max_iterations
  max_iter_exceeded = params$max_iter_exceeded

  write_time <- 0
  write_size <- 0
  for (id in 1:params$numDataPartners) {
    I.part = I[params$idx[[id]], , drop = FALSE]
    IDt.part = IDt[params$idx[[id]], , drop = FALSE]
    write_time <- write_time - proc.time()[3]
    save(I.part, IDt.part, file = file.path(params$write_path, paste0("update", id, ".rdata")))
    save(max_iter_exceeded, file = file.path(params$write_path, "maxiterexceeded.rdata"))
    write_size <- write_size + file.size(file.path(params$write_path, paste0("update", id, ".rdata"))) +
      file.size(file.path(params$write_path, "maxiterexceeded.rdata"))
    write_time <- write_time + proc.time()[3]
  }
  params <- add_to_log(params, "ComputeStWSCox.AC", read_time, read_size, write_time, write_size)
  return(params)
}


UpdateConvergeStatus.DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateConvergeStatus.DP\n\n")
  read_time <- 0
  read_size <- 0
  scale    = NULL
  converged = NULL
  max_iter_exceeded = NULL
  if (params$data_partner_id > 1) {
    read_time <- proc.time()[3]
    load(file.path(params$readPathDP[1], "converged.rdata"))
    read_size <- file.size(file.path(params$readPathDP[1], "converged.rdata"))
    read_time <- proc.time()[3] - read_time
    params$converged = converged
    params$scale     <- scale
  }
  read_time <- read_time - proc.time()[3]
  load(file.path(params$readPathAC, "maxiterexceeded.rdata"))
  read_size <- read_size + file.size(file.path(params$readPathAC, "maxiterexceeded.rdata"))
  read_time <- read_time + proc.time()[3]
  params$max_iter_exceeded = max_iter_exceeded
  params <- add_to_log(params, "UpdateConvergeStatus.DP", read_time, read_size, 0, 0)
}


#' @importFrom stats runif
UpdateBetasCox.DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetasCox.DP\n\n")
  I.part = NULL
  IDt.part = NULL
  tS.deltal.L = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, paste0("update", params$data_partner_id, ".rdata")))
  read_size <- file.size(file.path(params$readPathAC, paste0("update", params$data_partner_id, ".rdata")))
  if (params$data_partner_id > 1) {
    load(file.path(params$readPathDP[1], "tsdeltal.rdata"))
    read_size <- read_size + file.size(file.path(params$readPathDP[1], "tsdeltal.rdata"))
  }
  read_time <- proc.time()[3] - read_time

  if (params$data_partner_id == 1) {
    deltabeta = IDt.part + I.part %*% params$tS.deltal.L
    if (params$alg_iteration_counter == 1) {
      if (params$pReduct[1] == 0) {
        params$score = 0
      } else {
        idx <- 1:params$pReduct[1]
        params$score = 2 * t(IDt.part) %*% params$tS.deltal.L[idx, 1, drop = FALSE] +
          t(I.part %*% params$tS.deltal.L) %*% params$tS.deltal.L[idx, 1, drop = FALSE]
      }
    }
  } else {
    deltabeta = IDt.part + I.part %*% tS.deltal.L
    if (params$alg_iteration_counter == 1) {
      temp = sum(params$pReduct[1:(params$data_partner_id - 1)])
      idx <- (temp + 1):(temp + params$pReduct[params$data_partner_id])
      params$score = 2 * t(IDt.part) %*% tS.deltal.L[idx, 1, drop = FALSE] +
        t(I.part %*% tS.deltal.L) %*% tS.deltal.L[idx, 1, drop = FALSE]
    }
  }

  params$betas = params$betas + deltabeta - (1 - params$scale) * params$delta_beta_old
  params$delta_beta_old = deltabeta

  p = length(params$indicies[[params$data_partner_id]])
  if (p == 0) {
    u = 1
  } else {
    u = sum(runif(n = p, min = 1, max = 2) * abs(params$betas))
  }

  if (params$data_partner_id == 1) {
    loglikelihood <- params$loglikelihood
    nullloglikelihood <- params$nullLoglikelihood
  }
  write_time <- proc.time()[3]
  save(u, file = file.path(params$write_path, "u.rdata"))
  write_size <- file.size(file.path(params$write_path, "u.rdata"))
  if (params$converged) {
    betasnew  = params$betas
    scorePart = params$score
    if (params$data_partner_id == 1) {
      save(scorePart, nullLoglikelihood, loglikelihood, betasnew, file = file.path(params$write_path, "betas.rdata"))
    } else {
      save(scorePart, betasnew, file = file.path(params$write_path, "betas.rdata"))
    }
    write_size <- write_size + file.size(file.path(params$write_path, "betas.rdata"))
  }
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "UpdateBetasCox.DP", read_time, read_size, write_time, write_size)

  return(params)
}


survfit_cox_AC <- function(params, pred) {
  if (params$trace) cat(as.character(Sys.time()), "survfit_cox_AC\n\n")
  survival = params$survival
  surv = rep(1, length(survival$rank))
  for (i in 1:length(survival$strata)) {
    if (survival$strata[[i]]$J > 0) {
      start   = survival$strata[[i]]$start
      end     = survival$strata[[i]]$end
      risk    = exp(pred[start:end])
      dtime   = survival$rank[start:end]
      status  = survival$status[start:end]
      death   = status == 1                                  # times where death happened
      time    = sort(unique(dtime))                          # (A) get unique event times
      rcumsum <- function(x) rev(cumsum(rev(x)))
      nevent  = as.vector(rowsum(as.numeric(death), dtime))  # (A) Count the number of deaths at each event time
      ndeath  = rowsum(status, dtime)                        # (A) number of deaths at each unique event time
      nrisk   = rcumsum(rowsum(risk, dtime))                 # (A) rowsum = sum of risk at each time, then reverse cum sum, sorted by time
      erisk   = rowsum(risk * death, dtime)                  # (A) risk score sums of death at each unique event time
      n       = length(nevent)
      sum1 = double(n)   # a vector of 0's, length number of unique event times
      for (i in 1:n) {
        d = ndeath[i]
        if (d == 1) {
          sum1[i] = 1 / nrisk[i]
        } else if (d > 1) {
          for (j in 0:(d - 1)) {
            sum1[i] = sum1[i] + 1 / (d * nrisk[i] - erisk[i] * j)
          }
        }
      }
      temp = exp(-cumsum(nevent * sum1))
      for (i in start:end) {
        surv[i] = temp[which(time == survival$rank[i])]
      }
    }
  }
  return(surv)
}


#' @importFrom  stats pchisq pnorm qnorm
compute_results_cox_AC <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "compute_results_cox_AC\n\n")
  read_size          = 0
  betasnew          = NULL
  scorePart         = NULL
  loglikelihood     = NULL
  nullloglikelihood <- NULL
  score             = params$score
  read_time <- proc.time()[3]
  betas  <- matrix(0, nrow = 0, ncol = 1)
  for (id in 1:params$numDataPartners) {
    load(file.path(params$readPathDP[id], "betas.rdata"))
    betas = rbind(betas, betasnew)
    score = score + scorePart
    read_size <- read_size + file.size(file.path(params$readPathDP[id], "betas.rdata"))
  }
  read_time <- proc.time()[3] - read_time

  stats = params$stats
  stats$failed         = FALSE
  stats$converged      = params$converged
  names_old            = params$colnames[-(1:params$pStrata)]
  idx                  = params$fullindicies - params$pStrata
  stats$party          = params$party[-(1:params$pStrata)]
  stats$coefficients   = rep(NA, length(stats$party))
  stats$coefficients[idx] = betas / params$colrange
  stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
  stats$expncoef     = exp(-stats$coefficients)
  stats$var           <- matrix(0, length(names_old), length(names_old))
  stats$var[idx, idx] = params$I
  stats$secoef       = rep(NA, length(names_old))
  stats$secoef[idx]  = sqrt(diag(params$I)) / params$colrange  # se(coef)

  stats$zvals        = stats$coefficients / stats$secoef  # z values
  stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE)   # pvals
  stats$stars         <- matrix(sapply(stats$pvals, function(x) {
    if (is.na(x)) ""
    else if (x < 0.001) "***"
    else if (x < 0.01) "**"
    else if (x < 0.05) "*"
    else if (x < 0.1)  "."
    else " "
  }))
  stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  stats$loglik       = c(nullLoglikelihood, loglikelihood)
  stats$n            = params$n
  stats$nevent       = sum(params$survival$status)
  stats$df           = sum(params$pReduct)
  stats$iter         = params$alg_iteration_counter - 1
  stats$score        = c(score, 1 - pchisq(score, stats$df))
  stats$method       = "efron"
  stats$lrt          = 2 * (stats$loglik[2] - stats$loglik[1])
  stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      = c(1 - exp(-stats$lrt[1] / stats$n),
                         1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    = t(betas) %*% solve(params$I) %*% betas
  stats$wald.test    = c(stats$wald.test,
                         1 - pchisq(stats$wald.test, stats$df))
  pred = -params$sBeta[params$survival$sortedIdx]
  if (requireNamespace("survival", quietly = TRUE)) {
    if (!("package:survival" %in% search())) {
      attachNamespace("survival")
    }
    surv = survival::Surv(params$survival$rank, params$survival$status)
    strat = rep(0, length(surv))
    for (i in 1:length(params$survival$strata)) {
      strat[params$survival$strata[[i]]$start:params$survival$strata[[i]]$end] = i
    }
    results = survival::concordance(surv ~ pred + strata(strat))
    if (class(results$stats) == "matrix") {  # more than one strata
      stats$concordance = c(apply(results$count, 2, sum)[1:4], results$concordance, sqrt(results$var))
    } else {                                 # only one strata, so a numeric vector
      stats$concordance = c(results$count[1:4], results$concordance, sqrt(results$var))
    }
  } else {
    stats$concordance = c(NA, NA, NA, NA, NA, NA)
  }

  stats$survival <- data.frame(
    rank   = params$survival$rank,
    status = params$survival$status,
    sorted = params$survival$sortedIdx,
    surv   = survfit_cox_AC(params, pred)
  )
  stats$strata = as.data.frame(matrix(0, length(params$survival$strata), 3))
  stats$strata$label <- ""
  colnames(stats$strata) = c("start", "end", "events", "label")
  for (i in 1:length(params$survival$strata)) {
    stats$strata$start[i]  = params$survival$strata[[i]]$start
    stats$strata$end[i]    = params$survival$strata[[i]]$end
    stats$strata$events[i] = sum(params$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  = params$survival$strata[[i]]$label
  }

  names(stats$party)           = names_old
  names(stats$coefficients)    = names_old
  names(stats$expcoef)         = names_old
  names(stats$expncoef)        = names_old
  rownames(stats$var)          = names_old
  colnames(stats$var)          = names_old
  names(stats$secoef)          = names_old
  names(stats$zvals)           = names_old
  names(stats$pvals)           = names_old
  names(stats$stars)           = names_old
  names(stats$lower95)         = names_old
  names(stats$upper95)         = names_old
  names(stats$loglik)          = c("loglikelihood", "null loglikelihood")
  names(stats$score)           = c("score", "p-value")
  names(stats$lrt)             = c("likelihood ratio", "p-value")
  names(stats$rsquare)         = c("r-square", "max possible")
  names(stats$wald.test)       = c("wald", "p-value")
  names(stats$concordance)     = c("concordant", "discordant", "tied.risk", "tied.time",
                                   "concordance", "stderr")

  params$stats = stats
  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "compute_results_cox_AC", read_time, read_size, write_time, write_size)
  return(params)
}


get_results_cox_DP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "get_results_cox_DP\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$readPathAC, "stats.rdata"))
  read_size <- file.size(file.path(params$readPathAC, "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats

  params <- add_to_log(params, "get_results_cox_DP", read_time, read_size, 0, 0)
  return(params)
}


DoNothing.ACDP <- function(params) {
  if (params$trace) cat(as.character(Sys.time()), "DoNothing\n\n")
  params <- add_to_log(params, "--", 0, 0, 0, 0)
  return(params)
}

############################## PARENT FUNCTIONS ###############################


DataPartnerKCox <- function(data,
                           y_name           = NULL,
                           strata          = NULL,
                           mask            = TRUE,
                           numDataPartners = NULL,
                           data_partner_id   = NULL,
                           monitor_folder   = NULL,
                           sleep_time       = 10,
                           max_waiting_time  = 24 * 60 * 60,
                           popmednet       = TRUE,
                           trace           = FALSE,
                           verbose         = TRUE) {

  params <- prepare_params_kp("cox", data_partner_id, numDataPartners, ac = FALSE,
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

  data = prepare_data_cox_DP(params, data, y_name, strata, mask)
  params <- add_to_log(params, "prepare_data_cox_DP", 0, 0, 0, 0)

  if (data$failed) {
    params$error_message <- paste("Error processing data for data partner", params$data_partner_id)
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

  if (params$data_partner_id == 1) {
    params <- DoNothing.ACDP(params)
    params <- send_pause_continue_kp(params, filesAC = "empty.rdata", from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

    params <- check_strata_cox_DP(params, data)

    if (params$failed) {
      make_error_message(params$write_path, params$error_message)
      files <- "error_message.rdata"
    } else {
      files <- "empty.rdata"
    }
    params <- send_pause_continue_kp(params, filesDP = files, from = "AC",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

    params <- DoNothing.ACDP(params)

    if (params$failed) {
      warning(params$error_message)
      send_pause_quit_kp(params, filesAC = "error_message.rdata", sleep_time = sleep_time)
      return(params$stats)
    } else {
      params <- send_pause_continue_kp(params, filesDP = "empty.rdata", from = "DP",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    }
    params <- prepare_strata_cox_DP(params, data)
    data   = AddStrataToDataCox.DP(params, data)
    params <- add_to_log(params, "AddStrataToDataCox.DP", 0, 0, 0, 0)
  } else {
    params <- SendStrataNamesCox.DP(params, data)
    filesList = rep(list(list()), numDataPartners)
    filesList[[1]] = "strata_names.rdata"
    params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP1",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)


    possibleError = ReceivedError.kp(params, from = "DP1")
    if (possibleError$error) {
      params$error_message <- possibleError$message
      warning(possibleError$message)
      params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
      return(params$stats)
    }

    params <- send_strata_cox_DP(params, data)
    filesList = rep(list(list()), numDataPartners)
    filesList[[1]] = "strata.rdata"
    params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP1",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)
  }

  params <- prepare_params_cox_DP(params, data)
  params <- send_pause_continue_kp(params, filesDP = "p_scaler_seed.rdata", from = "DP",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  params <- PrepareSharesCox.DP(params, data)
  files <- c("products.rdata", "halfshare.rdata", "colstats.rdata")
  if (params$data_partner_id == 1) {
    files <- c(files, "survival.rdata")
  }
  params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

  possibleError = ReceivedError.kp(params, from = "AC")
  if (possibleError$error) {
    params$error_message <- possibleError$message
    warning(possibleError$message)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- update_params_cox_DP(params)
  data = update_data_cox_DP(params, data)
  params <- add_to_log(params, "update_data_cox_DP", 0, 0, 0, 0)

  params$alg_iteration_counter = 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- ComputeSBetaCox.DP(params, data)

    if (params$data_partner_id == 1) {
      files <- "sbeta.rdata"
      params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

      params <- compute_log_likelihood_cox_DP(params, data)

      files <- "sbeta.rdata"
      params <- send_pause_continue_kp(params, filesAC = files, from = "DP2",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

      params <- ComputeSDelLCox.DP(params, data)

      files <- c("tsdeltal.rdata", "scaledwsll.rdata", "converged.rdata")
      params <- send_pause_continue_kp(params, filesDP = files, from = "AC",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)
    } else if (params$data_partner_id == 2) {
      files <- "sbeta.rdata"
      params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

      params <- DoNothing.ACDP(params)

      filesList = rep(list(list()), numDataPartners)
      filesList[[1]] = "empty.rdata"
      params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP1",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)
    } else {
      files <- "sbeta.rdata"
      params <- send_pause_continue_kp(params, filesAC = files, from = "DP1",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    }

    params <- ComputeProductsCox.DP(params, data)
    if (params$data_partner_id == 1) {
      files <- c("products.rdata", "scaledwslr.rdata", "converged.rdata")
    } else {
      files <- c("products.rdata")
    }
    params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

    possibleError = ReceivedError.kp(params, from = "AC")
    if (possibleError$error) {
      params$error_message <- possibleError$message
      warning(possibleError$message)
      params <- send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
      return(params$stats)
    }

    params <- UpdateConvergeStatus.DP(params)
    params <- UpdateBetasCox.DP(params)

    if (params$converged || params$max_iter_exceeded) {
      files <- "betas.rdata"
    } else {
      files <- "u.rdata"
    }
    params <- send_pause_continue_kp(params, filesAC = files, from = "AC",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)
    EndingIteration(params)
    params$alg_iteration_counter = params$alg_iteration_counter + 1
  }
  params$lastIteration = TRUE
  params$completed = TRUE

  params <- get_results_cox_DP(params)
  send_pause_quit_kp(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


AnalysisCenterKCox <- function(numDataPartners = NULL,
                              monitor_folder   = NULL,
                              msreqid         = "v_default_0_000",
                              cutoff          = 1E-8,
                              max_iterations   = 25,
                              sleep_time       = 10,
                              max_waiting_time  = 24 * 60 * 60,
                              popmednet       = TRUE,
                              trace           = FALSE,
                              verbose         = TRUE) {

  filesList = rep(list(list()), numDataPartners)

  params <- prepare_params_kp("cox", 0, numDataPartners, msreqid, cutoff, max_iterations, ac = TRUE,
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
  params <- send_pause_continue_kp(params, filesDP = files, from = "DP1",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)


  params <- DoNothing.ACDP(params)
  filesList = rep(list(list()), numDataPartners)
  filesList[[1]] = "empty.rdata"
  params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  possibleError = ReceivedError.kp(params, from = "DP")
  if (possibleError$error) {
    params$error_message <- possibleError$message
    warning(possibleError$message)
    params <- send_pause_quit_kp(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.kp(params)
    return(params$stats)
  }

  params <- GetStrata.AC(params)
  params <- get_products_cox_AC(params)
  params <- check_colinearity_cox_AC(params)

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

  params$alg_iteration_counter = 1
  while (!params$converged && !params$max_iter_exceeded) {
    BeginningIteration(params)
    params <- ComputeUCox.AC(params)

    if (params$alg_iteration_counter == 1) {
      files <- c("indicies.rdata", "u.rdata")
    } else {
      files <- "u.rdata"
    }
    params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)

    params <- GetSBetaCox.AC(params)
    filesList = rep(list(list()), numDataPartners)
    filesList[[1]] = "sbeta.rdata"
    filesList[[2]] = "empty.rdata"
    params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP1",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)

    params <- ComputeSDelLCox.AC(params)
    filesList = rep(list(list()), numDataPartners)
    filesList[[1]] = "wsr1.rdata"
    params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time, waitForTurn = TRUE)

    params <- UpdateConvergeStatus.AC(params)
    params <- ComputeStWSCox.AC(params)

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

    filesList = rep(list(list()), numDataPartners)
    for (id in 1:params$numDataPartners) {
      filesList[[id]] = c(paste0("update", id, ".rdata"), "maxiterexceeded.rdata")
    }
    params <- send_pause_continue_kp(params, filesDP = filesList, from = "DP",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    EndingIteration(params)
    params$alg_iteration_counter = params$alg_iteration_counter + 1
  }
  params$lastIteration = TRUE
  params$completed = TRUE

  params <- compute_results_cox_AC(params)
  files <- "stats.rdata"
  params <- send_pause_continue_kp(params, filesDP = files, from = "DP",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)
  send_pause_quit_kp(params, sleep_time = sleep_time)
  SummarizeLog.kp(params)
  return(params$stats)
}
