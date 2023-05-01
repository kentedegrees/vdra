#################### DISTRIBUTED COX REGRESSION FUNCTIONS ####################

prepare_params_cox_a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_a3\n\n")
  params$n = nrow(data$X)
  params$p1 = ncol(data$X)
  params$p2 = 0
  params$colnames = colnames(data$X)

  params$survival_installed = requireNamespace("survival", quietly = TRUE)
  if (params$survival_installed & !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  pa           = list()
  pa$p1        = params$p1
  pa$n         = params$n
  pa$analysis  = params$analysis
  pa$colnames  = params$colnames
  pa$strata_from_a = data$strata$strata_from_a
  pa$strata_from_b = data$strata$strata_from_b
  pa$tags      = data$tags

  write_time <- proc.time()[3]
  save(pa, file = file.path(params$write_path, "pa.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pa.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_cox_a3", 0, 0, write_time, write_size)
  return(params)
}


prepare_params_cox_b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_b3\n\n")
  params$n = nrow(data$X)
  params$p1 = 0
  params$p2 = ncol(data$X)
  params$colnames = colnames(data$X)

  params$survival_installed = requireNamespace("survival", quietly = TRUE)
  if (params$survival_installed & !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  pb           = list()
  pb$p2        = params$p2
  pb$n         = params$n
  pb$analysis  = params$analysis
  pb$colnames  = params$colnames
  pb$strata_from_a = data$strata$strata_from_a
  pb$strata_from_b = data$strata$strata_from_b
  pb$tags      = data$tags

  write_time <- proc.time()[3]
  save(pb, file = file.path(params$write_path, "pb.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, "pb.rdata")))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_cox_b3", 0, 0, write_time, write_size)
  return(params)
}


prepare_params_cox_t3 = function(params, cutoff, max_iterations) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_params_cox_t3\n\n")
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

  params$survival_installed = requireNamespace("survival", quietly = TRUE)
  if (params$survival_installed & !("package:survival" %in% search())) {
    attachNamespace("survival")
  }

  params$n             = pa$n
  params$p1            = pa$p
  params$p2            = pb$p
  params$p1.old        = params$p1
  params$p2.old        = params$p2
  params$p             = pa$p + pb$p
  params$colnamesA     = pa$colnames
  params$colnamesB     = pb$colnames
  params$cutoff        = cutoff
  params$max_iterations = max_iterations

  params$Astrata_from_a  = pa$strata_from_a
  params$Astrata_from_b  = pa$strata_from_b
  params$Bstrata_from_a  = pb$strata_from_a
  params$Bstrata_from_b  = pb$strata_from_b

  params$Atags         = pa$tags
  params$Btags         = pb$tags

  write_time <- proc.time()[3]
  save(cutoff, max_iterations, file = file.path(params$write_path, "max_iterations.rdata"))
  write_size <- file.size(file.path(params$write_path, "max_iterations.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "prepare_params_cox_t3", read_time, read_size, write_time, write_size)
  return(params)
}


check_strata_cox_t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "check_strata_cox_t3\n\n")
  if (length(params$Astrata_from_a) == length(params$Bstrata_from_a) &&
      length(params$Astrata_from_b) == length(params$Bstrata_from_b) &&
      ifelse(length(params$Astrata_from_a) == 0, TRUE,
             order(params$Astrata_from_a) == order(params$Bstrata_from_a)) &&
      ifelse(length(params$Astrata_from_b) == 0, TRUE,
             order(params$Astrata_from_b) == order(params$Bstrata_from_b))) {
    params$strata_from_a = params$Astrata_from_a
    params$strata_from_b = params$Bstrata_from_b
    params$Astrata_from_a = params$Astrata_from_b =
      params$Bstrata_from_a = params$Bstrata_from_b = NULL
    params$getStrata = TRUE
  } else {
    params$getStrata = FALSE
    a_cap_b = intersect(params$Astrata_from_a, params$Bstrata_from_b)
    b_cap_a = intersect(params$Bstrata_from_a, params$Astrata_from_b)
    if (length(a_cap_b) > 0) {
      params$error_message <-
        paste("Party A and Party B have", length(a_cap_b), "variable(s) with the same name which are used in the strata.",
              "These variable(s) are <", paste0(a_cap_b, collapse = ", "), ">.",
              "Make sure the variables from each party have distinct names.")
    } else if (length(b_cap_a) > 0) {
      params$error_message <-
        paste("Party A and Party B have specified", length(b_cap_a), "variable(s) for the strata which are not found in the data.",
              "These variable(s) are <", paste0(b_cap_a, collapse = ", "), ">.",
              "Check the spelling of the variables names and / or remove them from the strata.")
    } else {
      params$error_message <-
        paste("Party A and Party B have specified different strata.",
              "Verify that both parties specify the same strata.")
    }
    params$failed <- TRUE
  }
  empty = 0
  save(empty, file = file.path(params$write_path, "empty.rdata"))
  params <- add_to_log(params, "check_strata_cox_t3", 0, 0, 0, 0)
  return(params)
}


send_strata_cox_a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "send_strata_cox_a3\n\n")
  Astrata = data$strata
  survival = data$survival
  write_time <- proc.time()[3]
  save(Astrata, survival, file = file.path(params$write_path, "Astrata.rdata"))
  write_size <- file.size(file.path(params$write_path, "Astrata.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "send_strata_cox_a3", 0, 0, write_time, write_size)
  return(params)
}


send_strata_cox_b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "send_strata_cox_b3\n\n")
  Bstrata = data$strata
  write_time <- proc.time()[3]
  save(Bstrata, file = file.path(params$write_path, "Bstrata.rdata"))
  write_size <- file.size(file.path(params$write_path, "Bstrata.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "send_strata_cox_b3", 0, 0, write_time, write_size)
  return(params)
}


prepare_strata_cox_t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "prepare_strata_cox_t3\n\n")
  Astrata  = NULL
  Bstrata  = NULL
  survival = NULL
  strata_temp   = list()

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "Astrata.rdata"))
  load(file.path(params$read_path[["B"]], "Bstrata.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "Astrata.rdata")) +
    file.size(file.path(params$read_path[["B"]], "Bstrata.rdata"))
  read_time <- proc.time()[3] - read_time


  if (length(params$strata_from_a) == 0 && length(params$strata_from_b) == 0) {
    strata_temp$X = data.frame(const__ = rep(1, params$n))
    strata_temp$legend = FALSE
  } else if (length(params$strata_from_a) == 0) {
    strata_temp$X = Bstrata$X
    strata_temp$legend = Bstrata$legend
  } else if (length(params$strata_from_b) == 0) {
    strata_temp$X = Astrata$X
    strata_temp$legend = Astrata$legend
  } else {
    strata_temp$X = cbind(Astrata$X, Bstrata$X)
    strata_temp$legend = c(Astrata$legend, Bstrata$legend)
  }

  sorted = do.call("order", cbind(strata_temp$X, survival$rank, survival$status))
  strata_temp$X = strata_temp$X[sorted, , drop = FALSE]
  survival$rank   = survival$rank[sorted]
  survival$status = survival$status[sorted]
  survival$sorted = sorted
  ranks = which(apply(abs(apply(strata_temp$X, 2, diff)), 1, sum) > 0)
  ranks = c(ranks, nrow(strata_temp$X))
  names(ranks) = NULL
  strata = rep(list(list()), length(ranks))
  if (length(ranks) == 1 && colnames(strata_temp$X)[1] == "const__") {
    strata[[1]]$start = 1
    strata[[1]]$end   = as.integer(length(survival$rank))
    strata[[1]]$label = ""
  } else {
    start = 1
    for (i in  1:length(ranks)) {
      strata[[i]]$start = start
      strata[[i]]$end   = as.integer(ranks[i])
      label = ""
      for (j in 1:ncol(strata_temp$X)) {
        temp = colnames(strata_temp$X)[j]
        label = paste0(label, temp, "=", strata_temp$legend[[temp]][strata_temp$X[start, j]])
        if (j < ncol(strata_temp$X)) {
          label = paste0(label, ", ")
        }
      }
      strata[[i]]$label = label
      start = as.numeric(ranks[i]) + 1
    }
  }
  for (i in 1:length(strata)) {
    idx = strata[[i]]$start:strata[[i]]$end
    temp  = table(survival$rank[idx])
    m = length(temp)   # number of unique observed times, including where no one fails
    # Count the number of 0's and 1's for each observed time
    temp0 = table(survival$rank[idx], survival$status[idx])
    # Check if there are all 1's or all 0's .  If so, add them into the table.
    if (ncol(temp0) == 1) {
      if (which(c(0, 1) %in% colnames(temp0)) == 1) {
        temp0 = cbind(temp0, 0)
        colnames(temp0) = c("0", "1")
      } else {
        temp0 = cbind(0, temp0)
        colnames(temp0) = c("0", "1")
      }
    }
    strata[[i]]$J = as.integer(sum(temp0[, 2] > 0))         # number of distinct failure times
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
  params$survival = survival

  write_time <- proc.time()[3]
  save(survival, file = file.path(params$write_path, "survival.rdata"))
  write_size <- file.size(file.path(params$write_path, "survival.rdata"))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "prepare_strata_cox_t3", read_time, read_size, write_time, write_size)

  return(params)
}


sort_data_cox_a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "sort_data_cox_a3\n\n")
  survival = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "survival.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "survival.rdata"))
  read_time <- proc.time()[3] - read_time
  data$X = data$X[survival$sorted, , drop = FALSE]
  data$survival = survival
  data$read_time <- read_time
  data$read_size <- read_size
  return(data)
}


get_z_cox_a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "get_z_cox_a3\n\n")
  write_time <- 0
  write_size <- 0

  num_blocks = params$blocks$num_blocks
  pbar <- MakeProgressBar1(num_blocks, "Z", params$verbose)
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
    Z = FindOrthogonalVectors(data$X[strt:stp, ], g)

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(Z), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$file_break_z || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "get_z_cox_a3", 0, 0, write_time, write_size)
  return(params)
}


sort_data_cox_b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "sort_data_cox_b3\n\n")
  survival = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "survival.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "survival.rdata"))
  read_time <- proc.time()[3] - read_time
  data$X = data$X[survival$sorted, , drop = FALSE]
  data$survival = survival
  data$read_time <- read_time
  data$read_size <- read_size
  return(data)
}


GetSXBCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetSXBCox.b3\n\n")
  S = matrix(0, nrow = params$n, ncol = length(data$survival$strata))
  for (i in 1:length(data$survival$strata)) {
    S[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
  }
  STXB = t(S) %*% data$X
  write_time <- proc.time()[3]
  save(STXB, file = file.path(params$write_path, "sxb.rdata"))
  write_size <- file.size(file.path(params$write_path, "sxb.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetSXBCox.b3", 0, 0, write_time, write_size)
  return(params)
}


GetWRCox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWRCox.a3\n\n")
  XATXA = t(data$X) %*% data$X
  write_time <- proc.time()[3]
  save(XATXA, file = file.path(params$write_path, "xatxa.rdata"))
  write_size <- file.size(file.path(params$write_path, "xatxa.rdata"))
  write_time <- proc.time()[3] - write_time

  p2 = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "p2.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "p2.rdata"))
  read_time <- proc.time()[3] - read_time
  params$p2 = p2

  num_blocks = params$blocks$num_blocks
  pbar <- MakeProgressBar1(num_blocks, "XA'(I - Z*Z')XB*R", params$verbose)

  container_ct_wR = 0
  container_ct_PR = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.WR) {
      container_ct_wR = container_ct_wR + 1
      filename1 <- paste0("cwr_", container_ct_wR, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
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
    WR = matrix(readBin(con = to_read, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    PR = t(data$X[strt:stp, ]) %*% WR
    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(PR), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.WR || i == num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.PR || i == num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }

    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "GetWRCox.a3", read_time, read_size, write_time, write_size)
  return(params)
}


GetSXACox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetSXACox.a3\n\n")
  S = matrix(0, nrow = params$n, ncol = length(data$survival$strata))
  for (i in 1:length(data$survival$strata)) {
    S[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
  }
  STXA = t(S) %*% data$X
  write_time <- proc.time()[3]
  save(STXA, file = file.path(params$write_path, "sxa.rdata"))
  write_size <- file.size(file.path(params$write_path, "sxa.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetSXACox.a3", 0, 0, write_time, write_size)
  return(params)
}


GetProductsCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsCox.t3\n\n")
  p1 = params$p1
  p2 = params$p2
  XATXA = 0
  xb_t_xb <- 0
  XATXB = 0
  STXA  = 0
  STXB  = 0

  num_blocks = params$blocks$num_blocks
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "xatxa.rdata"))
  load(file.path(params$read_path[["B"]], "xbtxb.rdata"))
  load(file.path(params$read_path[["A"]], "sxa.rdata"))
  load(file.path(params$read_path[["B"]], "sxb.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["B"]], "xbtxb.rdata")),
                 file.size(file.path(params$read_path[["A"]], "xatxa.rdata")),
                 file.size(file.path(params$read_path[["B"]], "sxb.rdata")),
                 file.size(file.path(params$read_path[["A"]], "sxa.rdata")))
  read_time <- proc.time()[3] - read_time

  pbar <- MakeProgressBar1(num_blocks, "X'X", params$verbose)

  container_ct_PR = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.PR) {
      container_ct_PR = container_ct_PR + 1
      filename1 <- paste0("cpr_", container_ct_PR, ".rdata")
      to_read <- file(file.path(params$read_path[["A"]], filename1), "rb")
      read_size <- read_size + file.size(file.path(params$read_path[["A"]], filename1))
    }

    filename1 <- paste0("r2_", i, ".rdata")
    read_time <- read_time - proc.time()[3]
    toRead1 = file(file.path(params$dp_local_path, filename1), "rb")
    R2 = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
                        endian = "little"), p2, p2)
    read_size <- read_size + file.size(file.path(params$dp_local_path, filename1))
    close(toRead1)
    PR = matrix(readBin(con = to_read, what = numeric(), n = p1 * p2,
                        endian = "little"), p1, p2)
    read_time <- read_time + proc.time()[3]
    XATXB = XATXB + PR %*% t(R2)
    if ((i + 1) %in% params$container$filebreak.PR || i == num_blocks) {
      close(to_read)
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }

  num = length(params$survival$strata)
  STS = matrix(0, nrow = num, ncol = num)
  for (i in 1:num) {
    STS[i, i] = params$survival$strata[[i]]$end - params$survival$strata[[i]]$start + 1
  }

  XTX = rbind(cbind(STS, STXA, STXB),
              cbind(t(STXA), XATXA, XATXB),
              cbind(t(STXB), t(XATXB), xb_t_xb))

  params$xtx = XTX

  params <- add_to_log(params, "GetProductsCox.t3", read_time, read_size, 0, 0)
  return(params)
}


check_colinearity_cox_t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "check_colinearity_cox_t3\n\n")
  xtx = params$xtx
  nrow = nrow(xtx)
  num_strata = length(params$survival$strata)
  indicies = 1:num_strata
  for (i in (1 + num_strata):nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }

  a_names = params$colnamesA
  b_names = params$colnamesB
  indicies = indicies[-(1:num_strata)] - num_strata# Get rid of the strata indicators
  a_index = which(indicies <= length(a_names))
  b_index = which(indicies > length(a_names))
  params$indicies      = indicies
  params$a_indicies_keep = indicies[a_index]
  params$b_indicies_keep = indicies[b_index] - length(a_names)
  AnamesKeep = a_names[params$a_indicies_keep]
  BnamesKeep = b_names[params$b_indicies_keep]
  params$colnamesA.old = params$colnamesA
  params$colnamesB.old = params$colnamesB
  params$colnamesA     = AnamesKeep
  params$colnamesB     = BnamesKeep
  params$p1.old        = params$p1
  params$p2.old        = params$p2
  params$p1            = length(AnamesKeep)
  params$p2            = length(BnamesKeep)
  params$p.old         = params$p1.old + params$p2.old
  params$p             = params$p1 + params$p2

  Aindicies = params$a_indicies_keep
  Bindicies = params$b_indicies_keep

  colnamesA.old = params$colnamesA.old
  p2 = params$p2

  write_time <- proc.time()[3]
  save(p2, Aindicies, file = file.path(params$write_path, "Aindicies.rdata"))
  save(colnamesA.old, Bindicies, file = file.path(params$write_path, "Bindicies.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("Aindicies.rdata",
                                                          "Bindicies.rdata"))))

  tags = params$Btags[Bindicies]

  if (length(unique(tags)) == 0) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has no covariates."
  }
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "check_colinearity_cox_t3", 0, 0, write_time, write_size)
  return(params)
}


ComputeInitialBetasCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasCox.t3\n\n")
  Abetas   = rep(0, params$p1)
  Bbetas   = rep(0, params$p2)
  betas    = c(Abetas, Bbetas)

  params$betas           = betas
  params$betasold        = betas
  params$Xbeta           = rep(0, params$n)
  params$algIterationCounter = 1
  params$deltabeta       = Inf
  params$loglikelihood   = -Inf
  params$converged       = FALSE
  params$maxIterExceeded = FALSE
  converged              = FALSE
  maxIterExceeded        = FALSE

  write_time <- proc.time()[3]
  save(Abetas, file = file.path(params$write_path, "betasA.rdata"))
  save(Bbetas, file = file.path(params$write_path, "betasB.rdata"))
  save(converged, maxIterExceeded,
       file = file.path(params$write_path, "converged.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("betasA.rdata",
                                                          "betasB.rdata",
                                                          "converged.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "ComputeInitialBetasCox.t3", 0, 0, write_time, write_size)
}


UpdateParamsCox.a3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.a3\n\n")
  Aindicies = NULL
  p2        = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Aindicies.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "Aindicies.rdata"))
  read_time <- proc.time()[3] - read_time
  params$p             = length(Aindicies)
  params$p2            = p2
  params$colnames.old  = params$colnames
  params$colnames      = params$colnames[Aindicies]
  params$a_indicies_keep = Aindicies
  params <- add_to_log(params, "UpdateParamsCox.a3, UpdateDataCox.a3", read_time, read_size, 0, 0)
  return(params)
}


UpdateParamsCox.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.b3\n\n")
  Bindicies     = NULL
  colnamesA.old = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Bindicies.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "Bindicies.rdata"))
  read_time <- proc.time()[3] - read_time
  params$p             = length(Bindicies)
  params$colnamesA.old = colnamesA.old
  params$colnames.old  = params$colnames
  params$colnames      = params$colnames[Bindicies]
  params$b_indicies_keep = Bindicies
  params <- add_to_log(params, "UpdateParamsCox.b3", read_time, read_size, 0, 0)
  return(params)
}


UpdateDataCox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.a3\n\n")
  data$X = as.matrix(data$X[, params$a_indicies_keep, drop = FALSE])
  return(data)
}


UpdateDataCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.b3\n\n")
  data$X = as.matrix(data$X[, params$b_indicies_keep, drop = FALSE])
  return(data)
}


GetBetaACox.a3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBeataACox.a3\n\n")
  converged = NULL
  maxIterExceeded = NULL
  Abetas = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasA.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]], c("converged.rdata",
                                                               "betasA.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  params$betas = Abetas
  params <- add_to_log(params, "GetBetaACox.a3", read_time, read_size, 0, 0)
  return(params)
}


GetBetaBCox.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaBCox.b3\n\n")
  converged       = NULL
  maxIterExceeded = NULL
  Bbetas          = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "converged.rdata"))
  load(file.path(params$read_path[["T"]], "betasB.rdata"))
  read_size <- sum(file.size(file.path(params$read_path[["T"]], c("converged.rdata",
                                                               "betasB.rdata"))))
  read_time <- proc.time()[3] - read_time
  params$converged = converged
  params$maxIterExceeded = maxIterExceeded
  params$betas = Bbetas
  params <- add_to_log(params, "GetBetaBCox.b3", read_time, read_size, 0, 0)
  return(params)
}


GetXABetaACox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXABetaACox.a3\n\n")
  XAbetaA = data$X %*% params$betas
  write_time <- proc.time()[3]
  save(XAbetaA, file = file.path(params$write_path, "xabetaa.rdata"))
  write_size <- file.size(file.path(params$write_path, "xabetaa.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetXABetaACox.a3", 0, 0, write_time, write_size)
  return(params)
}


GetXBBetaBCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBBetaBCox.b3\n\n")
  XBbetaB = data$X %*% params$betas
  write_time <- proc.time()[3]
  save(XBbetaB, file = file.path(params$write_path, "xbbetab.rdata"))
  write_size <- file.size(file.path(params$write_path, "xbbetab.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "GetXBBetaBCox.b3", 0, 0, write_time, write_size)
  return(params)
}


ComputeLogLikelihoodCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeLogLikelihoodCox.t3\n\n")
  n  = params$n
  XAbetaA = NULL
  XBbetaB = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "xabetaa.rdata"))
  load(file.path(params$read_path[["B"]], "xbbetab.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "xabetaa.rdata")) +
    file.size(file.path(params$read_path[["B"]], "xbbetab.rdata"))
  read_time <- proc.time()[3] - read_time
  params$Xbeta.old = params$Xbeta
  Xbeta = XAbetaA + XBbetaB
  params$Xbeta = Xbeta
  params$loglikelihood.old = params$loglikelihood
  stepSize = 1
  w = exp(Xbeta)
  while (max(w) == Inf) {
    Xbeta = (Xbeta + params$Xbetas.old) * 0.5
    stepSize = stepSize * 0.5
    w = exp(Xbeta)
  }
  computeLoglikelihood = TRUE

  while (computeLoglikelihood) {
    numEvents = sum(params$survival$status)
    stepCounter = 0
    pbar <- MakeProgressBar1(numEvents, "Loglikelihood", params$verbose)
    loglikelihood = 0
    for (i in 1:length(params$survival$strata)) {                    ##!
      if (params$survival$strata[[i]]$J > 0) {                       ##!
        for (j in 1:params$survival$strata[[i]]$J) {                 ##!
          nj = params$survival$strata[[i]]$nfails[j]                 ##!
          yIndex = params$survival$strata[[i]]$start0[j]:params$survival$strata[[i]]$end      ##!
          zIndex = params$survival$strata[[i]]$start1[j]:params$survival$strata[[i]]$stop1[j] ##!
          Aj1 = sum(w[yIndex])
          Aj2 = sum(w[zIndex]) / nj
          loglikelihood = loglikelihood + sum(log(w[zIndex]))
          for (r in 0:(nj - 1)) {
            Ajr = Aj1 - r * Aj2
            loglikelihood = loglikelihood - log(Ajr)
          }
          stepCounter = stepCounter + nj
          pbar <- MakeProgressBar2(stepCounter, pbar, params$verbose)
        }
      }
    }
    if (loglikelihood > params$loglikelihood.old || stepSize < 0.5^6) {
      computeLoglikelihood = FALSE
    } else {
      if (params$verbose) cat("Step Halving\n\n")
      Xbeta = (Xbeta + params$Xbeta.old) * 0.5
      stepSize = stepSize * 0.5
      w = exp(Xbeta)
    }
  }
  params$loglikelihoodold = params$loglikelihood
  params$loglikelihood = loglikelihood
  if (params$algIterationCounter == 1) {
    params$nullLoglikelihood = loglikelihood
  }
  params$Xbeta = Xbeta
  params$stepSize = stepSize
  write_time <- proc.time()[3]
  save(Xbeta, file = file.path(params$write_path, "Xbeta.rdata"))
  write_size <- file.size(file.path(params$write_path, "Xbeta.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeLogLikelihoodCox.t3", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeXBDeltaLCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeXBDeltaLCox.b3\n\n")
  Xbeta = NULL
  p2 = params$p
  n = params$n

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Xbeta.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "Xbeta.rdata"))
  read_time <- proc.time()[3] - read_time

  numEvents = sum(data$survival$status)

  w = exp(Xbeta)
  deltal = as.numeric(data$survival$status)
  deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  W.XB = matrix(0, n, p2)

  .Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XB,
        as.integer(n), as.integer(p2), as.integer(numEvents),
        as.integer(params$verbose))

  container_ct_RZ = 0
  container_ct_Cox = 0
  write_size <- 0
  write_time <- 0

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "R*(I-Z*Z')W*XB", params$verbose)
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.RZ) {
      container_ct_RZ = container_ct_RZ + 1
      filename1 <- paste0("crz_", container_ct_RZ, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.Cox) {
      container_ct_Cox = container_ct_Cox + 1
      filename2 <- paste0("cCox_", container_ct_Cox, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n2 <- stp - strt + 1
    g = params$blocks$g[i]

    read_time <- read_time - proc.time()[3]
    RZ = matrix(readBin(con = to_read, what = numeric(), n = n2 * n2,
                        endian = "little"), nrow = n2, ncol = n2)
    read_time <- read_time + proc.time()[3]

    IZ.tZ.W.XBtemp = RZ %*% W.XB[strt:stp, , drop = FALSE]

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(IZ.tZ.W.XBtemp), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.RZ || i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.Cox || i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size <- file.size(file.path(params$write_path, filename2))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }

  tXB.W.XB = t(data$X) %*% W.XB
  tXB.deltal = t(data$X) %*% deltal

  write_time <- write_time - proc.time()[3]
  save(tXB.deltal, tXB.W.XB, file = file.path(params$write_path, "tXB_w_XB.rdata"))
  write_size <- write_size + file.size(file.path(params$write_path, "tXB_w_XB.rdata"))
  write_time <- write_time + proc.time()[3]

  params <- add_to_log(params, "ComputeXBDeltaLCox.b3", read_time, read_size, write_time, write_size)
  return(params)
}


ComputeXADeltaLCox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeXADeltaLCox.a3\n\n")
  Xbeta = NULL
  p1 = params$p
  n = params$n

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "Xbeta.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "Xbeta.rdata"))
  read_time <- proc.time()[3] - read_time

  numEvents = sum(data$survival$status)

  w = exp(Xbeta)
  deltal = as.numeric(data$survival$status)
  deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
  # a pass by reference with the C call.
  W.XA = matrix(0, n, p1)

  .Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XA,
        as.integer(n), as.integer(p1), as.integer(numEvents),
        as.integer(params$verbose))

  tXA.W.XA = t(data$X) %*% W.XA
  tXA.deltal = t(data$X) %*% deltal

  write_time <- proc.time()[3]
  save(tXA.deltal, tXA.W.XA, file = file.path(params$write_path, "tXA_w_XA.rdata"))
  write_size <- file.size(file.path(params$write_path, "tXA_w_XA.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "ComputeXADeltaLCox.a3", read_time, read_size, write_time, write_size)
  return(params)
}


ProcessVCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessVCox.t3\n\n")
  write_time <- 0
  write_size <- 0
  read_time <- 0
  read_size <- 0
  p2 = params$p2

  num_blocks = params$blocks$num_blocks
  pbar <- MakeProgressBar1(num_blocks, "(I-Z*Z')W*XB*R", params$verbose)

  container_ct_RV = 0
  container_ct_VR = 0
  for (i in 1:num_blocks) {
    if (i %in% params$container$filebreak.RV) {
      container_ct_RV = container_ct_RV + 1
      filename2 <- paste0("cCox_", container_ct_RV, ".rdata")
      toRead2 = file(file.path(params$read_path[["B"]], filename2), "rb")
    }
    if (i %in% params$container$filebreak.VR) {
      container_ct_VR = container_ct_VR + 1
      filename3 = paste0("cvr_", container_ct_VR, ".rdata")
      to_write3 = file(file.path(params$write_path, filename3), "wb")
    }

    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1

    filename1 <- paste0("r1_", i, ".rdata")
    filename4 = paste0("r3_", i, ".rdata")

    read_time <- read_time - proc.time()[3]
    toRead1 = file(file.path(params$dp_local_path, filename1), "rb")
    R2 = matrix(readBin(con = toRead1, what = numeric(), n = n * n,
                        endian = "little"), nrow = n, ncol = n)
    read_size <- read_size + file.size(file.path(params$dp_local_path, filename1))
    close(toRead1)
    RV = matrix(readBin(con = toRead2, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]

    V = t(R2) %*% RV
    R3 = RandomOrthonomalMatrix(p2)
    VR = V %*% R3

    write_time <- write_time - proc.time()[3]
    to_write4 = file(file.path(params$dp_local_path, filename4), "wb")
    writeBin(as.vector(R3), con = to_write4, endian = "little")
    close(to_write4)
    write_size <- write_size + file.size(file.path(params$dp_local_path, filename4))
    writeBin(as.vector(VR), con = to_write3, endian = "little")
    write_time <- write_time + proc.time()[3]
    if ((i + 1) %in% params$container$filebreak.RV || i == num_blocks) {
      close(toRead2)
      read_size <- read_size + file.size(file.path(params$read_path[["B"]], filename2))
    }
    if ((i + 1) %in% params$container$filebreak.VR || i == num_blocks) {
      close(to_write3)
      write_size <- write_size + file.size(file.path(params$write_path, filename3))
    }

    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }
  params <- add_to_log(params, "ProcessVCox.t3", read_time, read_size, write_time, write_size)
  return(params)
}


GetXRCox.a3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXRCox.a3\n\n")
  write_time <- 0
  write_size <- 0
  read_time <- 0
  read_size <- 0
  p2 = params$p2
  container_ct_VR = 0
  container_ct_XR = 0
  pbar <- MakeProgressBar1(params$blocks$num_blocks, "XA'(I-Z*Z')W*XB*R", params$verbose)
  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.RV) {
      container_ct_VR = container_ct_VR + 1
      filename1 <- paste0("cvr_", container_ct_VR, ".rdata")
      to_read <- file(file.path(params$read_path[["T"]], filename1), "rb")
    }
    if (i %in% params$container$filebreak.XR) {
      container_ct_XR = container_ct_XR + 1
      filename2 <- paste0("cxr_", container_ct_XR, ".rdata")
      to_write <- file(file.path(params$write_path, filename2), "wb")
    }
    strt <- params$blocks$starts[i]
    stp <- params$blocks$stops[i]
    n = stp - strt + 1

    read_time <- read_time - proc.time()[3]
    VR = matrix(readBin(con = to_read, what = numeric(), n = n * p2,
                        endian = "little"), nrow = n, ncol = p2)
    read_time <- read_time + proc.time()[3]
    XR = t(data$X[strt:stp, , drop = FALSE]) %*% VR

    write_time <- write_time - proc.time()[3]
    writeBin(as.vector(XR), con = to_write, endian = "little")
    write_time <- write_time + proc.time()[3]

    if ((i + 1) %in% params$container$filebreak.VR || i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path[["T"]], filename1))
    }
    if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$num_blocks) {
      close(to_write)
      write_size <- write_size + file.size(file.path(params$write_path, filename2))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }

  params <- add_to_log(params, "GetXRCox.a3", read_time, read_size, write_time, write_size)
  return(params)
}


ProcessXtWXCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessXtWXCox.t3\n\n")
  p1 = params$p1
  p2 = params$p2

  tXA.deltal = NULL
  tXB.deltal = NULL
  tXA.W.XA   = NULL
  tXB.W.XB   = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["A"]], "tXA_w_XA.rdata"))
  load(file.path(params$read_path[["B"]], "tXB_w_XB.rdata"))
  read_size <- file.size(file.path(params$read_path[["A"]], "tXA_w_XA.rdata")) +
    file.size(file.path(params$read_path[["B"]], "tXB_w_XB.rdata"))
  read_time <- proc.time()[3] - read_time

  pbar <- MakeProgressBar1(params$blocks$num_blocks, "X'W*X", params$verbose)
  container_ct_XR = 0
  XATWXB = 0

  for (i in 1:params$blocks$num_blocks) {
    if (i %in% params$container$filebreak.XR) {
      container_ct_XR = container_ct_XR + 1
      filename1 <- paste0("cxr_", container_ct_XR, ".rdata")
      to_read <- file(file.path(params$read_path[["A"]], filename1), "rb")
    }

    filename2 <- paste0("r3_", i, ".rdata")
    read_time <- read_time - proc.time()[3]
    toRead1 = file(file.path(params$dp_local_path, filename2), "rb")
    R = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
                       endian = "little"), nrow = p2, ncol = p2)
    close(toRead1)
    XR = matrix(readBin(con = to_read, what = numeric(), n = p1 * p2,
                        endian = "little"), nrow = p1, ncol = p2)

    read_size <- read_size + file.size(file.path(params$dp_local_path, filename2))
    read_time <- read_time + proc.time()[3]

    XATWXB = XATWXB + XR %*% t(R)

    if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$num_blocks) {
      close(to_read)
      read_size <- read_size + file.size(file.path(params$read_path[["A"]], filename1))
    }
    pbar <- MakeProgressBar2(i, pbar, params$verbose)
  }


  xtwx = rbind(cbind(tXA.W.XA, XATWXB), cbind(t(XATWXB), tXB.W.XB))
  II = NULL
  tryCatch({
    II = solve(xtwx)
  },
  error = function(err) {
    II = NULL
  }
  )
  if (is.null(II)) {
    params$failed <- TRUE
    params$singular_matrix = TRUE
    params$error_message <-
      paste0("The matrix t(X)*W*X is not invertible.\n",
             "       This may be due to one of two possible problems.\n",
             "       1. Poor random initialization of the security vector.\n",
             "       2. Near multicollinearity in the data\n",
             "SOLUTIONS: \n",
             "       1. Rerun the data analysis.\n",
             "       2. If the problem persists, check the variables for\n",
             "          duplicates for both parties and / or reduce the\n",
             "          number of variables used. Once this is done,\n",
             "          rerun the data analysis.")
    params <- add_to_log(params, "ProcessXtWXCox.t3", read_time, read_size, 0, 0)
    return(params)
  }

  if (params$algIterationCounter == 1) {
    params$nullHessian = xtwx
    params$nullScore = rbind(tXA.deltal, tXB.deltal)
  }
  params$xtwx = xtwx
  deltaBeta = II %*% rbind(tXA.deltal, tXB.deltal)
  params$betas    = params$betasold + (params$betas - params$betasold) * params$stepSize
  params$betasold = params$betas
  params$betas    = params$betasold + deltaBeta
  converged       = abs(params$loglikelihood - params$loglikelihoodold) /
    (abs(params$loglikelihood) + 0.1) < params$cutoff
  maxIterExceeded = FALSE
  if (!converged) {
    maxIterExceeded = params$algIterationCounter >= params$max_iterations
  }

  params$converged = converged
  params$maxIterExceeded = maxIterExceeded

  Abetas = params$betas[1:p1]
  Bbetas = params$betas[(p1 + 1):(p1 + p2)]
  write_time <- proc.time()[3]
  save(Abetas, file = file.path(params$write_path, "betasA.rdata"))
  save(Bbetas, file = file.path(params$write_path, "betasB.rdata"))
  save(converged, maxIterExceeded,
       file = file.path(params$write_path, "converged.rdata"))
  write_size <- sum(file.size(file.path(params$write_path, c("betasA.rdata",
                                                          "betasB.rdata",
                                                          "converged.rdata"))))
  write_time <- proc.time()[3] - write_time

  params <- add_to_log(params, "ProcessXtWXCox.t3", read_time, read_size, write_time, write_size)

  return(params)
}


SurvFitCox.Bt3 = function(params, pred) {
  if (params$trace) cat(as.character(Sys.time()), "SurvFitCox.Bt3\n\n")
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
      rcumsum = function(x) rev(cumsum(rev(x)))
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
ComputeResultsCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComptueResultsCvox.t3\n\n")
  stats = params$stats
  stats$failed         = FALSE
  stats$converged      = params$converged
  names.new          = c(params$colnamesA, params$colnamesB)
  names.old          = c(params$colnamesA.old, params$colnamesB.old)
  idx                = params$indicies
  stats$party        = c(rep("dp1", params$p1.old), rep("dp2", params$p2.old))
  stats$coefficients = rep(NA, length(stats$party))
  stats$coefficients[idx] = params$betas
  stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
  stats$expncoef     = exp(-stats$coefficients)
  tempvar            = solve(params$xtwx)
  stats$var          = matrix(0, length(names.old), length(names.old))
  stats$var[idx, idx] = tempvar
  stats$secoef       = rep(NA, length(names.old))
  stats$secoef[idx]  = sqrt(diag(tempvar))  # se(coef)

  stats$zvals        = stats$coefficients / stats$secoef  # z values
  stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE)   # pvals
  stats$stars        = matrix(sapply(stats$pvals, function(x) {
    if (is.na(x))       ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  }))
  stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  stats$loglik       = c(params$nullLoglikelihood, params$loglikelihood)
  stats$n            = params$n
  stats$nevent       = sum(params$survival$status)
  stats$df           = params$p
  stats$iter         = params$algIterationCounter - 1
  stats$score        = t(params$nullScore) %*% solve(params$nullHessian) %*% params$nullScore
  stats$score        = c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       = "efron"
  stats$lrt          = 2 * (stats$loglik[2] - stats$loglik[1])
  stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      = c(1 - exp(-stats$lrt[1] / stats$n),
                         1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    = t(params$betas) %*% params$xtwx %*% params$betas
  stats$wald.test    = c(stats$wald.test,
                         1 - pchisq(stats$wald.test, stats$df))
  pred = -params$Xbeta
  if (params$survival_installed) {
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

  stats$survival = data.frame(
    rank   = params$survival$rank,
    status = params$survival$status,
    sorted = params$survival$sorted,
    surv   = SurvFitCox.Bt3(params, pred)
  )
  stats$strata = as.data.frame(matrix(0, length(params$survival$strata), 3))
  stats$strata$label = ""
  colnames(stats$strata) = c("start", "end", "events", "label")
  for (i in 1:length(params$survival$strata)) {
    stats$strata$start[i]  = params$survival$strata[[i]]$start
    stats$strata$end[i]    = params$survival$strata[[i]]$end
    stats$strata$events[i] = sum(params$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  = params$survival$strata[[i]]$label
  }

  names(stats$party)           = names.old
  names(stats$coefficients)    = names.old
  names(stats$expcoef)         = names.old
  names(stats$expncoef)        = names.old
  rownames(stats$var)          = names.old
  colnames(stats$var)          = names.old
  names(stats$secoef)          = names.old
  names(stats$zvals)           = names.old
  names(stats$pvals)           = names.old
  names(stats$stars)           = names.old
  names(stats$lower95)         = names.old
  names(stats$upper95)         = names.old
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

  params <- add_to_log(params, "ComputeResultsCox.t3", 0, 0, write_time, write_size)
  return(params)
}


GetResultsCox.a3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.a3\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsCox.a3", read_time, read_size, 0, 0)
  return(params)
}


GetResultsCox.b3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.b3\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats
  params <- add_to_log(params, "GetResultsCox.b3", read_time, read_size, 0, 0)
  return(params)
}


####################### REGRESSION BY B ONLY FUNCTIONS #######################

check_colinearity_cox_b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "check_colinearity_cox_b3\n\n")
  # Add in the strata here
  num_strata = length(data$survival$strata)
  X = cbind(matrix(0, nrow = nrow(data$X), ncol = num_strata), data$X)
  for (i in 1:num_strata) {
    X[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
  }

  XTX = t(X) %*% X
  nrow = nrow(XTX)
  indicies = 1:num_strata
  for (i in (num_strata + 1):nrow) {
    tempIndicies = c(indicies, i)
    if (rcond(XTX[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
      indicies = c(indicies, i)
    }
  }
  indicies = indicies[-(1:num_strata)] - num_strata  # Get rid of the Strata

  params$a_indicies_keep = c()
  params$colnamesA.old = c()
  params$colnamesA     = c()
  params$p1.old        = params$p1
  params$p1            = 0

  b_names               = params$colnames
  BnamesKeep           = b_names[indicies]
  params$b_indicies_keep = indicies
  params$colnames.old  = params$colnames
  params$colnames      = BnamesKeep
  params$p2.old        = params$p2
  params$p2            = length(BnamesKeep)
  params$p             = params$p1 + params$p2

  if (params$p2 == 0) {
    params$failed <- TRUE
    params$error_message <- "After removing colinear covariates, Party B has no covariates."
  }
  params <- add_to_log(params, "check_colinearity_cox_b3", 0, 0, 0, 0)

  return(params)
}

#' @importFrom stats as.formula
ComputeCoxFromSurvival.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCoxFromSurvival.b3\n\n")
  # We have loaded survival previously

  max_iterations = 25
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "max_iterations.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "max_iterations.rdata"))
  read_time <- proc.time()[3] - read_time
  strata = rep(0, nrow(data$X))
  for (i in 1:length(data$survival$strata)) {
    strata[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] = i
  }

  colnames(data$X) = paste0("V", 1:ncol(data$X))
  f = paste(c("Surv(rank, status) ~ strata(strata)", paste0("V", 1:ncol(data$X))), collapse = " + ")

  error = tryCatch(
    {fit = survival::coxph(as.formula(f),
                           data = data.frame(rank = data$survival$rank,
                                             status = data$survival$status,
                                             strata = strata,
                                             data$X),
                           iter.max = max_iterations)},
    error = function(e) {
      return(TRUE)
    },
    warning = function(e) {
      return(FALSE)
    }
  )

  if (class(error) == "logical" && error) {
    params$converged = FALSE
    params$failed <- TRUE
    params$error_message <- "Coxph in the survival package failed to converge."
  } else {
    params$converged = TRUE
    if (class(error) == "logical") {
      fit = suppressWarnings(survival::coxph(as.formula(f),
                                             data = data.frame(rank = data$survival$rank,
                                                               status = data$survival$status,
                                                               strata = strata,
                                                               data$X),
                                             iter.max = max_iterations))
      params$converged = FALSE
    }
    fit$linear.predictors = NULL
    fit$residuals = NULL
    fit$y = NULL
    fit$formula = NULL
    fit$call = NULL
    fit$assign = NULL
    fit$terms = NULL
    fit$means = NULL
    params$fit = fit
  }
  params <- add_to_log(params, "ComputeCoxFromSurvival.b3", read_time, read_size, 0, 0)
  return(params)
}


ComputeCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCox.b3\n\n")
  n           = params$n
  p2          = params$p
  params$algIterationCounter = 1
  X.betas.old = matrix(0, n, 1)
  X.betas     = matrix(0, n, 1)
  betasB      = matrix(0, p2, 1)
  betasBold   = betasB
  loglikelihood.old = -Inf
  max_iterations = 25
  cutoff        = 10^-8

  read_time <- proc.time()[3]
  load(file.path(params$read_path[["T"]], "max_iterations.rdata"))
  read_size <- file.size(file.path(params$read_path[["T"]], "max_iterations.rdata"))
  read_time <- proc.time()[3] - read_time


  while (params$algIterationCounter <= max_iterations && !params$converged) {
    BeginningIteration(params)
    loglikelihood = 0
    stepSize = 1
    w = exp(X.betas)
    while (max(w) == Inf) {
      if (params$verbose) cat("Step Halving\n\n")
      X.betas = (X.betas + X.betas.old) * 0.5
      stepSize = stepSize * 0.5
      w = exp(X.betas)
    }
    computeLoglikelihood = TRUE
    while (computeLoglikelihood) {
      numEvents = sum(data$survival$status)
      stepCounter = 0
      pbar <- MakeProgressBar1(numEvents, "Loglikelihood", params$verbose)
      loglikelihood = 0
      for (i in 1:length(data$survival$strata)) {                    ##!
        if (data$survival$strata[[i]]$J > 0) {                       ##!
          for (j in 1:data$survival$strata[[i]]$J) {                 ##!
            nj = data$survival$strata[[i]]$nfails[j]                 ##!
            yIndex = data$survival$strata[[i]]$start0[j]:data$survival$strata[[i]]$end      ##!
            zIndex = data$survival$strata[[i]]$start1[j]:data$survival$strata[[i]]$stop1[j] ##!
            Aj1 = sum(w[yIndex])
            Aj2 = sum(w[zIndex]) / nj
            loglikelihood = loglikelihood + sum(log(w[zIndex]))
            for (r in 0:(nj - 1)) {
              Ajr = Aj1 - r * Aj2
              loglikelihood = loglikelihood - log(Ajr)
            }
            stepCounter = stepCounter + nj
            pbar <- MakeProgressBar2(stepCounter, pbar, params$verbose)
          }
        }
      }
      if (loglikelihood > loglikelihood.old || stepSize < 0.5^6) {
        computeLoglikelihood = FALSE
      } else {
        if (params$verbose) cat("Step Halving\n\n")
        X.betas = (X.betas + X.betas.old) * 0.5
        stepSize = stepSize * 0.5
        w = exp(X.betas)
      }
    }
    numEvents = sum(data$survival$status)
    deltal = as.numeric(data$survival$status)
    deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
    # a pass by reference with the C call.
    W.XB = matrix(0, n, p2)

    .Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XB,
          as.integer(n), as.integer(p2), as.integer(numEvents),
          as.integer(params$verbose))

    m = t(data$X) %*% W.XB

    params$XtWX = m
    if (params$algIterationCounter == 1) {
      params$nullHessian = m
    }

    inv = NULL
    tryCatch({
      inv = solve(m)
    },
    error = function(err) {
      inv = NULL
    })
    m = inv
    if (is.null(m)) {
      params$failed <- TRUE
      params$error_message <- "The matrix t(X)WX is singular.  This is probably due to divergence of the coefficients."

      betas = rep(NA, length(params$Bcolnames.old))
      betas[params$b_indicies_keep] = betasB
      betas = data.frame(betas)
      rownames(betas) = params$Bcolnames.old
      # if (params$verbose) cat("Current Parameters:\n")
      # if (params$verbose) print(betas)
      # if (params$verbose) cat("\n")
      params <- add_to_log(params, "ComputeCox.b3", read_time, read_size, 0, 0)
      return(params)
    }

    deltaBeta = m %*% t(data$X) %*% deltal
    betasB    = betasBold + (betasB - betasBold) * stepSize
    betasBold = betasB
    betasB    = betasB + deltaBeta
    X.betas   = data$X %*% betasB

    converged = abs(loglikelihood - loglikelihood.old) /
      (abs(loglikelihood) + 0.1) < cutoff
    params$converged = converged

    if (params$algIterationCounter == 1) {
      params$nullScore         = t(data$X) %*% deltal
      params$nullLoglikelihood = loglikelihood
    }
    loglikelihood.old = loglikelihood
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }
  params$loglikelihood = loglikelihood
  params$betasB = betasB
  params$X.betas = X.betas
  params <- add_to_log(params, "ComputeCox.b3", read_time, read_size, 0, 0)
  return(params)
}


#' @importFrom  stats pchisq pnorm qnorm
ComputeResultsCox.b3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsCox.b3\n\n")
  stats = params$stats
  stats$converged = params$converged
  stats$party_name <- params$party_name
  stats$failed    = FALSE

  fitExists = !is.null(params$fit)
  names.old          = c(params$colnamesA.old, params$colnames.old)
  idxA               = params$a_indicies_keep
  idxB               = params$b_indicies_keep
  idx                = c(idxA, idxB + length(params$colnamesA.old))
  stats$party        = c(rep("dp1", length(params$colnamesA.old)),
                         rep("dp2", length(params$colnames.old)))
  stats$coefficients = rep(NA, length(stats$party))
  if (fitExists) {
    stats$coefficients[idx] = params$fit$coefficients
    tempvar          = params$fit$var
  } else {
    stats$coefficients[idx] = params$betasB
    tempvar            = solve(params$XtWX)
  }
  stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
  stats$expncoef     = exp(-stats$coefficients)
  stats$var          = matrix(0, length(names.old), length(names.old))
  stats$var[idx, idx] = tempvar
  stats$secoef       = rep(NA, length(names.old))
  stats$secoef[idx]  = sqrt(diag(tempvar))  # se(coef)
  stats$zvals        = stats$coefficients / stats$secoef  # z values
  stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE)   # pvals
  stats$stars        = matrix(sapply(stats$pvals, function(x) {
    if (is.na(x))       ""
    else if (x < 0.001) "***"
    else if (x < 0.01)  "**"
    else if (x < 0.05)  "*"
    else if (x < 0.1)   "."
    else " "
  }))
  stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
  stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
  if (fitExists) {
    stats$loglik     = params$fit$loglik
    stats$n          = params$fit$n
    stats$nevent     = params$fit$nevent
    stats$iter       = params$fit$iter
    stats$score      = params$fit$score
    stats$wald.test  = params$fit$wald.test
    stats$concordance = params$fit$concordance[c(2, 1, 3, 4, 6, 7)]
  } else {
    stats$loglik       = c(params$nullLoglikelihood, params$loglikelihood)
    stats$n            = params$n
    stats$nevent       = params$numEvents
    stats$iter         = params$algIterationCounter - 1
    stats$score        = t(params$nullScore) %*% solve(params$nullHessian) %*%
      params$nullScore
    stats$wald.test    = t(params$betasB) %*% params$XtWX %*% params$betasB
    stats$concordance = c(NA, NA, NA, NA, NA, NA)
  }
  stats$df           = params$p
  stats$score        = c(stats$score, 1 - pchisq(stats$score, stats$df))
  stats$method       = "efron"
  stats$lrt          = 2 * (stats$loglik[2] - stats$loglik[1])
  stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
  stats$rsquare      = c(1 - exp(-stats$lrt[1] / stats$n),
                         1 - exp(2 * stats$loglik[1] / stats$n))
  stats$wald.test    = c(stats$wald.test,
                         1 - pchisq(stats$wald.test, stats$df))

  params$survival = data$survival
  pred = data$X %*% stats$coefficients[idx]
  stats$survival = data.frame(
    rank   = data$survival$rank,
    status = data$survival$status,
    sorted = data$survival$sorted,
    surv   = SurvFitCox.Bt3(params, pred)
  )
  stats$strata = as.data.frame(matrix(0, length(data$survival$strata), 3))
  stats$strata$label = ""
  colnames(stats$strata) = c("start", "end", "events", "label")
  for (i in 1:length(data$survival$strata)) {
    stats$strata$start[i]  = data$survival$strata[[i]]$start
    stats$strata$end[i]    = data$survival$strata[[i]]$end
    stats$strata$events[i] = sum(data$survival$status[stats$strata$start[i]:stats$strata$end[i]])
    stats$strata$label[i]  = data$survival$strata[[i]]$label
  }

  names(stats$party)           = names.old
  names(stats$coefficients)    = names.old
  names(stats$expcoef)         = names.old
  names(stats$expncoef)        = names.old
  rownames(stats$var)          = names.old
  colnames(stats$var)          = names.old
  names(stats$secoef)          = names.old
  names(stats$zvals)           = names.old
  names(stats$pvals)           = names.old
  names(stats$stars)           = names.old
  names(stats$lower95)         = names.old
  names(stats$upper95)         = names.old
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
  params <- add_to_log(params, "ComputeResultsCox.b3", 0, 0, write_time, write_size)
  return(params)
}


TransferResultsCox.t3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "TransferResultsCox.t3\n\n")
  stats = NULL
  read_time <- proc.time()[3]
  load(file.path(params$read_path[["B"]], "stats.rdata"))
  read_size <- file.size(file.path(params$read_path[["B"]], "stats.rdata"))
  read_time <- proc.time()[3] - read_time
  params$stats = stats

  write_time <- proc.time()[3]
  save(stats, file = file.path(params$write_path, "stats.rdata"))
  write_size <- file.size(file.path(params$write_path, "stats.rdata"))
  write_time <- proc.time()[3] - write_time
  params <- add_to_log(params, "TransferResultsCox.t3", read_time, read_size, write_time, write_size)
  return(params)
}


############################## PARENT FUNCTIONS ##############################


PartyAProcess3Cox = function(data,
                             y_name                 = NULL,
                             strata                = NULL,
                             mask                  = TRUE,
                             monitor_folder         = NULL,
                             sleep_time             = 10,
                             max_waiting_time        = 24 * 60 * 60,
                             popmednet             = TRUE,
                             trace                 = FALSE,
                             verbose               = TRUE) {
  params <- prepare_params_3p("cox", "A",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)
  Header(params)
  params   = prepare_folder_linear_a3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = prepare_data_cox_23(params, data, y_name, strata, mask)
  params <- add_to_log(params, "prepare_data_cox_23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party A."
    make_error_message(params$write_path, message)
    files = c("error_message.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time,
                              job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_cox_a3(params, data)
  files = "pa.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }

  params <- send_strata_cox_a3(params, data)
  files = "Astrata.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "stats.rdata"))) {
    params$algIterationCounter = 1
    params <- GetResultsCox.a3(params)
    params$converged = params$stats$converged
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    return(params$stats)
  }

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }

  data = sort_data_cox_a3(params, data)
  params <- add_to_log(params, "sort_data_cox_a3", data$read_time, data$read_size, 0, 0)
  params <- prepare_blocks_linear_a3(params)
  params <- get_z_cox_a3(params, data)
  files = seq_zw("cz_", length(params$container$file_break_z))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- GetWRCox.a3(params, data)
  params <- GetSXACox.a3(params, data)
  files = c("sxa.rdata", "xatxa.rdata", seq_zw("cpr_", length(params$container$filebreak.PR)))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "stats.rdata"))) {
    params$algIterationCounter = 1
    params <- GetResultsCox.a3(params)
    params$converged = params$stats$converged
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    return(params$stats)
  }

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }


  params <- UpdateParamsCox.a3(params)
  data = UpdateDataCox.a3(params, data)

  params$algIterationCounter = 1
  repeat {
    params <- GetBetaACox.a3(params)
    if (params$converged || params$maxIterExceeded) break
    BeginningIteration(params)
    params <- GetXABetaACox.a3(params, data)

    files = c("xabetaa.rdata")
    params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                  waitForTurn = TRUE)

    params <- ComputeXADeltaLCox.a3(params, data)
    files = "tXA_w_XA.rdata"
    params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                  waitForTurn = TRUE)

    params <- GetXRCox.a3(params, data)
    files = seq_zw("cxr_", length(params$container$filebreak.XR))
    params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)

    if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["T"]]))
      params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                                waitForTurn = TRUE)
      return(params$stats)
    }
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- GetResultsCox.a3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}

PartyBProcess3Cox = function(data,
                             strata              = NULL,
                             mask                = TRUE,
                             monitor_folder       = NULL,
                             sleep_time           = 10,
                             max_waiting_time      = 24 * 60 * 60,
                             popmednet           = TRUE,
                             trace               = FALSE,
                             verbose             = TRUE) {
  params <- prepare_params_3p("cox", "B",
                            popmednet = popmednet, trace = trace, verbose = verbose)
  params <- initialize_log_3p(params)
  params <- initialize_time_stamps_3p(params)
  params <- initialize_tracking_table_3p(params)
  Header(params)
  params <- prepare_folder_linear_b3(params, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }
  data = prepare_data_cox_23(params, data, NULL, strata, mask)
  params <- add_to_log(params, "prepare_data_cox_23", 0, 0, 0, 0)

  if (data$failed) {
    message = "Error in processing the data for Party B."
    make_error_message(params$write_path, message)
    files = c("error_message.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time,
                              job_failed = TRUE, waitForTurn = TRUE)
    return(params$stats)
  }

  params <- prepare_params_cox_b3(params, data)
  files = "pb.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }

  params <- send_strata_cox_b3(params, data)

  files = "Bstrata.rdata"
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                waitForTurn = TRUE)

  if (file.exists(file.path(params$read_path[["T"]], "transferControl.rdata"))) {
    params$algIterationCounter = 1
    data = sort_data_cox_b3(params, data)
    params <- add_to_log(params, "sort_data_cox_b3", data$read_time, data$read_size, 0, 0)
    params <- check_colinearity_cox_b3(params, data)

    if (params$failed) {  # Happens if pB.new == 0
      params$complete = TRUE
      warning(params$error_message)
      make_error_message(params$write_path, params$error_message)
      files = c("error_message.rdata")
      params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time,
                                job_failed = TRUE)
      return(params$stats)
    }
    data = UpdateDataCox.b3(params, data)
    params <- add_to_log(params, "UpdateDataCox.b3", 0, 0, 0, 0)
    if (params$survival_installed) {
      params <- ComputeCoxFromSurvival.b3(params, data)
    } else {
      params <- ComputeCox.b3(params, data)
    }

    if (params$failed) {      # We could get a job_failed here from coefficient explosion
      params$complete = TRUE
      warning(params$error_message)
      make_error_message(params$write_path, params$error_message)
      files = c("error_message.rdata")
      params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time,
                                job_failed = TRUE)
      return(params$stats)
    }
    params <- ComputeResultsCox.b3(params, data)
    stats = params$stats
    save(stats, file = file.path(params$write_path, "stats.rdata"))
    files = c("stats.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time)
    return(params$stats)
  }


  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }


  data = sort_data_cox_b3(params, data)
  params <- add_to_log(params, "sort_data_cox_b3", data$read_time, data$read_size, 0, 0)
  params <- prepare_blocks_linear_b3(params)
  params <- GetRWLinear.b3(params, data)
  params <- GetSXBCox.b3(params, data)
  files = c("sxb.rdata", "xbtxb.rdata", seq_zw("crw_", length(params$container$filebreak.RW)))
  params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
    params$complete = TRUE
    warning(read_error_message(params$read_path[["T"]]))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                              waitForTurn = TRUE)
    return(params$stats)
  }

  if (file.exists(file.path(params$read_path[["T"]], "transferControl.rdata"))) {
    params$algIterationCounter = 1
    params <- UpdateParamsCox.b3(params)
    data = UpdateDataCox.b3(params, data)
    params <- add_to_log(params, "UpdateDataCox.b3", 0, 0, 0, 0)
    if (params$survival_installed) {
      params <- ComputeCoxFromSurvival.b3(params, data)
    } else {
      params <- ComputeCox.b3(params, data)
    }

    if (params$failed) {      # We could get a job_failed here from coefficient explosion
      params$complete = TRUE
      warning(params$error_message)
      make_error_message(params$write_path, params$error_message)
      files = c("error_message.rdata")
      params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time,
                                job_failed = TRUE)
      return(params$stats)
    }
    params <- ComputeResultsCox.b3(params, data)
    stats = params$stats
    save(stats, file = file.path(params$write_path, "stats.rdata"))
    files = c("stats.rdata")
    params <- send_pause_quit_3p(params, filesT = files, sleep_time = sleep_time)
    return(params$stats)
  }

  params <- UpdateParamsCox.b3(params)
  data = UpdateDataCox.b3(params, data)
  params$algIterationCounter = 1
  repeat {
    params <- GetBetaBCox.b3(params)
    if (params$converged || params$maxIterExceeded) break
    BeginningIteration(params)
    params <- GetXBBetaBCox.b3(params, data)

    files = "xbbetab.rdata"
    params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                  waitForTurn = TRUE)

    params <- ComputeXBDeltaLCox.b3(params, data)
    files = c("tXB_w_XB.rdata", seq_zw("cCox_", length(params$container$filebreak.Cox)))
    params <- send_pause_continue_3p(params, filesT = files, from = "T",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time,
                                  waitForTurn = TRUE)

    if (file.exists(file.path(params$read_path[["T"]], "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["T"]]))
      params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE,
                                waitForTurn = TRUE)
      return(params$stats)
    }
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- GetResultsCox.b3(params)
  params <- send_pause_quit_3p(params, sleep_time = sleep_time, waitForTurn = TRUE)
  return(params$stats)
}


PartyTProcess3Cox = function(monitor_folder         = NULL,
                             msreqid               = "v_default_0_000",
                             blocksize             = 500,
                             cutoff                = 1e-8,
                             max_iterations         = 25,
                             sleep_time             = 10,
                             max_waiting_time        = 24 * 60 * 60,
                             popmednet             = TRUE,
                             trace                 = FALSE,
                             verbose               = TRUE) {
  Tparams <- prepare_params_3p("cox", "T", msreqid = msreqid,
                             popmednet = popmednet, trace = trace, verbose = verbose)
  Tparams <- initialize_log_3p(Tparams)
  Tparams <- initialize_time_stamps_3p(Tparams)
  Tparams <- initialize_tracking_table_3p(Tparams)

  Header(Tparams)
  params   = prepare_folder_linear_t3(Tparams, monitor_folder)
  if (params$failed) {
    warning(params$error_message)
    return(invisible(NULL))
  }

  params <- PauseContinue.3p(params, from = c("A", "B"), max_waiting_time = max_waiting_time)

  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata")) &&
      file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
    warning(paste0(read_error_message(params$read_path[["A"]]), "\n",
                   read_error_message(params$read_path[["B"]])))
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }
  if (file.exists(file.path(params$read_path[["A"]], "error_message.rdata"))) {
    warning(read_error_message(params$read_path[["A"]]))
    file.copy(file.path(params$read_path[["A"]], "error_message.rdata"),
              file.path(params$write_path, "error_message.rdata"))
    files = "error_message.rdata"
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
    files = "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params   = prepare_params_cox_t3(params, cutoff, max_iterations)

  if (!params$failed) params <- check_strata_cox_t3(params)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files = "empty.rdata"
  params <- send_pause_continue_3p(params, filesA = files, filesB = files, from = c("A", "B"),
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- prepare_strata_cox_t3(params)

  if (params$p1 == 0) {
    params$algIterationCounter = 1
    MakeTransferMessage(params$write_path)
    files = c("transfercontrol.rdata", "max_iterations.rdata", "survival.rdata")
    params <- send_pause_continue_3p(params, filesB = files, from = "B",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["B"]]))
      file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
                file.path(params$write_path, "error_message.rdata"))
      files = "error_message.rdata"
      params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
      params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.3p(params)
      return(params$stats)
    }
    params <- TransferResultsCox.t3(params)
    params$converged = params$stats$converged
    files = "stats.rdata"
    params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params <- prepare_blocks_linear_t3(params, blocksize)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  files = c("survival.rdata", "blocks.rdata")
  params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- ProcessZLinear.t3(params)
  files = c("survival.rdata", "blocks.rdata", seq_zw("crz_", length(params$container$filebreak.RZ)))
  params <- send_pause_continue_3p(params, filesB = files, from = "B",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- ProcessWLinear.t3(params)
  files = c("p2.rdata", seq_zw("cwr_", length(params$container$filebreak.WR)))
  params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- GetProductsCox.t3(params)
  params <- check_colinearity_cox_t3(params)

  if (params$failed) {
    warning(params$error_message)
    make_error_message(params$write_path, params$error_message)
    files = "error_message.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                  from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  if (params$p1 == 0) {
    params$algIterationCounter = 1
    MakeTransferMessage(params$write_path)
    files = c("transfercontrol.rdata", "Bindicies.rdata", "max_iterations.rdata", "survival.rdata")
    params <- send_pause_continue_3p(params, filesB = files, from = "B",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    if (file.exists(file.path(params$read_path[["B"]], "error_message.rdata"))) {
      warning(read_error_message(params$read_path[["B"]]))
      file.copy(file.path(params$read_path[["B"]], "error_message.rdata"),
                file.path(params$write_path, "error_message.rdata"))
      files = "error_message.rdata"
      params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
      params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.3p(params)
      return(params$stats)
    }
    params <- TransferResultsCox.t3(params)
    params$converged = params$stats$converged
    files = "stats.rdata"
    params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    params <- send_pause_quit_3p(params, sleep_time = sleep_time)
    SummarizeLog.3p(params)
    return(params$stats)
  }

  params <- ComputeInitialBetasCox.t3(params)

  filesA = c("Aindicies.rdata", "betasA.rdata", "converged.rdata")
  filesB = c("Bindicies.rdata", "betasB.rdata", "converged.rdata")
  params <- send_pause_continue_3p(params, filesA = filesA, filesB = filesB,
                                from = c("A", "B"),
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  while (!params$converged && !params$maxIterExceeded) {
    BeginningIteration(params)
    if (params$algIterationCounter > 1) {
      filesA = c("converged.rdata", "betasA.rdata")
      filesB = c("converged.rdata", "betasB.rdata")
      params <- send_pause_continue_3p(params, filesA = filesA, filesB = filesB,
                                    from = c("A", "B"),
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
    }

    params <- ComputeLogLikelihoodCox.t3(params)
    files = "Xbeta.rdata"
    params <- send_pause_continue_3p(params, filesA = files, filesB = files, from = c("A", "B"),
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)

    params <- ProcessVCox.t3(params)

    files = seq_zw("cvr_", length(params$container$filebreak.VR))
    params <- send_pause_continue_3p(params, filesA = files, from = "A",
                                  sleep_time = sleep_time, max_waiting_time = max_waiting_time)

    params <- ProcessXtWXCox.t3(params)

    if (params$failed) {
      warning(params$error_message)
      make_error_message(params$write_path, params$error_message)
      files = "error_message.rdata"
      params <- send_pause_continue_3p(params, filesA = files, filesB = files,
                                    from = c("A", "B"),
                                    sleep_time = sleep_time, max_waiting_time = max_waiting_time)
      params <- send_pause_quit_3p(params, sleep_time = sleep_time, job_failed = TRUE)
      SummarizeLog.3p(params)
      return(params$stats)
    }
    EndingIteration(params)
    params$algIterationCounter = params$algIterationCounter + 1
  }

  params <- ComputeResultsCox.t3(params)

  filesA = c("converged.rdata", "betasA.rdata", "stats.rdata")
  filesB = c("converged.rdata", "betasB.rdata", "stats.rdata")
  params <- send_pause_continue_3p(params, filesA = filesA, filesB = filesB,
                                from = c("A", "B"),
                                sleep_time = sleep_time, max_waiting_time = max_waiting_time)

  params <- send_pause_quit_3p(params, sleep_time = sleep_time)
  SummarizeLog.3p(params)

  return(invisible(params$stats))
}
