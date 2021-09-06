################### DISTRIBUTED COX REGRESSION FUNCTIONS ##################

PrepareDataCox.DP = function(params, data, yname, strata, mask) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataCox.23\n\n")

  workdata = list()
  workdata$failed = CheckDataFormat(params, data)
  if (workdata$failed) {
    return(workdata)
  }
  data = data.frame(data) # convert to a clean data.frame
  workdata$strata = extractStrata(params, data, strata, mask)
  if (workdata$strata$failed) {
    workdata$failed = TRUE
    return(workdata)
  }
  strataIndex = workdata$strata$strataIndex
  responseIndex = numeric()
  if (params$dataPartnerID == 1) {
    responseIndex = CheckResponse(params, data, yname)
    if (is.null(responseIndex)) {
      workdata$failed = TRUE
      return(workdata)
    }
    workdata$survival        = list()
    workdata$survival$rank   = data[, responseIndex[1]]
    workdata$survival$status = data[, responseIndex[2]]
    if (length(intersect(strataIndex, responseIndex)) > 0) {
      warning("Response and strata share a variable.")
      workdata$failed = TRUE
      return(workdata)
    }
  }
  covariateIndex = setdiff(1:ncol(data), union(strataIndex, responseIndex))
  workdata$n = nrow(data)
  if (length(covariateIndex) == 0) {
    if (params$dataPartnerID == 1) {
      workdata$X = matrix(0, nrow = nrow(data), ncol = 0)
    } else {
      warning("After removing strata, data is empty.  Party B must supply at least one non-strata covariate.")
      workdata$failed = TRUE
      return(workdata)
    }
  } else {
    workdata$tags = CreateModelMatrixTags(data[, covariateIndex, drop = FALSE])
    # if (params$dataPartnerID != 1 & (ncol(data) < 2 | !("numeric" %in% names(workdata$tags)))) {
    #   warning("A data partner that does not have the response must have at least 2 covariates at least one of which must be numeric.")
    #   workdata$failed = TRUE
    #   return(workdata)
    # }
    workdata$X = model.matrix(~ ., data[, covariateIndex, drop = FALSE])
    workdata$X = workdata$X[, -1, drop = FALSE]
    workdata$colmin   = apply(workdata$X, 2, min)
    workdata$colmax   = apply(workdata$X, 2, max)
    workdata$colsum   = apply(workdata$X, 2, sum)
    workdata$colrange = workdata$colmax - workdata$colmin
    for (i in 1:ncol(workdata$X)) {
      if (workdata$colmin[i] == workdata$colmax[i]) {
        workdata$colmin[i] = 0
        workdata$colrange[i] = 1
      }
      workdata$X[, i] = (workdata$X[, i] - workdata$colmin[i]) / workdata$colrange[i]
    }
  }

  return(workdata)
}


SendStrataNamesCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataNamesCox.DP\n\n")
  strataNames = list()
	strataNames$strataFromMe = data$strata$strataFromMe
	strataNames$strataFromOthers = data$strata$strataFromOthers
	writeTime = proc.time()[3]
	save(strataNames, file = file.path(params$writePath, "strata_names.rdata"))
	writeSize = file.size(file.path(params$writePath, "strata_names.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendStrataNamesCox.DP", 0, 0, writeTime, writeSize)
	return(params)
}


CheckStrataCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckStrataCox.DP\n\n")
  readTime = 0
	readSize = 0
	strataNames          = NULL
	strataClaimed        = rep(list(list()), params$numDataPartners)
	strataUnclaimed      = rep(list(list()), params$numDataPartners)
	if (length(data$strata$strataFromMe) > 0) {
		strataClaimed[[1]]   = data$strata$strataFromMe
	}
	if (length(data$strata$strataFromOthers) > 0) {
		strataUnclaimed[[1]] = data$strata$strataFromOthers
	}
	for (i in 2:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[i], "strata_names.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[i], "strata_names.rdata"))
		readTime = readTime + proc.time()[3]
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
		params$failed = TRUE
		params$errorMessage = "Data partners specified different strata:\n"
		for (i in 1:params$numDataPartners) {
			temp = NULL
			if (length(specified[[i]] > 0)) {
				temp = paste0(specified[[i]], collapse = ", ")
			}
			params$errorMessage = paste0(params$errorMessage,
																	 paste("Data Partner", i, "specified strata:", temp, "\n"))
		}
		params = AddToLog(params, "CheckStrataCox.DP", readTime, readSize, 0, 0)
		return(params)
	}
	passed = TRUE
	if (length(claimed) > 0) {
		tab = table(claimed)
		passed = all(tab == 1)
	}

	if (!passed) {
		params$failed = TRUE
		params$errorMessage = "The following strata are claimed by two or more data partners: "
		params$errorMessage = paste0(params$errorMessage, paste0(names(tab[which(tab > 1)]), collapse = ", "), "\n")
		params$errorMessage = paste0(params$errorMessage, "Make Sure that strata covariate names are unique to each data partner.")
		params = AddToLog(params, "CheckStrataCox.DP", readTime, readSize, 0, 0)
		return(params)
	}

	passed = TRUE
	claimed1 = sort(unique(claimed))
	passed = length(claimed1) == length(unclaimed) && all(claimed1 == unclaimed)

	if (!passed) {
		params$failed = TRUE
		params$errorMessage = "No data partner has the following specified strata: "
		params$errorMessage = paste0(params$errorMessage, paste0(unclaimed[which(!(unclaimed %in% claimed1))], collapse = ", "), ".")
		params = AddToLog(params, "CheckStrataCox.DP", readTime, readSize, 0, 0)
		return(params)
	}

	params = AddToLog(params, "CheckStrataCox.DP", readTime, readSize, 0, 0)
	return(params)
}


SendStrataCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataCox.DP\n\n")
  strata = list()
	strata$X = data$strata$X
	strata$legend = data$strata$legend
	writeTime = proc.time()[3]
	save(strata, file = file.path(params$writePath, "strata.rdata"))
	writeSize = file.size(file.path(params$writePath, "strata.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendStrataCox.DP", 0, 0, writeTime, writeSize)
	return(params)
}


PrepareStrataCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareStrataCox.DP\n\n")
  strata = NULL
	survival = data$survival
	strataTemp   = list()
	totalStrata = 0
	readSize = 0
	readTime = 0
	for (id in 1:params$numDataPartners) {
		totalStrata = totalStrata + length(params$strataClaimed[[id]])
	}

	if (totalStrata == 0) {
		strataTemp$X = data.frame(const__ = rep(1, params$n))
		strataTemp$legend = FALSE
	} else {
		first = TRUE
		strataTemp$legend = list()
		for (id in 1:params$numDataPartners) {
			if (length(params$strataClaimed[[id]]) > 0) {
				if (id == 1) {
					strataTemp$X = data$strata$X
					strataTemp$legend = data$strata$legend
					first = FALSE
				} else {
					readTime = readTime - proc.time()[3]
					load(file.path(params$readPathDP[[id]], "strata.rdata"))
					readSize = readSize + file.size(file.path(params$readPathDP[[id]], "strata.rdata"))
					readTime = readTime + proc.time()[3]
					if (first) {
						strataTemp$X = strata$X
						strataTemp$legend = strata$legend
						first = FALSE
					} else {
						strataTemp$X = cbind(strataTemp$X, strata$X)
						strataTemp$legend = c(strataTemp$legend, strata$legend)
					}
				}
			}
		}
	}

	sorted = do.call("order", cbind(strataTemp$X, survival$rank, survival$status))
	unsort = order(sorted)
	strataTemp$X = strataTemp$X[sorted, , drop = FALSE]
	survival$rank   = survival$rank[sorted]
	survival$status = survival$status[sorted]
	survival$sortedIdx = sorted
	survival$unsortIdx = unsort
	ranks = which(apply(abs(apply(strataTemp$X, 2, diff)), 1, sum) > 0)
	ranks = c(ranks, nrow(strataTemp$X))
	names(ranks) = NULL
	strata = rep(list(list()), length(ranks))
	if (length(ranks) == 1 && colnames(strataTemp$X)[1] == "const__") {
		strata[[1]]$start = 1
		strata[[1]]$end   = as.integer(length(survival$rank))
		strata[[1]]$label = ""
	} else {
		start = 1
		for (i in  1:length(ranks)) {
			strata[[i]]$start = start
			strata[[i]]$end   = as.integer(ranks[i])
			label = ""
			for (j in 1:ncol(strataTemp$X)) {
				temp = colnames(strataTemp$X)[j]
				label = paste0(label, temp, "=", strataTemp$legend[[temp]][strataTemp$X[start, j]])
				if (j < ncol(strataTemp$X)) {
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
		M = length(temp)   # number of unique observed times, including where no one fails
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
		strata[[i]]$J = as.integer(sum(temp0[, 2] > 0))    # number of distinct failure times
		# The number of failures at each rank which has a failure
		strata[[i]]$nfails = as.numeric(temp0[which(temp0[, 2] > 0), 2])
		# The first index of the ranks for which the number of failures is > 0.
		strata[[i]]$start0 = c(1, (cumsum(temp)[1:(M - 1)] + 1))[which(temp0[, 2] > 0)]
		# The first index of a failure for each rank which has a failure
		strata[[i]]$start1 = strata[[i]]$start0 + temp0[which(temp0[, 2] > 0), 1]
		# The last index of a failure for each rank which has a failure
		strata[[i]]$stop1  = as.numeric(strata[[i]]$start1 + strata[[i]]$nfails  - 1)
		strata[[i]]$start0 = as.numeric(strata[[i]]$start0 + strata[[i]]$start - 1)
		strata[[i]]$start1 = as.numeric(strata[[i]]$start1 + strata[[i]]$start - 1)
		strata[[i]]$stop1  = as.numeric(strata[[i]]$stop1  + strata[[i]]$start - 1)
	}

	survival$strata = strata
	survival$numEvents = sum(survival$status)
	params$survival = survival

	if (totalStrata == 0) {
		params$strata = as.matrix(strataTemp$X)
	} else {
		ptemp = 1
		for (i in 1:ncol(strataTemp$X)) {
			ptemp = ptemp + max(strataTemp$X[, i] - min(strataTemp$X[, i]))
		}
		params$strata = matrix(0, nrow = params$n, ncol = ptemp)
		colnames(params$strata) = paste0("strata:", 1:ptemp)
		params$strata[, 1] = 1
		idx = 2
		for (i in 1:ncol(strataTemp$X)) {
			min1 = min(strataTemp$X[, i])
			max1 = max(strataTemp$X[, i])
			if (min1 < max1) {
				for (j in (min1 + 1):max1) {
					params$strata[which(strataTemp$X[, i] == j), idx] = 1
					idx = idx + 1
				}
			}
		}
	}

	pStrata = ncol(params$strata)
	params$pStrata = pStrata

	writeTime = proc.time()[3]
	save(pStrata, survival, file = file.path(params$writePath, "survival.rdata"))
	writeSize = file.size(file.path(params$writePath, "survival.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareStrataCox.DP", readTime, readSize, writeTime, writeSize)

	return(params)
}


AddStrataToDataCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "AddStrataToDataCox.DP\n\n")
  data$X = cbind(params$strata, data$X)

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


PrepareParamsCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.DP\n\n")
  params$strata     = NULL
	params$n          = nrow(data$X)
	params$p          = ncol(data$X)

	temp = as.numeric(Sys.time())
	set.seed((temp - trunc(temp)) * .Machine$integer.max)
	params$seed       = floor(runif(1) * .Machine$integer.max)
	params$scaler     = 1 + runif(1)

	p = params$p
	seed = params$seed
	scaler = params$scaler

	writeTime = proc.time()[3]
	save(p, scaler, seed, file = file.path(params$writePath, "p_scaler_seed.rdata"))
	writeSize = file.size(file.path(params$writePath, "p_scaler_seed.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareParamsCox.DP", 0, 0, writeTime, writeSize)
	return(params)
}

PrepareSharesCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareSharesCox.DP\n\n")
  readTime = 0
	readSize = 0
	n = params$n
	p = params$p
	scaler = NULL
	seed   = NULL
	set.seed(params$seed, kind = "Mersenne-Twister")
	halfshare.L  = matrix(rnorm(n * p, sd = 20), nrow = n, ncol = p)

	halfshare.R  = data$X - halfshare.L

	products = rep(list(list()), params$numDataPartners)

	params$ps = c()
	params$scalers = c()
	params$seeds = c()

	for (id in 1:params$numDataPartners) {
		if (id == params$dataPartnerID) {
			products[[id]] = t(data$X) %*% data$X
			params$ps      = c(params$ps, params$p)
			params$scalers = c(params$scalers, params$scaler)
			params$seeds   = c(params$seeds, params$seed)
			next
		}
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "p_scaler_seed.rdata"))
		readTime = readTime + proc.time()[3]
		params$ps      = c(params$ps, p)
		params$scalers = c(params$scalers, scaler)
		params$seeds   = c(params$seeds, seed)

		set.seed(seed, kind = "Mersenne-Twister")
		halfShare2 = matrix(rnorm(params$n * p, sd = 20), nrow = params$n, ncol = p)

		if (id < params$dataPartnerID) {
			products[[id]] = t(halfShare2) %*% (data$X - scaler / (scaler + params$scaler) * halfshare.L)
		}

		if (id > params$dataPartnerID) {
			products[[id]] = t(data$X - scaler / (scaler + params$scaler) * halfshare.L) %*% halfShare2
		}
	}

	colmin    = data$colmin
	colrange  = data$colrange
	colsum    = data$colsum
	colnames  = colnames(data$X)
  tags      = data$tags

	writeTime = proc.time()[3]
	save(products, file = file.path(params$writePath, "products.rdata"))
	save(halfshare.R, file = file.path(params$writePath, "halfshare.rdata"))
	save(colmin, colrange, colsum, colnames, tags, file = file.path(params$writePath, "colstats.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("products.rdata",
																													"halfshare.rdata",
																													"colstats.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareSharesCox.DP", readTime, readSize, writeTime, writeSize)

	return(params)
}


GetStrata.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetStrata.AC\n\n")
  survival = NULL
  pStrata  = NULL
  readTime = proc.time()[3]
	load(file.path(params$readPathDP[1], "survival.rdata"))
	readSize = file.size(file.path(params$readPathDP[1], "survival.rdata"))
	readTime = proc.time()[3] - readTime
	params$survival = survival
	params$pStrata = pStrata
	params = AddToLog(params, "GetStrata.AC", readTime, readSize, 0, 0)
	return(params)
}


GetProductsCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsCox.AC\n\n")
  readTime = 0
	readSize = 0
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
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "products.rdata"))
		load(file.path(params$readPathDP[id], "halfshare.rdata"))
		load(file.path(params$readPathDP[id], "colstats.rdata"))
		readSize = readSize + sum(file.size(file.path(params$readPathDP[id],
																									c("products.rdata",
																										"halfshare.rdata",
																										"colstats.rdata"))))
		readTime = readTime + proc.time()[3]

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

	M = matrix(0, p, p)
	colnames(M) = allcolnames
	rownames(M) = allcolnames
	offset1 = 1
	params$pi = rep(0, params$numDataPartners)
	for (id1 in 1:params$numDataPartners) {
		p1 = ncol(allhalfshare[[id1]])
		params$pi[id1] = p1
		offset2 = offset1
		for (id2 in id1:params$numDataPartners) {
			p2 = ncol(allhalfshare[[id2]])
			if (id1 == id2) {
				M[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = allproducts[[id1]][[id2]]
			} else {
				temp = allproducts[[id1]][[id2]] + allproducts[[id2]][[id1]] +
					t(allhalfshare[[id1]]) %*% allhalfshare[[id2]]
				M[offset1:(offset1 + p1 - 1), offset2:(offset2 + p2 - 1)] = temp
				M[offset2:(offset2 + p2 - 1), offset1:(offset1 + p1 - 1)] = t(temp)
			}
			offset2 = offset2 + p2
		}
		offset1 = offset1 + p1
	}

	params$sts          = M
	params$n            = n
	params$p            = p
	params$colmin       = allcolmin
	params$colrange     = allcolrange
	params$colsum       = allcolsum
	params$colnames     = allcolnames
	params$party        = party
  params$tags         = alltags

	params$halfshare    = allhalfshare

	params = AddToLog(params, "GetProductsCox.AC", readTime, readSize, 0, 0)
	return(params)
}

CheckColinearityCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityCox.AC\n\n")
  sts = params$sts

	nrow = nrow(sts)
	indicies = c(1)
	for (i in 2:nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(sts[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
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
		idx = which(min <= indicies & indicies <= max)
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

  params$errorMessage = ""
  numeric_found = FALSE
  for (id in 2:params$numDataPartners) {
    if (length(unique(tags[[id]])) == 0) {
      params$failed = TRUE
      params$errorMessage = paste0(params$errorMessage,
                                   paste("After removing colinear covariates, Data Partner", id, "has no covariates."))
    }
  # else {
  #     numeric_found = numeric_found | "numeric" %in% names(tags[[id]])
  #   }
  # }
  # if (!numeric_found) {
  #   params$failed = TRUE
  #   params$errorMessage = paste0(params$errorMessage,
  #                                paste("After removing colinear covariates, no Data Partner > DP1 has a numeric covariate.\n"))
  }

	# if (params$numDataPartners == 2) {
	# 	if (length(params$indicies[[2]]) == 1) {
	# 		params$failed = TRUE
	# 		errorMessage = "Data Partner 2 has only one covariate. This is not secure.\n"
	# 		params$errorMessage = paste0(params$errorMessage, errorMessage)
	# 	} else if (length(params$indicies[[2]]) == 0) {
	# 		params$failed = TRUE
	# 		errorMessage = "All of Data Partner 2's covariates are colinear with Data Partner 1's covariates or the strata.\n"
	# 		params$errorMessage = paste0(params$errorMessage, errorMessage)
	# 	}
	# } else {
	# 	errorid = c()
	# 	for (id in 2:params$numDataPartners) {
	# 		if (length(params$indicies[[id]]) == 0) {
	# 			params$failed = TRUE
	# 			errorid = c(errorid, id)
	# 		}
	# 	}
	# 	if (params$failed) {
	# 		for (id in 1:length(errorid)) {
	# 			params$errorMessage = paste0(params$errorMessage,
	# 																	 paste("All of Data Partner", errorid[id], "covariates are colinear with covariates from other data partners or the strata.\n"))
	# 		}
	# 	}
	# }
	if (params$failed) {
		params = AddToLog(params, "CheckColinearityLogistic.AC", 0, 0, 0, 0)
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
	params$halfshare = do.call(cbind, params$halfshare)

	params$halfsharecolsum = apply(params$halfshare, 2, sum)

	cutoff = params$cutoff
	pReduct = params$pReduct

	writeTime = proc.time()[3]
	save(cutoff, idx, indicies, pReduct, file = file.path(params$writePath, "indicies.rdata"))
	writeSize = file.size(file.path(params$writePath, "indicies.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "CheckColinearityLogistic.AC", 0, 0, writeTime, writeSize)
	return(params)
}

ComputeUCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeUCox.AC\n\n")
  readTime = 0
	readSize = 0
	if (params$algIterationCounter == 1) {
		u = 1
	} else {
		uTemp = 0
		for (id in 1:params$numDataPartners) {
			readTime = readTime - proc.time()[3]
			load(file.path(params$readPathDP[id], "u.rdata"))
			readSize = readSize + file.size(file.path(params$readPathDP[id], "u.rdata"))
			readTime = readTime + proc.time()[3]
			uTemp = uTemp + u
		}
		u = uTemp
	}
	params$u = u
	writeTime = proc.time()[3]
	save(u, file = file.path(params$writePath, "u.rdata"))
	writeSize = file.size(file.path(params$writePath, "u.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeUCox.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}

UpdateParamsCox.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.DP\n\n")
  indicies = NULL
	idx      = NULL
	u        = NULL
	cutoff   = NULL
	pReduct  = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "indicies.rdata"))
	load(file.path(params$readPathAC, "u.rdata"))
	readSize = file.size(file.path(params$readPathAC, "indicies.rdata")) +
		file.size(file.path(params$readPathAC, "u.rdata"))
	readTime = proc.time()[3] - readTime
	betas = matrix(0, nrow = length(indicies[[params$dataPartnerID]]), ncol = 1)
	params$u             = u
	params$idx           = idx
	params$betas         = betas
	params$deltabeta.old = betas
	params$indicies      = indicies
	params$pReduct       = pReduct
	if (params$dataPartnerID == 1) {
		params$loglikelihood = -Inf
		params$sBeta.old = rep(0, params$n)
		params$cutoff    = cutoff
	}
	params = AddToLog(params, "UpdateParamsCox.DP", readTime, readSize, 0, 0)
	return(params)
}

UpdateDataCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.DP\n\n")
  idx = params$indicies[[params$dataPartnerID]]
	data$X = data$X[, idx, drop = FALSE]
	X = data$X
	data$colmin = data$colmin[idx]
	data$colmax = data$colmax[idx]
	data$colsum = data$colsum[idx]
	data$colrange = data$colrange[idx]

	if (params$dataPartnerID == 1) {
		halfshare = rep(list(list()), params$numDataPartners)
		for (id in 1:params$numDataPartners) {
			set.seed(params$seeds[id], kind = "Mersenne-Twister")
			halfshare[[id]] = matrix(rnorm(params$n * params$ps[id], sd = 20), nrow = params$n, ncol = params$ps[id])
			halfshare[[id]] = halfshare[[id]][, params$indicies[[id]], drop = FALSE]
		}
		data$halfshare = do.call(cbind, halfshare)
	}
	return(data)
}

ComputeSBetaCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSBetaCox.DP\n\n")
  u = NULL
	n = params$n
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "u.rdata"))
	readSize = file.size(file.path(params$readPathAC, "u.rdata"))
	readTime = proc.time()[3] - readTime

	sBetaPart = (data$X %*% params$betas + u) / (2 * u)
	V = 0
	for (id in 1:params$numDataPartners) {
		set.seed(params$seeds[id] + params$algIterationCounter)
		v = rnorm(n, mean = runif(n = 1, min = -1, max = 1), sd = 20)
		V = V + v
		if (id == params$dataPartnerID) {
			sBetaPart = sBetaPart + v
		}
	}

	sBetaPart = sBetaPart - params$scalers[params$dataPartnerID] / sum(params$scalers) * V

	writeTime = proc.time()[3]
	save(sBetaPart, file = file.path(params$writePath, "sbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "sbeta.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeSBetaCox.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}

GetSBetaCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetSBetaCox.AC\n\n")
  sBeta = 0
	sBetaPart = NULL
	readTime = 0
	readSize = 0
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "sbeta.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
		readTime = readTime + proc.time()[3]
		sBeta = sBeta + sBetaPart
	}
	sBeta = 2 * params$u * sBeta - params$u * params$numDataPartners
	params$sBeta = sBeta
	writeTime = proc.time()[3]
	save(sBeta, file = file.path(params$writePath, "sbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "sbeta.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetBetaCox.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}

ComputeLogLikelihoodCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeLogLikelihoodCox.DP\n\n")

	sBeta = 0
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "sbeta.rdata"))
	readSize = file.size(file.path(params$readPathAC, "sbeta.rdata"))
	readTime = proc.time()[3] - readTime
	sBeta = sBeta[params$survival$sortedIdx]
	loglikelihood.old = params$loglikelihood
	stephalving = TRUE
	scale = 1
	while (stephalving) {
		w = exp(sBeta)
		loglikelihood = 0
		stepCounter = 0
		pbar = MakeProgressBar1(params$survival$numEvents, "loglikelihood", params$verbose)
		for (i in 1:length(params$survival$strata)) {
			if (params$survival$strata[[i]]$J > 0) {
				for (j in 1:params$survival$strata[[i]]$J) {
					nj = params$survival$strata[[i]]$nfails[j]
					yIndex = params$survival$strata[[i]]$start0[j]:params$survival$strata[[i]]$end
					zIndex = params$survival$strata[[i]]$start1[j]:params$survival$strata[[i]]$stop1[j]
					Aj1 = sum(w[yIndex])
					Aj2 = sum(w[zIndex]) / nj
					loglikelihood = loglikelihood + sum(log(w[zIndex]))
					for (r in 0:(nj - 1)) {
						Ajr = Aj1 - r * Aj2
						loglikelihood = loglikelihood - log(Ajr)
					}
					stepCounter = stepCounter + nj
					pbar = MakeProgressBar2(stepCounter, pbar, params$verbose)
				}
			}
		}
		if (loglikelihood < loglikelihood.old && scale > 0.5^6) {
			sBeta = 0.5 * (sBeta + params$sBeta)
			scale = scale / 2
			if (params$verbose) cat("Step Halving\n\n")
		} else {
			stephalving = FALSE
		}
	}

	if (params$algIterationCounter == 1) {
		params$nullLoglikelihood = loglikelihood
	}

	converged = abs(loglikelihood - loglikelihood.old) / (abs(loglikelihood) + 0.1) < params$cutoff
	params$converged = converged
	params$scale     = scale
	params$sBeta     = sBeta
	params$loglikelihood = loglikelihood

	writeTime = proc.time()[3]
	save(sBeta, file = file.path(params$writePath, "sbeta.rdata"))
	save(scale, converged, file = file.path(params$writePath, "converged.rdata"))
	writeSize = file.size(file.path(params$writePath, "converged.rdata")) +
	  file.size(file.path(params$writePath, "sbeta.rdata"))
	writeTime = proc.time()[3] - proc.time()[3]

	params = AddToLog(params, "ComputeLogLikelihoodCox.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeSDelLCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSDeltaCox.AC\n\n")
  sBeta = 0
  readTime = proc.time()[3]
  load(file.path(params$readPathDP[1], "sbeta.rdata"))
  readSize = file.size(file.path(params$readPathDP[1], "sbeta.rdata"))
  readTime = proc.time()[3] - readTime
	halfshare = params$halfshare[params$survival$sortedIdx, , drop = FALSE]
	p = ncol(halfshare)
	n = params$n
	w = exp(sBeta)

	deltal = as.numeric(params$survival$status)
	deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
	                       # a pass by reference with the C call.
	W.S.R = matrix(0, n, p)
	numEvents = params$survival$numEvents

	.Call("ComputeCox", params$survival$strata, halfshare, w, deltal, W.S.R,
				as.integer(n), as.integer(p), as.integer(numEvents),
				as.integer(params$verbose))

	W.S.R = W.S.R[params$survival$unsortIdx, , drop = FALSE]
	tS.deltal.R = t(halfshare) %*% deltal
	W.S.R.1 = W.S.R[, params$idx[[1]], drop = FALSE]
	params$W.S.R = W.S.R
	params$tS.deltal.R = tS.deltal.R

	writeTime = proc.time()[3]
	save(W.S.R.1, file = file.path(params$writePath, "wsr1.rdata"))
	writeSize = file.size(file.path(params$writePath, "wsr1.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeSDelLCox.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeSDelLCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSDelLCox.DP\n\n")
  halfshare = data$halfshare[params$survival$sortedIdx, , drop = FALSE]
	p = ncol(halfshare)
	n = params$n
	w = exp(params$sBeta)

	deltal = as.numeric(params$survival$status)
	deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
	# a pass by reference with the C call.
	W.S.L = matrix(0, n, p)
	numEvents = params$survival$numEvents

	.Call("ComputeCox", params$survival$strata, halfshare, w, deltal, W.S.L,
				as.integer(n), as.integer(p), as.integer(numEvents),
				as.integer(params$verbose))

	W.S.L = W.S.L[params$survival$unsortIdx, , drop = FALSE]
	tS.deltal.L = t(halfshare) %*% deltal
	params$W.S.L = W.S.L
	params$tS.deltal.L = tS.deltal.L

	colmin.W.S.L   = apply(W.S.L, 2, min)
	colmax.W.S.L   = apply(W.S.L, 2, max)
	colrange.W.S.L = colmax.W.S.L - colmin.W.S.L

	for (i in 1:ncol(W.S.L)) {
		if (colmin.W.S.L[i] == colmax.W.S.L[i]) {
			colmin.W.S.L[i] = 0
		}
		if (colrange.W.S.L[i] == 0) {
			colrange.W.S.L[i] = 1
		}
	}

	scaled.W.S.L = W.S.L

	for (i in 1:ncol(scaled.W.S.L)) {
		scaled.W.S.L[, i] = (scaled.W.S.L[, i] - colmin.W.S.L[i]) / colrange.W.S.L[i]
	}

	temp = as.numeric(Sys.time())
	set.seed((temp - trunc(temp)) * .Machine$integer.max)
	scaled.W.S.L.L = matrix(rnorm(n * p, sd = 20), nrow = n, ncol = p)
	scaled.W.S.L.R = scaled.W.S.L - scaled.W.S.L.L

	params$colmin.W.S.L = colmin.W.S.L
	params$colrange.W.S.L = colrange.W.S.L
	params$scaled.W.S.L.L = scaled.W.S.L.L

	writeTime = proc.time()[3]
	save(colmin.W.S.L, colrange.W.S.L, scaled.W.S.L.R,
			 file = file.path(params$writePath, "scaledwslr.rdata"))
	save(tS.deltal.L, file = file.path(params$writePath, "tsdeltal.rdata"))
	save(scaled.W.S.L.L, file = file.path(params$writePath, "scaledwsll.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, "scaledwslr.rdata"),
														file.path(params$writePath, "tsdeltal.rdata"),
														file.path(params$writePath, "scaledwsll.rdata")))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeSDelLCox.DP", 0, 0, writeTime, writeSize)

	return(params)
}


ComputeProductsCox.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeProductsCox.DP\n\n")
  if (params$dataPartnerID == 1) {
		W.S.R.1 = NULL
		readTime = proc.time()[3]
		load(file.path(params$readPathAC, "wsr1.rdata"))
		readSize = file.size(file.path(params$readPathAC, "wsr1.rdata"))
		readTime = proc.time()[3] - readTime
		E = rep(list(list()), params$numDataPartners)
		E1 = rep(list(matrix(0, 0, 0)), params$numDataPartners)
		p = length(params$idx[[1]])
		if (p == 0) {
			for (id2 in 1:params$numDataPartners) {
				E1[[id2]] = matrix(0, nrow = 0, ncol = length(params$idx[[id2]]))
			}
		} else {
			E1[[1]] = t(data$X) %*% (params$W.S.L[, params$idx[[1]], drop = FALSE] + W.S.R.1)
			D = diag(params$colrange.W.S.L[params$idx[[1]]], ncol = p, nrow = p)
			for (id2 in 2:params$numDataPartners) {
				set.seed(params$seeds[id2], kind = "Mersenne-Twister")
				halfshare.RL = matrix(rnorm(params$n * params$ps[id2], sd = 20),
															nrow = params$n, ncol = params$ps[id2])  # needed to get randomization to right spot
				halfshare.RL = matrix(rnorm(params$n * params$ps[id2], sd = 20),
															nrow = params$n, ncol = params$ps[id2])[, params$indicies[[id2]], drop = FALSE]
				E1[[id2]] = solve(D) %*% t(data$X) %*% params$W.S.L[, params$idx[[id2]], drop = FALSE] +
					params$scaler / (params$scaler + params$scalers[id2]) *
					t(params$scaled.W.S.L.L[, params$idx[[1]], drop = FALSE]) %*% halfshare.RL
				}
		}
		E[[1]] = E1
		for (id1 in 2:params$numDataPartners) {
			E1 = rep(list(matrix(0, 0, 0)), params$numDataPartners)
			p = length(params$idx[[id1]])
			D = diag(params$colrange.W.S.L[params$idx[[id1]]], ncol = p, nrow = p)
			set.seed(params$seeds[id1], kind = "Mersenne-Twister")
			halfshare.L = matrix(rnorm(params$n * params$ps[id1], sd = 20),
													 nrow = params$n, ncol = params$ps[id1])[, params$indicies[[id1]], drop = FALSE]
			for (id2 in 2:params$numDataPartners) {
				set.seed(params$seeds[id2], kind = "Mersenne-Twister")
				halfshare.RL = matrix(rnorm(params$n * params$ps[id2], sd = 20),
											 			 nrow = params$n, ncol = params$ps[id2])
				halfshare.RL = matrix(rnorm(params$n * params$ps[id2], sd = 20),
															nrow = params$n, ncol = params$ps[id2])[, params$indicies[[id2]], drop = FALSE]
				E1[[id2]] = 0.5 * solve(D) %*% t(halfshare.L) %*% params$W.S.L[, params$idx[[id2]], drop = FALSE] +
					params$scaler / (params$scaler + params$scalers[id2]) *
					t(params$scaled.W.S.L.L[, params$idx[[id1]], drop = FALSE]) %*% halfshare.RL
			}
			E[[id1]] = E1
		}
		writeTime = proc.time()[3]
		save(E, file = file.path(params$writePath, "products.rdata"))
		writeSize = file.size(file.path(params$writePath, "products.rdata"))
		writeTime = proc.time()[3] - writeTime
	} else {
		scaled.W.S.L.L = NULL
		readTime = proc.time()[3]
		load(file.path(params$readPathDP[1], "scaledwsll.rdata"))
		readSize = file.size(file.path(params$readPathDP[1], "scaledwsll.rdata"))
		readTime = proc.time()[3] - readTime

		F = rep(list(list()), params$numDataPartners)
		set.seed(params$seed, kind = "Mersenne-Twister")
		halfshare.L = matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)[, params$indicies[[params$dataPartnerID]], drop = FALSE]
		halfshare.R = data$X - halfshare.L
		halfshare.R.L = matrix(rnorm(params$n * params$p, sd = 20), nrow = params$n, ncol = params$p)[, params$indicies[[params$dataPartnerID]], drop = FALSE]
		halfshare.R.R = halfshare.R - halfshare.R.L  # This is halfshare.R.L and halfshare.R.R
		for (id in 1:params$numDataPartners) {
			F[[id]] = params$scaler / (params$scalers[1] + params$scaler) * t(scaled.W.S.L.L[, params$idx[[id]], drop = FALSE]) %*% halfshare.R.L +
				t(scaled.W.S.L.L[, params$idx[[id]], drop = FALSE]) %*% halfshare.R.R
		}
		writeTime = proc.time()[3]
		save(F, file = file.path(params$writePath, "products.rdata"))
		writeSize = file.size(file.path(params$writePath, "products.rdata"))
		writeTime = proc.time()[3] - writeTime
	}
	params = AddToLog(params, "ComputeProductsCox.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateConvergeStatus.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateConvergeStatus.AC\n\n")
  converged = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathDP[1], "converged.rdata"))
	readSize = file.size(file.path(params$readPathDP[1], "converged.rdata"))
	readTime = proc.time()[3] - readTime
	params$converged = converged
	params = AddToLog(params, "UpdateConvergeStatus.AC", readTime, readSize, 0, 0)
}


ComputeStWSCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSCox.AC\n\n")
  E = NULL
	F = NULL
	scaled.W.S.L.R = NULL
	colmin.W.S.L = NULL
	colrange.W.S.L = NULL
	F1 = rep(list(list()), params$numDataPartners)
	readTime = proc.time()[3]
	load(file.path(params$readPathDP[1], "scaledwslr.rdata"))
	load(file.path(params$readPathDP[1], "products.rdata"))
	readSize = file.size(file.path(params$readPathDP[1], "scaledwslr.rdata")) +
		file.size(file.path(params$readPathDP[1], "products.rdata"))
	readTime = proc.time()[3] - readTime
	for (id in 2:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "products.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "products.rdata"))
		readTime = readTime + proc.time()[3]
		F1[[id]] = F
	}
	p = 0
	for (id in 1:params$numDataPartners) {
		p = p + length(params$idx[[id]])
	}
	M = matrix(0, p, p)
	startrow = 1
	for (id1 in 1:params$numDataPartners) {
		p1 = length(params$idx[[id1]])
		if (p1 == 0) next
		endrow = startrow + p1 - 1
		startcol = startrow
		D = diag(x = colrange.W.S.L[params$idx[[id1]]], nrow = p1, ncol = p1)
		for (id2  in id1:params$numDataPartners) {
			p2 = length(params$idx[[id2]])
		  endcol = startcol + p2 - 1
			if (p2 == 0) next
			if (id1 == 1 && id2 == 1) {
				M[startrow:endrow, startcol:endcol] = E[[1]][[1]]
			} else if (id1 == id2) {
				idx = params$idx[[id1]]
				G = D %*% (E[[id1]][[id1]] + F1[[id1]][[id1]] + t(scaled.W.S.L.R[, idx, drop = FALSE]) %*%
									 	params$halfshare[, idx, drop = FALSE]) +
									 	outer(colmin.W.S.L[idx], params$halfsharecolsum[idx])
				M[startrow:endrow, startcol:endcol] = G + t(G) +
					t(params$halfshare[, idx, drop = FALSE]) %*% params$W.S.R[, idx, drop = FALSE]
			} else {
				idx1 = params$idx[[id1]]
				idx2 = params$idx[[id2]]
				if (id1 == 1) {
					temp = D %*% (E[[id1]][[id2]] + F1[[id2]][[id1]] + t(scaled.W.S.L.R[, idx1, drop = FALSE]) %*%
													params$halfshare[, idx2, drop = FALSE]) +
						t(params$halfshare[, idx1, drop = FALSE]) %*% params$W.S.R[, idx2, drop = FALSE] +
						outer(colmin.W.S.L[idx1], params$halfsharecolsum[idx2])
					M[startrow:endrow, startcol:endcol] = temp
					M[startcol:endcol, startrow:endrow] = t(temp)
				} else {
					p1 = length(params$idx[[id2]])
					D2 = diag(x = colrange.W.S.L[params$idx[[id2]]], nrow = p2, ncol = p2)
					G23 = D %*% (E[[id1]][[id2]] + F1[[id2]][[id1]] + t(scaled.W.S.L.R[, idx1, drop = FALSE]) %*%
											 	params$halfshare[, idx2, drop = FALSE]) +
						outer(colmin.W.S.L[idx1], params$halfsharecolsum[idx2])
					G32 = D2 %*% (E[[id2]][[id1]] + F1[[id1]][[id2]] + t(scaled.W.S.L.R[, idx2, drop = FALSE]) %*%
											 	params$halfshare[, idx1, drop = FALSE]) +
						outer(colmin.W.S.L[idx2], params$halfsharecolsum[idx1])
					temp = G23 + t(G32) + t(params$halfshare[, idx1, drop = FALSE]) %*% params$W.S.R[, idx2, drop = FALSE]
					M[startrow:endrow, startcol:endcol] = temp
					M[startcol:endcol, startrow:endrow] = t(temp)
				}
			}
		  startcol = endcol + 1
		}
		startrow = endrow + 1
	}

	I = NULL
	tryCatch({I = solve(M)},
					 error = function(err) { I = NULL }
	)
	if (is.null(I)) {
		params$failed = TRUE
		params$singularMatrix = TRUE
		params$errorMessage =
			paste0("The matrix t(S)*W*S is not invertible.\n",
						 "       This may be due to one of two possible problems.\n",
						 "       1. Poor random initialization of the security halfshares.\n",
						 "       2. Near multicollinearity in the data\n",
						 "SOLUTIONS: \n",
						 "       1. Rerun the data analysis.\n",
						 "       2. If the problem persists, check the variables for\n",
						 "          duplicates for both parties and / or reduce the\n",
						 "          number of variables used. Once this is done,\n",
						 "          rerun the data analysis.")
		params = AddToLog(params, "computeStWSCox.AC", readTime, readSize, 0, 0)
		return(params)
	}

	params$I = I

	IDt = I %*% params$tS.deltal.R

	if (params$algIterationCounter == 1) {
		params$score = t(params$tS.deltal.R) %*% I %*% params$tS.deltal.R
	}
	params$maxIterExceeded = params$algIterationCounter > params$maxIterations
	maxIterExceeded = params$maxIterExceeded

	writeTime = 0
	writeSize = 0
	for (id in 1:params$numDataPartners) {
		I.part = I[params$idx[[id]], , drop = FALSE]
		IDt.part = IDt[params$idx[[id]], , drop = FALSE]
		writeTime = writeTime - proc.time()[3]
		save(I.part, IDt.part, file = file.path(params$writePath, paste0("update", id, ".rdata")))
		save(maxIterExceeded, file = file.path(params$writePath, "maxiterexceeded.rdata"))
		writeSize = writeSize + file.size(file.path(params$writePath, paste0("update", id, ".rdata"))) +
			file.size(file.path(params$writePath, "maxiterexceeded.rdata"))
		writeTime = writeTime + proc.time()[3]
	}
	params = AddToLog(params, "ComputeStWSCox.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateConvergeStatus.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateConvergeStatus.DP\n\n")
  readTime = 0
	readSize = 0
	scale    = NULL
	converged = NULL
	maxIterExceeded = NULL
	if (params$dataPartnerID > 1) {
		readTime = proc.time()[3]
		load(file.path(params$readPathDP[1], "converged.rdata"))
		readSize = file.size(file.path(params$readPathDP[1], "converged.rdata"))
		readTime = proc.time()[3] - readTime
		params$converged = converged
		params$scale     = scale
	}
	readTime = readTime - proc.time()[3]
	load(file.path(params$readPathAC, "maxiterexceeded.rdata"))
	readSize = readSize + file.size(file.path(params$readPathAC, "maxiterexceeded.rdata"))
	readTime = readTime + proc.time()[3]
	params$maxIterExceeded = maxIterExceeded
	params = AddToLog(params, "UpdateConvergeStatus.DP", readTime, readSize, 0, 0)
}


UpdateBetasCox.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetasCox.DP\n\n")
  I.part = NULL
	IDt.part = NULL
	tS.deltal.L = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, paste0("update", params$dataPartnerID, ".rdata")))
	readSize = file.size(file.path(params$readPathAC, paste0("update", params$dataPartnerID, ".rdata")))
	if (params$dataPartnerID > 1) {
		load(file.path(params$readPathDP[1], "tsdeltal.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[1], "tsdeltal.rdata"))
	}
	readTime = proc.time()[3] - readTime

	if (params$dataPartnerID == 1) {
		deltabeta = IDt.part + I.part %*% params$tS.deltal.L
		if (params$algIterationCounter == 1) {
			if (params$pReduct[1] == 0) {
				params$score = 0
			} else {
				idx = 1:params$pReduct[1]
				params$score = 2 * t(IDt.part) %*% params$tS.deltal.L[idx, 1, drop = FALSE] +
					t(I.part %*% params$tS.deltal.L) %*% params$tS.deltal.L[idx, 1, drop = FALSE]
			}
		}
	} else {
		deltabeta = IDt.part + I.part %*% tS.deltal.L
		if (params$algIterationCounter == 1) {
			temp = sum(params$pReduct[1:(params$dataPartnerID - 1)])
			idx = (temp + 1):(temp + params$pReduct[params$dataPartnerID])
			params$score = 2 * t(IDt.part) %*% tS.deltal.L[idx, 1, drop = FALSE] +
				t(I.part %*% tS.deltal.L) %*% tS.deltal.L[idx, 1, drop = FALSE]
		}
	}

	params$betas = params$betas + deltabeta - (1 - params$scale) * params$deltabeta.old
	params$deltabeta.old = deltabeta

	p = length(params$indicies[[params$dataPartnerID]])
	if (p == 0) {
		u = 1
	} else {
		u = sum(runif(n = p, min = 1, max = 2) * abs(params$betas))
	}

	if (params$dataPartnerID == 1) {
		loglikelihood = params$loglikelihood
		nullLoglikelihood = params$nullLoglikelihood
	}
	writeTime = proc.time()[3]
	save(u, file = file.path(params$writePath, "u.rdata"))
	writeSize = file.size(file.path(params$writePath, "u.rdata"))
	if (params$converged) {
		betasnew  = params$betas
		scorePart = params$score
		if (params$dataPartnerID == 1) {
			save(scorePart, nullLoglikelihood, loglikelihood, betasnew, file = file.path(params$writePath, "betas.rdata"))
		} else {
			save(scorePart, betasnew, file = file.path(params$writePath, "betas.rdata"))
		}
		writeSize = writeSize + file.size(file.path(params$writePath, "betas.rdata"))
	}
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "UpdateBetasCox.DP", readTime, readSize, writeTime, writeSize)

	return(params)
}


SurvFitCox.AC = function(params, pred) {
  if (params$trace) cat(as.character(Sys.time()), "SurvFitCox.AC\n\n")
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
				d = ndeath[i];
				if (d == 1) {
					sum1[i] = 1 / nrisk[i];
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


ComputeResultsCox.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsCox.AC\n\n")
  readSize          = 0
	betasnew          = NULL
	scorePart         = NULL
	loglikelihood     = NULL
	nullLoglikelihood = NULL
	score             = params$score
	readTime = proc.time()[3]
	betas = matrix(0, nrow = 0, ncol = 1)
	for (id in 1:params$numDataPartners) {
		load(file.path(params$readPathDP[id], "betas.rdata"))
		betas = rbind(betas, betasnew)
		score = score + scorePart
		readSize = readSize + file.size(file.path(params$readPathDP[id], "betas.rdata"))
	}
	readTime = proc.time()[3] - readTime

	stats = params$stats
	stats$failed         = FALSE
	stats$converged      = params$converged
	names.old            = params$colnames[-(1:params$pStrata)]
	idx                  = params$fullindicies - params$pStrata
	stats$party          = params$party[-(1:params$pStrata)]
	stats$coefficients   = rep(NA, length(stats$party))
	stats$coefficients[idx] = betas / params$colrange
	stats$expcoef      = exp(stats$coefficients)  # exp(coef) = hazard ratios
	stats$expncoef     = exp(-stats$coefficients)
	stats$var          = matrix(0, length(names.old), length(names.old))
	stats$var[idx, idx] = params$I
	stats$secoef       = rep(NA, length(names.old))
	stats$secoef[idx]  = sqrt(diag(params$I)) / params$colrange  # se(coef)

	stats$zvals        = stats$coefficients / stats$secoef  # z values
	stats$pvals        = 2 * pnorm(abs(stats$zvals), lower.tail = FALSE )   # pvals
	stats$stars        = matrix(sapply(stats$pvals, function(x) {
		if (is.na(x)) ''
		else if (x < 0.001) '***'
		else if (x < 0.01) '**'
		else if (x < 0.05) '*'
		else if (x < 0.1) '.'
		else ' '
	}))
	stats$lower95      = exp(stats$coefficients - qnorm(0.975) * stats$secoef)
	stats$upper95      = exp(stats$coefficients + qnorm(0.975) * stats$secoef)
	stats$loglik       = c(nullLoglikelihood, loglikelihood)
	stats$n            = params$n
	stats$nevent       = sum(params$survival$status)
	stats$df           = sum(params$pReduct)
	stats$iter         = params$algIterationCounter - 1
	stats$score        = c(score, 1 - pchisq(score, stats$df))
	stats$method       = "efron"
	stats$lrt          = 2*(stats$loglik[2] - stats$loglik[1])
	stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
	stats$rsquare      = c(1 - exp(-stats$lrt[1]/stats$n),
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
		results = survival::concordance(surv~pred + strata(strat))
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
		sorted = params$survival$sortedIdx,
		surv   = SurvFitCox.AC(params, pred)
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
	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeResultsCox.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetResultsCox.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.DP\n\n")
  stats = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "stats.rdata"))
	readSize = file.size(file.path(params$readPathAC, "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats

	params = AddToLog(params, "GetResultsCox.DP", readTime, readSize, 0, 0)
	return(params)
}


DoNothing.ACDP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "DoNothing\n\n")
  params = AddToLog(params, "--", 0, 0, 0, 0)
  return(params)
}

############################## PARENT FUNCTIONS ###############################


DataPartnerKCox = function(data,
													 yname           = NULL,
													 strata          = NULL,
													 mask            = TRUE,
													 numDataPartners = NULL,
													 dataPartnerID   = NULL,
													 monitorFolder   = NULL,
													 sleepTime       = 10,
													 maxWaitingTime  = 24 * 60 * 60,
													 popmednet       = TRUE,
													 trace           = FALSE,
													 verbose         = TRUE) {

	params = PrepareParams.kp("cox", dataPartnerID, numDataPartners, ac = FALSE,
	                          popmednet = popmednet, trace = trace, verbose = verbose)
	if (params$failed) {
	  warning(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		warning(params$errorMessage)
		return(invisible(NULL))
	}

	data = PrepareDataCox.DP(params, data, yname, strata, mask)
	params = AddToLog(params, "PrepareDataCox.DP", 0, 0, 0, 0)

	if (data$failed) {
		params$errorMessage = paste("Error processing data for data partner", params$dataPartnerID)
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		params$errorMessage = ReadErrorMessage(params$readPathAC)
		warning(params$errorMessage)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = SendBasicInfo.DP(params, data)
	files = "n_analysis.rdata"
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	possibleError = ReceivedError.kp(params, from = "AC")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	if (params$dataPartnerID == 1) {
    params = DoNothing.ACDP(params)
		params = SendPauseContinue.kp(params, filesAC = "empty.rdata", from = "DP",
											sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = CheckStrataCox.DP(params, data)

		if (params$failed) {
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = "errorMessage.rdata"
		} else {
			files = "empty.rdata"
		}
		params = SendPauseContinue.kp(params, filesDP = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = DoNothing.ACDP(params)

		if (params$failed) {
			warning(params$errorMessage)
			SendPauseQuit.kp(params, filesAC = "errorMessage.rdata", sleepTime = sleepTime)
			return(params$stats)
		} else {
		  params = SendPauseContinue.kp(params, filesDP = "empty.rdata", from = "DP",
		                                sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		}
		params = PrepareStrataCox.DP(params, data)
		data   = AddStrataToDataCox.DP(params, data)
		params = AddToLog(params, "AddStrataToDataCox.DP", 0, 0, 0, 0)
	} else {
		params = SendStrataNamesCox.DP(params, data)
		filesList = rep(list(list()), numDataPartners)
		filesList[[1]] = "strata_names.rdata"
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)


		possibleError = ReceivedError.kp(params, from = "DP1")
		if (possibleError$error) {
			params$errorMessage = possibleError$message
			warning(possibleError$message)
			params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
			return(params$stats)
		}

		params = SendStrataCox.DP(params, data)
		filesList = rep(list(list()), numDataPartners)
		filesList[[1]] = "strata.rdata"
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
	}

	params = PrepareParamsCox.DP(params, data)
	params = SendPauseContinue.kp(params, filesDP = "p_scaler_seed.rdata", from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	params = PrepareSharesCox.DP(params, data)
	files = c("products.rdata", "halfshare.rdata", "colstats.rdata")
	if (params$dataPartnerID == 1) {
		files = c(files, "survival.rdata")
	}
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	possibleError = ReceivedError.kp(params, from = "AC")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = UpdateParamsCox.DP(params)
	data = UpdateDataCox.DP(params, data)
	params = AddToLog(params, "UpdateDataCox.DP", 0, 0, 0, 0)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = ComputeSBetaCox.DP(params, data)

		if (params$dataPartnerID == 1) {
			files = "sbeta.rdata"
			params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

			params = ComputeLogLikelihoodCox.DP(params, data)

			files = "sbeta.rdata"
			params = SendPauseContinue.kp(params, filesAC = files, from = "DP2",
			                           sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

			params = ComputeSDelLCox.DP(params, data)

			files = c("tsdeltal.rdata", "scaledwsll.rdata", "converged.rdata")
			params = SendPauseContinue.kp(params, filesDP = files, from = "AC",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		} else if (params$dataPartnerID == 2) {
		  files = "sbeta.rdata"
		  params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
		                             sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		  params = DoNothing.ACDP(params)

		  filesList = rep(list(list()), numDataPartners)
			filesList[[1]] = "empty.rdata"
			params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		} else {
		  files = "sbeta.rdata"
		  params = SendPauseContinue.kp(params, filesAC = files, from = "DP1",
		                             sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		}

		params = ComputeProductsCox.DP(params, data)
		if (params$dataPartnerID == 1) {
			files = c("products.rdata", "scaledwslr.rdata", "converged.rdata")
		} else {
			files = c("products.rdata")
		}
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		possibleError = ReceivedError.kp(params, from = "AC")
		if (possibleError$error) {
			params$errorMessage = possibleError$message
			warning(possibleError$message)
			params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
			return(params$stats)
		}

		params = UpdateConvergeStatus.DP(params)
		params = UpdateBetasCox.DP(params)

		if (params$converged || params$maxIterExceeded) {
			files = "betas.rdata"
		} else {
			files = "u.rdata"
		}
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}
	params$lastIteration = TRUE
	params$completed = TRUE

	params = GetResultsCox.DP(params)
	SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}


AnalysisCenterKCox = function(numDataPartners = NULL,
															monitorFolder   = NULL,
															msreqid         = "v_default_0_000",
															cutoff          = 1E-8,
															maxIterations   = 25,
															sleepTime       = 10,
															maxWaitingTime  = 24 * 60 * 60,
															popmednet       = TRUE,
															trace           = FALSE,
															verbose         = TRUE) {

	filesList = rep(list(list()), numDataPartners)

	params = PrepareParams.kp("cox", 0, numDataPartners, msreqid, cutoff, maxIterations, ac = TRUE,
	                          popmednet = popmednet, trace = trace, verbose = verbose)
	if (params$failed) {
	  warning(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		warning(params$errorMessage)
		return(invisible(NULL))
	}

	params = PauseContinue.kp(params, from = "DP", maxWaitingTime = maxWaitingTime)

	possibleError = ReceivedError.kp(params, from = "DP")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		MakeErrorMessage(params$writePath, possibleError$message)
		files = "errorMessage.rdata"
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	params = CheckAgreement.AC(params)

	if (params$failed) {
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		warning(params$errorMessage)
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	files = "empty.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP1",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)


	params = DoNothing.ACDP(params)
	filesList = rep(list(list()), numDataPartners)
	filesList[[1]] = "empty.rdata"
	params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	possibleError = ReceivedError.kp(params, from = "DP")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		warning(possibleError$message)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	params = GetStrata.AC(params)
	params = GetProductsCox.AC(params)
	params = CheckColinearityCox.AC(params)

	if (params$failed) {
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		warning(params$errorMessage)
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = ComputeUCox.AC(params)

		if (params$algIterationCounter == 1) {
			files = c("indicies.rdata", "u.rdata")
		} else {
			files = "u.rdata"
		}
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = GetSBetaCox.AC(params)
		filesList = rep(list(list()), numDataPartners)
		filesList[[1]] = "sbeta.rdata"
		filesList[[2]] = "empty.rdata"
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ComputeSDelLCox.AC(params)
		filesList = rep(list(list()), numDataPartners)
		filesList[[1]] = "wsr1.rdata"
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = UpdateConvergeStatus.AC(params)
		params = ComputeStWSCox.AC(params)

		if (params$failed) {
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = "errorMessage.rdata"
			warning(params$errorMessage)
			params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
			params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
			SummarizeLog.kp(params)
			return(params$stats)
		}

		filesList = rep(list(list()), numDataPartners)
		for (id in 1:params$numDataPartners) {
			filesList[[id]] = c(paste0("update", id, ".rdata"), "maxiterexceeded.rdata")
		}
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}
	params$lastIteration = TRUE
	params$completed = TRUE

	params = ComputeResultsCox.AC(params)
	files = "stats.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
	SendPauseQuit.kp(params, sleepTime = sleepTime)
	SummarizeLog.kp(params)
	return(params$stats)
}
