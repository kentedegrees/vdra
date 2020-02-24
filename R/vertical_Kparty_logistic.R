################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ###################

PrepareDataLogistic.DP = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.DP\n\n")
  workdata = list()
	workdata$failed = FALSE

	if (class(data) %in% c("character", "double", "integer", "logical",
												 "numeric", "single", "factor")) {
		data = as.data.frame(data)
	}
	if (class(data) != "data.frame" && class(data) != "matrix") {
		cat("Error: data is not a matrix or a data frame.\n\n")
		workdata$failed = TRUE
		return(workdata)
	}
	if (nrow(data) < 1 || ncol(data) < 1) {
		cat("Error: the data is empty.\n\n")
		workdata$failed = TRUE
		return(workdata)
	}
	badValue = rep(FALSE, nrow(data))
	for (i in 1:ncol(data)) {
		if (class(data[, i]) %in% c("integer", "single", "double", "numeric")) {
			badValue = badValue | !is.finite(data[, i])
		} else {
			badValue = badValue | is.na(data[, i])
		}
	}
	idx = data.frame(which(badValue))
	colnames(idx) = "Observations with invalid entries"
	if (nrow(idx) > 0) {
		cat("Error: Some observations contain invalid values: NA, NaN, or InF.",
				"A list of all such observations has been outputted to",
				file.path(params$writePath, "invalidEntries.csv"),
				". Terinating program.\n\n")
		write.csv(idx, file.path(params$writePath, "invalidEntries.csv"))
		workdata$failed = TRUE
		return(workdata)
	}
	if ("(Intercept)" %in% colnames(data)) {
		cat("Error: \"(Intercept)\" is not a valid covariate name.  Please change it.  Terminiating Program.\n\n")
		workdata$failed = TRUE
		return(workdata)
	}
	if (max(table(colnames(data))) > 1) {
		cat("Error: duplicate column names found.\n\n")
		workdata$failed = TRUE
		return(workdata)
	}
	if (params$dataPartnerID == 1) {
		if (!is.null(yname) && class(yname) != "character") {
			cat("Error: response label is not a string.\n\n")
			workdata$failed = TRUE
			return(workdata)
		}
		if (length(yname) > 1) {
			cat("Error: More than one name for response variable provided.\n\n")
			workdata$failed = TRUE
			return(workdata)
		}
		if (is.null(colnames(data))) {
			cat("Warning: variables are not named. Assuming variable 1 is response.  Assigning labels to the rest of the columns.\n\n")
			if (is.null(yname)) {
				yname = c("y")
			}
			if (ncol(data) == 1) {
				colnames(data) = yname
			} else {
				colnames(data) = c(yname, paste0("A", 1:(ncol(data) - 1)))
			}
		} else {
			if (is.null(yname)) {
				yname = colnames(data)[1]
			}
		}
		if (!(yname %in% colnames(data))) {
			cat("Error: response variable", paste0("'", yname, "'"), "not found.\n\n")
			workdata$failed = TRUE
			return(workdata)
		}
		responseColIndex = which(colnames(data) %in% yname)
		if (class(data[1, responseColIndex]) != "numeric" & class(data[1, responseColIndex]) != "integer") {
			cat("Error: response variable", paste0("'", yname, "'"), "is not numeric.\n\n")
			workdata$failed = TRUE
			return(workdata)
		}
		workdata$X = matrix(data[, responseColIndex], ncol = 1)
		workdata$X = cbind(workdata$X, 1)
		if (sum(!(workdata$X[, 1] %in% c(0, 1))) > 0) {
			cat("Error: response variable is not binary.  It should only be 0's and 1's.  Terminating program.\n\n")
			workdata$failed = TRUE
			return(workdata)
		}
		data = data[, -responseColIndex, drop = FALSE]
		colnames(workdata$X) = c(yname, "(Intercept)")
	} else {
		workdata$X = matrix(0, nrow = nrow(data), ncol = 0)
	}
	if (class(data) == "matrix") {
		workdata$X = cbind(workdata$X, data)
	} else if (ncol(data) >= 1) { # class(data) == "data.frame"
		varNames = c()
		covariateNames = setdiff(colnames(data), yname)
		for (name in covariateNames) {
			if (class(data[[name]]) %in% c("integer", "numeric", "single", "double")) {
				varNames = c(varNames, name)
			} else {
				if (class(data[[name]]) %in% c("character", "logical")) {
					data[[name]] = as.factor(data[[name]])
				}
				if (class(data[[name]]) != "factor") {
					cat("Error: variable", name, "is not numeric and cannot be converted to a factor.\n\n")
					workdata$failed = TRUE
					return(workdata)
				}
				levelNames = levels(data[[name]])
				if (length(levelNames) == 1) {
					data[[name]] = rep(1, nrow(data))
					varNames = c(varNames, name)
				} else {
					for (lname in levelNames[-1]) {
						newName = paste0(name, ":", lname)
						if (newName %in% varNames) {
							cat("Error: variable", newName, "already exists.\n\n")
							workdata$failed = TRUE
							return(workdata)
						}
						data[[newName]] = ifelse(data[[name]] == lname, 1, 0)
						varNames = c(varNames, newName)
					}
					data[[name]] = NULL
				}
			}
		}
		index = match(varNames, colnames(data))
		workdata$X = cbind(workdata$X, as.matrix(data[index]))
	}

	workdata$n        = nrow(workdata$X)
	workdata$colmin   = apply(workdata$X, 2, min)
	workdata$colmax   = apply(workdata$X, 2, max)
	workdata$colsum   = apply(workdata$X, 2, sum)
	workdata$colrange = workdata$colmax - workdata$colmin
	for (i in 1:ncol(workdata$X)) {
		if (workdata$colmin[i] == workdata$colmax[i]) {
			workdata$colmin[i] = 0
			workdata$colrange[i] = workdata$colmax[i]
			if (workdata$colrange[i] == 0) {
				workdata$colrange[i] = 1
			}
		}
	}
	for (i in 1:ncol(workdata$X)) {
		workdata$X[, i] = (workdata$X[, i] - workdata$colmin[i]) / workdata$colrange[i]
	}
	return(workdata)
}


GetProductsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLogistic.AC\n\n")
	readTime = 0
	readSize = 0
	p = 0
	n = 0
	pi = c()

	allproducts = rep(list(list()), params$numDataPartners)
	allhalfshare = rep(list(list()), params$numDataPartners)
	products = NULL
	halfshare = NULL
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
		allhalfshare[[id]] = halfshare
		allcolmin          = c(allcolmin, colmin)
		allcolrange        = c(allcolrange, colrange)
		allcolsum          = c(allcolsum, colsum)
		allcolnames        = c(allcolnames, colnames)
		party              = c(party, rep(paste0("dp", id), length(colnames)))
		p = p + ncol(halfshare)
		pi = c(pi, ncol(halfshare))
		if (id == 1) n = nrow(halfshare)
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


	params$halfshare    = allhalfshare
	params$sts          = M[2:p, 2:p, drop = FALSE]
	params$sty          = M[2:p, 1, drop = FALSE]
	params$yty          = M[1, 1]
	params$meansy       = allcolsum[1] / n
	params$means        = allcolsum[-1] / n
	params$n            = n
	params$p            = p
	params$pi           = pi
	params$colmin       = allcolmin[-1]
	params$colrange     = allcolrange[-1]
	params$colsum       = allcolsum[-1]
	params$colnames     = allcolnames[-1]
	params$party        = party[-1]

	params = AddToLog(params, "GetProductsLogistic.AC", readTime, readSize, 0, 0)
	return(params)
}


CheckColinearityLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.AC\n\n")
  sts = params$sts
	sty = params$sty

	nrow = nrow(sts)
	indicies = c(1)
	for (i in 2:nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(sts[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}

	sts = sts[indicies, indicies, drop = FALSE]
	sty = sty[indicies, drop = FALSE]

	params$sts = sts
	params$sty = sty

	# Extract the indicies to keep for each party and check for errors.

	params$colmin       = params$colmin[indicies]
	params$colrange     = params$colrange[indicies]
	params$colsum       = params$colsum[indicies]
	params$fullindicies = indicies       # To be used when computing stats
	params$p            = params$p - 1   # Get rid of the response from the count

	indicies = indicies + 1  # take into account that pi still counts sty, which we removed earlier.

	params$indicies = rep(list(list()), params$numDataPartners)

	min = 1
	for (id in 1:params$numDataPartners) {
		max = min + params$pi[id] - 1
		params$indicies[[id]] =  indicies[which(min <= indicies & indicies <= max)] - min + 1
		min = max + 1
	}

	if (params$numDataPartners == 2) {
		if (length(params$indicies[[1]]) == 2) {
			params$failed = TRUE
			params$errorMessage = "Data Partner 1 has only one covaraite. This is not secure.\n"
		}
		if (length(params$indicies[[2]]) == 1) {
			params$failed = TRUE
			errorMessage = "Data Partner 2 has only one covariate. This is not secure.\n"
			params$errorMessage = paste0(params$errorMessage, errorMessage)
		} else if (length(params$indicies[[2]]) == 0) {
			params$failed = TRUE
			errorMessage = "All of Data Partner 2's covariates are coliner with Data Partner 1's covariates.\n"
			params$errorMessage = paste0(params$errorMessage, errorMessage)
		}
	} else {
		errorid = c()
		for (id in 1:params$numDataPartners) {
			if (length(params$indicies[[id]]) == 0) {
				params$failed = TRUE
				errorid = c(errorid, id)
			}
		}
		if (params$failed) {
			for (id in 1:length(errorid)) {
				params$errorMessage = paste0(params$errorMessage,
						 paste("All of Data Partner", errorid[id], "covariates are colinear with covariates from other data partners.\n"))
			}
		}
	}
	indicies = params$indicies

	params$pReduct = c()
	for (id in 1:params$numDataPartners) {
		params$pReduct = c(params$pReduct, length(indicies[[id]]))
	}

	for (id in 1:params$numDataPartners) {
		params$halfshare[[id]] = params$halfshare[[id]][, indicies[[id]], drop = FALSE]
	}

	writeTime = proc.time()[3]
	save(indicies, file = file.path(params$writePath, "indicies.rdata"))
	writeSize = file.size(file.path(params$writePath, "indicies.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "CheckColinearityLogistic.AC", 0, 0, writeTime, writeSize)

	return(params)
}


ComputeInitialBetasLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.AC\n\n")
  writeTime = 0
	writeSize = 0
	colsumS = (params$colsum - params$n * params$colmin) / params$colran
	beta = 4 * solve(params$sts) %*% (params$sty - 0.5 * colsumS)

	u = sum(runif(length(beta), min = 1, max = 5) * abs(beta))
	params$u = u
	start = 1
	for (id in 1:params$numDataPartners) {
		end = start + length(params$indicies[[id]]) - 1
		betas = beta[start:end]

		writeTime = writeTime - proc.time()[3]
		save(u, betas, file = file.path(params$writePath, paste0("u_beta_", id, ".rdata")))
		writeSize = writeSize + file.size(file.path(params$writePath, paste0("u_beta_", id, ".rdata")))
		writeTime = writeTime + proc.time()[3]
		start = end + 1
	}
	params = AddToLog(params, "ComputeInitialBetasLogistic.AC", 0, 0, writeTime, writeSize)
	return(params)
}


UpdateParamsLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.DP\n\n")
  indicies = NULL
	u = NULL
	betas = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "indicies.rdata"))
	filename = paste0("u_beta_", params$dataPartnerID, ".rdata")
	load(file.path(params$readPathAC, filename))
	readSize = file.size(file.path(params$readPathAC, "indicies.rdata")) +
		file.size(file.path(params$readPathAC, filename))
	readTime = proc.time()[3] - readTime
	params$u = u
	params$betas = betas
	params$indicies = indicies
	params = AddToLog(params, "UpdateParamsLogistic.DP", readTime, readSize, 0, 0)
	return(params)
}


UpdateDataLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.DP\n\n")
  if (params$dataPartnerID == 1) {
		data$Y = data$X[, 1, drop = FALSE]
	}
	idx = params$indicies[[params$dataPartnerID]]
	data$X = data$X[, idx, drop = FALSE]
	data$colmin = data$colmin[idx]
	data$colmax = data$colmax[idx]
	data$colsum = data$colsum[idx]
	data$colrange = data$colrange[idx]
	return(data)
}


ComputeSbetaLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeSbetaLogistic.DP\n\n")
  set.seed(params$seed + params$algIterationCounter, kind = "Mersenne-Twister")
	V = matrix(rnorm(params$n, mean = runif(n = 1, min = -1, max = 1), sd = 10), ncol = 1)
	Vsum = 0
	for (id in 1:params$numDataPartners) {
		set.seed(params$seeds[id] + params$algIterationCounter, kind = "Mersenne-Twister")
		Vsum = Vsum + matrix(rnorm(params$n, mean = runif(n = 1, min = -1, max = 1), sd = 10), ncol = 1)
	}

	Sbeta = (data$X %*% params$betas + params$u) / (2 * params$u) + V - params$scaler / sum(params$scalers) * Vsum

	writeTime = proc.time()[3]
	save(Sbeta, file = file.path(params$writePath, "sbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "sbeta.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeSbetaLogistic.DP", 0, 0, writeTime, writeSize)
	return(params)
}


ComputeWeightsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeWeightsLogistic.AC\n\n")
  Sbeta = 0
	readTime = 0
	readSize = 0
	sbeta = 0
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "sbeta.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
		readTime = readTime + proc.time()[3]
		sbeta = sbeta + Sbeta
	}
	sbeta = 2 * params$u * sbeta - params$numDataPartners * params$u
	pi_ = 1 / (1 + exp(-sbeta))
	params$pi_ = pi_

	writeTime = proc.time()[3]
	save(pi_, file = file.path(params$writePath, "pi.rdata"))
	writeSize = file.size(file.path(params$writePath, "pi.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComptueWeightsLogistic.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeStWSLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSLogistic.DP\n\n")
  pi_ = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "pi.rdata"))
	readSize = file.size(file.path(params$readPathAC, "pi.rdata"))
	readTime = proc.time()[3] - readTime
	params$pi_ = pi_

	W = pi_ * (1 - pi_)
	C = rep(list(list()), params$numDataPartners)

	idx = params$indicies[[params$dataPartnerID]]
	set.seed(params$seed, kind = "Mersenne-Twister")
	halfshare = matrix(rnorm(params$n * params$p, sd = 20),
										 nrow = params$n, ncol = params$p)[, idx, drop = FALSE]

	for (id in 1:params$numDataPartners) {
		if (id < params$dataPartnerID) {
			set.seed(params$seeds[id], kind = "Mersenne-Twister")
			idx = params$indicies[[id]]
			halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
								  				 nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
			C[[id]] = params$scaler / (params$scaler + params$scalers[id]) *
				t(halfshareDP) %*% MultiplyDiagonalWTimesX(W, halfshare) +
				t(halfshareDP) %*% MultiplyDiagonalWTimesX(W, data$X - halfshare)
		} else if (id == params$dataPartnerID) {
			C[[id]] = t(data$X) %*% MultiplyDiagonalWTimesX(W, data$X)
		} else { # id > params$dataPartnerID
			set.seed(params$seeds[id], kind = "Mersenne-Twister")
			idx = params$indicies[[id]]
			halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
													 nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
			C[[id]] = params$scaler / (params$scaler + params$scalers[id]) *
				t(halfshare) %*% MultiplyDiagonalWTimesX(W, halfshareDP) +
				t(data$X - halfshare) %*% MultiplyDiagonalWTimesX(W, halfshareDP)
		}
	}

	writeTime = proc.time()[3]
	save(C, file = file.path(params$writePath, "stwsshare.rdata"))
	writeSize = file.size(file.path(params$writePath, "stwsshare.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeStWSLogistic.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeStWSLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeStWSLogistic.AC\n\n")
  readTime = 0
	readSize = 0
	W = params$pi_ * (1 - params$pi_)

	StWS = matrix(0, sum(params$pReduct), sum(params$pReduct))

	for (id1 in 1:params$numDataPartners) {
		end = sum(params$pReduct[1:id1])
		start = end - params$pReduct[id1] + 1
		idx1 = start:end
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id1], "stwsshare.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id1], "stwsshare.rdata"))
		readTime = readTime + proc.time()[3]
		for (id2 in 1:params$numDataPartners) {
			end = sum(params$pReduct[1:id2])
			start = end - params$pReduct[id2] + 1
			idx2 = start:end
			if (id1 < id2) {
				StWS[idx1, idx2] = StWS[idx1, idx2] + C[[id2]]
				StWS[idx2, idx1] = StWS[idx2, idx1] + t(C[[id2]])
			} else if (id1 == id2) {
				StWS[idx1, idx1] = C[[id1]]
			} else {
				StWS[idx2, idx1] = StWS[idx2, idx1] + C[[id2]]
				StWS[idx1, idx2] = StWS[idx1, idx2] + t(C[[id2]])
			}
		}
		if (id1 < params$numDataPartners) {
			for (id2 in (id1 + 1):params$numDataPartners) {
				end = sum(params$pReduct[1:id2])
				start = end - params$pReduct[id2] + 1
				idx2 = start:end
				temp = t(params$halfshare[[id1]]) %*% MultiplyDiagonalWTimesX(W, params$halfshare[[id2]])
				StWS[idx1, idx2] = StWS[idx1, idx2] + temp
				StWS[idx2, idx1] = StWS[idx2, idx1] + t(temp)
			}
		}
	}

	I = NULL
	tryCatch({I = solve(StWS)},
					 error = function(err) { I = NULL }
	)
	if (is.null(I)) {
		params$failed = TRUE
		params$singularMatrix = TRUE
		params$errorMessage =
			paste0("ERROR: The matrix t(X)*W*X is not invertible.\n",
						 "       This may be due to one of two possible problems.\n",
						 "       1. Poor random initilization of the security matricies.\n",
						 "       2. Near multicolinearity in the data\n",
						 "SOLUTIONS: \n",
						 "       1. Rerun the data analaysis.\n",
						 "       2. If the problem persists, check the variables for\n",
						 "          duplicates for both parties and / or reduce the\n",
						 "          number of variables used. Once this is done,\n",
						 "          rerun the data analysis.\n\n")
		params = AddToLog(params, "ComputeStWSLogistic.AC", readTime, readSize, 0, 0)
		return(params)
	}
	params$I = I
	halfshare = params$halfshare[[1]]
	for (id in 2:params$numDataPartners) {
		halfshare = cbind(halfshare, params$halfshare[[id]])
	}
	IDt = I %*% (params$sty - t(halfshare) %*% params$pi_)
	Itemp = I
	IDttemp = IDt

	writeTime = 0
	writeSize = 0
	start = 1
	stop  = params$pReduct[1]
	for (id in 1:params$numDataPartners) {
		I = Itemp[start:stop, , drop = FALSE]
		IDt = IDttemp[start:stop, , drop = FALSE]
    writeTime = writeTime - proc.time()[3]
    save(I, IDt, file = file.path(params$writePath, paste0("ID", id, ".rdata")))
    writeSize = writeSize + file.size(file.path(params$writePath, paste0("ID", id, ".rdata")))
    writeTime = writeTime + proc.time()[3]
		start = stop + 1
		stop = stop + params$pReduct[id + 1]
	}

	params = AddToLog(params, "ComputeStWSLogistic.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateBetaLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.DP\n\n")
  I = IDt = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, paste0("ID", params$dataPartnerID, ".rdata")))
	readSize = file.size(file.path(params$readPathAC, paste0("ID", params$dataPartnerID, ".rdata")))
	readTime = proc.time()[3] - readTime

	id = 1
	set.seed(params$seeds[id], kind = "Mersenne-Twister")
	idx = params$indicies[[id]]
	halfshareDP = matrix(rnorm(params$n * params$ps[id], sd = 20),
											 nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE]
	for (id in 2:params$numDataPartners) {
		set.seed(params$seeds[id], kind = "Mersenne-Twister")
		idx = params$indicies[[id]]
		halfshareDP = cbind(halfshareDP,
												matrix(rnorm(params$n * params$ps[id], sd = 20),
															 nrow = params$n, ncol = params$ps[id])[, idx, drop = FALSE])
	}

	D0 = t(halfshareDP) %*% params$pi_
	deltaBeta = IDt - I %*% D0
	params$betas = params$betas + deltaBeta
	maxdifference = max(abs(deltaBeta) / (abs(params$betas) + .1))
	utemp = sum(runif(length(deltaBeta), min = 1, max = 5) * abs(params$betas))

	writeTime = proc.time()[3]
	save(utemp, maxdifference, file = file.path(params$writePath, "u_converge.rdata"))
	writeSize = file.size(file.path(params$writePath, "u_converge.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "UpdateBetaLogistic.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeConvergeStatusLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeConvergeStatusLogistic.AC\n\n")
  readTime = 0
	readSize = 0
	u = 0
	converged = TRUE
	utemp = NULL
	maxdifference = NULL
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "u_converge.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "u_converge.rdata"))
		readTime = readTime + proc.time()[3]
		u = u + utemp
		converged = converged && (maxdifference < params$cutoff)
	}
  maxIterExceeded = params$algIterationCounter >= params$maxIterations
	params$maxIterExceeded = maxIterExceeded
	params$u = u
	params$converged = converged
	writeTime = proc.time()[3]
	save(u, converged, maxIterExceeded, file = file.path(params$writePath, "u_converge.rdata"))
	writeSize = file.size(file.path(params$writePath, "u_converge.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeConvergeStatusLogistic.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetConvergeStatusLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetconvergeStatusLogistic.DP\n\n")
  u = converge = maxIterExceeded = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "u_converge.rdata"))
	readSize = file.size(file.path(params$readPathAC, "u_converge.rdata"))
	readTime = proc.time()[3] - readTime
	params$u = u
	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	params = AddToLog(params, "GetConvergeStatusLogistic.DP", readTime, readSize, 0, 0)
	return(params)
}

SendFinalBetasLogistic.DP = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "SendFinalBetasLogistic.DP\n\n")
  betas = params$betas
	writeTime = proc.time()[3]
	save(betas, file = file.path(params$writePath, "finalbetas.rdata"))
	writeSize = file.size(file.path(params$writePath, "finalbetas.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendFinalBetasLogistic.DP", 0, 0, writeTime, writeSize)
	return(params)
}

ComputeFinalSBetaLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeFinalSBetaLogistic.AC\n\n")
  Sbeta = 0
	readTime = 0
	readSize = 0
	sbeta = 0
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "sbeta.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "sbeta.rdata"))
		readTime = readTime + proc.time()[3]
		sbeta = sbeta + Sbeta
	}
	sbeta = 2 * params$u * sbeta - params$numDataPartners * params$u

	writeTime = proc.time()[3]
	save(sbeta, file = file.path(params$writePath, "sbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "sbeta.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeFinalSBetaLogistic.AC", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeResultsLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.DP\n\n")
  sbeta = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "sbeta.rdata"))
	readSize = file.size(file.path(params$readPathAC, "sbeta.rdata"))
	readTime = proc.time()[3] - readTime

	n       = params$n
	ct      = sum(data$Y)
	params$FinalFitted = sbeta
	resdev  = -2 * (sum(data$Y * sbeta) - sum(log(1 + exp(sbeta))))
	nulldev = -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))

	hoslem  = HoslemInternal(params, data)
	ROC     = RocInternal(params, data)

	writeTime = proc.time()[3]
	save(resdev, nulldev, hoslem, ROC, file = file.path(params$writePath, "logisticstats.rdata"))
	writeSize = file.size(file.path(params$writePath, "logisticstats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeResultsLogistic.DP", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeResultsLogistic.AC = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.AC\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPathDP[1], "logisticstats.rdata"))
	readSize = file.size(file.path(params$readPathDP[1], "logisticstats.rdata"))
	readTime = proc.time()[3] - readTime
	coefficients = c()
	p            = 0
	betas = NULL
	for (id in 1:params$numDataPartners) {
		readTime = readTime - proc.time()[3]
		load(file.path(params$readPathDP[id], "finalbetas.rdata"))
		readSize = readSize + file.size(file.path(params$readPathDP[id], "finalbetas.rdata"))
		readTime = readTime + proc.time()[3]
		coefficients = c(coefficients, betas)
		p            = p + length(params$indicies[[id]])
	}

	coefficients[2:p] = coefficients[2:p] / params$colran[2:p]
	coefficients[1] = coefficients[1] - sum(coefficients[2:p] * params$colmin[2:p])

	serror = rep(0, p)
	serror[2:p] = sqrt(diag(params$I)[2:p]) / params$colran[2:p]
	d1 = diag(c(1, params$colmin[-1] / params$colran[-1]))
	temp = d1 %*% params$I %*% d1
	serror[1] = sqrt(temp[1, 1] - 2 * sum(temp[1, 2:p]) + sum(temp[2:p, 2:p]))

	stats = params$stats
	stats$failed         = FALSE
	stats$converged      = params$converged


	# If xtwx were singular, it would have been caught in GetII.A2(), so we may
	# assume that xtwx is NOT singular and so we do not have to do a check.
	stats$party = params$party
	stats$coefficients = rep(NA, params$p)
	stats$secoef = rep(NA, params$p)
	stats$tvals  = rep(NA, params$p)
	stats$pvals  = rep(NA, params$p)
	stats$n  = params$n
	stats$nulldev = nulldev
	stats$resdev = resdev
	stats$aic = resdev + 2 * sum(params$pReduct)
	stats$bic = resdev + sum(params$pReduct) * log(params$n)
	stats$nulldev_df = params$n - 1
	stats$resdev_df = params$n - sum(params$pReduct)
	stats$coefficients[params$fullindicies] = coefficients
	stats$secoef[params$fullindicies] = serror
	tvals = coefficients / serror
	pvals = 2 * pnorm(abs(tvals), lower.tail = FALSE)
	stats$tvals[params$fullindicies] = tvals
	stats$pvals[params$fullindicies] = pvals
	stats$hoslem  = hoslem
	stats$ROC     = ROC
	stats$iter    = params$algIterationCounter - 1
	names(stats$coefficients) = params$colnames
	names(stats$party) = params$colnames
	names(stats$secoef) = params$colnames
	names(stats$tvals) = params$colnames
	names(stats$pvals) = params$colnames

	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params$stats      = stats

	params = AddToLog(params, "ComputeResultsLogistic.AC", readTime, readSize, writeTime, writeSize)
	return(params)}


GetResultsLogistic.DP = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.DP\n\n")
  stats = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPathAC, "stats.rdata"))
	readSize = file.size(file.path(params$readPathAC, "stats.rdata"))
	readTime = proc.time()[3] - readTime
	if (params$dataPartnerID == 1) {
		stats$Y           = data$Y # For Hoslem and ROC
		stats$FinalFitted = params$FinalFitted
	}
	params$stats      = stats
	params = AddToLog(params, "GetResultsLogistic.DP", readTime, readSize, 0, 0)
	return(params)
}

############################## PARENT FUNCTIONS ###############################

DataPartnerKLogistic = function(data,
																yname           = NULL,
																numDataPartners = NULL,
																dataPartnerID   = NULL,
																monitorFolder   = NULL,
																sleepTime       = 10,
																maxWaitingTime  = 24 * 60 * 60,
																popmednet       = TRUE,
																trace           = FALSE) {

	params = PrepareParams.kp("logistic", dataPartnerID, numDataPartners, ac = FALSE,
	                          popmednet = popmednet, trace = trace)
	if (params$failed) {
	  cat(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}

	data = PrepareDataLogistic.DP(params, data, yname)
	params = AddToLog(params, "PrepareParamsLinear.DP", 0, 0, 0, 0)

	if (data$failed) {
		params$errorMessage = paste("Error processing data for data partner", params$dataPartnerID, "\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)
		params$errorMessage = ReadErrorMessage(params$readPathAC)
		cat(params$errorMessage, "\n")
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
		cat(possibleError$message, "\n")
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = PrepareParamsLinear.DP(params, data)
	files = "p_scaler_seed.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	params = PrepareSharesLinear.DP(params, data)
	files = c("products.rdata", "halfshare.rdata", "colstats.rdata")
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	possibleError = ReceivedError.kp(params, from = "AC")
	if (possibleError$error) {
	  params$errorMessage = possibleError$message
		cat(possibleError$message, "\n")
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
		return(params$stats)
	}

	params = UpdateParamsLogistic.DP(params)

	data = UpdateDataLogistic.DP(params, data)
	params = AddToLog(params, "UpdateDataLogistic.DP", 0, 0, 0, 0)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)

		params = ComputeSbetaLogistic.DP(params, data)
		files = "Sbeta.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = ComputeStWSLogistic.DP(params, data)
		files = "stwsshare.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		possibleError = ReceivedError.kp(params, from = "AC")
		if (possibleError$error) {
		  params$errorMessage = possibleError$message
			cat(possibleError$message, "\n")
			params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
			return(params$stats)
		}

		params = UpdateBetaLogistic.DP(params)
		files = "u_converge.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)


		params = GetConvergeStatusLogistic.DP(params)

		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = ComputeSbetaLogistic.DP(params, data)

	params = SendFinalBetasLogistic.DP(params)

	files = c("sbeta.rdata", "finalbetas.rdata")
	params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	if (dataPartnerID == 1) {
		params = ComputeResultsLogistic.DP(params, data)
		files = "logisticstats.rdata"
		params = SendPauseContinue.kp(params, filesAC = files, from = "AC",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
	}

	params = GetResultsLogistic.DP(params, data)
	params = SendPauseQuit.kp(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}


AnalysisCenterKLogistic = function(numDataPartners = NULL,
																	 monitorFolder   = NULL,
																	 msreqid         = "v_default_0_000",
																	 cutoff          = 1E-8,
																	 maxIterations   = 25,
																	 sleepTime       = 10,
																	 maxWaitingTime  = 24 * 60 * 60,
																	 popmednet       = TRUE,
																	 trace           = FALSE) {

	filesList = rep(list(list()), numDataPartners)

	params = PrepareParams.kp("logistic", 0, numDataPartners, msreqid, cutoff, maxIterations, ac = TRUE,
	                          popmednet = popmednet, trace = trace)
	if (params$failed) {
	  cat(params$errorMessage)
	  return(invisible(NULL))
	}
	params = InitializeLog.kp(params)
	params = InitializeStamps.kp(params)
	params = InitializeTrackingTable.kp(params)
	Header(params)

	params   = PrepareFolder.ACDP(params, monitorFolder)

	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}

	params = PauseContinue.kp(params, from = "DP", maxWaitingTime = maxWaitingTime)

	possibleError = ReceivedError.kp(params, from = "DP")
	if (possibleError$error) {
		params$errorMessage = possibleError$message
		cat(possibleError$message, "\n")
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
		cat(params$errorMessage, "\n")
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	files = "empty.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetProductsLogistic.AC(params)


	params = CheckColinearityLogistic.AC(params)

	if (params$failed) {
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		cat(params$errorMessage, "\n")
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.kp(params)
		return(params$stats)
	}

	params = ComputeInitialBetasLogistic.AC(params)

	for (id in 1:numDataPartners) {
		filesList[[id]] = c(paste0("u_beta_", id, ".rdata"), "indicies.rdata")
	}

	params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = ComputeWeightsLogistic.AC(params)
		files = "pi.rdata"
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ComputeStWSLogistic.AC(params)

		if (params$failed) {
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = "errorMessage.rdata"
			cat(params$errorMessage, "\n")
			params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
			params = SendPauseQuit.kp(params, sleepTime = sleepTime, job_failed = TRUE)
			SummarizeLog.kp(params)
			return(params$stats)
		}

		for (id in 1:numDataPartners) {
			filesList[[id]] = paste0("id", id, ".rdata")
		}
		params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ComputeConvergeStatusLogistic.AC(params)
		files = "u_converge.rdata"
		params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = ComputeFinalSBetaLogistic.AC(params)
	filesList = rep(list(list()), numDataPartners)
	filesList[[1]] = "sbeta.rdata"
	params = SendPauseContinue.kp(params, filesDP = filesList, from = "DP1",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = ComputeResultsLogistic.AC(params)
	files = "stats.rdata"
	params = SendPauseContinue.kp(params, filesDP = files, from = "DP",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
	params = SendPauseQuit.kp(params, sleepTime = sleepTime)
	SummarizeLog.kp(params)
	return(params$stats)
}
