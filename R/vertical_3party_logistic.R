################### DISTRIBUTED LOGISTIC REGRESSION FUNCTIONS ##################

PrepareDataLogistic.A3 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.A3\n\n")
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
      colnames(data) = c(yname, paste0("dp1:", 1:(ncol(data) - 1)))
    }
  } else {
    if (is.null(yname)) {
      yname = colnames(data)[1]
    }
  }
  if ("(Intercept)" %in% colnames(data)) {
    cat("Error: \"(Intercept)\" is not a valid covariate name.  Please change it.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (max(table(colnames(data))) > 1) {
    cat("Error: duplicate column names found.\n\n")
    workdata$failed = TRUE
    return(workdata)
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
  workdata$Y = matrix(data[, responseColIndex], ncol = 1)
  if (sum(!(workdata$Y %in% c(0, 1))) > 0) {
    cat("Error: response variable is not binary.  It should only be 0's and 1's.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (class(data) == "matrix") {
    workdata$X = cbind(1, data[, -responseColIndex, drop = FALSE])
  } else if (ncol(data) > 1) { # class(data) == "data.frame"
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
    workdata$X = cbind(1, as.matrix(data[index]))
  } else {
    workdata$X = matrix(1, nrow = nrow(data), ncol = 1)
  }
  colnames(workdata$X)[1] = "(Intercept)"

  workdata$yty    = t(workdata$Y) %*% workdata$Y
  colnames(workdata$Y) = yname

  if (ncol(workdata$X) == 1) { # Only the constant column 1
    workdata$means = 1         # Should I set this to 0 since I don't actually transform it?
    workdata$sd    = 1         # Should be 0, but this makes future computations easier
  } else {
    workdata$means = apply(workdata$X, 2, mean)
    workdata$sd    = apply(workdata$X, 2, sd)
    workdata$sd       = sapply(workdata$sd, function(x) { if (x > 0) x else 1})
    for (i in 2:ncol(workdata$X)) {
      workdata$X[, i] = matrix((workdata$X[, i] - workdata$means[i]) / workdata$sd[i], ncol = 1)
    }
  }

  return(workdata)
}

PrepareDataLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLogistic.B3\n\n")
  workdata = list()
  workdata$failed = FALSE
  if (class(data) %in% c("character", "double", "integer", "logical",
                         "numeric", "single", "factor")) {
    data = as.data.frame(data)
  }
  if (class(data) != "data.frame" & class(data) != "matrix") {
    cat("Error: data is not a matrix or a data frame.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (ncol(data) < 1 | nrow(data) < 1) {
    cat("Error: data is empty.\n\n")
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
  if (is.null(colnames(data))) {
    cat("Warning: variables are not named.  Assigning labels.\n\n")
    colnames(data) = paste0("dp2:", 1:ncol(data))
  }
  if (max(table(colnames(data))) > 1) {
    cat("Error: duplicate variable names found.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (class(data) == "matrix") {
    if (class(data[1, 1]) != "numeric" && class(data[1, 1]) != "integer") {
      data = as.data.frame(data)
    } else {
      workdata$X = data
    }
  }
  if (class(data) == "data.frame") {
    varNames = c()
    for (name in colnames(data)) {
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
              return(invisible(NULL))
            }
            data[[newName]] = ifelse(data[[name]] == lname, 1, 0)
            varNames = c(varNames, newName)
          }
          data[[name]] = NULL
        }
      }
    }
    index = match(varNames, colnames(data))
    workdata$X = as.matrix(data[index])
  }

  cnames = colnames(workdata$X)

  if (ncol(workdata$X) == 1) {
    workdata$means = mean(workdata$X)
    workdata$sd    = sd(workdata$X)
    if (workdata$sd == 0) workdata$sd = 1
    workdata$X = matrix((workdata$X - workdata$means) / workdata$sd, ncol = 1)
  } else {
    workdata$means = apply(workdata$X, 2, mean)
    workdata$sd    = apply(workdata$X, 2, sd)
    workdata$sd    = sapply(workdata$sd, function(x) { if (x > 0) x else 1})
    for (i in 1:ncol(workdata$X)) {
      workdata$X[, i] = matrix((workdata$X[, i] - workdata$means[i]) / workdata$sd[i], ncol = 1)
    }
  }

  colnames(workdata$X) = cnames

  return(workdata)
}


CheckColinearityLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityLogistic.T3\n\n")
  xtx = params$xtx
	xty = params$xty

	nrow = nrow(xtx)
	indicies = c(1)
	for (i in 2:nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}

	xtx = xtx[indicies, indicies, drop = FALSE]
	xty = xty[indicies, drop = FALSE]

	Anames               = params$colnamesA
	Bnames               = params$colnamesB
	Aindex               = which(indicies <= length(Anames))
	params$IndiciesKeep  = indicies
	params$AIndiciesKeep = indicies[Aindex]
	params$BIndiciesKeep = indicies[-Aindex] - length(Anames)
	AnamesKeep           = Anames[params$AIndiciesKeep]
	BnamesKeep           = Bnames[params$BIndiciesKeep]
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
	params$meansA        = params$meansA[params$AIndiciesKeep]
	params$meansB        = params$meansB[params$BIndiciesKeep]
	params$sdA           = params$sdA[params$AIndiciesKeep]
	params$sdB           = params$sdB[params$BIndiciesKeep]
	params$xtx           = xtx
	params$xty           = xty

	Aindicies = params$AIndiciesKeep
	Bindicies = params$BIndiciesKeep

	writeTime = proc.time()[3]
	save(Aindicies, file = file.path(params$writePath, "Aindicies.rdata"))
	save(Bindicies, file = file.path(params$writePath, "Bindicies.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("Aindicies.rdata",
																													"Bindicies.rdata"))))
	writeTime = proc.time()[3] - writeTime

	if (params$p2 == 0) {
		params$failed = TRUE
		params$errorMessage = "All of party B's covariates are either linear or are colienar with Party A's covariates."
	}
	params = AddToLog(params, "CheckColinearityLogistic.T3", 0, 0, writeTime, writeSize)
	return(params)
}


ComputeInitialBetasLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasLogistic.T3\n\n")
  # de-standardize xty
	p1     = params$p1
	p2     = params$p2
	xty    = params$xty
	xtx    = params$xtx

	betas = 4 * solve(xtx) %*% xty
	betas[1] = betas[1] - 2

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

	params$algIterationCounter      = 1
	params$deltabeta = Inf
	params$converged = FALSE
	converged = FALSE
	maxIterExceeded = FALSE

	writeTime = proc.time()[3]
	save(Abetas, file = file.path(params$writePath, "betasA.rdata"))
	save(p2, Axty,   file = file.path(params$writePath, "Axty.rdata"))
	save(Bbetas, file = file.path(params$writePath, "betasB.rdata"))
	save(Bxty,   file = file.path(params$writePath, "Bxty.rdata"))
	save(converged, maxIterExceeded,
			 file = file.path(params$writePath, "converged.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("betasA.rdata",
																													"betasB.rdata",
																													"Axty.rdata",
																													"Bxty.rdata",
																													"converged.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeInitialBetasLogistic.T3", 0, 0, writeTime, writeSize)

	return(params)}


UpdateParamsLogistic.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Aindicies.rdata"))
	load(file.path(params$readPath[["T"]], "Axty.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("Aindicies.rdata",
															                                 "Axty.rdata"))))

	readTime = proc.time()[3] - readTime
	params$colnamesA.old = params$colnamesA
	params$colnamesA     = params$colnamesA.old[Aindicies]
	params$p.old         = params$p
	params$p             = length(Aindicies)
	params$p2            = p2
	params$AIndiciesKeep = Aindicies
	params$means         = params$means[Aindicies]
	params$sd            = params$sd[Aindicies]
	params$Axty          = Axty
	params = AddToLog(params, "UpdateParamsLogistic.A3, UpdateDataLogistic.A3", readTime, readSize, 0, 0)
	return(params)
}


UpdateParamsLogistic.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsLogistic.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Bindicies.rdata"))
	load(file.path(params$readPath[["T"]], "Bxty.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("Bindicies.rdata",
																												"Bxty.rdata"))))

	readTime = proc.time()[3] - readTime
	params$colnamesB.old = params$colnamesB
	params$colnamesB     = params$colnamesB.old[Bindicies]
	params$p.old         = params$p
	params$p             = length(Bindicies)
	params$BIndiciesKeep = Bindicies
	params$means         = params$means[Bindicies]
	params$sd            = params$sd[Bindicies]
	params$Bxty          = Bxty
	params = AddToLog(params, "UpdateParamsLogistic.B3, UpdateDataLogistic.B3", readTime, readSize, 0, 0)
	return(params)
}


UpdateDataLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.A3\n\n")
  data$X = as.matrix(data$X[, params$AIndiciesKeep, drop = FALSE])
	return(data)
}


UpdateDataLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataLogistic.B3\n\n")
  data$X = as.matrix(data$X[, params$BIndiciesKeep, drop = FALSE])
	return(data)
}


GetBetaALogistic.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "converged.rdata"))
	load(file.path(params$readPath[["T"]], "betasA.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("converged.rdata",
																															 "betasA.rdata"))))
	readTime = proc.time()[3] - readTime
	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	params$betas = Abetas
	params = AddToLog(params, "GetBetaALogistic.A3", readTime, readSize, 0, 0)
	return(params)
}


GetBetaBLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaLogistic.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "converged.rdata"))
	load(file.path(params$readPath[["T"]], "betasB.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("converged.rdata",
																															 "betasB.rdata"))))
	readTime = proc.time()[3] - readTime
	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	params$betas = Bbetas
	params = AddToLog(params, "GetBetaBLogistic.B3", readTime, readSize, 0, 0)
	return(params)
}


GetXAbetaALogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaALogistic.A3\n\n")
  XAbeta = data$X %*% params$betas

	writeTime = proc.time()[3]
	save(XAbeta, file = file.path(params$writePath, "xabeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "xabeta.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "GetXAbetaALogistic.A3", 0, 0, writeTime, writeSize)
	return(params)
}


GetXBbetaBLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBetaBLogistic.B3\n\n")
  XBbeta = data$X %*% params$betas

	writeTime = proc.time()[3]
	save(XBbeta, file = file.path(params$writePath, "xbbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "xbbeta.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "GetXBbetaBLogistic.B3", 0, 0, writeTime, writeSize)
	return(params)
}


GetWeightsLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetWeightsLogistic.T3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "xabeta.rdata"))  # Load XbetaB
	load(file.path(params$readPath[["B"]], "xbbeta.rdata"))  # Load XbetaB
	readSize = file.size(file.path(params$readPath[["A"]], "xabeta.rdata")) +
		file.size(file.path(params$readPath[["B"]], "xbbeta.rdata"))
	readTime = proc.time()[3] - readTime

	Xbeta = XAbeta + XBbeta
	pi_ = (1 + exp(-Xbeta))^(-1)
	params$pi_ = pi_

	writeTime = proc.time()[3]
	save(pi_, file = file.path(params$writePath, "pi.rdata"))
	writeSize = file.size(file.path(params$writePath, "pi.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetWeightsLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetRVLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetRVLogistic.B3\n\n")
  writeTime = 0
	writeSize = 0
	readTime  = proc.time()[3]
	load(file.path(params$readPath[["T"]], "pi.rdata"))
	readSize  = file.size(file.path(params$readPath[["T"]], "pi.rdata"))
	readTime  = proc.time()[3] - readTime

	params$pi_ = pi_
	W = pi_ * (1 - params$pi_)
	XBTWXB = 0
	pbar = MakeProgressBar1(params$blocks$numBlocks, "R(I-Z*Z')W*XB")
	containerCt.RZ = 0
	containerCt.RV = 0
	for (i in 1:params$blocks$numBlocks) {
		if (i %in% params$container$filebreak.RZ) {
			containerCt.RZ = containerCt.RZ + 1
			filename1 = paste0("crz_", containerCt.RZ, ".rdata")
			toRead = file(file.path(params$readPath[["T"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.RV) {
			containerCt.RV = containerCt.RV + 1
			filename2 = paste0("crv_", containerCt.RV, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}
		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1

		Xblock  = data$X[strt:stp, , drop = FALSE]
		Wblock  = W[strt:stp]
		WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

		readTime = readTime - proc.time()[3]
		RZ = matrix(readBin(con = toRead, what = numeric(), n = n * n,
		    								endian = "little"), nrow = n, ncol = n)
		readTime = readTime + proc.time()[3]

		RV = RZ %*% WXblock

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(RV), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]

		XBTWXB = XBTWXB + t(Xblock) %*% WXblock

		if ((i + 1) %in% params$container$filebreak.RZ || i == params$blocks$numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["T"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.RV || i == params$blocks$numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename2))
		}
		pbar = MakeProgressBar2(i, pbar)
	}

	writeTime = writeTime - proc.time()[3]
	save(XBTWXB, file = file.path(params$writePath, "xbtwxb.rdata"))
	writeSize = writeSize + sum(file.size(c(file.path(params$writePath, "xbtwxb.rdata"))))
	writeTime = writeTime + proc.time()[3]

	params = AddToLog(params, "GetRVLogistic.B3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ProcessVLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessVLogistic.T3\n\n")
  writeTime = 0
	writeSize = 0
	p2 = params$p2
	readTime = proc.time()[3]
	load(file.path(params$readPath[["B"]], "xbtwxb.rdata"))
	readSize = file.size(file.path(params$readPath[["B"]], "xbtwxb.rdata"))
	readTime = proc.time()[3] - readTime

	params$xbtwxb = XBTWXB

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "(I-Z*Z')W*XB*R")

	containerCt.RV = 0
	containerCt.VR = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.RV) {
			containerCt.RV = containerCt.RV + 1
			filename2 = paste0("crv_", containerCt.RV, ".rdata")
			toRead2 = file(file.path(params$readPath[["B"]], filename2), "rb")
		}
		if (i %in% params$container$filebreak.VR) {
			containerCt.VR = containerCt.VR + 1
			filename3 = paste0("cvr_", containerCt.VR, ".rdata")
			toWrite3 = file(file.path(params$writePath, filename3), "wb")
		}

		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1

		filename1 = paste0("r1_", i, ".rdata")
		filename4 = paste0("r3_", i, ".rdata")

		readTime = readTime - proc.time()[3]
		toRead1 = file(file.path(params$dplocalPath, filename1), "rb")
		R2 = matrix(readBin(con = toRead1, what = numeric(), n = n * n,
												endian = "little"), nrow = n, ncol = n)
		readSize = readSize + file.size(file.path(params$dplocalPath, filename1))
		close(toRead1)
		RV = matrix(readBin(con = toRead2, what = numeric(), n = n * p2,
												endian = "little"), nrow = n, ncol = p2)
		readTime = readTime + proc.time()[3]

		V = t(R2) %*% RV
		R3 = RandomOrthonomalMatrix(p2)
		VR = V %*% R3

		writeTime = writeTime - proc.time()[3]
		toWrite4 = file(file.path(params$dplocalPath, filename4), "wb")
		writeBin(as.vector(R3), con = toWrite4, endian = "little")
		close(toWrite4)
		writeSize = writeSize + file.size(file.path(params$dplocalPath, filename4))
		writeBin(as.vector(VR), con = toWrite3, endian = "little")
		writeTime = writeTime + proc.time()[3]
		if ((i + 1) %in% params$container$filebreak.RV || i == numBlocks) {
			close(toRead2)
			readSize = readSize + file.size(file.path(params$dplocalPath, filename1))
		}
		if ((i + 1) %in% params$container$filebreak.VR || i == numBlocks) {
			close(toWrite3)
			writeSize = writeSize + file.size(file.path(params$writePath, filename3))
		}

		pbar = MakeProgressBar2(i, pbar)
	}
	params = AddToLog(params, "ProcessVLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetXRLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXRLogistic.A3\n\n")
  p2 = params$p2
	writeTime = 0
	writeSize = 0
	readTime  = proc.time()[3]
	load(file.path(params$readPath[["T"]], "pi.rdata"))
	readSize  = file.size(file.path(params$readPath[["T"]], "pi.rdata"))
	readTime  = proc.time()[3] - readTime

	params$pi_ = pi_
	W = pi_ * (1 - params$pi_)
	XATWXA = 0

	pbar = MakeProgressBar1(params$blocks$numBlocks, "XA'(I-Z*Z')W*XB*R")

	containerCt.VR = 0
	containerCt.XR = 0
	for (i in 1:params$blocks$numBlocks) {
		if (i %in% params$container$filebreak.RV) {
			containerCt.VR = containerCt.VR + 1
			filename1 = paste0("cvr_", containerCt.VR, ".rdata")
			toRead = file(file.path(params$readPath[["T"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.XR) {
			containerCt.XR = containerCt.XR + 1
			filename2 = paste0("cxr_", containerCt.XR, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}
		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1

		Xblock  = data$X[strt:stp, , drop = FALSE]
		Wblock  = W[strt:stp]
		WXblock = MultiplyDiagonalWTimesX(Wblock, Xblock)

		readTime = readTime - proc.time()[3]
		VR = matrix(readBin(con = toRead, what = numeric(), n = n * p2,
												endian = "little"), nrow = n, ncol = p2)
		readTime = readTime + proc.time()[3]

		XR = t(Xblock) %*% VR

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(XR), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]

		XATWXA = XATWXA + t(Xblock) %*% WXblock

		if ((i + 1) %in% params$container$filebreak.VR || i == params$blocks$numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["T"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename2))
		}
		pbar = MakeProgressBar2(i, pbar)
	}

	writeTime = writeTime - proc.time()[3]
	save(XATWXA, file = file.path(params$writePath, "xatwxa.rdata"))
	writeSize = writeSize + sum(file.size(c(file.path(params$writePath, "xatwxa.rdata"))))
	writeTime = writeTime + proc.time()[3]

	params = AddToLog(params, "GetXRLogistic.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ProcessXtWXLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessXtWXLogistic.T3\n\n")
  p1 = params$p1
	p2 = params$p2

	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "xatwxa.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "xatwxa.rdata"))
	readTime = proc.time()[3] - readTime

	params$xatwxa = XATWXA

	pbar = MakeProgressBar1(params$blocks$numBlocks, "X'W*X")
	containerCt.XR = 0
	XATWXB = 0

	for (i in 1:params$blocks$numBlocks) {
		if (i %in% params$container$filebreak.XR) {
			containerCt.XR = containerCt.XR + 1
			filename1 = paste0("cxr_", containerCt.XR, ".rdata")
			toRead = file(file.path(params$readPath[["A"]], filename1), "rb")
		}

		filename2 = paste0("r3_", i, ".rdata")
		readTime = readTime - proc.time()[3]
		toRead1 = file(file.path(params$dplocalPath, filename2), "rb")
		R = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
		 									 endian = "little"), nrow = p2, ncol = p2)
		close(toRead1)
		XR = matrix(readBin(con = toRead, what = numeric(), n = p1 * p2,
											  endian = "little"), nrow = p1, ncol = p2)
		readSize = readSize + file.size(file.path(params$dplocalPath, filename2))
		readTime = readTime + proc.time()[3]

		XATWXB = XATWXB + XR %*% t(R)

		if ((i + 1) %in% params$container$filebreak.XR || i == params$blocks$numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["A"]], filename1))
		}
		pbar = MakeProgressBar2(i, pbar)
	}

	xtwx = rbind(cbind(params$xatwxa, XATWXB), cbind(t(XATWXB), params$xbtwxb))
	params$xtwx = xtwx

	II = NULL
	tryCatch({II = solve(xtwx)},
					 error = function(err) { II = NULL }
	)
	if (is.null(II)) {
		params$failed = TRUE
		params$singularMatrix = TRUE
		params$errorMessage =
			paste0("ERROR: The matrix t(X)*W*X is not invertible.\n",
				 		 "       This may be due to one of two possible problems.\n",
				     "       1. Poor random initilization of the security vector.\n",
				     "       2. Near multicolinearity in the data\n",
				     "SOLUTIONS: \n",
				     "       1. Rerun the data analaysis.\n",
				     "       2. If the problem persists, check the variables for\n",
				     "          duplicates for both parties and / or reduce the\n",
				     "          number of variables used. Once this is done,\n",
				     "          rerun the data analysis.\n\n")
		return(params)
	}
	params$II = II
	IIA = II[, 1:p1, drop = FALSE]
	IIB = II[, (p1 + 1):(p1 + p2), drop = FALSE]
	writeTime = proc.time()[3]
	save(IIA, file = file.path(params$writePath, "IIA.rdata"))
	save(IIB, file = file.path(params$writePath, "IIB.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("IIA.rdata", "IIB.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ProcessXtWXLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateBetaLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "IIA.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "IIA.rdata"))
	readTime = proc.time()[3] - readTime

	IA = params$Axty - t(data$X) %*% params$pi_
	AI = IIA %*% IA


	writeTime = proc.time()[3]
	save(AI, file = file.path(params$writePath, "AI.rdata"))
	writeSize = file.size(file.path(params$writePath, "AI.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "UpdateBetaLogistic.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateBetaLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "IIB.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "IIB.rdata"))
	readTime = proc.time()[3] - readTime

	IB = params$Bxty - t(data$X) %*% params$pi_
	BI = IIB %*% IB

	writeTime = proc.time()[3]
	save(BI, file = file.path(params$writePath, "BI.rdata"))
	writeSize = file.size(file.path(params$writePath, "BI.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "UpdateBetaLogistic.B3", readTime, readSize, writeTime, writeSize)
	return(params)
}


UpdateBetaLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateBetaLogistic.T3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "AI.rdata"))
	load(file.path(params$readPath[["B"]], "BI.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "AI.rdata")) +
		file.size(file.path(params$readPath[["B"]], "BI.rdata"))
	readTime = proc.time()[3] - readTime

	delta = AI + BI
	betas = params$betas + delta
	params$betas = betas
	converged = all(abs(delta) / (abs(betas) + .1) < params$cutoff)
	maxIterExceeded = (params$algIterationCounter >= params$maxIterations) && !converged

	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	Abetas = betas[1:params$p1]
	Bbetas = betas[(params$p1 + 1):(params$p1 + params$p2)]

	writeTime = proc.time()[3]
	save(converged, maxIterExceeded,
			 file = file.path(params$writePath, "converged.rdata"))
	save(Abetas, file = file.path(params$writePath, "betasA.rdata"))
	save(Bbetas, file = file.path(params$writePath, "betasB.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("betasA.rdata",
																											"betasB.rdata",
																											"converged.rdata"))))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "UpdateBetaLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetFinalBetaLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.A3\n\n")
  betas = params$betas / params$sd
	offsetA = sum(betas[-1] * params$means[-1])
	AFinalFitted = t(params$sd * t(data$X) + params$means) %*% betas -
		             t(params$sd[1] * t(data$X[, 1]) + params$means[1]) %*% betas[1]
	writeTime = proc.time()[3]
	save(offsetA, AFinalFitted, file = file.path(params$writePath, "Afinalfitted.rdata"))
	writeSize = file.size(file.path(params$writePath, "Afinalfitted.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetFinalBetaLogistic.A3", 0, 0, writeTime, writeSize)
	return(params)
}


GetFinalBetaLogistic.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "getFinalBetaLogistic.B3\n\n")
  betas = params$betas / params$sd
	offsetB = sum(betas * params$means)
	BFinalFitted = t(params$sd * t(data$X) + params$means) %*% betas

	writeTime = proc.time()[3]
	save(offsetB, BFinalFitted, file = file.path(params$writePath, "Bfinalfitted.rdata"))
	writeSize = file.size(file.path(params$writePath, "Bfinalfitted.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetFinalBetaLogistic.B3", 0, 0, writeTime, writeSize)
	return(params)
}


GetFinalFittedLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetFinalFittedLogistic.T3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "Afinalfitted.rdata"))
	load(file.path(params$readPath[["B"]], "Bfinalfitted.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "Afinalfitted.rdata")) +
		file.size(file.path(params$readPath[["B"]], "Bfinalfitted.rdata"))
	readTime = proc.time()[3] - readTime

	betas = params$betas / c(params$sdA, params$sdB)
	betas[1] = betas[1] - offsetA - offsetB

	finalFitted = AFinalFitted + BFinalFitted + betas[1]

	params$betas = betas
	params$finalFitted = finalFitted

	writeTime = proc.time()[3]
	save(finalFitted, file = file.path(params$writePath, "finalFitted.rdata"))
	writeSize = file.size(file.path(params$writePath, "finalFitted.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "GetFinalBetaLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}

ComputeResultsLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "finalFitted.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "finalFitted.rdata"))
	readTime = proc.time()[3] - readTime

	n = params$n
	ct      = sum(data$Y)
	params$FinalFitted = finalFitted
	resdev  = -2 * (sum(data$Y * finalFitted) - sum(log(1 + exp(finalFitted))))
	nulldev = -2 * (ct * log(ct / n) + (n - ct) * log(1 - ct / n))
	hoslem  = HoslemInternal(params, data)
	ROC     = RocInternal(params, data)

	writeTime = proc.time()[3]
	save(resdev, nulldev, hoslem, ROC, file = file.path(params$writePath, "logisticstats.rdata"))
	writeSize = file.size(file.path(params$writePath, "logisticstats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeResultsLogistic.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeResultsLogistic.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLogistic.T3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "logisticstats.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "logisticstats.rdata"))
	readTime = proc.time()[3] - readTime

	stats = params$stats
	stats$failed         = FALSE
	stats$converged      = params$converged

	n      = params$n
	p1     = params$p1
	p2     = params$p2
	sdA    = params$sdA
	sdB    = params$sdB
	meansA = params$meansA
	meansB = params$meansB
	Anames = params$colnamesA.old
	Bnames = params$colnamesB.old
	p1.old = params$p1.old
	p2.old = params$p2.old
	p.old  = params$p.old
	indicies = params$IndiciesKeep

	# If xtwx were singular, it would have been caught in GetII.A2(), so we may
	# assume that xtwx is NOT singular and so we do not have to do a check.
	cov1 = solve(params$xtwx)
	secoef = sqrt(diag(cov1)) / c(sdA, sdB)
	tmp = matrix(c(1, (-meansA / sdA)[-1], -meansB / sdB), ncol = 1)
	secoef[1] = sqrt(t(tmp) %*% cov1 %*% tmp)

	stats$party = c(rep("dp1", p1.old), rep("dp2", p2.old))
	stats$coefficients = rep(NA, p.old)
	stats$secoef = rep(NA, p.old)
	stats$tvals  = rep(NA, p.old)
	stats$pvals  = rep(NA, p.old)
	stats$n  = n
	stats$nulldev = nulldev
	stats$resdev = resdev
	stats$aic = resdev + 2 * (p1 + p2)
	stats$bic = resdev + (p1 + p2) * log(n)
	stats$nulldev_df = n - 1
	stats$resdev_df = n - (p1 + p2)
	stats$coefficients[indicies] = params$betas
	stats$secoef[indicies] = secoef
	tvals = params$betas / secoef
	pvals = 2 * pnorm(abs(tvals), lower.tail = FALSE)
	stats$tvals[indicies] = tvals
	stats$pvals[indicies] = pvals
	stats$hoslem  = hoslem
	stats$ROC     = ROC
	stats$iter    = params$algIterationCounter - 1
	names.old = c(Anames, Bnames)
	names(stats$coefficients) = names.old
	names(stats$party) = names.old
	names(stats$secoef) = names.old
	names(stats$tvals) = names.old
	names(stats$pvals) = names.old

	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params$stats      = stats

	params = AddToLog(params, "ComputeResultsLogistic.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetResultsLogistic.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	stats$Y           = data$Y # For Hoslem and ROC
	stats$FinalFitted = params$FinalFitted
	params$stats      = stats
	params = AddToLog(params, "GetResultsLogistic.A3", readTime, readSize, 0, 0)
	return(params)
}


GetResultsLogistic.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLogistic.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats
	params = AddToLog(params, "GetResultsLogistic.B3", readTime, readSize, 0, 0)
	return(params)
}


############################### PARENT FUNCTIONS ###############################


PartyAProcess3Logistic = function(data,
																	yname          = NULL,
																	monitorFolder  = NULL,
																	sleepTime      = 10,
																	maxWaitingTime = 24 * 60 * 60,
																	popmednet      = TRUE,
																	trace          = FALSE) {

	params = PrepareParams.3p("logistic", "A",
	                          popmednet = popmednet, trace = trace)
	params = InitializeLog.3p(params)
	params = InitializeStamps.3p(params)
	params = InitializeTrackingTable.3p(params)
	Header(params)

	params   = PrepareFolderLinear.A3(params, monitorFolder)
	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}
	data = PrepareDataLogistic.A3(params, data, yname)
	params = AddToLog(params, "PrepareDataLogistic.A3", 0, 0, 0, 0)

	if (data$failed) {
		message = "Error in processing the data for Party A."
		MakeErrorMessage(params$writePath, message)
		files = c("errorMessage.rdata")
		params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = PrepareParamsLinear.A3(params, data)
	files = "pa.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params$algIterationCounter = 1
	params = PrepareBlocksLinear.A3(params)

	params = GetZLinear.A3(params, data)
	files = SeqZW("cz_", length(params$container$filebreak.Z))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetWRLinear.A3(params, data)
	files = c("xatxa.rdata", SeqZW("cpr_", length(params$container$filebreak.PR)))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = UpdateParamsLogistic.A3(params)
	data = UpdateDataLogistic.A3(params, data)
	params = GetBetaALogistic.A3(params)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = GetXAbetaALogistic.A3(params, data)
		files = c("xabeta.rdata")
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = GetXRLogistic.A3(params, data)
		files = c("xatwxa.rdata", SeqZW("cxr_", length(params$container$filebreak.XR)))
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
			cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
			params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
			return(params$stats)
		}

		params = UpdateBetaLogistic.A3(params, data)
		files = c("ai.rdata")
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = GetBetaALogistic.A3(params)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = GetFinalBetaLogistic.A3(params, data)
	files = "Afinalfitted.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	params = ComputeResultsLogistic.A3(params, data)
	files = c("logisticstats.rdata")
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetResultsLogistic.A3(params, data)

	params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}


PartyBProcess3Logistic = function(data,
																	monitorFolder  = NULL,
																	sleepTime      = 10,
																	maxWaitingTime = 24 * 60 * 60,
																	popmednet      = TRUE,
																	trace          = FALSE) {
	params = PrepareParams.3p("logistic", "B",
	                          popmednet = popmednet, trace = trace)
	params = InitializeLog.3p(params)
	params = InitializeStamps.3p(params)
	params = InitializeTrackingTable.3p(params)

	Header(params)
	params   = PrepareFolderLinear.B3(params, monitorFolder)
	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}

	data = PrepareDataLogistic.B3(params, data)
	params = AddToLog(params, "PrepareDataLogistic.B3", 0, 0, 0, 0)

	if (data$failed) {
		message = "Error in processing the data for Party B."
		MakeErrorMessage(params$writePath, message)
		files = c("errorMessage.rdata")
		params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = PrepareParamsLinear.B3(params, data)
	files = "pb.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params$algIterationCounter = 1
	params = PrepareBlocksLinear.B3(params)

	params = GetRWLinear.B3(params, data)
	files = c("xbtxb.rdata", SeqZW("crw_", length(params$container$filebreak.RW)))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = UpdateParamsLogistic.B3(params)
	data = UpdateDataLogistic.B3(params, data)
	params = GetBetaBLogistic.B3(params, data)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = GetXBbetaBLogistic.B3(params, data)
		files = c("xbbeta.rdata")
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = GetRVLogistic.B3(params, data)
		files = c("xbtwxb.rdata", SeqZW("crv_", length(params$container$filebreak.RV)))
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
			cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
			params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
			return(params$stats)
		}

		params = UpdateBetaLogistic.B3(params, data)
		files = c("bi.rdata")
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = GetBetaBLogistic.B3(params)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = GetFinalBetaLogistic.B3(params, data)
	files = "Bfinalfitted.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	params = GetResultsLogistic.B3(params)

	params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}


PartyTProcess3Logistic = function(monitorFolder         = NULL,
																	msreqid               = "v_default_0_000",
																	blocksize             = 500,
																	cutoff                = 1e-8,
																	maxIterations         = 25,
																	sleepTime             = 10,
																	maxWaitingTime        = 24 * 60 * 60,
																	popmednet             = TRUE,
																	trace                 = FALSE) {
	params = PrepareParams.3p("logistic", "T", msreqid = msreqid,
	                          popmednet = popmednet, trace = trace)
	params = InitializeLog.3p(params)
	params = InitializeStamps.3p(params)
	params = InitializeTrackingTable.3p(params)

	Header(params)
	params   = PrepareFolderLinear.T3(params, monitorFolder)
	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}

	params = PauseContinue.3p(params, from = c("A", "B"), maxWaitingTime = maxWaitingTime)

	if (file.exists(file.path(params$readPath[["A"]], "errorMessage.rdata")) &&
			file.exists(file.path(params$readPath[["B"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["A"]]), "\n\n")
		cat("Error:", ReadErrorMessage(params$readPath[["B"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}
	if (file.exists(file.path(params$readPath[["A"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["A"]]), "\n\n")
		file.copy(file.path(params$readPath[["A"]], "errorMessage.rdata"),
							file.path(params$writePath, "errorMessage.rdata"))
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesB = files, from = "B",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}
	if (file.exists(file.path(params$readPath[["B"]], "errorMessage.rdata"))) {
		cat("Error:", ReadErrorMessage(params$readPath[["B"]]), "\n\n")
		file.copy(file.path(params$readPath[["B"]], "errorMessage.rdata"),
							file.path(params$writePath, "errorMessage.rdata"))
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, from = "A",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	params   = PrepareParamsLinear.T3(params, cutoff, maxIterations)

	if (!params$failed) params = PrepareBlocksLinear.T3(params, blocksize)

	if (params$failed) {
		cat("Error:", params$errorMessage, "\n\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files,
															 from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	files = "blocks.rdata"
	params = SendPauseContinue.3p(params, filesA = files, from = "A",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params$algIterationCounter = 1
	params = ProcessZLinear.T3(params)
	files = c("blocks.rdata", SeqZW("crz_", length(params$container$filebreak.RZ)))
	params = SendPauseContinue.3p(params, filesB = files, from  = "B",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = ProcessWLinear.T3(params)
	files = c("p2.rdata", SeqZW("cwr_", length(params$container$filebreak.WR)))
	params = SendPauseContinue.3p(params, filesA = files, from  = "A",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetProductsLinear.T3(params)

	params = CheckColinearityLogistic.T3(params)

	if (params$failed) {
		cat("Error:", params$errorMessage, "\n\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files,
															 from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	params = ComputeInitialBetasLogistic.T3(params)
	filesA = c("Aindicies.rdata", "betasA.rdata", "Axty.rdata", "converged.rdata")
	filesB = c("Bindicies.rdata", "betasB.rdata", "Bxty.rdata", "converged.rdata")
	params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params$algIterationCounter = 1
	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		params = GetWeightsLogistic.T3(params)
		files = "pi.rdata"
		params = SendPauseContinue.3p(params, filesB = files, from  = "B",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ProcessVLogistic.T3(params)
		files = c("pi.rdata", SeqZW("cvr_", length(params$container$filebreak.RV)))
		params = SendPauseContinue.3p(params, filesA = files, from  = "A",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ProcessXtWXLogistic.T3(params)

		if (params$failed) {
			cat("Error:", params$errorMessage, "\n\n")
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = "errorMessage.rdata"
			params = SendPauseContinue.3p(params, filesA = files, filesB = files,
																 from = c("A", "B"),
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
			params = SendPauseQuit.3p(params, sleepTime = sleepTime)
			SummarizeLog.3p(params)
			return(params$stats)
		}
		filesA = c("IIA.rdata")
		filesB = c("IIB.rdata")
		params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = UpdateBetaLogistic.T3(params)
		filesA = c("betasA.rdata", "converged.rdata")
		filesB = c("betasB.rdata", "converged.rdata")
		params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB, from  = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = GetFinalFittedLogistic.T3(params)
	filesA = "finalfitted.rdata"
	params = SendPauseContinue.3p(params, filesA = filesA, from  = "A",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = ComputeResultsLogistic.T3(params)
	files = "stats.rdata"
	params = SendPauseContinue.3p(params, filesA = files, filesB = files, from  = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = SendPauseQuit.3p(params, sleepTime = sleepTime)
	SummarizeLog.3p(params)
	return(params$stats)
}
