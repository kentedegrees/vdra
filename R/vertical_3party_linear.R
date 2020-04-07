################### DISTRIBUTED LINEAR REGRESSION FUNCTIONS ###################

PrepareFolderLinear.A3 = function(params, monitorFolder = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.A3\n\n")
  if (is.null(monitorFolder)) {
		cat("monitorFolder must be specified.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
  }
	if (class(monitorFolder) != "character") {
		cat("monitorFolder directory is not valid.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
	}
  while (!dir.exists(monitorFolder)) {
    Sys.sleep(1)
  }

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "msoc")
	params$readPath      = c(file.path(monitorFolder, "inputfiles"),
													 NA,
													 file.path(monitorFolder, "msoc2"))
	names(params$readPath) = c("T", "A", "B")

	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["T"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc2")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["B"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.\n\n")
	}

	Sys.sleep(1)
	DeleteTrigger("files_done.ok", params$readPath[1])
	DeleteTrigger("files_done.ok", params$readPath[3])

	params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)

	return(params)
}


PrepareFolderLinear.B3 = function(params, monitorFolder = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.B3\n\n")
  if (is.null(monitorFolder)) {
		cat("monitorFolder must be specified.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
	}
	if (class(monitorFolder) != "character") {
		cat("monitorFolder directory is not valid.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
	}
  while (!dir.exists(monitorFolder)) {
	  Sys.sleep(1)
	}

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "msoc")
	params$readPath      = c(file.path(monitorFolder, "inputfiles"),
													 file.path(monitorFolder, "msoc1"),
													 NA)
	names(params$readPath) = c("T", "A", "B")

	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "inputfiles")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["T"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc1")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["A"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.\n\n")
	}

	Sys.sleep(1)
	DeleteTrigger("files_done.ok", params$readPath[1])
	DeleteTrigger("files_done.ok", params$readPath[2])

	params = AddToLog(params, "PrepareFolderLinear.B3", 0, 0, 0, 0)

	return(params)
}


PrepareFolderLinear.T3 = function(params, monitorFolder = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareFolderLinear.T3\n\n")
  if (is.null(monitorFolder)) {
		cat("monitorFolder must be specified.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
	}
	if (class(monitorFolder) != "character") {
		cat("monitorFolder directory is not valid.  Please use the same monitorFolder as the Datamart Client.\n")
		params$failed = TRUE
		params = AddToLog(params, "PrepareFolderLinear.A3", 0, 0, 0, 0)
		return(params)
	}
  while (!dir.exists(monitorFolder)) {
    Sys.sleep(1)
  }

	params$dplocalPath   = file.path(monitorFolder, "dplocal")
	params$rprogramsPath = file.path(monitorFolder, "rprograms")
	params$macrosPath    = file.path(monitorFolder, "macros")
	params$writePath     = file.path(monitorFolder, "inputfiles")
	params$readPath      = c(NA,
													 file.path(monitorFolder, "msoc1"),
													 file.path(monitorFolder, "msoc2"))
	names(params$readPath) = c("T", "A", "B")

	if (!CreateIOLocation(monitorFolder, "dplocal")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$dplocalPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "rprograms")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$rprogramsPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "macros")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$macrosPath, "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc1")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["A"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc2")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$readPath[["B"]], "."),
																"Check the path and restart the program.\n\n")
	}
	if (!CreateIOLocation(monitorFolder, "msoc")) {
		params$failed = TRUE
		params$errorMessage = paste(params$errorMessage,
																"Error: Could not create directory",
																paste0(params$writePath, "."),
																"Check the path and restart the program.\n\n")
	}

	params = AddToLog(params, "PrepareFolderLinear.T3", 0, 0, 0, 0)

	return(params)
}


PrepareDataLinear.A3 = function(params, data, yname = NULL) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.A3\n\n")
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
  if (class(data) == "matrix") {
    workdata$X = cbind(1, data[, -responseColIndex])
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
  workdata$meansy = mean(workdata$Y)
  workdata$sdy    = sd(workdata$Y)
  if (workdata$sdy == 0) {
    workdata$sdy = 1           # Makes future computations easier
  }
  workdata$Y      = matrix((workdata$Y - workdata$meansy) / workdata$sdy, ncol = 1)
  colnames(workdata$Y) = yname

  if (ncol(workdata$X) == 1) { # Only the constant column 1
    workdata$means = 1         # Should this be 0, since I don't actually transform it?
    workdata$sd    = 1         # Should be 0, but this makes future computations easier
  } else {
    workdata$means = apply(workdata$X, 2, mean)
    workdata$sd    = apply(workdata$X, 2, sd)
    workdata$sd    = sapply(workdata$sd, function(x) { if (x > 0) x else 1})
    for (i in 2:ncol(workdata$X)) {
      workdata$X[, i] = matrix((workdata$X[, i] - workdata$means[i]) / workdata$sd[i], ncol = 1)
    }
  }

  return(workdata)
}

PrepareDataLinear.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataLinear.B3\n\n")
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
    workdata$X = as.matrix(data[, index])
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


PrepareParamsLinear.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.A3\n\n")
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
	pa$yname    = data$yname
	pa$colnames = colnames(data$X)

	writeTime = proc.time()[3]
	save(pa, file = file.path(params$writePath, "pa.rdata"))
	writeSize = file.size(file.path(params$writePath, "pa.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "PrepareParamsLinear.A3", 0, 0, writeTime, writeSize)
	return(params)
}


PrepareParamsLinear.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.B3\n\n")
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

  writeTime = proc.time()[3]
  save(pb, file = file.path(params$writePath, "pb.rdata"))
  writeSize = file.size(file.path(params$writePath, "pb.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareParamsLinear.B3", 0, 0, writeTime, writeSize)
  return(params)
}


PrepareParamsLinear.T3 = function(params, cutoff = 1e-8, maxIterations = 25) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsLinear.T3\n\n")
  pa = NULL
	pb = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "pa.rdata"))
	load(file.path(params$readPath[["B"]], "pb.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "pa.rdata")) +
		         file.size(file.path(params$readPath[["B"]], "pb.rdata"))
	readTime = proc.time()[3] - readTime
	if (length(table(c(pa$analysis, pb$analysis, params$analysis))) > 1) {
		params$failed = TRUE
		params$errorMessage = paste("Party A specified", pa$analysis, "regression, ",
                                "Party B specified", pb$analysis, "regression, ",
																"and Party T specified", params$analysis, "regression. ")
	}
	if (pa$n != pb$n) {
		params$failed = TRUE
		params$errorMessage = paste0(params$errorMessage,
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
	params$yname         = pa$yname
	params$cutoff        = cutoff
	params$maxIterations = maxIterations

	params = AddToLog(params, "PrepareParamsLinear.T3", readTime, readSize, 0, 0)
	return(params)
}


PrepareBlocksLinear.T3 = function(params, blocksize) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.T3\n\n")
  # For now, assuming that p1 > 0 and p2 > 0
  n  = params$n
  p1 = params$p1
  p2 = params$p2


  minimumBlocksize = GetBlockSize(p1, p2)
  if (n < minimumBlocksize) {
    maxACovariates = trunc(sqrt(p2 * n) - p2 - 1)

    params$errorMessage =
    paste("The minimum secure blocksize of", minimumBlocksize,
    			"is larger than the number of observations", paste0(n, ".\n"),
    		  "Your options are:\n",
          "Increase the number of observations to at least",
           paste0(minimumBlocksize, ".\n"),
          "Decrease the number of A covariates to", maxACovariates, "or less.")
    b = n - 2 * p1 - 2
    discrim = b^2 - 4 * (p1 + 1)^2
    if (discrim >= 0) {
      minBCovariates = trunc(1 + (b - sqrt(discrim)) / 2)
      maxBCovariates = trunc((b + sqrt(discrim)) / 2)
      params$errorMessage = paste0(params$errorMessage,
      														 "\nSet the number of B covariates to be between ", minBCovariates, "and",
          paste0(maxBCovariates, "."))
    }
    params$failed = TRUE
    params = AddToLog(params, "PrepareBlocksLinear.T3", 0, 0, 0, 0)
    return(params)
  }

  if (is.null(blocksize)) {
    blocksize = minimumBlocksize
  }
  if (blocksize < minimumBlocksize) {
    cat("Warning:  Block size of", blocksize,
        "is too small. Proceeding with minimum blocksize of",
        paste0(minimumBlocksize, ".\n\n"))
    blocksize = minimumBlocksize
  } else if (n < blocksize) {
    cat("Warning: Block size of", blocksize,
        "is larger than size of data.  Proceeding with blocksize of",
        paste0(n, ".\n\n"))
  }

  params$blocks    = CreateBlocks(p1, p2, n, blocksize)
  params$container = CreateContainers(p1, p2, params$blocks)
  blocks     = params$blocks
  containers = params$container
  writeTime = proc.time()[3]
  save(blocks, containers, file = file.path(params$writePath, "blocks.rdata"))
  writeSize = file.size(file.path(params$writePath, "blocks.rdata"))
  writeTime = proc.time()[3] - writeTime
  params = AddToLog(params, "PrepareBlocksLinear.T3", 0, 0, writeTime, writeSize)
  return(params)
}


PrepareBlocksLinear.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "blocks.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "blocks.rdata"))
	readTime = proc.time()[3] - readTime

	params$blocks = blocks
	params$containers = containers
	params = AddToLog(params, "PrepareBlocksLinear.A3", readTime, readSize, 0, 0)
	return(params)
}


GetZLinear.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZLinear.A3\n\n")
  writeTime = 0
	writeSize = 0

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "Z")
	containerCt.Z = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.Z) {
			containerCt.Z = containerCt.Z + 1
			filename = paste0("cz_", containerCt.Z, ".rdata")
			toWrite = file(file.path(params$writePath, filename), "wb")
		}
		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1
		g = params$blocks$g[i]
		Z = FindOrthogonalVectors(cbind(data$Y[strt:stp, ], data$X[strt:stp, ]), g)

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(Z), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]
		if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename))
		}
		pbar = MakeProgressBar2(i, pbar)
	}
	params = AddToLog(params, "GetZLinear.A3", 0, 0, writeTime, writeSize)
	return(params)
}


ProcessZLinear.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessZLinear.T3\n\n")
  readTime = 0
	readSize = 0
	writeTime = 0
	writeSize = 0
	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "R(I-Z*Z')")
	containerCt.Z = 0
	containerCt.RZ = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.Z) {
			containerCt.Z = containerCt.Z + 1
			filename1 = paste0("cz_", containerCt.Z, ".rdata")
			toRead = file(file.path(params$readPath[["A"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.RZ) {
			containerCt.RZ = containerCt.RZ + 1
			filename2 = paste0("crz_", containerCt.RZ, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}
		filename3 = paste0("r1_", i, ".rdata")

		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1
		g = params$blocks$g[i]

		readTime = readTime - proc.time()[3]
		Z = matrix(readBin(con = toRead, what = numeric(), n = n * g,
											endian = "little"), nrow = n, ncol = g)
		readTime = readTime + proc.time()[3]
		R = RandomOrthonomalMatrix(n)
		RZ = R - (R %*% Z) %*% t(Z)

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(RZ), con = toWrite, endian = "little")
		toWrite2 = file(file.path(params$dplocalPath, filename3), "wb")
		writeBin(as.vector(R), con = toWrite2, endian = "little")
		close(toWrite2)
		writeSize = writeSize + file.size(file.path(params$dplocalPath, filename3))
		writeTime = writeTime + proc.time()[3]

		if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["A"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.RZ || i == numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename2))
		}

		pbar = MakeProgressBar2(i, pbar)
	}

	params = AddToLog(params, "ProcessZLinear.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


PrepareBlocksLinear.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareBlocksLinear.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "blocks.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "blocks.rdata"))
	readTime = proc.time()[3] - readTime

	params$blocks = blocks
	params$containers = containers
	params = AddToLog(params, "PrepareBlocksLinear.B3", readTime, readSize, 0, 0)
	return(params)
}


GetRWLinear.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetRWLinear.B3\n\n")

	readTime = 0
	readSize = 0

	numBlocks = params$blocks$numBlocks

	XBTXB = t(data$X) %*% data$X
	writeTime = proc.time()[3]
	save(XBTXB, file = file.path(params$writePath, "xbtxb.rdata"))
	writeSize = file.size(file.path(params$writePath, "xbtxb.rdata"))
	writeTime = proc.time()[3] - writeTime

	pbar = MakeProgressBar1(numBlocks, "R(I-Z*Z')XB")
	containerCt.RZ = 0
	containerCt.RW = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.RZ) {
			containerCt.RZ = containerCt.RZ + 1
			filename1 = paste0("crz_", containerCt.RZ, ".rdata")
			toRead = file(file.path(params$readPath[["T"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.RW) {
			containerCt.RW = containerCt.RW + 1
			filename2 = paste0("crw_", containerCt.RW, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}
		strt = params$blocks$starts[i]
		stp  = params$blocks$stops[i]
		n    = stp - strt + 1
		g = params$blocks$g[i]

		XB = data$X[strt:stp, , drop = FALSE]
		readTime = readTime - proc.time()[3]
		RZ = matrix(readBin(con = toRead, what = numeric(), n = n * n,
											 endian = "little"), nrow = n, ncol = n)
		readTime = readTime + proc.time()[3]

		RW = RZ %*% XB

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(RW), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]
		if ((i + 1) %in% params$container$filebreak.RZ || i == numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["T"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.RW || i == numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename2))
		}

		pbar = MakeProgressBar2(i, pbar)
	}

	params = AddToLog(params, "GetRWLinear.B3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ProcessWLinear.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessWLinear.T3\n\n")
  readTime = 0
	readSize = 0
	writeTime = 0
	writeSize = 0

	p2 = params$p2

	writeTime = proc.time()[3]
	save(p2, file = file.path(params$writePath, "p2.rdata"))
	writeSize = file.size(file.path(params$writePath, "p2.rdata"))
	writeTime = proc.time()[3] - writeTime

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "(I-Z*Z')XB*R")

	containerCt.RW = 0
	containerCt.WR = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.RW) {
			containerCt.RW = containerCt.RW + 1
			filename2 = paste0("crw_", containerCt.RW, ".rdata")
			toRead2 = file(file.path(params$readPath[["B"]], filename2), "rb")
		}
		if (i %in% params$container$filebreak.WR) {
			containerCt.WR = containerCt.WR + 1
			filename3 = paste0("cwr_", containerCt.WR, ".rdata")
			toWrite3 = file(file.path(params$writePath, filename3), "wb")
		}

		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1

		filename1 = paste0("r1_", i, ".rdata")
		filename4 = paste0("r2_", i, ".rdata")

		readTime = readTime - proc.time()[3]
		toRead1 = file(file.path(params$dplocalPath, filename1), "rb")
		R1 = matrix(readBin(con = toRead1, what = numeric(), n = n * n,
		  								  endian = "little"), nrow = n, ncol = n)
		readSize = readSize + file.size(file.path(params$dplocalPath, filename1))
		close(toRead1)
		RW = matrix(readBin(con = toRead2, what = numeric(), n = n * p2,
											 endian = "little"), nrow = n, ncol = p2)
		readTime = readTime + proc.time()[3]

		W = t(R1) %*% RW
		R2 = RandomOrthonomalMatrix(p2)
		WR2 = W %*% R2

		writeTime = writeTime - proc.time()[3]
		toWrite4 = file(file.path(params$dplocalPath, filename4), "wb")
		writeBin(as.vector(R2), con = toWrite4, endian = "little")
		close(toWrite4)
		writeSize = writeSize + file.size(file.path(params$dplocalPath, filename4))
		writeBin(as.vector(WR2), con = toWrite3, endian = "little")
		writeTime = writeTime + proc.time()[3]
		if ((i + 1) %in% params$container$filebreak.RW || i == numBlocks) {
			close(toRead2)
			readSize = readSize + file.size(file.path(params$readPath[["B"]], filename2))
		}
		if ((i + 1) %in% params$container$filebreak.WR || i == numBlocks) {
			close(toWrite3)
			writeSize = writeSize + file.size(file.path(params$writePath, filename3))
		}

		pbar = MakeProgressBar2(i, pbar)
	}

	params = AddToLog(params, "ProcessWLinear.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetWRLinear.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWRLinear.A3\n\n")
  XATXA = t(data$X) %*% data$X
	XATY  = t(data$X) %*% data$Y
	writeTime = proc.time()[3]
	save(XATXA, XATY, file = file.path(params$writePath, "xatxa.rdata"))
	writeSize = file.size(file.path(params$writePath, "xatxa.rdata"))
	writeTime = proc.time()[3] - writeTime

	p2 = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "p2.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "p2.rdata"))
	readTime = proc.time()[3] - readTime

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "XA'(I-Z*Z')XB*R")

	containerCt.WR = 0
	containerCt.PR = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.WR) {
			containerCt.WR = containerCt.WR + 1
			filename1 = paste0("cwr_", containerCt.WR, ".rdata")
			toRead = file(file.path(params$readPath[["T"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.PR) {
			containerCt.PR = containerCt.PR + 1
			filename2 = paste0("cpr_", containerCt.PR, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}

		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n = stp - strt + 1

		readTime = readTime - proc.time()[3]
		WR = matrix(readBin(con = toRead, what = numeric(), n = n * p2,
			  								 endian = "little"), nrow = n, ncol = p2)
		readTime = readTime + proc.time()[3]

		YXA = cbind(data$Y[strt:stp, ], data$X[strt:stp, ])
		PR = t(YXA) %*% WR
		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(PR), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]

		if ((i + 1) %in% params$container$filebreak.WR || i == numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["T"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.PR || i == numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename2))
		}

		pbar = MakeProgressBar2(i, pbar)
	}
	params = AddToLog(params, "GetWRLinear.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetProductsLinear.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsLinear.T3\n\n")
  n  = params$n
	p1 = params$p1
	p2 = params$p2
	XATXA = 0
	XBTXB = 0
	XATY  = 0
	YXATXB = 0

	numBlocks = params$blocks$numBlocks
	readTime = proc.time()[3]
	load(file.path(params$readPath[["B"]], "xbtxb.rdata"))
	load(file.path(params$readPath[["A"]], "xatxa.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["B"]], "xbtxb.rdata")),
								 file.size(file.path(params$readPath[["A"]], "xatxa.rdata")))
	readTime = proc.time()[3] - readTime

	pbar = MakeProgressBar1(numBlocks, "X'X")

	containerCt.PR = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.PR) {
			containerCt.PR = containerCt.PR + 1
			filename1 = paste0("cpr_", containerCt.PR, ".rdata")
			toRead = file(file.path(params$readPath[["A"]], filename1), "rb")
			readSize = readSize + file.size(file.path(params$readPath[["A"]], filename1))
		}

		filename1 = paste0("r2_", i, ".rdata")

		readTime = readTime - proc.time()[3]
		toRead1 = file(file.path(params$dplocalPath, filename1), "rb")
		R2 = matrix(readBin(con = toRead1, what = numeric(), n = p2 * p2,
												endian = "little"), p2, p2)
		readSize = readSize + file.size(file.path(params$dplocalPath, filename1))
		close(toRead1)
		PR = matrix(readBin(con = toRead, what = numeric(), n = (p1 + 1) * p2,
			 								  endian = "little"), p1 + 1, p2)
		readTime = readTime + proc.time()[3]

		YXATXB = YXATXB + PR %*% t(R2)

		if ((i + 1) %in% params$container$filebreak.PR || i == numBlocks) {
			close(toRead)
		}

		pbar = MakeProgressBar2(i, pbar)
	}

	YTXB = YXATXB[1, , drop = FALSE]
	XATXB = YXATXB[-1, , drop = FALSE]
	XTX = rbind(cbind(XATXA, XATXB), cbind(t(XATXB), XBTXB))
	XTY = rbind(XATY, t(YTXB))

	XTXLasso = XTX / (n - 1)
	XTYLasso = params$sdy * XTY / sqrt(n - 1)

	params$xtx = XTX
	params$xty = XTY
	params$xtxLasso = XTXLasso
	params$xtyLasso = XTYLasso

	params$converged = TRUE

	params = AddToLog(params, "GetProductsLinear.T3", readTime, readSize, 0, 0)
	return(params)
}


ComputeResultsLinear.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsLinear.T3\n\n")
  stats    = params$stats
	stats$converged = params$converged
	stats$failed    = FALSE
	Anames   = params$colnamesA
	Bnames   = params$colnamesB
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
		if (rcond(xtx[tempIndicies, tempIndicies]) > 10^3 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}

	Aindex        = which(indicies <= length(Anames))
	AIndiciesKeep = indicies[Aindex]
	BIndiciesKeep = indicies[-Aindex] - length(Anames)
	names.old     = c(Anames, Bnames)
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

	names(stats$party)           = names.old
	names(stats$coefficients)    = names.old
	names(stats$secoef)          = names.old
	names(stats$tvals)           = names.old
	names(stats$pvals)           = names.old
	colnames(stats$xtx)          = names.old
	rownames(stats$xtx)          = names.old
	colnames(stats$xty)          = colnames(params$xty)
	rownames(stats$xty)          = names.old

	params$stats = stats
	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeResultsLinear.T3", 0, 0, writeTime, writeSize)
	return(params)
}


GetResultsLinear.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLinear.A3\n\n")
  params$converged = TRUE
	stats = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats

	params = AddToLog(params, "GetResultsLinear.A3", readTime, readSize, 0, 0)
	return(params)
}


GetResultsLinear.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsLinear.B3\n\n")
  params$converged = TRUE
	stats = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats

	params = AddToLog(params, "GetResultsLinear.B3", readTime, readSize, 0, 0)
	return(params)
}


############################## PARENT FUNCTIONS ###############################


PartyAProcess3Linear = function(data,
                                yname          = NULL,
																monitorFolder  = NULL,
                                sleepTime      = 10,
                                maxWaitingTime = 24 * 60 * 60,
																popmednet      = TRUE,
																trace          = FALSE) {

  params = PrepareParams.3p("linear", "A",
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
  data = PrepareDataLinear.A3(params, data, yname)
  params = AddToLog(params, "PrepareDataLinear.A3", 0, 0, 0, 0)

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


  params$algIterationCounter = 2
  params = GetResultsLinear.A3(params)
  params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
  return(params$stats)
}


PartyBProcess3Linear = function(data,
																monitorFolder  = NULL,
																sleepTime      = 10,
																maxWaitingTime = 24 * 60 * 60,
																popmednet      = TRUE,
																trace          = FALSE) {
  params = PrepareParams.3p("linear", "B",
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

  data = PrepareDataLinear.B3(params, data)
  params = AddToLog(params, "PrepareDataLinear.B3", 0, 0, 0, 0)

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

  params$algIterationCounter = 2
  params = GetResultsLinear.B3(params)
  params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
  return(params$stats)
}


PartyTProcess3Linear = function(monitorFolder         = NULL,
																msreqid               = "v_default_0_000",
																blocksize             = 500,
																sleepTime             = 10,
																maxWaitingTime        = 24 * 60 * 60,
																popmednet             = TRUE,
																trace                 = FALSE) {
	params = PrepareParams.3p("linear", "T", msreqid = msreqid,
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

	params   = PrepareParamsLinear.T3(params)
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

	params$algIterationCounter = 2
	params = GetProductsLinear.T3(params)
	params = ComputeResultsLinear.T3(params)
	files = "stats.rdata"
	params = SendPauseContinue.3p(params, filesA = files, filesB = files, from  = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = SendPauseQuit.3p(params, sleepTime = sleepTime)
	SummarizeLog.3p(params)
	return(params$stats)
}
