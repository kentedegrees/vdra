#################### DISTRIBUTED COX REGRESSION FUNCTIONS ####################

PrepareDataCox.A3 = function(params, data, yname, strata, mask) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataCox.A3\n\n")
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
  if (ncol(data) == 1) {
    cat("Error: There is only one variable.  Need at least two variables for time and censoring.\n ")
    workdata$failed = TRUE
    return(workdata)
  }
  if (!is.null(yname) && class(yname) != "character") {
    cat("Error: time and censor labels are not strings.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (length(yname) == 1) {
    cat("Error: Only one name for time and censor variables provided.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (length(yname) > 2) {
    cat("Error: More than two names for time and censor variables provided.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (is.null(colnames(data))) {
    cat("Warning: variables are not named. Assuming variable 1 is time and variable 2 is censoring.  Assigning labels to the rest of the columns.\n\n")
    if (is.null(yname)) {
      yname = c("time", "censor")
    }
    if (ncol(data) == 2) {
      colnames(data) = yname
    } else {
      colnames(data) = c(yname, paste0("dp1:", 1:(ncol(data) - 2)))
    }
  } else {
    if (is.null(yname)) {
      yname = colnames(data)[1:2]
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
  if (!(yname[1] %in% colnames(data))) {
    cat("Error: time variable", paste0("'", yname[1], "'"), "not found.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (!(yname[2] %in% colnames(data))) {
    cat("Error: censoring variable", paste0("'", yname[2], "'"), "not found.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  timeColIndex = which(colnames(data) %in% yname[1])
  if (class(data[1, timeColIndex]) != "numeric" & class(data[1, timeColIndex]) != "integer") {
    cat("Error: time variable", paste0("'", yname[1], "'"), "is not numeric.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  censorColIndex = which(colnames(data) %in% yname[2])
  if (class(data[1, censorColIndex]) != "numeric" & class(data[1, censorColIndex]) != "integer") {
    cat("Error: censoring variable", paste0("'", yname[2], "'"), "is not numeric.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  if (sum(data[, censorColIndex] %in% c(0, 1)) < nrow(data)) {
    cat("Error: censoring data should be only 0's and 1's.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  workdata$survival        = list()
  workdata$survival$rank   = data[, timeColIndex]
  workdata$survival$status = data[, censorColIndex]
  data = data[, -c(timeColIndex, censorColIndex), drop = FALSE]

  # Extract the strata First
  workdata$strata = list()
  if (!is.null(strata)) {
    if (class(strata) != "character") {
      cat("Error: strata is not a valid variable name.\n\n")
      workdata$failed = TRUE
      return(workdata)
    }
    if (length(strata) > 0) {
      idx = which(strata %in% colnames(data))
      if (length(idx) > 0) {
        workdata$strata$strataFromA = strata[idx]
        workdata$strata$strataFromB = strata[-idx]
      } else {
        workdata$strata$strataFromA = c()
        workdata$strata$strataFromB = strata
      }
      idx = which(colnames(data) %in% workdata$strata$strataFromA)
      if (length(idx) > 0) {
        workdata$strata$X = as.data.frame(data[, idx, drop = FALSE])
        colnames(workdata$strata$X) = colnames(data)[idx]
        workdata$strata$legend = list()
        data              = data[, -idx, drop = FALSE]
        #randomize values for levels for each strata
        for (i in 1:ncol(workdata$strata$X)) {
          levels = levels(as.factor(workdata$strata$X[, i]))
          workdata$strata$legend[[colnames(workdata$strata$X)[i]]] = levels
          if (mask) {
            levels = sample(levels, length(levels))
            workdata$strata$legend[[colnames(workdata$strata$X)[i]]] = rep("NA", length(levels))
          }
          workdata$strata$X[, i] = sapply(workdata$strata$X[, i], function(x) { which(levels %in% x)})
        }
      }
    }
  }

  # Now covert what is left to a matrix and convert all non-numeric data to indicators

  if (class(data) == "matrix") {
    if (ncol(data) == 0) {
      workdata$X = data
    } else {
      workdata$X = scale(data, center = TRUE, scale = FALSE)
    }
  } else if (ncol(data) > 0) { # class(data) == "data.frame"
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
    workdata$X = scale(as.matrix(data[index]), center = TRUE, scale = FALSE)
  } else {
    workdata$X = as.matrix(data)
  }
  return(workdata)
}


PrepareDataCox.B3 = function(params, data, strata, mask) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareDataCox.B3\n\n")
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

  # Extract the strata First
  workdata$strata = list()
  if (!is.null(strata)) {
    if (class(strata) != "character") {
      cat("Error: strata is not a valid variable name.\n\n")
      workdata$failed = TRUE
      return(workdata)
    }
    if (length(strata) > 0) {
      idx = which(strata %in% colnames(data))
      if (length(idx) > 0) {
        workdata$strata$strataFromB = strata[idx]
        workdata$strata$strataFromA = strata[-idx]
      } else {
        workdata$strata$strataFromB = c()
        workdata$strata$strataFromA = strata
      }
      idx = which(colnames(data) %in% workdata$strata$strataFromB)
      if (length(idx) > 0) {
        workdata$strata$X = as.data.frame(data[, idx])
        colnames(workdata$strata$X) = colnames(data)[idx]
        workdata$strata$legend = list()
        data              = data[, -idx, drop = FALSE]
        #randomize values for levels for each strata
        for (i in 1:ncol(workdata$strata$X)) {
          levels = levels(as.factor(workdata$strata$X[, i]))
          workdata$strata$legend[[colnames(workdata$strata$X)[i]]] = levels
          if (mask) {
            levels = sample(levels, length(levels))
            workdata$strata$legend[[colnames(workdata$strata$X)[i]]] = rep("NA", length(levels))
          }
          outcome = sapply(workdata$strata$X[, i], function(x) { which(levels %in% x)})
          workdata$strata$X[, i] = outcome
        }
      }
    }
  }
  if (ncol(data) < 1) {
    cat("Error: After removing strata, data is empty.  Party B must supply at least one non-strata covariate.\n\n")
    workdata$failed = TRUE
    return(workdata)
  }
  # Now covert what is left to a matrix and convert all non-numeric data to indicators
  if (class(data) == "matrix") {
    if (class(data[1, 1]) != "numeric" && class(data[1, 1]) != "integer") {
      data = as.data.frame(data)
    } else {
      workdata$X = scale(data, center = TRUE, scale = FALSE)
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
    workdata$X = scale(as.matrix(data[index]), center = TRUE, scale = FALSE)
  }
  return(workdata)
}


PrepareParamsCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.A3\n\n")
  params$n = nrow(data$X)
	params$p1 = ncol(data$X)
	params$p2 = 0
	params$colnames = colnames(data$X)

	if (requireNamespace("survival", quietly = TRUE)) {
	  library(survival)
	  params$survivalInstalled = TRUE
	} else {
	  params$survivalInstalled = FALSE
	}

	pa           = list()
	pa$p1        = params$p1
	pa$n         = params$n
	pa$analysis  = params$analysis
	pa$colnames  = params$colnames
	pa$strataFromA = data$strata$strataFromA
	pa$strataFromB = data$strata$strataFromB

	writeTime = proc.time()[3]
	save(pa, file = file.path(params$writePath, "pa.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, "pa.rdata")))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "PrepareParamsCox.A3", 0, 0, writeTime, writeSize)
	return(params)
}


PrepareParamsCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.B3\n\n")
  params$n = nrow(data$X)
	params$p1 = 0
	params$p2 = ncol(data$X)
	params$colnames = colnames(data$X)

	if (requireNamespace("survival", quietly = TRUE)) {
	  library(survival)
	  params$survivalInstalled = TRUE
	} else {
	  params$survivalInstalled = FALSE
	}

	pb           = list()
	pb$p2        = params$p2
	pb$n         = params$n
	pb$analysis  = params$analysis
	pb$colnames  = params$colnames
	pb$strataFromA = data$strata$strataFromA
	pb$strataFromB = data$strata$strataFromB

	writeTime = proc.time()[3]
	save(pb, file = file.path(params$writePath, "pb.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, "pb.rdata")))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "PrepareParamsCox.B3", 0, 0, writeTime, writeSize)
	return(params)
}


PrepareParamsCox.T3 = function(params, cutoff, maxIterations) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareParamsCox.T3\n\n")
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

	if (requireNamespace("survival", quietly = TRUE)) {
	  library(survival)
	  params$survivalInstalled = TRUE
	} else {
	  params$survivalInstalled = FALSE
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
	params$maxIterations = maxIterations

	params$AstrataFromA  = pa$strataFromA
	params$AstrataFromB  = pa$strataFromB
	params$BstrataFromA  = pb$strataFromA
	params$BstrataFromB  = pb$strataFromB

	writeTime = proc.time()[3]
	save(cutoff, maxIterations, file = file.path(params$writePath, "maxiterations.rdata"))
	writeSize = file.size(file.path(params$writePath, "maxiterations.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "PrepareParamsCox.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


CheckStrataCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckStrataCox.T3\n\n")
  if (length(params$AstrataFromA) == length(params$BstrataFromA) &&
			length(params$AstrataFromB) == length(params$BstrataFromB) &&
			ifelse(length(params$AstrataFromA) == 0, TRUE,
						 order(params$AstrataFromA) == order(params$BstrataFromA)) &&
			ifelse(length(params$AstrataFromB) == 0, TRUE,
						 order(params$AstrataFromB) == order(params$BstrataFromB))) {
		params$strataFromA = params$AstrataFromA
		params$strataFromB = params$BstrataFromB
		params$AstrataFromA = params$AstrataFromB =
			params$BstrataFromA = params$BstrataFromB = NULL
		params$getStrata = TRUE
	} else {
		params$getStrata = FALSE
		AcapB = intersect(params$AstrataFromA, params$BstrataFromB)
		BcapA = intersect(params$BstrataFromA, params$AstrataFromB)
		if (length(AcapB) > 0) {
			params$errorMessage =
				paste("Party A and Party B have", length(AcapB), "variable(s) with the same name which are used in the strata.",
							"These variable(s) are <", paste0(AcapB, collapse = ", "), ">.",
							"Make sure the variables from each party have distinct names.")
		}
		else if (length(BcapA) > 0) {
			params$errorMessage =
				paste("Party A and Party B have specified", length(BcapA), "variable(s) for the strata which are not found in the data.",
							"These variable(s) are <", paste0(BcapA, collapse = ", "), ">.",
							"Check the spelling of the variables names and / or remove them from the strata.")
		} else {
			params$errorMessage =
				paste("Party A and Party B have specified different stata.",
							"Verify that both parties specify the same strata.")
		}
		params$failed = TRUE
	}
	empty = 0
	save(empty, file = file.path(params$writePath, "empty.rdata"))
	params = AddToLog(params, "CheckStrataCox.T3", 0, 0, 0, 0)
	return(params)
}


SendStrataCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataCox.A3\n\n")
  Astrata = data$strata
	survival = data$survival
	writeTime = proc.time()[3]
	save(Astrata, survival, file = file.path(params$writePath, "Astrata.rdata"))
	writeSize = file.size(file.path(params$writePath, "Astrata.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendStrataCox.A3", 0, 0, writeTime, writeSize)
	return(params)
}


SendStrataCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SendStrataCox.B3\n\n")
  Bstrata = data$strata
	writeTime = proc.time()[3]
	save(Bstrata, file = file.path(params$writePath, "Bstrata.rdata"))
	writeSize = file.size(file.path(params$writePath, "Bstrata.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "SendStrataCox.B3", 0, 0, writeTime, writeSize)
	return(params)
}


PrepareStrataCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "PrepareStrataCox.T3\n\n")
  Astrata  = NULL
	Bstrata  = NULL
	survival = NULL
	strataTemp   = list()

	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "Astrata.rdata"))
	load(file.path(params$readPath[["B"]], "Bstrata.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "Astrata.rdata")) +
		file.size(file.path(params$readPath[["B"]], "Bstrata.rdata"))
	readTime = proc.time()[3] - readTime


	if (length(params$strataFromA) == 0 && length(params$strataFromB) == 0) {
		strataTemp$X = data.frame(const__ = rep(1, params$n))
		strataTemp$legend = FALSE
	} else if (length(params$strataFromA) == 0) {
		strataTemp$X = Bstrata$X
		strataTemp$legend = Bstrata$legend
	} else if (length(params$strataFromB) == 0) {
		strataTemp$X = Astrata$X
		strataTemp$legend = Astrata$legend
	} else {
		strataTemp$X = cbind(Astrata$X, Bstrata$X)
		strataTemp$legend = c(Astrata$legend, Bstrata$legend)
	}

	sorted = do.call("order", cbind(strataTemp$X, survival$rank, survival$status))
	strataTemp$X = strataTemp$X[sorted, , drop = FALSE]
	survival$rank   = survival$rank[sorted]
	survival$status = survival$status[sorted]
	survival$sorted = sorted
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
		strata[[i]]$J = as.integer(sum(temp0[, 2] > 0))         # number of distinct failure times
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
	params$survival = survival

	writeTime = proc.time()[3]
	save(survival, file = file.path(params$writePath, "survival.rdata"))
	writeSize = file.size(file.path(params$writePath, "survival.rdata"))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "PrepareStrataCox.T3", readTime, readSize, writeTime, writeSize)

	return(params)
}


SortDataCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SortDataCox.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "survival.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "survival.rdata"))
	readTime = proc.time()[3] - readTime
	data$X = data$X[survival$sorted, , drop = FALSE]
	data$survival = survival
	data$readTime = readTime
	data$readSize = readSize
	return(data)
}


GetZCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetZCox.A3\n\n")
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
		Z = FindOrthogonalVectors(data$X[strt:stp, ], g)

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(Z), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]
		if ((i + 1) %in% params$container$filebreak.Z || i == numBlocks) {
			close(toWrite)
			writeSize = writeSize + file.size(file.path(params$writePath, filename))
		}
		pbar = MakeProgressBar2(i, pbar)
	}
	params = AddToLog(params, "GetZCox.A3", 0, 0, writeTime, writeSize)
	return(params)
}


SortDataCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "SortDataCox.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "survival.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "survival.rdata"))
	readTime = proc.time()[3] - readTime
	data$X = data$X[survival$sorted, , drop = FALSE]
	data$survival = survival
	data$readTime = readTime
	data$readSize = readSize
	return(data)
}


GetSXBCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetSXBCox.B3\n\n")
  S = matrix(0, nrow = params$n, ncol = length(data$survival$strata))
	for (i in 1:length(data$survival$strata)) {
		S[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
	}
	STXB = t(S) %*% data$X
	writeTime = proc.time()[3]
	save(STXB, file = file.path(params$writePath, "sxb.rdata"))
	writeSize = file.size(file.path(params$writePath, "sxb.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetSXBCox.B3", 0, 0, writeTime, writeSize)
	return(params)
}


GetWRCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetWRCox.A3\n\n")
  XATXA = t(data$X) %*% data$X
	writeTime = proc.time()[3]
	save(XATXA, file = file.path(params$writePath, "xatxa.rdata"))
	writeSize = file.size(file.path(params$writePath, "xatxa.rdata"))
	writeTime = proc.time()[3] - writeTime

	p2 = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "p2.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "p2.rdata"))
	readTime = proc.time()[3] - readTime
	params$p2 = p2

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "XA'(I - Z*Z')XB*R")

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

		PR = t(data$X[strt:stp, ]) %*% WR
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
	params = AddToLog(params, "GetWRCox.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetSXACox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetSXACox.A3\n\n")
  S = matrix(0, nrow = params$n, ncol = length(data$survival$strata))
	for (i in 1:length(data$survival$strata)) {
		S[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
	}
	STXA = t(S) %*% data$X
	writeTime = proc.time()[3]
	save(STXA, file = file.path(params$writePath, "sxa.rdata"))
	writeSize = file.size(file.path(params$writePath, "sxa.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetSXACox.A3", 0, 0, writeTime, writeSize)
	return(params)
}


GetProductsCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetProductsCox.T3\n\n")
  p1 = params$p1
	p2 = params$p2
	XATXA = 0
	XBTXB = 0
	XATXB = 0
	STXA  = 0
	STXB  = 0

	numBlocks = params$blocks$numBlocks
	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "xatxa.rdata"))
	load(file.path(params$readPath[["B"]], "xbtxb.rdata"))
	load(file.path(params$readPath[["A"]], "sxa.rdata"))
	load(file.path(params$readPath[["B"]], "sxb.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["B"]], "xbtxb.rdata")),
								 file.size(file.path(params$readPath[["A"]], "xatxa.rdata")),
								 file.size(file.path(params$readPath[["B"]], "sxb.rdata")),
								 file.size(file.path(params$readPath[["A"]], "sxa.rdata")))
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
		PR = matrix(readBin(con = toRead, what = numeric(), n = p1 * p2,
												endian = "little"), p1, p2)
		readTime = readTime + proc.time()[3]
		XATXB = XATXB + PR %*% t(R2)
		if ((i + 1) %in% params$container$filebreak.PR || i == numBlocks) {
			close(toRead)
		}
		pbar = MakeProgressBar2(i, pbar)
	}

	num = length(params$survival$strata)
	STS = matrix(0, nrow = num, ncol = num)
	for (i in 1:num) {
		STS[i, i] = params$survival$strata[[i]]$end - params$survival$strata[[i]]$start + 1
	}

	XTX = rbind(cbind(STS, STXA, STXB),
							cbind(t(STXA), XATXA, XATXB),
							cbind(t(STXB), t(XATXB), XBTXB))

	params$xtx = XTX

	params = AddToLog(params, "GetProductsCox.T3", readTime, readSize, 0, 0)
	return(params)
}


CheckColinearityCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityCox.T3\n\n")
  xtx = params$xtx
	nrow = nrow(xtx)
	numStrata = length(params$survival$strata)
	indicies = 1:numStrata
	for (i in (1 + numStrata):nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(xtx[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}

	Anames = params$colnamesA
	Bnames = params$colnamesB
	indicies = indicies[-(1:numStrata)] - numStrata# Get rid of the strata indicators
	Aindex = which(indicies <= length(Anames))
	Bindex = which(indicies > length(Anames))
	params$indicies      = indicies
	params$AIndiciesKeep = indicies[Aindex]
	params$BIndiciesKeep = indicies[Bindex] - length(Anames)
	AnamesKeep = Anames[params$AIndiciesKeep]
	BnamesKeep = Bnames[params$BIndiciesKeep]
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

	Aindicies = params$AIndiciesKeep - numStrata
	Bindicies = params$BIndiciesKeep

	Aindicies = params$AIndiciesKeep
	Bindicies = params$BIndiciesKeep

	colnamesA.old = params$colnamesA.old
	p2 = params$p2

	writeTime = proc.time()[3]
	save(p2, Aindicies, file = file.path(params$writePath, "Aindicies.rdata"))
	save(colnamesA.old, Bindicies, file = file.path(params$writePath, "Bindicies.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("Aindicies.rdata",
																													"Bindicies.rdata"))))
	writeTime = proc.time()[3] - writeTime

	if (params$p2 == 0) {
		params$failed = TRUE
		params$errorMessage = "All of party B's covariates are either linear or are colienar with Party A's covariates."
	}
	params = AddToLog(params, "CheckColinearityCox.T3", 0, 0, writeTime, writeSize)
	return(params)
}


ComputeInitialBetasCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeInitialBetasCox.T3\n\n")
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

	writeTime = proc.time()[3]
	save(Abetas, file = file.path(params$writePath, "betasA.rdata"))
	save(Bbetas, file = file.path(params$writePath, "betasB.rdata"))
	save(converged, maxIterExceeded,
			 file = file.path(params$writePath, "converged.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("betasA.rdata",
																													"betasB.rdata",
																													"converged.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ComputeInitialBetasCox.T3", 0, 0, writeTime, writeSize)
}


UpdateParamsCox.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Aindicies.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "Aindicies.rdata"))
	readTime = proc.time()[3] - readTime
	params$p             = length(Aindicies)
	params$p2            = p2
	params$colnames.old  = params$colnames
	params$colnames      = params$colnames[Aindicies]
	params$AIndiciesKeep = Aindicies
	params = AddToLog(params, "UpdateParamsCox.A3, UpdateDataCox.A3", readTime, readSize, 0, 0)
	return(params)
}


UpdateParamsCox.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateParamsCox.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Bindicies.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "Bindicies.rdata"))
	readTime = proc.time()[3] - readTime
	params$p             = length(Bindicies)
	params$colnamesA.old = colnamesA.old
	params$colnames.old  = params$colnames
	params$colnames      = params$colnames[Bindicies]
	params$BIndiciesKeep = Bindicies
	params = AddToLog(params, "UpdateParamsCox.B3", readTime, readSize, 0, 0)
	return(params)
}


UpdateDataCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.A3\n\n")
  data$X = as.matrix(data$X[, params$AIndiciesKeep, drop = FALSE])
	return(data)
}


UpdateDataCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "UpdateDataCox.B3\n\n")
  data$X = as.matrix(data$X[, params$BIndiciesKeep, drop = FALSE])
	return(data)
}


GetBetaACox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetBeataACox.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "converged.rdata"))
	load(file.path(params$readPath[["T"]], "betasA.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("converged.rdata",
																															 "betasA.rdata"))))
	readTime = proc.time()[3] - readTime
	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	params$betas = Abetas
	params = AddToLog(params, "GetBetaACox.A3", readTime, readSize, 0, 0)
	return(params)
}


GetBetaBCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetBetaBCox.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "converged.rdata"))
	load(file.path(params$readPath[["T"]], "betasB.rdata"))
	readSize = sum(file.size(file.path(params$readPath[["T"]], c("converged.rdata",
																															 "betasB.rdata"))))
	readTime = proc.time()[3] - readTime
	params$converged = converged
	params$maxIterExceeded = maxIterExceeded
	params$betas = Bbetas
	params = AddToLog(params, "GetBetaBCox.B3", readTime, readSize, 0, 0)
	return(params)
}


GetXABetaACox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXABetaACox.A3\n\n")
  XAbetaA = data$X %*% params$betas
	writeTime = proc.time()[3]
	save(XAbetaA, file = file.path(params$writePath, "xabetaa.rdata"))
	writeSize = file.size(file.path(params$writePath, "xabetaa.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetXABetaACox.A3", 0, 0, writeTime, writeSize)
	return(params)
}


GetXBBetaBCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXBBetaBCox.B3\n\n")
  XBbetaB = data$X %*% params$betas
	writeTime = proc.time()[3]
	save(XBbetaB, file = file.path(params$writePath, "xbbetab.rdata"))
	writeSize = file.size(file.path(params$writePath, "xbbetab.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "GetXBBetaBCox.B3", 0, 0, writeTime, writeSize)
	return(params)
}


ComputeLogLikelihoodCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeLogLikelihoodCox.T3\n\n")
  n  = params$n
	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "xabetaa.rdata"))
	load(file.path(params$readPath[["B"]], "xbbetab.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "xabetaa.rdata")) +
		         file.size(file.path(params$readPath[["B"]], "xbbetab.rdata"))
	readTime = proc.time()[3] - readTime
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
		pbar = MakeProgressBar1(numEvents, "Loglikelihood")
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
					pbar = MakeProgressBar2(stepCounter, pbar)
				}
			}
		}
		if (loglikelihood > params$loglikelihood.old || stepSize < 0.5^6) {
			computeLoglikelihood = FALSE
		} else {
			cat("Step Halving\n\n")
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
	writeTime = proc.time()[3]
	save(Xbeta, file = file.path(params$writePath, "Xbeta.rdata"))
	writeSize = file.size(file.path(params$writePath, "Xbeta.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeLogLikelihoodCox.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeXBDeltaLCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeXBDeltaLCox.B3\n\n")
  Xbeta = NULL
	p2 = params$p
	n = params$n

	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Xbeta.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "Xbeta.rdata"))
	readTime = proc.time()[3] - readTime

	numEvents = sum(data$survival$status)

	w = exp(Xbeta)
	deltal = as.numeric(data$survival$status)
	deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
	# a pass by reference with the C call.
	W.XB = matrix(0, n, p2)

	.Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XB,
				as.integer(n), as.integer(p2), as.integer(numEvents))

  containerCt.RZ = 0
	containerCt.Cox = 0
	writeSize = 0
	writeTime = 0

	pbar = MakeProgressBar1(params$blocks$numBlocks, "R*(I-Z*Z')W*XB")
	for (i in 1:params$blocks$numBlocks) {
		if (i %in% params$container$filebreak.RZ) {
			containerCt.RZ = containerCt.RZ + 1
			filename1 = paste0("crz_", containerCt.RZ, ".rdata")
			toRead = file(file.path(params$readPath[["T"]], filename1), "rb")
		}
		if (i %in% params$container$filebreak.Cox) {
			containerCt.Cox = containerCt.Cox + 1
			filename2 = paste0("cCox_", containerCt.Cox, ".rdata")
			toWrite = file(file.path(params$writePath, filename2), "wb")
		}
		strt = params$blocks$starts[i]
		stp = params$blocks$stops[i]
		n2 = stp - strt + 1
		g = params$blocks$g[i]

		readTime = readTime - proc.time()[3]
		RZ = matrix(readBin(con = toRead, what = numeric(), n = n2 * n2,
											 endian = "little"), nrow = n2, ncol = n2)
		readTime = readTime + proc.time()[3]

		IZ.tZ.W.XBtemp = RZ %*% W.XB[strt:stp, , drop = FALSE]

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(IZ.tZ.W.XBtemp), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]

		if ((i + 1) %in% params$container$filebreak.RZ || i == params$blocks$numBlocks) {
			close(toRead)
			readSize = readSize + file.size(file.path(params$readPath[["T"]], filename1))
		}
		if ((i + 1) %in% params$container$filebreak.Cox || i == params$blocks$numBlocks) {
			close(toWrite)
			writeSize = writeSize = file.size(file.path(params$writePath, filename2))
		}
		pbar = MakeProgressBar2(i, pbar)
	}

	tXB.W.XB = t(data$X) %*% W.XB
	tXB.deltal = t(data$X) %*% deltal

	writeTime = writeTime - proc.time()[3]
	save(tXB.deltal, tXB.W.XB, file = file.path(params$writePath, "tXB_W_XB.rdata"))
	writeSize = writeSize + file.size(file.path(params$writePath, "tXB_W_XB.rdata"))
	writeTime = writeTime + proc.time()[3]

	params = AddToLog(params, "ComputeXBDeltaLCox.B3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ComputeXADeltaLCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeXADeltaLCox.A3\n\n")
  Xbeta = NULL
	p1 = params$p
	n = params$n

	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "Xbeta.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "Xbeta.rdata"))
	readTime = proc.time()[3] - readTime

	numEvents = sum(data$survival$status)

	w = exp(Xbeta)
	deltal = as.numeric(data$survival$status)
	deltal[1] = deltal[1]  # This is to force R to make a copy since we are exploiting
	# a pass by reference with the C call.
	W.XA = matrix(0, n, p1)

	.Call("ComputeCox", data$survival$strata, data$X, w, deltal, W.XA,
				as.integer(n), as.integer(p1), as.integer(numEvents))

	tXA.W.XA = t(data$X) %*% W.XA
	tXA.deltal = t(data$X) %*% deltal

	writeTime = proc.time()[3]
	save(tXA.deltal, tXA.W.XA, file = file.path(params$writePath, "tXA_W_XA.rdata"))
	writeSize = file.size(file.path(params$writePath, "tXA_W_XA.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeXADeltaLCox.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ProcessVCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessVCox.T3\n\n")
  writeTime = 0
	writeSize = 0
	readTime  = 0
	readSize  = 0
	p2 = params$p2

	numBlocks = params$blocks$numBlocks
	pbar = MakeProgressBar1(numBlocks, "(I-Z*Z')W*XB*R")

	containerCt.RV = 0
	containerCt.VR = 0
	for (i in 1:numBlocks) {
		if (i %in% params$container$filebreak.RV) {
			containerCt.RV = containerCt.RV + 1
			filename2 = paste0("cCox_", containerCt.RV, ".rdata")
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
			readSize = readSize + file.size(file.path(params$readPath[["B"]], filename2))
		}
		if ((i + 1) %in% params$container$filebreak.VR || i == numBlocks) {
			close(toWrite3)
			writeSize = writeSize + file.size(file.path(params$writePath, filename3))
		}

		pbar = MakeProgressBar2(i, pbar)
	}
	params = AddToLog(params, "ProcessVCox.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


GetXRCox.A3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "GetXRCox.A3\n\n")
  writeTime = 0
	writeSize = 0
	readTime  = 0
	readSize  = 0
	p2 = params$p2
	containerCt.VR = 0
	containerCt.XR = 0
	pbar = MakeProgressBar1(params$blocks$numBlocks, "XA'(I-Z*Z')W*XB*R")
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

		readTime = readTime - proc.time()[3]
		VR = matrix(readBin(con = toRead, what = numeric(), n = n * p2,
												endian = "little"), nrow = n, ncol = p2)
		readTime = readTime + proc.time()[3]
		XR = t(data$X[strt:stp, , drop = FALSE]) %*% VR

		writeTime = writeTime - proc.time()[3]
		writeBin(as.vector(XR), con = toWrite, endian = "little")
		writeTime = writeTime + proc.time()[3]

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

	params = AddToLog(params, "GetXRCox.A3", readTime, readSize, writeTime, writeSize)
	return(params)
}


ProcessXtWXCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ProcessXtWXCox.T3\n\n")
  p1 = params$p1
	p2 = params$p2

	tXA.deltal = NULL
	tXB.deltal = NULL
	tXA.W.XA   = NULL
	tXB.W.XB   = NULL
	readTime = proc.time()[3]
	load(file.path(params$readPath[["A"]], "tXA_W_XA.rdata"))
	load(file.path(params$readPath[["B"]], "tXB_W_XB.rdata"))
	readSize = file.size(file.path(params$readPath[["A"]], "tXA_W_XA.rdata")) +
		file.size(file.path(params$readPath[["B"]], "tXB_W_XB.rdata"))
	readTime = proc.time()[3] - readTime

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


	xtwx = rbind(cbind(tXA.W.XA, XATWXB), cbind(t(XATWXB), tXB.W.XB))
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
		params = AddToLog(params, "ProcessXtWXCox.T3", readTime, readSize, 0, 0)
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
		maxIterExceeded = params$algIterationCounter >= params$maxIterations
	}

	params$converged = converged
	params$maxIterExceeded = maxIterExceeded

	Abetas = params$betas[1:p1]
	Bbetas = params$betas[(p1 + 1):(p1 + p2)]
	writeTime = proc.time()[3]
	save(Abetas, file = file.path(params$writePath, "betasA.rdata"))
	save(Bbetas, file = file.path(params$writePath, "betasB.rdata"))
	save(converged, maxIterExceeded,
			 file = file.path(params$writePath, "converged.rdata"))
	writeSize = sum(file.size(file.path(params$writePath, c("betasA.rdata",
																													"betasB.rdata",
																													"converged.rdata"))))
	writeTime = proc.time()[3] - writeTime

	params = AddToLog(params, "ProcessXtWXCox.T3", readTime, readSize, writeTime, writeSize)

	return(params)
}


SurvFitCox.BT3 = function(params, pred) {
  if (params$trace) cat(as.character(Sys.time()), "SurvFitCox.BT3\n\n")
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


ComputeResultsCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "ComptueResultsCvox.T3\n\n")
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
	stats$loglik       = c(params$nullLoglikelihood, params$loglikelihood)
	stats$n            = params$n
	stats$nevent       = sum(params$survival$status)
	stats$df           = params$p
	stats$iter         = params$algIterationCounter - 1
	stats$score        = t(params$nullScore) %*% solve(params$nullHessian) %*% params$nullScore
	stats$score        = c(stats$score, 1 - pchisq(stats$score, stats$df))
	stats$method       = "efron"
	stats$lrt          = 2*(stats$loglik[2] - stats$loglik[1])
	stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
	stats$rsquare      = c(1 - exp(-stats$lrt[1]/stats$n),
												 1 - exp(2 * stats$loglik[1] / stats$n))
	stats$wald.test    = t(params$betas) %*% params$xtwx %*% params$betas
	stats$wald.test    = c(stats$wald.test,
												 1 - pchisq(stats$wald.test, stats$df))
	pred = -params$Xbeta
	if (params$survivalInstalled) {
		surv = Surv(params$survival$rank, params$survival$status)
		strat = rep(0, length(surv))
		for (i in 1:length(params$survival$strata)) {
			strat[params$survival$strata[[i]]$start:params$survival$strata[[i]]$end] = i
		}
		results = concordance(surv~pred + strata(strat))
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
		surv   = SurvFitCox.BT3(params, pred)
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

	params = AddToLog(params, "ComputeResultsCox.T3", 0, 0, writeTime, writeSize)
	return(params)
}


GetResultsCox.A3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.A3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats      = stats
	params = AddToLog(params, "GetResultsCox.A3", readTime, readSize, 0, 0)
	return(params)
}


GetResultsCox.B3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "GetResultsCox.B3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats
	params = AddToLog(params, "GetResultsCox.B3", readTime, readSize, 0, 0)
	return(params)
}


####################### REGRESSION BY B ONLY FUNCTIONS #######################

CheckColinearityCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "CheckColinearityCox.B3\n\n")
  # Add in the strata here
	numStrata = length(data$survival$strata)
	X = cbind(matrix(0, nrow = nrow(data$X), ncol = numStrata), data$X)
	for (i in 1:numStrata) {
		X[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end, i] = 1
	}

	XTX = t(X) %*% X
	nrow = nrow(XTX)
	indicies = 1:numStrata
	for (i in (numStrata + 1):nrow) {
		tempIndicies = c(indicies, i)
		if (rcond(XTX[tempIndicies, tempIndicies]) > 10^8 * .Machine$double.eps) {
			indicies = c(indicies, i)
		}
	}
	indicies = indicies[-(1:numStrata)] - numStrata  # Get rid of the Strata

	params$AIndiciesKeep = c()
	params$colnamesA.old = c()
	params$colnamesA     = c()
	params$p1.old        = params$p1
	params$p1            = 0

	Bnames               = params$colnames
	BnamesKeep           = Bnames[indicies]
	params$BIndiciesKeep = indicies
	params$colnames.old  = params$colnames
	params$colnames      = BnamesKeep
	params$p2.old        = params$p2
	params$p2            = length(BnamesKeep)
	params$p             = params$p1 + params$p2

	if (params$p2 == 0) {
		params$failed = TRUE
		params$errorMessage = "Party A has no covariates and all of Party B's covariates are linear."
	}
	params = AddToLog(params, "CheckColinearityCox.B3", 0, 0, 0, 0)

	return(params)
}


ComputeCoxFromSurvival.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCoxFromSurvival.B3\n\n")
  # We have loaded survival previously

	maxIterations = 25
	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "maxiterations.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "maxiterations.rdata"))
	readTime = proc.time()[3] - readTime
	strata = rep(0, nrow(data$X))
	for (i in 1:length(data$survival$strata)) {
		strata[data$survival$strata[[i]]$start:data$survival$strata[[i]]$end] = i
	}

	colnames(data$X) = paste0("V", 1:ncol(data$X))
	f = paste(c("Surv(rank, status) ~ strata(strata)", paste0("V", 1:ncol(data$X))), collapse = " + ")

	error = tryCatch(
		{fit = coxph(as.formula(f),
								 data = data.frame(rank = data$survival$rank,
								 									status = data$survival$status,
								 									strata = strata,
								 									data$X),
								 iter.max = maxIterations)},
		error = function(e) { return(TRUE)},
		warning = function(e) { return(FALSE)}
	)

	if (class(error) == "logical" && error) {
		params$converged = FALSE
		params$failed = TRUE
		params$errorMessage = "Coxph in the survival package failed to converge."
	} else {
		params$converged = TRUE
		if (class(error) == "logical") {
			fit = suppressWarnings(coxph(as.formula(f),
									data = data.frame(rank = data$survival$rank,
																		status = data$survival$status,
																		strata = strata,
																		data$X),
									iter.max = maxIterations))
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
	params = AddToLog(params, "ComputeCoxFromSurvival.B3", readTime, readSize, 0, 0)
	return(params)
}


ComputeCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeCox.B3\n\n")
  n           = params$n
	p2          = params$p
	params$algIterationCounter = 1
	X.betas.old = matrix(0, n, 1)
	X.betas     = matrix(0, n, 1)
	betasB      = matrix(0, p2, 1)
	betasBold   = betasB
	loglikelihood.old = -Inf
	maxIterations = 25
	cutoff        = 10^-8

	readTime = proc.time()[3]
	load(file.path(params$readPath[["T"]], "maxiterations.rdata"))
	readSize = file.size(file.path(params$readPath[["T"]], "maxiterations.rdata"))
	readTime = proc.time()[3] - readTime


	while (params$algIterationCounter <= maxIterations && !params$converged) {
		BeginningIteration(params)
		loglikelihood = 0
		stepSize = 1
		w = exp(X.betas)
		while (max(w) == Inf) {
			cat("Step Halving\n\n")
			X.betas = (X.betas + X.betas.old) * 0.5
			stepSize = stepSize * 0.5
			w = exp(X.betas)
		}
		computeLoglikelihood = TRUE
		while (computeLoglikelihood) {
		  numEvents = sum(data$survival$status)
		  stepCounter = 0
		  pbar = MakeProgressBar1(numEvents, "Loglikelihood")
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
		        pbar = MakeProgressBar2(stepCounter, pbar)
		      }
		    }
		  }
		  if (loglikelihood > loglikelihood.old || stepSize < 0.5^6) {
		    computeLoglikelihood = FALSE
		  } else {
		    cat("Step Halving\n\n")
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
					as.integer(n), as.integer(p2), as.integer(numEvents))

		M = t(data$X) %*% W.XB

		params$XtWX = M
		if (params$algIterationCounter == 1) {
			params$nullHessian = M
		}

		inv = NULL
		tryCatch({inv = solve(M)},
						 error = function(err) { inv = NULL } )
		M = inv
		if (is.null(M)) {
			params$failed = TRUE
			params$errorMessage = "The matrix t(X)WX is singular.  This is probably due to divergence of the coefficients."

			betas = rep(NA, length(params$Bcolnames.old))
			betas[params$BIndiciesKeep] = betasB
			betas = data.frame(betas)
			rownames(betas) = params$Bcolnames.old
			cat("Current Parameters:\n")
			print(betas)
			cat("\n")
			params = AddToLog(params, "ComputeCox.B3", readTime, readSize, 0, 0)
			return(params)
		}

		deltaBeta = M %*% t(data$X) %*% deltal
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
	params = AddToLog(params, "ComputeCox.B3", readTime, readSize, 0, 0)
	return(params)
}


ComputeResultsCox.B3 = function(params, data) {
  if (params$trace) cat(as.character(Sys.time()), "ComputeResultsCox.B3\n\n")
  stats = params$stats
	stats$converged = params$converged
	stats$partyName = params$partyName
	stats$failed    = FALSE

	fitExists = !is.null(params$fit)
	names.old          = c(params$colnamesA.old, params$colnames.old)
	idxA               = params$AIndiciesKeep
	idxB               = params$BIndiciesKeep
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
	stats$lrt          = 2*(stats$loglik[2] - stats$loglik[1])
	stats$lrt          = c(stats$lrt, 1 - pchisq(stats$lrt, stats$df))
	stats$rsquare      = c(1 - exp(-stats$lrt[1]/stats$n),
												 1 - exp(2 * stats$loglik[1] / stats$n))
	stats$wald.test    = c(stats$wald.test,
												 1 - pchisq(stats$wald.test, stats$df))

	params$survival = data$survival
	pred = data$X %*% stats$coefficients[idx]
	stats$survival = data.frame(
		rank   = data$survival$rank,
		status = data$survival$status,
		sorted = data$survival$sorted,
		surv   = SurvFitCox.BT3(params, pred)
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
	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "ComputeResultsCox.B3", 0, 0, writeTime, writeSize)
	return(params)
}


TransferResultsCox.T3 = function(params) {
  if (params$trace) cat(as.character(Sys.time()), "TransferResultsCox.T3\n\n")
  readTime = proc.time()[3]
	load(file.path(params$readPath[["B"]], "stats.rdata"))
	readSize = file.size(file.path(params$readPath[["B"]], "stats.rdata"))
	readTime = proc.time()[3] - readTime
	params$stats = stats

	writeTime = proc.time()[3]
	save(stats, file = file.path(params$writePath, "stats.rdata"))
	writeSize = file.size(file.path(params$writePath, "stats.rdata"))
	writeTime = proc.time()[3] - writeTime
	params = AddToLog(params, "TransferResultsCox.T3", readTime, readSize, writeTime, writeSize)
	return(params)
}


############################## PARENT FUNCTIONS ##############################


PartyAProcess3Cox = function(data,
                             yname                 = NULL,
                             strata                = NULL,
                             mask                  = TRUE,
														 monitorFolder         = NULL,
                             sleepTime             = 10,
                             maxWaitingTime        = 24 * 60 * 60,
														 popmednet             = TRUE,
														 trace                 = FALSE) {
  params = PrepareParams.3p("cox", "A",
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
  data = PrepareDataCox.A3(params, data, yname, strata, mask)
  params = AddToLog(params, "PrepareDataCox.A3", 0, 0, 0, 0)

  if (data$failed) {
  	message = "Error in processing the data for Party A."
  	MakeErrorMessage(params$writePath, message)
  	files = c("errorMessage.rdata")
  	params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
  	return(params$stats)
  }

  params = PrepareParamsCox.A3(params, data)
  files = "pa.rdata"
  params = SendPauseContinue.3p(params, filesT = files, from = "T",
  													 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
  	params$complete = TRUE
  	cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
  	params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
  	return(params$stats)
  }

  params = SendStrataCox.A3(params, data)
  files = "Astrata.rdata"
  params = SendPauseContinue.3p(params, filesT = files, from = "T",
  													 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

  if (file.exists(file.path(params$readPath[["T"]], "stats.rdata"))) {
  	params$algIterationCounter = 1
  	params = GetResultsCox.A3(params)
  	params$converged = params$stats$converged
  	params = SendPauseQuit.3p(params, sleepTime = sleepTime)
  	return(params$stats)
  }

  if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
  	params$complete = TRUE
  	cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
  	params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
  	return(params$stats)
  }

  data = SortDataCox.A3(params, data)
  params = AddToLog(params, "SortDataCox.A3", data$readTime, data$readSize, 0, 0)
  params = PrepareBlocksLinear.A3(params)
  params = GetZCox.A3(params, data)
	files = SeqZW("cz_", length(params$container$filebreak.Z))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetWRCox.A3(params, data)
	params = GetSXACox.A3(params, data)
	files = c("sxa.rdata", "xatxa.rdata", SeqZW("cpr_", length(params$container$filebreak.PR)))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	if (file.exists(file.path(params$readPath[["T"]], "stats.rdata"))) {
		params$algIterationCounter = 1
		params = GetResultsCox.A3(params)
		params$converged = params$stats$converged
		params = SendPauseQuit.3p(params, sleepTime = sleepTime)
		return(params$stats)
	}

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		params$complete = TRUE
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}


	params = UpdateParamsCox.A3(params)
	data = UpdateDataCox.A3(params, data)

	params$algIterationCounter = 1
	repeat {
		params = GetBetaACox.A3(params, data)
		if (params$converged || params$maxIterExceeded) break
		BeginningIteration(params)
		params = GetXABetaACox.A3(params, data)

		files = c("xabetaa.rdata")
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = ComputeXADeltaLCox.A3(params, data)
		files = "tXA_W_XA.rdata"
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = GetXRCox.A3(params, data)
		files = SeqZW("cxr_", length(params$container$filebreak.XR))
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
			cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
			params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
			return(params$stats)
		}
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = GetResultsCox.A3(params)
	params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}

PartyBProcess3Cox = function(data,
														 strata              = NULL,
														 mask                = TRUE,
														 monitorFolder       = NULL,
														 sleepTime           = 10,
														 maxWaitingTime      = 24 * 60 * 60,
														 popmednet           = TRUE,
														 trace               = FALSE) {
	params = PrepareParams.3p("cox", "B",
	                          popmednet = popmednet, trace = trace)
	params = InitializeLog.3p(params)
	params = InitializeStamps.3p(params)
	params = InitializeTrackingTable.3p(params)
	Header(params)
	params = PrepareFolderLinear.B3(params, monitorFolder)
	if (params$failed) {
		cat(params$errorMessage)
		return(invisible(NULL))
	}
	data = PrepareDataCox.B3(params, data, strata, mask)
	params = AddToLog(params, "PrepareDataCox.B3", 0, 0, 0, 0)

	if (data$failed) {
		message = "Error in processing the data for Party B."
		MakeErrorMessage(params$writePath, message)
		files = c("errorMessage.rdata")
		params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = PrepareParamsCox.B3(params, data)
	files = "pb.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		params$complete = TRUE
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	params = SendStrataCox.B3(params, data)

	files = "Bstrata.rdata"
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

	if (file.exists(file.path(params$readPath[["T"]], "transferControl.rdata"))) {
		params$algIterationCounter = 1
		data = SortDataCox.B3(params, data)
		params = AddToLog(params, "SortDataCox.B3", data$readTime, data$readSize, 0, 0)
		params = CheckColinearityCox.B3(params, data)

		if (params$failed) {  # Happens if pB.new == 0
			params$complete = TRUE
			cat("Error:", params$errorMessage, "\n\n")
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = c("errorMessage.rdata")
			params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE)
			return(params$stats)
		}
		data = UpdateDataCox.B3(params, data)
		params = AddToLog(params, "UpdateDataCox.B3", 0, 0, 0, 0)
		if (params$survivalInstalled) {
			params = ComputeCoxFromSurvival.B3(params, data)
		} else {
			params = ComputeCox.B3(params, data)
		}

		if (params$failed) {      # We could get a job_failed here from coefficient explosion
			params$complete = TRUE
			cat("Error:", params$errorMessage, "\n\n")
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = c("errorMessage.rdata")
			params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE)
			return(params$stats)
		}
		params = ComputeResultsCox.B3(params, data)
		stats = params$stats
		save(stats, file = file.path(params$writePath, "stats.rdata"))
		files = c("stats.rdata")
		params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime)
		return(params$stats)
	}


	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		params$complete = TRUE
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}


	data = SortDataCox.B3(params, data)
	params = AddToLog(params, "SortDataCox.B3", data$readTime, data$readSize, 0, 0)
	params = PrepareBlocksLinear.B3(params)
	params = GetRWLinear.B3(params, data)
	params = GetSXBCox.B3(params, data)
	files = c("sxb.rdata", "xbtxb.rdata", SeqZW("crw_", length(params$container$filebreak.RW)))
	params = SendPauseContinue.3p(params, filesT = files, from = "T",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
		params$complete = TRUE
		cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
		return(params$stats)
	}

	if (file.exists(file.path(params$readPath[["T"]], "transferControl.rdata"))) {
		params$algIterationCounter = 1
		params = UpdateParamsCox.B3(params)
		data = UpdateDataCox.B3(params, data)
		params = AddToLog(params, "UpdateDataCox.B3", 0, 0, 0, 0)
		if (params$survivalInstalled) {
			params = ComputeCoxFromSurvival.B3(params, data)
		} else {
			params = ComputeCox.B3(params, data)
		}

		if (params$failed) {      # We could get a job_failed here from coefficient explosion
			params$complete = TRUE
			cat("Error:", params$errorMessage, "\n\n")
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = c("errorMessage.rdata")
			params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime, job_failed = TRUE)
			return(params$stats)
		}
		params = ComputeResultsCox.B3(params, data)
		stats = params$stats
		save(stats, file = file.path(params$writePath, "stats.rdata"))
		files = c("stats.rdata")
		params = SendPauseQuit.3p(params, filesT = files, sleepTime = sleepTime)
		return(params$stats)
	}

	params = UpdateParamsCox.B3(params)
	data = UpdateDataCox.B3(params, data)
	params$algIterationCounter = 1
	repeat {
		params = GetBetaBCox.B3(params)
		if (params$converged || params$maxIterExceeded) break
		BeginningIteration(params)
		params = GetXBBetaBCox.B3(params, data)

		files = "xbbetab.rdata"
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		params = ComputeXBDeltaLCox.B3(params, data)
		files = c("tXB_W_XB.rdata", SeqZW("cCox_", length(params$container$filebreak.Cox)))
		params = SendPauseContinue.3p(params, filesT = files, from = "T",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime, waitForTurn = TRUE)

		if (file.exists(file.path(params$readPath[["T"]], "errorMessage.rdata"))) {
			cat("Error:", ReadErrorMessage(params$readPath[["T"]]), "\n\n")
			params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE, waitForTurn = TRUE)
			return(params$stats)
		}
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = GetResultsCox.B3(params)
	params = SendPauseQuit.3p(params, sleepTime = sleepTime, waitForTurn = TRUE)
	return(params$stats)
}


PartyTProcess3Cox = function(monitorFolder         = NULL,
														 msreqid               = "v_default_0_000",
														 blocksize             = 500,
														 cutoff                = 1e-8,
														 maxIterations         = 25,
														 sleepTime             = 10,
														 maxWaitingTime        = 24 * 60 * 60,
														 popmednet             = TRUE,
														 trace                 = FALSE) {
	Tparams = PrepareParams.3p("cox", "T", msreqid = msreqid,
	                           popmednet = popmednet, trace = trace)
	Tparams = InitializeLog.3p(Tparams)
	Tparams = InitializeStamps.3p(Tparams)
	Tparams = InitializeTrackingTable.3p(Tparams)

	Header(Tparams)
	params   = PrepareFolderLinear.T3(Tparams, monitorFolder)
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

	params   = PrepareParamsCox.T3(params, cutoff, maxIterations)

	if (!params$failed) params = CheckStrataCox.T3(params)

	if (params$failed) {
		cat("Error:", params$errorMessage, "\n\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files,
															 from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	files = "empty.rdata"
	params = SendPauseContinue.3p(params, filesA = files, filesB = files, from = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = PrepareStrataCox.T3(params)

	if (params$p1 == 0) {
		params$algIterationCounter = 1
		MakeTransferMessage(params$writePath)
		files = c("transfercontrol.rdata", "maxiterations.rdata", "survival.rdata")
		params = SendPauseContinue.3p(params, filesB = files, from = "B",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
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
		params = TransferResultsCox.T3(params)
		params$converged = params$stats$converged
		files = "stats.rdata"
		params = SendPauseContinue.3p(params, filesA = files, from = "A",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	params = PrepareBlocksLinear.T3(params, blocksize)

	if (params$failed) {
		cat("Error:", params$errorMessage, "\n\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files,
															 from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	files = c("survival.rdata", "blocks.rdata")
	params = SendPauseContinue.3p(params, filesA = files, from = "A",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = ProcessZLinear.T3(params)
	files = c("survival.rdata", "blocks.rdata", SeqZW("crz_", length(params$container$filebreak.RZ)))
	params = SendPauseContinue.3p(params, filesB = files, from = "B",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = ProcessWLinear.T3(params)
	files = c("p2.rdata", SeqZW("cwr_", length(params$container$filebreak.WR)))
	params = SendPauseContinue.3p(params, filesA = files, from = "A",
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = GetProductsCox.T3(params)
	params = CheckColinearityCox.T3(params)

	if (params$failed) {
		cat("Error:", params$errorMessage, "\n\n")
		MakeErrorMessage(params$writePath, params$errorMessage)
		files = "errorMessage.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files,
															 from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	if (params$p1 == 0) {
		params$algIterationCounter = 1
		MakeTransferMessage(params$writePath)
		files = c("transfercontrol.rdata", "Bindicies.rdata", "maxiterations.rdata", "survival.rdata")
		params = SendPauseContinue.3p(params, filesB = files, from = "B",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
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
		params = TransferResultsCox.T3(params)
		params$converged = params$stats$converged
		files = "stats.rdata"
		params = SendPauseContinue.3p(params, filesA = files, from = "A",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		params = SendPauseQuit.3p(params, sleepTime = sleepTime)
		SummarizeLog.3p(params)
		return(params$stats)
	}

	params = ComputeInitialBetasCox.T3(params)

	filesA = c("Aindicies.rdata", "betasA.rdata", "converged.rdata")
	filesB = c("Bindicies.rdata", "betasB.rdata", "converged.rdata")
	params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB,
														 from = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	while (!params$converged && !params$maxIterExceeded) {
		BeginningIteration(params)
		if (params$algIterationCounter > 1) {
			filesA = c("converged.rdata", "betasA.rdata")
			filesB = c("converged.rdata", "betasB.rdata")
			params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB,
																 from = c("A", "B"),
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
		}

		params = ComputeLogLikelihoodCox.T3(params)
		files = "Xbeta.rdata"
		params = SendPauseContinue.3p(params, filesA = files, filesB = files, from = c("A", "B"),
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ProcessVCox.T3(params)

		files = SeqZW("cvr_", length(params$container$filebreak.VR))
		params = SendPauseContinue.3p(params, filesA = files, from = "A",
															 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

		params = ProcessXtWXCox.T3(params)

		if (params$failed) {
			cat("Error:", params$errorMessage, "\n\n")
			MakeErrorMessage(params$writePath, params$errorMessage)
			files = "errorMessage.rdata"
			params = SendPauseContinue.3p(params, filesA = files, filesB = files,
																 from = c("A", "B"),
																 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)
			params = SendPauseQuit.3p(params, sleepTime = sleepTime, job_failed = TRUE)
			SummarizeLog.3p(params)
			return(params$stats)
		}
		EndingIteration(params)
		params$algIterationCounter = params$algIterationCounter + 1
	}

	params = ComputeResultsCox.T3(params)

	filesA = c("converged.rdata", "betasA.rdata", "stats.rdata")
	filesB = c("converged.rdata", "betasB.rdata", "stats.rdata")
	params = SendPauseContinue.3p(params, filesA = filesA, filesB = filesB,
														 from = c("A", "B"),
														 sleepTime = sleepTime, maxWaitingTime = maxWaitingTime)

	params = SendPauseQuit.3p(params, sleepTime = sleepTime)
	SummarizeLog.3p(params)

	return(invisible(params$stats))
}
