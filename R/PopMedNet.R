pmn = function(numParty, directory = getwd()) {
  sleepTime = 0.5

  hms = function(time) {
    time = round(time, 1)
    ss = time %% 60
    mm = time %/% 60 %% 60
    hh = time %/% 3600
    paste0(formatC(hh, width = 2, flag = "0"), ":",
           formatC(mm, width = 2, flag = "0"), ":",
           formatC(ss, width = 4, digits = 1, flag = "0", format = "f"))
  }

  startTime = proc.time()[3]


  if (missing(numParty) ||
      !is.numeric(numParty) ||
      (round(numParty, 0) != numParty) ||
      (numParty <= 1)) {
    cat("Parameter numParty not specified or not valid.\n")
    cat("Usage: pmn(numParty, directory)\n")
    cat("  numParty: number of parties being simulated, at least 2\n")
    cat("  directory: the directory where the data partner directories will be located\n")
    return(invisible(NULL))
  }
  partyName = paste0("dp", 0:(numParty - 1))
  paths = file.path(directory, partyName)

  if (!dir.exists(directory)) {
    cat("Directory", directory, "does not exist.\n")
    cat("Usage: pmn(numParty, directory)\n")
    cat("  numParty: number of parties being simulated, at least 2\n")
    cat("  directory: the directory where the data partner directories will be located\n")
    return(invisible(NULL))
  }

  cat("Running PMN simulation in directory", directory, "\n")
  for (p in paths) {
    if (!dir.exists(p)) {
      dir.create(p)
    }
  }

  writeDirectory = rep("", numParty)
  readDirectory  = matrix("", numParty, numParty)
  names(writeDirectory) = partyName
  colnames(readDirectory) = partyName
  rownames(readDirectory) = partyName
  for (i in 1:numParty) {
    writeDirectory[i] = file.path(paths[i], ifelse(i == 1, "inputfiles", "msoc"))
    if (!dir.exists(writeDirectory[i])) dir.create(writeDirectory[i])
    for (j in 1:numParty) {
      readDirectory[i, j] = file.path(paths[j], ifelse(i == 1, "inputfiles", paste0("msoc", i - 1)))
      if (!dir.exists(readDirectory[i, j])) dir.create(readDirectory[i, j])
    }
    readDirectory[i, i] = NA
  }

  source = c(FALSE, rep(TRUE, numParty - 1))
  quit = FALSE
  copyit = 1
  repeat {
    # Gather all files to copy

    # look for files_done.ok -> remove the file and add everything to the copy list
    # job_started.ok -> just remove the file and add nothing to the copy list

    filesToSend = NULL
    cat("\nWaiting for", partyName[source], "-", hms(proc.time()[3] - startTime), "\n")
    while (sum(source) > 0) {  # We are waiting parties to write
      for (i in 1:numParty) {
        if (source[i]) {
          # if (file.exists(file.path(writeDirectory[i], "job_started.ok"))) {
          # 	Sys.sleep(sleepTime)
          # 	file.remove(file.path(writeDirectory[i], "job_started.ok"))
          # 	cat("  Party", partyName[i], "- job_started.ok\n")
          # 	source[i] = FALSE
          # 	if (sum(source) > 0) cat("Waiting for", partyName[source], "\n")
          # }
          if (file.exists(file.path(writeDirectory[i], "files_done.ok"))) {
            Sys.sleep(sleepTime)
            files = read.csv(file.path(writeDirectory[i], "file_list.csv"))
            cat("\n")
            print(files)
            cat("\n")
            if (i == 1) {
              jobdone = files$file_nm == "job_done.ok" &
                files$dp_cd_list == 10 &
                files$transfer_to_site_in == 10
              jobfail = files$file_nm == "job_fail.ok" &
                files$dp_cd_list == 10 &
                files$transfer_to_site_in == 10
              if (sum(jobdone | jobfail) > 0) quit = TRUE

            }
            idx = which(files$transfer_to_site_in == 1)
            if (length(idx) > 0) {
              files = files[idx, c(1, 3)]
              files$dp_cd_list = files$dp_cd_list + 1
              files$source = i
              if (is.null(filesToSend)) {
                filesToSend = files
              } else {
                filesToSend = rbind(filesToSend, files)
              }
            }
            # Sys.sleep(sleepTime)
            file.remove(file.path(writeDirectory[i], "files_done.ok"))
            cat("  Party", partyName[i], "- files_done.ok\n")
            source[i] = FALSE
            if (sum(source) > 0) cat("Waiting for", partyName[source], "-", hms(proc.time()[3] - startTime), "\n")
          }
        }
      }
      Sys.sleep(sleepTime)
    }

    if (quit) break

    # cat("Simulating PMN Delay\n")
    # Sys.sleep(runif(n = 1, min = 0, max = 15))

    sink   = rep(FALSE, numParty)

    if (is.null(filesToSend) || nrow(filesToSend) == 0) {
      cat("Error: No files to transfer and job_done.ok not specified.\n")
      quit = TRUE
    } else {
      cat("\nCOPYING", paste0("(", copyit, ") -"),  hms(proc.time()[3] - startTime), "\n")
      copyit = copyit + 1
      mark = matrix(FALSE, numParty, numParty)
      for (i in 1:nrow(filesToSend)) {
        origin      = filesToSend$source[i]
        destination = filesToSend$dp_cd_list[i]
        fn          = filesToSend$file_nm[i]
        exists = file.exists(file.path(writeDirectory[origin], fn))
        if (exists) {
          cat("  ")
        } else {
          cat("X ")
        }
        size = format(file.size(file.path(writeDirectory[origin], fn)), big.mark = ",", scientific = FALSE)
        space = paste0(rep(" ", 33 - nchar(size) - nchar(as.character(fn))), collapse = "")

        cat(partyName[origin], "->", partyName[destination], ":", as.character(fn), space, size, "\n")
        file.copy(file.path(writeDirectory[origin], fn),
                  file.path(readDirectory[origin, destination], fn),
                  overwrite = TRUE)
        mark[origin, destination] = TRUE
        sink[destination] = TRUE
      }
      for (origin in 1:numParty) {
        for (destination in 1:numParty) {
          if (mark[origin, destination]) {
            save(sink, file = file.path(readDirectory[origin, destination], "files_done.ok"))
          }
        }
      }
    }

    source = sink
  }

  cat("\nFinished -", hms(proc.time()[3] - startTime), "\n")

  return(invisible(NULL))
}
