#' @title PopMedNet Simulator
#' @description This function is intended to act as a proxy for PopMedNet when
#'   developing code to run on PopMedNet, or when testing out the distributed
#'   regression programs provided with this package.  When used, it is expected
#'   that this function, the analysis center, and the data partner(s) will each
#'   be run in their own instances of R.  The analysis center and the data
#'   partner(s) will share a common directory where the subdirectories dp0, dp1,
#'   ... will be stored.  The directory dp0 is the monitor folder for the
#'   analysis center.  With the exception of 2-party regression, it is assumed
#'   that the data partner using dp1 is also the data partner which holds the
#'   response.  In the case of 2-party regression, the analysis center holds the
#'   response.
#' @param num_party The number of parties (analysis center + data partners)
#'   involved in the multiple regression.  If a data partner is also acting as
#'   the analysis center, then that data partner is only counted once.
#' @param directory The directory where the directories dp0, dp1, ... are
#'   located, which are used by the analysis center and data partner(s) to save
#'   data and receive data from each other.
#' @param verbose logical value.  If \code{TRUE}, prints out information to
#'   document the progression of the computation.
#' @return \code{NULL}.
#' @seealso \code{\link{analysis_center_2_party}}
#'   \code{\link{AnalysisCenter.3Party}} \code{\link{AnalysisCenter.KParty}}
#' @importFrom utils read.csv
#' @export
pmn <- function(num_party, directory = NULL, verbose = TRUE) {
  sleep_time <- 0.5

  hms <- function(time) {
    time <- round(time, 1)
    ss <- time %% 60
    mm <- time %/% 60 %% 60
    hh <- time %/% 3600
    paste0(formatC(hh, width = 2, flag = "0"), ":",
           formatC(mm, width = 2, flag = "0"), ":",
           formatC(ss, width = 4, digits = 1, flag = "0", format = "f"))
  }

  start_time <- proc.time()[3]


  if (missing(num_party) ||
      !is.numeric(num_party) ||
      (round(num_party, 0) != num_party) ||
      (num_party <= 1)) {
    warning(
      paste0("Parameter num_party not specified or not valid.\n",
             "Usage: pmn(num_party, directory)\n",
             "  num_party: number of parties being simulated, at least 2\n",
             "  directory: the directory where the data partner directories",
             " will be located\n")
    )
    return(invisible(NULL))
  }
  party_name <- paste0("dp", 0:(num_party - 1))
  paths <- file.path(directory, party_name)

  if (is.null(directory) || !dir.exists(directory)) {
    warning(
      paste0("Directory ", directory, " not specified or does not exist.\n",
             "Usage: pmn(num_party, directory)\n",
             "  num_party: number of parties being simulated, at least 2\n",
             "  directory: the directory where the data partner directories",
             "will be located\n")
    )
    return(invisible(NULL))
  }

  if (verbose) cat("Running PMN simulation in directory", directory, "\n")
  for (p in paths) {
    if (!dir.exists(p)) {
      dir.create(p)
    }
  }

  write_directory <- rep("", num_party)
  read_directory <- matrix("", num_party, num_party)
  names(write_directory) <- party_name
  colnames(read_directory) <- party_name
  rownames(read_directory) <- party_name
  for (i in 1:num_party) {
    write_directory[i] <- file.path(paths[i],
                                    ifelse(i == 1, "inputfiles", "msoc"))
    if (!dir.exists(write_directory[i])) dir.create(write_directory[i])
    for (j in 1:num_party) {
      read_directory[i, j] <- file.path(paths[j],
                                        ifelse(i == 1, "inputfiles",
                                               paste0("msoc", i - 1)))
      if (!dir.exists(read_directory[i, j])) dir.create(read_directory[i, j])
    }
    read_directory[i, i] <- NA
  }

  source <- c(FALSE, rep(TRUE, num_party - 1))
  quit <- FALSE
  copy_it <- 1
  repeat {
    # Gather all files to copy

    # look for files_done.ok -> remove the file and add everything to the copy
    # list job_started.ok -> just remove the file and add nothing to the copy
    # list

    files_to_send <- NULL
    if (verbose) cat("\nWaiting for", party_name[source], "-",
                     hms(proc.time()[3] - start_time), "\n")
    while (sum(source) > 0) {  # We are waiting parties to write
      for (i in 1:num_party) {
        if (source[i]) {
          if (file.exists(file.path(write_directory[i], "files_done.ok"))) {
            Sys.sleep(sleep_time)
            files <- read.csv(file.path(write_directory[i], "file_list.csv"))
            if (verbose) cat("\n")
            if (verbose) print(files)
            if (verbose) cat("\n")
            if (i == 1) {
              job_done <- files$file_nm == "job_done.ok" &
                files$dp_cd_list == 10 &
                files$transfer_to_site_in == 10
              job_fail <- files$file_nm == "job_fail.ok" &
                files$dp_cd_list == 10 &
                files$transfer_to_site_in == 10
              if (sum(job_done | job_fail) > 0) quit <- TRUE

            }
            idx <- which(files$transfer_to_site_in == 1)
            if (length(idx) > 0) {
              files <- files[idx, c(1, 3)]
              files$dp_cd_list <- files$dp_cd_list + 1
              files$source <- i
              if (is.null(files_to_send)) {
                files_to_send <- files
              } else {
                files_to_send <- rbind(files_to_send, files)
              }
            }
            file.remove(file.path(write_directory[i], "files_done.ok"))
            if (verbose) cat("  Party", party_name[i], "- files_done.ok\n")
            source[i] <- FALSE
            if (sum(source) > 0) {
              if (verbose) cat("Waiting for", party_name[source], "-",
                               hms(proc.time()[3] - start_time), "\n")
            }
          }
        }
      }
      Sys.sleep(sleep_time)
    }

    if (quit) break

    sink   <- rep(FALSE, num_party)

    if (is.null(files_to_send) || nrow(files_to_send) == 0) {
      warning("No files to transfer and job_done.ok not specified.")
      quit <- TRUE
    } else {
      if (verbose) cat("\nCOPYING", paste0("(", copy_it, ") -"),
                       hms(proc.time()[3] - start_time), "\n")
      copy_it <- copy_it + 1
      mark <- matrix(FALSE, num_party, num_party)
      for (i in seq_len(nrow(files_to_send))) {
        origin      <- files_to_send$source[i]
        destination <- files_to_send$dp_cd_list[i]
        fn          <- files_to_send$file_nm[i]
        exists <- file.exists(file.path(write_directory[origin], fn))
        if (exists) {
          if (verbose) cat("  ")
        } else {
          if (verbose) cat("x ")
        }
        size  <- format(file.size(file.path(write_directory[origin], fn)),
                        big.mark = ",", scientific = FALSE)
        space <- paste0(rep(" ", 33 - nchar(size) - nchar(as.character(fn))),
                        collapse = "")

        if (verbose) cat(party_name[origin], "->", party_name[destination], ":",
                         as.character(fn), space, size, "\n")
        file.copy(file.path(write_directory[origin], fn),
                  file.path(read_directory[origin, destination], fn),
                  overwrite = TRUE)
        mark[origin, destination] <- TRUE
        sink[destination] <- TRUE
      }
      for (origin in 1:num_party) {
        for (destination in 1:num_party) {
          if (mark[origin, destination]) {
            save(sink, file = file.path(read_directory[origin, destination],
                                        "files_done.ok"))
          }
        }
      }
    }

    source <- sink
  }

  if (verbose) cat("\nFinished -", hms(proc.time()[3] - start_time), "\n")

  return(invisible(NULL))
}
