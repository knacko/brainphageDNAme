#' Starts an internal stopwatch
#' @details Save the current time to later use for split/lap and overall times
#' @return NULL
#' @export
start_time <- function() {
  assign("time.all", proc.time()["elapsed"], envir=timer.env)
  assign("time.split", proc.time()["elapsed"], envir=timer.env)
  invisible(NULL)
}

#' Outputs the split/lap/iteration time
#' @details Gets the stored elapsed \\code{\link{proc.time}} from either the initial
#' \code{\link{start_time}} or the previous \code{split_time}
#' @return Returns formatted elapsed time since \code{\link{start_time}} or last \code{\link{split_time}}
#' @export
split_time <- function() {
  time <- get("time.split", envir=timer.env)
  if (!is.numeric(time)) {
    warning("start_time() not set. Starting from now.")
    start_time()
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-time
  assign("time.split", proc.time()["elapsed"], envir=timer.env)
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}

#' Stops an internal stopwatch and outputs overall time
#' @details Gets the stored elapsed \code{proc.time()} from initial \code{\link{start_time}} to calculate
#' overall runtime
#' @return Returns formatted elapsed time since \code{\link{start_time}}
#' @export
stop_time <- function() {
  time <- get("time.all", envir=timer.env)
  if (!is.numeric(time)) {
    #warning("start_time() not set")
    return("[unknown time]")
  }
  time <- proc.time()["elapsed"]-get("time.all", envir=timer.env)
  assign("time.split", NA, envir=timer.env)
  assign("time.all", NA, envir=timer.env)
  return(paste0(sprintf(time[[1]], fmt = '%#.2f'),"s"))
}
