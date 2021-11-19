mkdirs <- function (...) {

  dirs <- list(...)
  for (i in 1:length(dirs)) {dir.create(dirs[[i]], showWarnings = FALSE)}
}

get_sample_name = function(s) {
  if (!is.character(s)) stop("Must be a string file path")
  return(tools::file_path_sans_ext(basename(s)))
}

#' A faster version of cbind when trying to combine lists of data.tables
#' @param ... A list of data.tables with identical # of rows
#' @return data.table; the cbinded output 
#' @export
colbind = function(...) {
  setDT(
    unlist(..., recursive = FALSE),
    check.names = TRUE
  )[]
}




corner <- function(mtx, n=30) {
  mtx <- mtx[1:n,1:n]
  colnames(mtx) <- NULL
  row.names(mtx) <- NULL
  mtx
}

headless <- function(mtx, n=5) {
  mtx <- mtx[1:n,]
  colnames(mtx) <- NULL
  row.names(mtx) <- NULL
  mtx
}


#------------------------------------------------------------------------------------------------------------
#' Loads reference CpGs from files. Will download and save them if the files are not present
#' @param dir string; the \code{\link{file.path}} of the directory to save/read the files
#' @param genome string; the desired genome. This should be available from BSgenome
#' @return data.table; the reference CpGs in bedgraph format
#' @export
load_ref_cpgs = function(dir, genome = c("BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10")) {  # BSgenome::available.genomes()) {
  
  genome = .validateArg(genome,load_ref_cpgs)
  
  name = gsub(".*\\.(.*)$","\\1",genome)
  
  cpg_file = paste0(dir,name,"_cpgs.rds")
  if (!file.exists(cpg_file)) {
    saveRDS(scMethrix::extract_CpGs(ref_genome = genome), file = cpg_file)
  }
  
  return(readRDS(cpg_file))
}

#' Starts an internal stopwatch
#' @details Save the current time to later use for split/lap and overall times
#' @return NULL
#' @export
start_time <- function() {
  
  return(NULL)
  
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
  
  return("5")
  
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
  
  return("5")
  
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
