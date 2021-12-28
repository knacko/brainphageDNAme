#--- .validateArg -------------------------------------------------------------------------------------------
#' Validates arguments. Allows partial matching.
#' @details Check the parent function input arguments to see whether the inputted value is part of the set. Will return a formatted error message with the incorrect variable name and all the acceptable inputs.
#' 
#' To use, it can be called as such:
#' 
#' genericFunc <- function(values = c("apple","orange","banana")) {
#' 
#'    values <- .validateArg(values)
#' 
#' }
#' 
#' If the argument for values is acceptable (e.g. "apple" or "ban"), it will return the matched string.
#' 
#' This can be used with pipes, but will give erroneous names for the input variable
#' 
#' @param parent closure; the parent function in which to check the input
#' @param arg variable; the variable in which to check
#' @param ignore.case boolean; ignores case of the choices
#' @return arg, if the value is in the function definition.
.validateArg <- function(arg, parent = NULL, ignore.case = T, partial.match = T) {
  
  #.validateType(ignore.case,"boolean")
  
  #- Function code -----------------------------------------------------------------------------
  if (is.null(parent)) {
    parent <- deparse(sys.calls()[[sys.nframe()-1]])
    parent <- unlist(strsplit(parent, "[(]"))[[1]]
  }
  
  name = substitute(arg)
  choices <- eval(formals(parent)[[name]])
  
  arg <- head(arg,1)
  
  if (partial.match) {
    if (ignore.case) {
      m <- grepl(tolower(arg), tolower(choices), fixed = TRUE)
    } else {
      m <- grepl(arg, choices, fixed = TRUE)
    }
  } else {
    if (ignore.case) {
      m <- choices %in% arg
    } else {
      m <- tolower(choices) %in% tolower(arg)
    }
  }
  
  if (sum(m) != 1) stop(paste0("Invalid arg input for '",paste(name),"'. Found: '",arg,"'; Must match one of: '",
                               paste0(eval(formals(parent)[[name]]), collapse="', '"),"'"), call. = FALSE)
  
  return(choices[which(m)])
}

#--- .validateAssay -----------------------------------------------------------------------------------------
#' Validates an assay is in the object. Allows partial pattern.
#' @details Check the assays in an scMethrix object and partial matches 
#' @param scm scMethrix; the experiment object
#' @param assay string; the name of the assay
#' @param check.absent boolean; Checks if the assay is present
#' @return string or boolean; if \code{check.absent == T}, the name of the matched assay or error if it doesn't exist. If \code{check_assay == F}, the boolean value for if the assay exists in the experiment
.validateAssay <- function(scm = NULL, assay = NULL, check.absent = F) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  .validateType(assay,"string")
  .validateType(check.absent, "boolean")
  
  #- Function code -----------------------------------------------------------------------------  
  if (!check.absent) {
    assay <- tryCatch(
      match.arg(arg = assay, choices = SummarizedExperiment::assayNames(scm)),
      error=function(cond) 
        stop(paste0("Invalid assay. No assay named '",assay,"' found in the experiment"), call. = FALSE)
    )
    return(invisible(assay))
  } else {
    return(!(assay %in% SummarizedExperiment::assayNames(scm)))
  }
}

#--- .validateType ------------------------------------------------------------------------------------------
.validateType <- function(input = NULL, type=c("Integer","Numeric","Character","String","Boolean","Logical","Vector",
                                               "List","File","Directory","GRanges","GenomicRanges","Function","Null",
                                               "NA","Dataframe","DF","S4","Distance","Chain","Soft"), throws=T, recursive_sub = NULL) {
  
  #- Input Validation --------------------------------------------------------------------------
  if (length(type) == length(eval(formals(.validateType)[["type"]]))) {
    stop("No valid type specified.")
  }
  
  types <- sapply(type,function(type) .validateArg(type,.validateType))
  
  if (is.null(recursive_sub)) recursive_sub = gsub('"', "'", deparse(substitute(input)))
  
  #- Function code -----------------------------------------------------------------------------
  valid <- F

  # Check for list structures
  if (length(input) > 1) {
    for (type in types) {
      if(type == "List"){
        valid <- is.list(input)
      } else if(type == "GRanges" | type == "GenomicRanges"){
        valid <- is(input, "GRanges")
      } else if (type == "Dataframe" | type == "DF") {
        valid <- is.data.frame(input)
      } else if (type == "Distance") {
        valid <- is(input,"dist")
      } else if (type == "Chain") {
        valid <- is(input,"Chain")
      }
    }
    if (valid) return(invisible(TRUE))
  }
  
  if (length(input) <= 1) {
    for (type in types) {
      if (type == "Null") {
        valid <- is.null(input)
      }  else if (type == "NA") {
        valid <- is.na(input)
      } else if (type == "Integer") {
        if(is.numeric(input)) valid <- (input == round(input))
      } else if (type == "Numeric") {
        valid = is.numeric(input)
      } else if (type == "Character") {
        valid = is.character(input) && nchar(input)==1
      } else if (type == "String") {
        valid = is.character(input)
      } else if(type == "Boolean" | type == "Logical"){
        valid <- is.logical(input)
      } else if(type == "List"){
        valid <- is.list(input)
      } else if (type == "File") {
        valid <- file_test("-f", input)
      } else if (type == "Directory") {
        valid <- file_test("-d", input)
      } else if(type == "GRanges" | type == "GenomicRanges"){
        valid <- is(input, "GRanges")
      } else if (type == "Function") {
        valid <- is.function(input)
      } else if (type == "Dataframe" | type == "DF") {
        valid <- is.data.frame(input)
      } else if (type == "S4") {
        valid <- isS4(input)
      } else if (type == "Distance") {
        valid <- is(input,"dist")
      } else if (type == "Soft") {
        valid <- (class(input) == "GSE")
      } else {
        stop("Invalid type with '",type,"'. This type is not supported for validation.", call. = FALSE)
      }
      
      if (valid) {
        break
      } else if (type == types[length(types)]) {
        if (throws) {
          stop("Invalid type input for '",recursive_sub,"'. Must be of type: '",
               paste0(types, collapse="', '"),"'", call. = FALSE)
        } else {
          return(invisible(FALSE)) 
        }
      }
    } 
  } else {valid = any(sapply(input, .validateType, type = types, throws = throws, recursive_sub = recursive_sub))}
  
  return(invisible(valid))
}

#--- .validateExp -------------------------------------------------------------------------------------------
#' Validates to see if object is a proper scMethrix object
#' @param scm scMethrix; the experiment object to test. Can be either an in-memory object, a file path, or a directory path.
#' @param throws boolean; whether to throw an error on a missing experiment. Will return FALSE on missing otherwise.
#' @param verbose boolean; be chatty
#' @return invisible(TRUE), if the object is valid. Error or FALSE if not.
.validateExp <- function(scm = NULL, throws = TRUE, verbose = FALSE) {

  if (!is(scm, "scMethrix")) {
  
    if (!.validateType(scm,c("file","directory"), throws = F)) {
      
      if (!throws)  return(invisible(FALSE))
        
      stop(paste0("Invalid scMethrix object supplied for '",substitute(scm),
                    "'. Input must either be type 'scMethrix', or a 'file', or 'directory' for",
                    " an experiment save with save_scMethrix()."), call. = FALSE)
    }
  
    tryCatch(
      expr = {
        scm <- load_scMethrix(scm, verbose = verbose)
      },
      error = function(e){ 
        if (throws) {stop("Invalid scMethrix object found at ",scm, call. = FALSE)
        } else {return(invisible(FALSE))}
      }
    )
  }
   
  return(invisible(TRUE)) 
}

#--- .validateValue -----------------------------------------------------------------------------------------
#' Validates numeric values based on some experession
#' @param value numeric; the value to test
#' @param ... string; the expressions to test
#' @return invisible(TRUE), if the object is valid. Error if not.
.validateValue <- function(value,...) {
  
  if (!is.null(value) && !is.na(value)) {
    
    if (!.validateType(value,c("numeric","integer"),throws=F))
      stop("Invalid value for '",substitute(value),"'. Must be of 'numeric' or 'integer' type.")
    
    for (condition in list(...)) {
      if (!(eval(parse(text=paste0(value,condition))))) {
        stop ("Invalid value: '",substitute(value)," = ",value,"'. Must fit condition: ",substitute(value)," ",condition, call. = FALSE)
      }
    }
  }
  return(invisible(TRUE))
}

#--- .validateThreads ---------------------------------------------------------------------------------------
#' Validates the number of threads for the session. Windows can only support one thread
#' @param n_threads numeric; the number of threads
#' @return integer; 1 if windows, or some number of threads between 1 and parallel::detectCores
.validateThreads <- function(n_threads) {
  
  .validateType(n_threads,"integer")
  
  # if (grepl("Windows", Sys.getenv("OS"))) {
  #   if (n_threads > 1) warning("Invalid threads. Parallel processing is not enabled for non-POSIX system (e.g., Windows). ")
  #   return(0)
  # } 
  
  return(max(min(parallel::detectCores(),n_threads),1))
}