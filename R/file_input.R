# raw.idat.to.scMethrix -------------------------------------------------------------------------------------------
raw.idat.to.scMethrix <- function(GEO, raw_dir, colData, regex, array = c("IlluminaHumanMethylation450k",
                                                                          "IlluminaHumanMethylationEPIC",
                                                                          "IlluminaHumanMethylation27k"), 
                                  target_array = "IlluminaHumanMethylation450k", verbose = TRUE) {
  #- Input Validation --------------------------------------------------------------------------
  .validateType(GEO,"string")
  .validateType(raw_dir,"directory")
  .validateType(colData,"dataframe")
  .validateType(verbose,"boolean")
  
  array <- .validateArg(array,raw.idat.to.scMethrix)
  annotation <- ifelse(array == "IlluminaHumanMethylationEPIC", "ilm10b4.hg19", "ilmn12.hg19")
  
  #- Function code -----------------------------------------------------------------------------
  
  if (verbose) message("Convert .idat to scMethrix  (",start_time(),")")
  
  files <- list.files (raw_dir,full.names = TRUE)
  files <- files[grepl(".*idat.*$", files,ignore.case = TRUE)]

  if (length(files) == 0) {
    
    supp_file = paste0(raw_dir,GEO,"_RAW.tar")
    
    if (!.validateType(supp_file,"file",throw=F)) {
      supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, 
                                             baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
      supp_file <- rownames(supp_file)
    }
    
    supp_files <- untar(tarfile = supp_file, list=TRUE)
    supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
    untar(tarfile = supp_file, exdir = raw_dir, files = supp_files)
    supp_files <- paste0(raw_dir,supp_files)
    #file.remove(supp_file)
    sapply(supp_files, GEOquery::gunzip, overwrite=TRUE)
    files <- list.files(raw_dir, full.names=TRUE, pattern = ".*idat$",ignore.case = T)
    
    row_idx <- !apply(outer(files, row.names(colData), str_detect),1,any)
    # 
    # row_idx <- rep(TRUE,length(files))
    # 
    # for (gsm in row.names(colData)) {
    #   row_idx[which(str_detect(files,gsm))] <- FALSE
    # }
    
    file.remove(files[row_idx])
    files <- files[!row_idx]
    renamer <- str_replace(files,"(.*?GSM.*?)_.*?(_(Red|Grn)\\.idat)","\\1\\2")
    file.rename(files, renamer)
  }
  
  RGset <- minfi::read.metharray.exp(raw_dir, force=TRUE)
  stopifnot(length(setdiff(colnames(RGset),row.names(colData)))==0)
  
  Mset <- preprocessNoob(RGset)
  Rset <- minfi::ratioConvert(Mset)
  minfi::annotation(Rset) = c(array = array, annotation = annotation)
  GRset <-  minfi::mapToGenome(Rset)
  GRset <- convertArray(GRset, outType = target_array, verbose = verbose)
  
  if (verbose) message("Data extracted  (",stop_time(),")")
  
  return(scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData, verbose = verbose))
}

tar.bed.to.scMethrix <- function(GEO, raw_dir, colData, array = c("IlluminaHumanMethylation450k",
                                                                          "IlluminaHumanMethylationEPIC",
                                                                          "IlluminaHumanMethylation27k"), 
                                  target_array = "IlluminaHumanMethylation450k", verbose = TRUE) {
  #- Input Validation --------------------------------------------------------------------------
  .validateType(GEO,"string")
  .validateType(raw_dir,"directory")
  .validateType(colData,"dataframe")
  .validateType(verbose,"boolean")
  
  array <- .validateArg(array,raw.idat.to.scMethrix)
  annotation <- ifelse(array == "IlluminaHumanMethylationEPIC", "ilm10b4.hg19", "ilmn12.hg19")
  
  #- Function code -----------------------------------------------------------------------------
  
  if (verbose) message("Convert .bed to scMethrix  (",start_time(),")")
  
  supp_file = paste0(raw_dir,GEO,"_RAW.tar")
  
  if (!.validateType(supp_file,"file",throw=F)) {
    supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, 
                                             baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
    supp_file <- rownames(supp_file)
  }
    
  files <- untar(tarfile = supp_file, list=TRUE)
  #supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
  untar(tarfile = supp_file, exdir = raw_dir, files = files)
  files <- paste0(raw_dir,files)
  #file.remove(supp_file)
  #sapply(supp_files, GEOquery::gunzip, overwrite=TRUE)
  
  row_idx <- !apply(outer(files, row.names(colData), str_detect),1,any)
  
  file.remove(files[row_idx])
  files <- files[!row_idx]
  renamer <- paste0(raw_dir,str_replace(files,".*?(GSM.*?)_.*(\\.gz$)","\\1\\.bed\\.\\2"))
  file.rename(files, renamer)
  files <- renamer

  sapply(files, function(file) {
    GEOquery::gunzip(files,destname = tools::file_path_sans_ext(file), overwrite=TRUE)
  })

  files <- tools::file_path_sans_ext(files)

  scm <- read_beds(files, ref_cpgs = ref.gen.hg38, colData = colData, chr_idx = 1, start_idx = 2, M_idx = 3, cov_idx=4, h5 = T, h5_dir = exp_dir)
  
  if (verbose) message("Data extracted  (",stop_time(),")")
  
  return(scm)
}


#------------------------------------------------------------------------------------------------------------
#' If data is already normalized and stored within dataTables for each respective sample, this function will parse the downloaded soft file from the GEO and generate bedgraphs from it.
#' 
#' For example, in this repo (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166844), 
#' 
#' 
#' 
#' the ratio values are desired for each respective sample. If you look at a particular sample (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM861635), the data table lists these values. 
#' @param soft string; the soft file
#' @param colData string; the \code{\link{file.path}} for the bedgraphs to be output into
#' @param exclude_id string; which IDs to exclude
#' @param array string; the illumina array
#' @param annotation string; the genome build
#' @param verbose boolean; whether to be chatty
#' @return data.table; the reference CpGs in bedgraph format
#' @export
var.proc.to.scMethrix <- function(GEO, raw_dir, colData, proc_file, array = c("IlluminaHumanMethylation450k",
                                                                              "IlluminaHumanMethylationEPIC",
                                                                              "IlluminaHumanMethylation27k"), 
                                  target_array = "IlluminaHumanMethylation450k", verbose = TRUE) {
  #- Input Validation --------------------------------------------------------------------------
  .validateType(GEO,"string")
  .validateType(raw_dir,"directory")
  .validateType(colData,"dataframe")
  .validateType(verbose,"boolean")
  
  array <- .validateArg(array,raw.idat.to.scMethrix)
  annotation <- ifelse(array == "IlluminaHumanMethylationEPIC", "ilm10b4.hg19", "ilmn12.hg19")
  
  #- Function code -----------------------------------------------------------------------------

  if (verbose) message("Convert var.proc to scMethrix  (",start_time(),")")

  proc_file <- paste0(raw_dir,proc_file)
  
  if (!.validateType(proc_file,"file",throw=F)) {
    proc_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, 
                                           baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = basename(proc_file))
    proc_file <- rownames(proc_file)
  }
  
  header <- fread(file=proc_file,nrows=1,header = F,sep=",")
  cols <- which(!str_detect(header,"Pval"))
  
  meths <- as.data.frame(fread(file = proc_file,select=as.integer(cols)))
  row.names(meths) <- meths[,1]
  meths$V1 <- NULL
  
  colMatch <- match(colnames(meths),colData$ID)
  keep_idx <- !is.na(colMatch)
  
  stopifnot(identical(colnames(meths)[1],colData$ID[colMatch[1]]))
  
  colnames(meths) <- row.names(colData)[colMatch]
  meths = meths[!is.na(colMatch)]
  
  meths <- meths[order(colnames(meths))]
  
  colData$ID <- NULL
  
  Rset = RatioSet(Beta = meths)
  minfi::annotation(Rset) = c(array = array, annotation = annotation)
  GRset <-  minfi::mapToGenome(Rset)
  GRset <- convertArray(GRset, outType = target_array, verbose = verbose)
  
  if (verbose) message("Data extracted  (",stop_time(),")")
  
  return(scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData, verbose = verbose))
}

#------------------------------------------------------------------------------------------------------------
#' If data is already normalized and stored within dataTables for each respective sample, this function will parse the downloaded soft file from the GEO and generate bedgraphs from it.
#' 
#' For example, in this repo (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069), only the ratio values are desired for each respective sample. If you look at a particular sample (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM861635), the data table lists these values. 
#' @param soft string; the soft file
#' @param colData string; the \code{\link{file.path}} for the bedgraphs to be output into
#' @param exclude_id string; which IDs to exclude
#' @param array string; the illumina array
#' @param annotation string; the genome build
#' @param verbose boolean; whether to be chatty
#' @return data.table; the reference CpGs in bedgraph format
#' @export
soft.to.scMethrix <- function(soft, colData = NULL,  array = c("IlluminaHumanMethylation450k",
                                                               "IlluminaHumanMethylationEPIC",
                                                               "IlluminaHumanMethylation27k"), 
                              target_array = "IlluminaHumanMethylation450k", verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(soft,"soft")
  .validateType(colData,"dataframe")
  .validateType(verbose,"boolean")
  
  array <- .validateArg(array,raw.idat.to.scMethrix)
  annotation <- ifelse(array == "IlluminaHumanMethylationEPIC", "ilm10b4.hg19", "ilmn12.hg19")
  
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Converting soft to bed...",start_time())
  
  gsmids <- sort(intersect(row.names(colData),names(soft@gsms)))
  meths = NULL
  
  for(gsmid in gsmids) {
    meth <- soft@gsms[[gsmid]]@dataTable@table
    meth <- subset(meth, select = c("ID_REF","VALUE"))
    colnames(meth) <- c("ID_REF",gsmid)
    if (is.null(meths)) meths <- meth else meths <- merge(meths,meth,by="ID_REF")
  }
  
  ids <- meths$ID_REF
  meths$ID_REF <- NULL
  
  GRset <- makeGenomicRatioSetFromMatrix(as.matrix(meths),rownames = ids, pData = colData, array = array, 
                                         annotation = annotation,what="Beta")
  
  GRset <- convertArray(GRset, outType = target_array, verbose = verbose)
  
  if (verbose) message("Data extracted  (",stop_time(),")")

  return(scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData, verbose = verbose))
}

