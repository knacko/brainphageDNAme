#---- standardize.scMethrix --------------------------------------------------------------------------------------
standardize.scMethrix <- function(scm, GEO, src_genome=NULL, out_genome=NULL, regions = NULL) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  #.validateType(chain,c("chain","null"))
  .validateType(regions,c("granges","null"))
  
  #- Function code -----------------------------------------------------------------------------
  colData(scm)$Experiment <- GEO
  metadata(scm)["Experiment"] <- GEO
  # if (str_detect(metadata(scm)$genome,"hg19") && !is.null(chain)) 
  
  if (!is.null(src_genome) && !is.null(out_genome)) {
    if (src_genome == "hg19" && out_genome == "hg38") {
      scm <- liftover_CpGs(scm,get_chains()[["hg19ToHg38"]],target_genome = "hg38")
    } else if (src_genome == "hg38" && out_genome == "hg19" ) {
      scm <- liftover_CpGs(scm,get_chains()[["hg38ToHg19"]],target_genome = "hg19")
    }
  }
  
  if (!is.null(regions)) {
    scm <- subset_scMethrix(scm,regions = regions)
    scm <- convert_scMethrix(scm,type="memory")
    scm <- bin_scMethrix(scm,regions = regions, fill=T, batch_size = 10)
  }
  
  #names(rowRanges(scm)) <- paste0("rid_",1:nrow(scm))
  #if (is_h5(scm)) scm <- convert_HDF5_scMethrix(scm)
  #if (nrow(scm) != length(probe.set)) {message("Wrong rows");browser()}
  return (scm) 
}

#---- is.empty ---------------------------------------------------------------------------------------
#' Checks to see if an object is NULL, NA, zero length, whitespace
#' Adapted from https://github.com/cran/rapportools/blob/master/R/utils.R
#' @param x an object to check its emptiness
#' @param trim trim whitespace? (\code{TRUE} by default)
#' @param ... additional arguments for \code{\link{sapply}}
#' @return integer; 1 if windows, or some number of threads between 1 and parallel::detectCores
is.empty <- function(input, trim = TRUE, ...) {
  if (length(input) <= 1) {
    if (is.null(input))
      return (TRUE)
    if (length(input) == 0)
      return (TRUE)
    if (is.na(input) || is.nan(input))
      return (TRUE)
    if (is.character(input) && nchar(ifelse(trim, trim.space(input), input)) == 0)
      return (TRUE)
    if (is.logical(input) && !isTRUE(input))
      return (TRUE)
    if (is.numeric(input) && input == 0)
      return (TRUE)
    return (FALSE)
  } else
    any(sapply(input, is.empty, trim = trim, ...))
}

#---- mkdirs ----------------------------------------------------------------------------------------------------------
#' Creates multiple directories at the same time
#' @param ... string; a list of directories to create
#' @return nothing
#' @export
#' @examples
mkdirs <- function (...) {
  
  dirs <- list(...)
  for (i in 1:length(dirs)) {dir.create(dirs[[i]], showWarnings = FALSE)}
}

#---- get_sample_name -------------------------------------------------------------------------------------------------
#' Strips the path and extension from a file path, so only file name remains
#' @param s string; a file path for a file
#' @return string; the name of the file without the path or extension
#' @export
#' @examples
get_sample_name = function(s) {
  if (!is.character(s)) stop("Must be a string file path")
  return(tools::file_path_sans_ext(basename(s)))
}

#---- get_timestamp ----------------------------------------------------------------------------------------------------
#' Formats current time into a string
#' @return string; the current time
#' @export
get_timestamp <- function() {
  ts <- Sys.time()
  ts <- str_replace_all(ts,":","-")
  ts <- str_replace_all(ts," ","_")
}


#---- remove_col ------------------------------------------------------------------------------------------------------
#' Remove a sample from an scMethrix object
#' @param scm scMethrix; an scMethrix object
#' @param col string; the name of a sample in the scMethrix object
#' @return scMethrix; the output object with the specified column removed
#' @export
#' @examples
remove_col <- function(scm,col) {
  scm[,!(colnames(scm) %in% col)]
}

#---- getcells --------------------------------------------------------------------------------------------------------
#' Get table of cells for an scMethrix object
#' @param scm scMethrix; the object to the table of cell types
#' @return table; count of cells in the scMethrix object
#' @export
#' @examples
getcells <- function(scm) {
  return(table(colData(scm)$Cell))
}

#---- colbind ---------------------------------------------------------------------------------------------------------
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

#---- generate_heatmap ------------------------------------------------------------------------------------------------
#' Generates a ComplexHeatmap from an scMethrix object
#' @param scm scMethrix; the object to plot
#' @param assay string; the assay in the scMethrix object
#' @param grouping string; the column in colData in which to group into
#' @param n_cpg string; the number of CpGs to plot
#' @param ... Additional arguments to Heatmap()
#' @return
#' @export
#' @examples
generate_heatmap <- function(scm,assay = "score",grouping = NULL, n_cpg = NULL,...) {
  
  if (!is.null(n_cpg)) scm <- reduce_scMethrix(scm,assay=assay,n_cpg = n_cpg,var="rand")

  mat <- as.matrix(get_matrix(scm,assay))
  type <- colData(scm)$Cell
  ha = HeatmapAnnotation(
    df = data.frame(type = type),
    annotation_height = unit(4, "mm")
  )
  
  Heatmap(mat, name = "expression", km = 5, top_annotation = ha, cluster_columns = agnes(t(mat)), 
          show_row_names = FALSE, column_dend_side = "bottom", column_names_side = "top", show_column_names = FALSE,...) 
  
}

#---- %allin% ----------------------------------------------------------------------------------------------------
#' Checks for presence of all elements in another vector#'
#' @param x vector; the elements in which to search for  
#' @param y vector; the vector to search in
#' @return boolean; TRUE if all elements in x are in y, FALSE if not
#' @export
#' @examples
'%allin%' <- function(x,y) {length(setdiff(x,y))==0}

#---- %arein% ----------------------------------------------------------------------------------------------------
#' Returns all elements of x that are in y
#' @param x vector; the elements in which to search for  
#' @param y vector; the vector to search in
#' @return vector; all elements of x that are in y
#' @export
#' @examples
'%arein%' <- function(x,y) {x[!is.na(match(x,y))]}

#---- left_match ------------------------------------------------------------------------------------------------------
#' Checks whether elements of x are in y
#' @param x vector; the elements in which to search for  
#' @param y vector; the vector to search in
#' @return vector; all of the elements in x that appear in y
#' @export
#'
#' @examples
left_match <- function(x,y) {
  !is.na(match(x,y))
}

#---- headless --------------------------------------------------------------------------------------------------------
#' Displays n rows of a matrix with the column and row names removed. Basically head() with only data
#' @param mtx matrix; the matrix to display
#' @param n integer; the number of rows to display
#' @return matrix; the formatted matrix
#' @export
#' @examples
headless <- function(mtx, n=5) {
  mtx <- head(mtx,n)
  colnames(mtx) <- NULL
  row.names(mtx) <- NULL
  mtx
}

#---- headless --------------------------------------------------------------------------------------------------------
#' Displays n rows and n cols of a matrix with the column and row names removed.
#' @param mtx matrix; the matrix to display
#' @param n integer; the number of rows and columns to display
#' @return matrix; the formatted matrix
#' @export
#' @examples
corner <- function(mtx, n=30) {
  mtx <- mtx[1:n,1:n]
  colnames(mtx) <- NULL
  row.names(mtx) <- NULL
  mtx
}

#---- load_ref_cpgs ---------------------------------------------------------------------------------------------------
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

#---- remove_bad_probes -------------------------------------------------------------------------------------------
#' Removes cross-reactive or polymorphic CpG-matching probes. 
#' @param probes character, list, GRanges; Input probes
#' @param array string; Which array to use
#' @return (invisible) List of non-problematic probes 
#' @export
#' @examples
#' probes <- c("cg14817997", "cg00764514", "cg26928153", "cg16269199","cg00756362")
#' remove_bad_probes(probes)
#' exportAnnoToBed(p450k,tempfile())
remove_bad_probes = function(probes, array = c("450K","EPIC")) {
  
  array <- .validateArg(array,remove_bad_probes)
  bad.probes <- unlist(maxprobes::xreactive_probes(array_type = array))
  
  if (.validateType(probes,"list",throws = F)) {
    ids = probes$ID
    ids <- ids %in% bad.probes
    message("Found ", sum(ids), " bad probes.")
    probes = probes[!ids,]
  } else {#if (.validateType(probes,"character",throws = F)){
    ids <- probes %in% bad.probes
    message("Found ", sum(ids), " bad probes.")
    probes <- probes[!ids]  
  }
  
  return(invisible(probes))
}

#---- exportAnnoToBed -------------------------------------------------------------------------------------------------
#' Exports an annotation to bedgraph
#' TODO: Add liftover support directly
#' @param anno The probe annotation package
#' @param file The file name to export with
#' @param ... Additional options for fwrite
#' @return invisible(TRUE)
#' @export
#' @examples
#' p450k <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
#' library(p450k)
#' exportAnnoToBed(p450k,tempfile())
exportAnnoToBed <- function(anno,file,...) {
  probes <- setDT(data.frame(get(data(Locations, list = character(), package = anno))), keep.rownames = TRUE)[]
  colnames(probes) <- c("ID","chr","start","strand")
  setkeyv(probes,c("chr","start"))
  probes[, end := start + 1] 
  setcolorder(probes, c("chr", "start", "end", "strand", "ID"))
  probes[order(chr, start)]
  fwrite(probes,file=file,sep="\t",...)
  invisible(return(TRUE))
}

#---- liftover_beds ---------------------------------------------------------------------------------------------------
#' @param files The files to convert
#' @param chain The chain to use
#' @param func Function to rename each file
#' @param verbose Be chatty
#' @return 
#' @export
#'
#' @examples
#' library(AnnotationHub)
#' <- AnnotationHub()
#' ins <- query(ah , c("hg19","hg38", "chainfile"))
#' 9ToHg38 <- chains[['AH14150']]
#' 
#' func <- function(str) stringr::str_replace(str,"hg19","hg38)
liftover_beds <- function (files, chain, func = NULL, verbose = T) {
  
  if (verbose) message("Starting liftover...")#,start_time())
  
  for(file in files) {
    
    bed <- data.table::fread(file,header = T,key=c("chr","start"))
    bed <- makeGRangesFromDataFrame(bed,keep.extra.columns = T)
    bed <- unlist(rtracklayer::liftOver(bed,chain))
    bed <- subset(as.data.frame(bed), select = -c(width))
    
    if (!is.null(func)) file <- func(file)
    
    data.table::fwrite(bed, file, append = FALSE, sep = "\t", row.names = FALSE, 
                       col.names = TRUE, quote = FALSE)
    
    if (verbose) message("   ",basename(file), " liftover'd in ")#,split_time())
    
  }
  invisible(return(TRUE))
}

#---- makeGRfromArrayProbes ----------------------------------------------------------------------------------
#' Creates a \code{\link{GenomicRanges}} object from Illumina array probes
#' @details Converts probe IDs from array annotation data into genomic positions. 
#' @param array_anno Either the array data package off Bioconductor, or a custom data.frame consisting of probes in bedgraph format
#' @param probes string; A list of desired array probe IDs. If null, it will return all available probes in array_anno
#' @return \code{\link{GenomicRanges}}; Array CpG probes 
#' @examples
#' #From Illumina annotation
#' library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
#' probes <- c("cg14817997", "cg26928153", "cg16269199")
#' makeGRfromArrayProbes(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, probes)
#' 
#' #From a data.frame
#' array  <- list(chr = c("chr1", "chr1", "chr1", "chr1", "chr1"), 
#'                start = c(10524L, 10847L, 10849L, 15864L, 18826L), 
#'                end = c(10526L, 10849L, 10851L, 15866L, 18828L), 
#'                strand = c("-", "+", "+", "-", "-"), 
#'                ID = c("cg14817997", "cg26928153", "cg16269199", "cg13869341", "cg14008030"))
#' 
#' makeGRfromArrayProbes(data.frame(array), probes, id_col="ID")
#' 
#' @export
makeGRfromArrayProbes <- function(array_anno, probes = NULL, id_col = "ID") {
  
  if (typeof(array_anno) == "S4") {
    probe_anno <- minfi::getAnnotation(array_anno)
    probe_anno <- data.frame(chr = probe_anno$chr, start = probe_anno$pos, end = probe_anno$pos+1, 
                             strand = probe_anno$strand, ID = row.names(probe_anno))
    id_col = "ID"
  } else {
    if (!all(c("chr","start","end") %in% colnames(array_anno)))
      stop("Input data.frame must have columns 'chr','start','end'", call. = FALSE)
    
    probe_anno = array_anno
  }
  
  probe_anno <- probe_anno[!is.na(probe_anno$chr),] # Remove invalid values
  probe_anno <- probe_anno[!is.na(probe_anno$start),]
  probe_anno <- probe_anno[!is.na(probe_anno$end),]
  
  if (!is.null(probes)) probe_anno <- probe_anno[which(probe_anno[,id_col] %in% probes),]
  
  probe_anno <- sort(GenomicRanges::makeGRangesFromDataFrame(probe_anno, keep.extra.columns=T))
  
  return(probe_anno)
}

#---- disjointWindow -----------------------------------------------------------------------------------------
#' Creates disjoint (non-overlapping) windows for all regions
#' @details Creates windows for each input range. If multiple ranges overlap, the window will be shifted to fit the overlapping sites with a minimum of windows.  
#' 
#' @param gr GenomicRanges; The ranges to window
#' @param window integer; > 0; The size of the window in bp
#' @param neg.rm boolean; Should negative start positions be switched to 0? Ranges that are truncated will have an output size less than the window
#' @return GenomicRanges; All ranges will be disjoint and width equal to window size
#' @examples
#' gr <- GRanges(
#'    seqnames = Rle(c("chr1"), c(1)),
#'    ranges = IRanges(1, end = 4600))
#' gr
#' 
#' gr <- glm_probes_mani
#' gr <- disjointWindow(gr,window=1000)
#' gr <- disjointWindow(gr,reduce=T)
#' gr <- disjointWindow(gr,reduce=T,keep_col="ID")
#' gr <- reduceWithMcols(gr,keep_col="ID")
#' 
#' @export
disjointWindow <- function(gr, window = 1000, neg.rm = T) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(gr,"GRanges")
  .validateType(window,"integer")
  .validateValue(window,"> 0")
  
  #- Function code -----------------------------------------------------------------------------
  mcol <- data.frame(1:length(gr))
  names(mcol) <- id_col <- rawToChar(as.raw(sample(c(65:90), 8, replace=T)))
  mcols(gr) <- append(mcols(gr),mcol)
  orig <- gr
  n_probes = length(orig)
  
  #gr <- gr+(window/2)
  gr <- resize(gr,window,fix="center")
  
  gr <- reduce(gr)
  
  row_idx <- which(width(gr) > window)
  new_size <- floor(width(gr)/window)*window
  
  for (size in unique(new_size)) {
    gr[which(new_size == size)] <- resize(gr[which(new_size == size)], size ,fix="center")}
  
  gr <- unlist(tile(gr,width=window))
  overlaps <- findOverlaps(orig,gr,type="any")
  orig <- orig[queryHits(overlaps)]
  gr <- gr[subjectHits(overlaps)]
  mcols(gr) <- mcols(orig)
  gr <- gr[!duplicated(mcols(gr)[,id_col])]
  stopifnot(n_probes == length(gr))
  mcols(gr) <- mcols(gr)[, !(names(mcols(gr)) %in% id_col),drop=FALSE]
  
  if (neg.rm) start(gr)[start(gr) < 0] <- 0
  
  return(gr)
}


#' Reduces a GenomicRanges object, but keeps specified mcols() data. The data with be appended and added as a column in the reduced object. Wraps keepMcols().
#' @param gr GenomicRanges; the object to reduce
#' @param keep_col string; the column name in mcols() to retain during reducing
#' @return GenomicRanges; the reduced object with columns for the kept data
#' @export
#' @examples
reduceWithMcols <- function(gr, keep_col = "ID") {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateType(gr,"GRanges")
  .validateType(keep_col,c("null","string"))
  
  if (.validateType(keep_col,"string",throws=F)) {
    if(!keep_col %in% colnames(mcols(gr))) {
      stop(paste0("Column '",keep_col,"' for reduction not in mcols(): ",paste(colnames(mcols(gr)),collapse=", ")))
    }
  }
  
  #- Function code -----------------------------------------------------------------------------
  return(keepMcols(gr,GenomicRanges::reduce,col=keep_col))
}  


#---- keepMcols ---------------------------------------------------------------------------------------------------
#' Does an operation on a GenomicRanges object but keeps mcols()
#' @param gr GenomicRanges; the object to operate on
#' @param op function; the function to do
#' @param col string; the column name in mcols() to retain
#' @return GenomicRanges; the transformed GR object with the retained mcols() data
#' @export
#' @examples
keepMcols <- function(gr,op,col = NULL) {
  
  out <- op(gr)
  
  if (!is.null(col)) {
    overlaps <- as.data.table(findOverlaps(gr,out))
    
    if (length(gr) > length(out)) {
      overlaps[, `:=` (mcols = unlist(mcols(gr)[col]))]
      overlaps <- overlaps[, .(mcols = paste(mcols,collapse=",")), by = .(subjectHits)]
      mcols(out)$ID <- unlist(overlaps[["mcols"]])
      colnames(mcols(out)) = col
    } else {
      mcols(out)[col] <- mcols(gr)[col][overlaps$queryHits,]
    }
  }
  
  return(out)
}





getCellCount <- function(GEO, raw_dir, colData, cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"), 
                         array = c("IlluminaHumanMethylation450k",
                                   "IlluminaHumanMethylationEPIC",
                                   "IlluminaHumanMethylation27k"), verbose = TRUE,...) {
  #- Input Validation --------------------------------------------------------------------------
  .validateType(GEO,"string")
  .validateType(raw_dir,"directory")
  .validateType(colData,"dataframe")
  .validateType(verbose,"boolean")
  
  immune <- c("Bcell","CD4T","CD8T","Eos","Gran","Mono","Neu","NK")
  cellTypes <- sapply(immune, like, vector = cellTypes)
  cellTypes <- immune[colSums(cellTypes) > 0]
  
  if (length(cellTypes == 0)) return(invisible(TRUE))
  
  #- Function code -----------------------------------------------------------------------------
  
  if (verbose) message("Convert .idat to scMethrix  (",start_time(),")")
  
  files <- list.files (raw_dir,full.names = TRUE)
  files <- files[grepl(".*idat.*$", files,ignore.case = TRUE)]
  
  RGset <- minfi::read.metharray.exp(raw_dir, force=TRUE)
  
  est <- estimateCellCounts(RGset, cellTypes = cellTypes,
                            referencePlatform = array,
                            returnAll = FALSE, meanPlot = FALSE, verbose = verbose, ...)
  return(est)
}

#---- start_time ------------------------------------------------------------------------------------------------
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


#---- split_time -------------------------------------------------------------------------------------------------
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


#---- stop_time --------------------------------------------------------------------------------------------------
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
