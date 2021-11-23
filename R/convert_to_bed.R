#--- remove_bad_probes -----------------------------------------------------------------------------------------------
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

#--- exportAnnoToBed -------------------------------------------------------------------------------------------------
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

# liftover_beds ---------------------------------------------------------------------------------------------------
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

#--- makeGRfromArrayProbes ----------------------------------------------------------------------------------
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

#--- disjointWindow -----------------------------------------------------------------------------------------
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

# disjointTile <- function(gr, n_tiles = NULL, width = NULL) {
#   
#   #- Input Validation --------------------------------------------------------------------------
#   .validateType(gr,"GRanges")
#   .validateType(n_tiles,c("integer","null"))
#   .validateType(width,c("integer","null"))
#   
#   if (!xor(is.null(n_tiles),is.null(width))) {
#     stop("Either n_tiles or width must be filled.")
#   }
#   
#   #- Function code -----------------------------------------------------------------------------
#   
#   if (!is.null(width)) {
#   
#   row_idx <- which(width(gr) > width)
#   new_size <- ceiling(width(gr)/width)*width
#   
#   
#   
#   
#   } else {
#     
#     
#     
#     
#   }
#   
#   
# }
  
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

standardize.scMethrix <- function(scm, chain = NULL, probe.set = NULL) {
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  .validateType(raw_dir,"directory")
  .validateType(chain,"chain")
  .validateType(probe.set,"granges")
  
  #- Function code -----------------------------------------------------------------------------

  if (str_detect(metadata(scm)$genome,"hg19") && !is.null(chain)) scm <- liftover_CpGs(scm,chain,target_genome = "hg38")
  scm <- bin_scMethrix(scm,regions = probe.set, fill=T)
  
  if (nrow(scm) != length(probe.set)) stop("After standardizing, probes are missing in the experiment")

  return (scm) 
}


# raw.idat.to.scMethrix -------------------------------------------------------------------------------------------
getCellCount <- function(GEO, raw_dir, colData, cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"), array = c("IlluminaHumanMethylation450k",
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
  
#' 
#' #------------------------------------------------------------------------------------------------------------
#' #' Converts green and red idat files into bedgraph files via minfi
#' #' @param in_dir string; the \code{\link{file.path}} of the directory containing the .IDAT files
#' #' @param out_dir string; the \code{\link{file.path}} for the bedgraphs to be output into
#' #' @param regex string; the regex to format the samples names with
#' #' @param verbose boolean; whether to be chatty
#' #' @return data.table; the reference CpGs in bedgraph format
#' #' @export
#' idat.to.scMethrix <- function(dir, colData, regex = "(.*)", array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", verbose = TRUE) {
#'   
#'   if (verbose) message("Converting idat to bed...",start_time())
#'   
#'   RGSet <- minfi::read.metharray.exp(dir, force=TRUE)
#'   MSet <- preprocessNoob(RGSet)
#'  # MSet <- minfi::preprocessIllumina(RGSet, bg.correct = TRUE, normalize = "controls") 
#'   
#'   RSet <- minfi::ratioConvert(MSet)
#'   minfi::annotation(Rset) = c(array = array, annotation = annotation)
#'   GRset <- minfi::mapToGenome(Rset)
#'   
#'   if (verbose) message("Data extracted  (",stop_time(),")")
#'   
#'   return(scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData, verbose = verbose))
#' }
#' 
#' raw.idat.to.scMethrix <- function(GEO, raw_dir, colData, array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", verbose = TRUE) {
#'   
#'   supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
#'   supp_file <- rownames(supp_file)
#'   supp_files <- untar(tarfile = supp_file, list=TRUE)
#'   supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
#'   supp_files <- untar(tarfile = supp_file, exdir = raw_dir, files = supp_files)
#'   file.remove(supp_file)
#'   #sapply(supp_files, function(file) GEOquery::gunzip(file, remove=TRUE))
#'   #mix_files <- list.files(raw_dir, full.names=TRUE)
#'   #file.remove(mix_files[which(rowSums(sapply(remove_mix, like, vector = mix_files)) == 1)])
#'   
#'   RGSet <- minfi::read.metharray.exp(raw_dir, force=TRUE)
#'   MSet <- preprocessNoob(RGSet)
#'   
#'   RSet <- minfi::ratioConvert(MSet)
#'   minfi::annotation(Rset) = c(array = array, annotation = annotation)
#'   
#'   if (verbose) message("Data extracted  (",stop_time(),")")
#'   
#'   return(scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData, verbose = verbose))
#' }




#------------------------------------------------------------------------------------------------------------
#' If data is stored within dataTables for each respective sample, this function will parse the downloaded soft file from the GEO and generate bedgraphs from it.
#' 
#' For example, in this repo (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069), only the ratio values are desired for each respective sample. If you look at a particular sample (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM861635), the data table lists these values. 
#' 
#' Adapted from https://support.bioconductor.org/p/73941/
#' @param soft string; the \code{\link{file.path}} of the .SOFT file downloaded from the GSE
#' @param out_dir string; the \code{\link{file.path}} for the bedgraphs to be output into
#' @param regex string; the regex to format the samples names with
#' @param verbose boolean; whether to be chatty
#' @return data.table; the reference CpGs in bedgraph format
#' @export
# soft.to.bed.old <- function(soft, out_dir, regex = "(.*)", verbose = TRUE) {
#   
#   if (verbose) message("Converting soft to bed...",start_time())
#   
#   soft <- GEOquery::getGEO(filename=in_file)
#   gsmids <- names(soft@gsms)
#   
#   meths = NULL
#   
#   for(gsmid in gsmids) {
#     
#     meth <- soft@gsms[[gsmid]]@dataTable@table[,1:2]
#     colnames(meth) <- c("ID_REF",gsmid)
#     if (is.null(meths)) meths <- meth else meths <- merge(meths,meth,by="ID_REF")
#     
#   }
#   
#   if (verbose) message("Data extracted  (",split_time(),")")
#   
#   row.names(meths) <- meths[,"ID_REF"]
#   meths$ID_REF <- NULL
#   
#   RSet = RatioSet(Beta = meths)
# 
# 
#   return(ratioset.to.bed(RSet, out_dir = out_dir, regex = regex, verbose = verbose))
# }



# ## 
# signal.to.bed <- function(in_dir, out_dir, regex = "(.*)", verbose = TRUE) {
#   
#   if (verbose) message("Converting signal intensities to bed...",start_time())
#   
#   files = list.files(in_dir, full.names=TRUE)
#   
#   meths = NULL
#   
#   for(file in files) {
#     meth <- data.table::fread(file,sep="\t",header=T,skip="ID",col.names=c("cpg","beta"), key = "cpg", select = c(1:2))
#     if (is.null(meths)) meths <- meth else meths <- merge(meths,meth,by="cpg")
#   }
#   
# 
#   
#   RatioSet(Beta = betaM)
# }

#------------------------------------------------------------------------------------------------------------
#' Helper function for idat.to.bed and soft.to.bed
#' 
#' Converts a minfi::ratioset into bedgraph files
#' 
#' @param RSet RatioSet; the RatioSet to convert into bedgraph files
#' @inheritParams idat.to.bed
#' @return data.table; the reference CpGs in bedgraph format
#' @export
ratioset.to.bed <- function(RSet, out_dir, regex = "(.*)", verbose = TRUE) {
  
  #if (verbose) message("Converting ratioset to bed...",start_time())
  
  
  
  rrng <- data.table::as.data.table(rowRanges(GRSet))
  rrng[,c("width"):=NULL]
  
  if (verbose) message("Mapped to genome (",split_time(),")")
  
  for (i in 1:nrow(colData(GRSet))) {
    
    name <- gsub(regex,"\\1",row.names(colData(GRSet))[i])
    
    out <- cbind(rrng,as.data.table(getBeta(GRSet)[,i,drop=FALSE])[, lapply(.SD, round, 2)])
    data.table::fwrite(out, paste0(out_dir,name,".bedgraph"), append = FALSE, sep = "\t", row.names = FALSE, 
                       col.names = FALSE, quote = FALSE)
    
    if (verbose) message("Exported #",i,": ",name," (",split_time(),")")
    
  }
  
  if (verbose) message("All bed files exported (",split_time(),")")
  
  return(invisible(TRUE))
}



# 
# 
# array_to_granges <- function(anno, probes = NULL, chain = NULL, window = 0) {
# 
#   anno <- minfi::getAnnotation(anno)
#   
#   if (!is.null(probes)) {
#     row_idx <- which(row.names(anno) %in% probes)
#     if (length(row_idx) != length(probes)) message(length(probes)-length(row_idx), " probes not found in the annotation.")
#     anno <- anno[row_idx,]
#   }
#   
#   anno <- data.table(chr = anno$chr, start = anno$pos, end = anno$pos+1, strand = anno$strand, probeID = row.names(anno))
#   anno <- GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns=T)
# 
#   if (!is.null(chain)) anno <- unlist(rtracklayer::liftOver(anno,chain)) # This is a one-to-many conversion. Flank if using this.
#   if (window != 0) anno <- reduceKeepMcols(anno + (window/2))
# 
#   return(anno)
# }
# 
# manifest_to_granges <- function(mani, probes = NULL, chain = NULL, window = 0) {
# 
#   if (!is.null(probes)) {
#     row_idx <- which(row.names(anno) %in% probes)
#     if (length(row_idx) != length(probes)) message(length(probes)-length(row_idx), " probes not found in the annotation.")
#     anno <- anno[row_idx,]
#   }
#   
#   anno <- data.table(chr = anno$chr, start = anno$pos, end = anno$pos+1, strand = anno$strand, probeID = row.names(anno))
#   anno <- GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns=T)
#   
#   if (!is.null(chain)) anno <- unlist(rtracklayer::liftOver(anno,chain)) # This is a one-to-many conversion. Flank if using this.
#   if (window != 0) anno <- reduceKeepMcols(anno + (window/2))
#   
#   return(anno)
# }



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
