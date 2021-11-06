list.of.packages <- c("DMRcate","minfi","data.table","scMethrix","minfiData","minfiDataEPIC","GEOquery")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)


#------------------------------------------------------------------------------------------------------------
#' Converts green and red idat files into bedgraph files via minfi
#' @param in_dir string; the \code{\link{file.path}} of the directory containing the .IDAT files
#' @param out_dir string; the \code{\link{file.path}} for the bedgraphs to be output into
#' @param regex string; the regex to format the samples names with
#' @param verbose boolean; whether to be chatty
#' @return data.table; the reference CpGs in bedgraph format
#' @export
idat.to.bed <- function(in_dir, out_dir, regex = "(.*)", verbose = TRUE) {

  if (verbose) message("Converting idat to bed...",start_time())
  
  RGSet <- minfi::read.metharray.exp(in_dir, force=TRUE)
  MSet <- minfi::preprocessIllumina(RGSet, bg.correct = TRUE, normalize = "controls") 
  
  if (verbose) message("Data extracted  (",split_time(),")")
  
  RSet <- minfi::ratioConvert(MSet)
  
  if (verbose) message("Beta values calculated (",split_time(),")")

  return(ratioset.to.bed(RSet, out_dir = out_dir, regex = regex, verbose = verbose))
}


#------------------------------------------------------------------------------------------------------------
#' If data is stored within dataTables for each respective sample, this function will parse the downloaded soft file from the GEO and generate bedgraphs from it.
#' 
#' For example, in this repo (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069), only the ratio values are desired for each respective sample. If you look at a particular sample (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM861635), the data table lists these values. 
#' 
#' Adapted from https://support.bioconductor.org/p/73941/
#' @param in_file string; the \code{\link{file.path}} of the .SOFT file downloaded from the GSE
#' @param out_dir string; the \code{\link{file.path}} for the bedgraphs to be output into
#' @param regex string; the regex to format the samples names with
#' @param verbose boolean; whether to be chatty
#' @return data.table; the reference CpGs in bedgraph format
#' @export
soft.to.bed <- function(in_file, out_dir, regex = "(.*)", verbose = TRUE) {
  
  if (verbose) message("Converting soft to bed...",start_time())
  
  soft <- GEOquery::getGEO(filename=in_file)
  gsmids <- names(soft@gsms)
  
  meths = NULL
  
  for(gsmid in gsmids) {
    
    meth <- soft@gsms[[gsmid]]@dataTable@table[,1:2]
    colnames(meth) <- c("ID_REF",gsmid)
    if (is.null(meths)) meths <- meth else meths <- merge(meths,meth,by="ID_REF")
  
  }
  
  if (verbose) message("Data extracted  (",split_time(),")")
  
  row.names(meths) <- meths[,"ID_REF"]
  meths$ID_REF <- NULL
  
  RSet = RatioSet(Beta = meths)
  annotation(RSet) = annotation(minfiData::MsetEx)

  return(ratioset.to.bed(RSet, out_dir = out_dir, regex = regex, verbose = verbose))
}

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
  
  GRSet <- minfi::mapToGenome(RSet)
  
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



chain <- "D:/Git/thesis_data/chains/hg19ToHg38.over.chain"

liftover_beds <- function (files, chain) {
  
  for(file in files) {
    
    bed <- data.table::fread(file,col.names = c("chr","start","end","strand","beta"),key=c("chr","start"))
    bed <- makeGRangesFromDataFrame(bed,keep.extra.columns = T)
    
    ch = rtracklayer::import.chain(chain)
    bed <- unlist(rtracklayer::liftOver(bed,ch))
    bed <- subset(as.data.frame(bed), select = -c(width))
    
    data.table::fwrite(bed, paste0(file), append = FALSE, sep = "\t", row.names = FALSE, 
                       col.names = FALSE, quote = FALSE)
    
  }
}