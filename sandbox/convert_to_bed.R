#---- Get hg19 -> hg38 chain
library(AnnotationHub)
ah <- AnnotationHub()
chains <- query(ah , c("hg19","hg38", "chainfile"))
hg19ToHg38 <- chains[['AH14150']]

#--- Get Illumina probes
library("IlluminaHumanMethylation27kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

array_to_granges <- function(anno, probes = NULL, chain = NULL, window = 0) {

  anno <- minfi::getAnnotation(anno)
  
  if (!is.null(probes)) {
    row_idx <- which(row.names(anno) %in% probes)
    if (length(row_idx) != length(probes)) message(length(probes)-length(row_idx), " probes not found in the annotation.")
    anno <- anno[row_idx,]
  }
  
  anno <- data.table(chr = anno$chr, start = anno$pos, end = anno$pos+1, strand = anno$strand, probeID = row.names(anno))
  anno <- GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns=T)

  if (!is.null(chain)) anno <- unlist(rtracklayer::liftOver(anno,chain)) # This is a one-to-many conversion. Flank if using this.
  if (window != 0) anno <- reduceKeepMcols(anno + (window/2))

  return(anno)
}

manifest_to_granges <- function(mani, probes = NULL, chain = NULL, window = 0) {

  if (!is.null(probes)) {
    row_idx <- which(row.names(anno) %in% probes)
    if (length(row_idx) != length(probes)) message(length(probes)-length(row_idx), " probes not found in the annotation.")
    anno <- anno[row_idx,]
  }
  
  anno <- data.table(chr = anno$chr, start = anno$pos, end = anno$pos+1, strand = anno$strand, probeID = row.names(anno))
  anno <- GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns=T)
  
  if (!is.null(chain)) anno <- unlist(rtracklayer::liftOver(anno,chain)) # This is a one-to-many conversion. Flank if using this.
  if (window != 0) anno <- reduceKeepMcols(anno + (window/2))
  
  return(anno)
}


window <- 0

# anno27k  <- array_to_granges(anno = IlluminaHumanMethylation27kanno.ilmn12.hg19,   hg19ToHg38,window)
# anno450k <- array_to_granges(anno = IlluminaHumanMethylation450kanno.ilmn12.hg19,  hg19ToHg38,window)
# annoEPIC <- array_to_granges(anno = IlluminaHumanMethylationEPICanno.ilm10b4.hg19, hg19ToHg38,window)

#--- Get glioma probes
library(openxlsx)

file <- "https://api.gdc.cancer.gov/data/d9027b0c-8d24-47ff-98fb-4066852e3ab3"
filename <- "PanGlioma_MethylationSignatures.xlsx"
if(!file.exists(filename)) downloader::download(file,filename)

glm_probelist <- openxlsx::read.xlsx(filename,sheet=1,startRow = 2,
                              colNames = FALSE)[,1]
                              
glm_probes <- array_to_granges(anno = IlluminaHumanMethylation450kanno.ilmn12.hg19, probes = glm_probelist, chain = hg19ToHg38, window = 0)
glm_probes_anno <- glm_probes

## Or

window = 1000
 
ill_glm_probelist <- fread("D:/Git/thesis_data/methSignatures/TCGA_HG38.csv",sep=",",header=T)
colnames(ill_glm_probelist) <- c("chr","start","end","strand","ID")
ill_glm_probelist <- ill_glm_probelist[which(ill_glm_probelist$ID %in% glm_probelist)]
ill_glm_probelist$strand <- "*"

glm_probes_mani <- sort(makeGRangesFromDataFrame(ill_glm_probelist, keep.extra.columns=T))
glm_probes_mani <- resize(glm_probes_mani,1000,fix="center")
glm_probes_mani <- keepMcols(glm_probes_mani,reduce,col="ID")

row_idx <- which(width(glm_probes_mani) > window)

expandwindow <- function(idx,glm_probes_mani,window){
  resize(glm_probes_mani[idx], floor(width(glm_probes_mani[idx])/window)*window,fix="center")}

glm_probes_mani[row_idx] <- unlist(GRangesList(lapply(row_idx,function(idx) expandwindow(idx,glm_probes_mani,window))))
glm_probes_mani <- keepMcols(glm_probes_mani,
  function(x) {unlist(slidingWindows(x,width=window,step=window))},col="ID")

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