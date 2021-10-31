list.of.packages <- c("DMRcate","minfi","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages,new.packages)

dir = "D:/Git/sampleData/ArrayData/GSE121483"
idat.to.bed(dir,".*_.*_(.*)")

idat.to.bed <- function(dir, regex = "(.*)", verbose = TRUE) {
  
  if (verbose) message("Extracting data from IDATs...",start_time(),appendLF = F)
  
  RGSet <- minfi::read.metharray.exp(dir, force=TRUE)
  MSet <- minfi::preprocessIllumina(RGSet, bg.correct = TRUE, normalize = "controls") 
  
  if (verbose) message(split_time(),"\nCalculating beta values... ", appendLF = F)
  
  RSet <- minfi::ratioConvert(MSet)
  
  if (verbose) message(split_time(),"\nMapping to genome... ", appendLF = F)
  
  GRSet <- minfi::mapToGenome(RSet)
  
  rrng <- data.table::as.data.table(rowRanges(GRSet))
  rrng[,c("width"):=NULL]
  
  for (i in 1:ncol(colData(GRSet))) {
    
    name <- gsub(regex,"\\1",row.names(colData(GRSet))[i])
    
    if (verbose) message("Exporing #",i,": ",name," ", appendLF = F)
    
    out <- cbind(rrng,as.data.table(getBeta(GRSet)[,i,drop=FALSE])[, lapply(.SD, round, 2)])
    data.table::fwrite(out, paste0(dir,"/",name,".bedgraph"), append = FALSE, sep = "\t", row.names = FALSE, 
                       col.names = FALSE, quote = FALSE)
    
    if (verbose) message(split_time())
    
  }
  
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

