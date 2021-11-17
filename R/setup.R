list.of.packages <- c("DMRcate","minfi","data.table","scMethrix","minfiData","minfiDataEPIC","GEOquery","testthat","stringr","IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylation27kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
status <- lapply(list.of.packages, require, character.only = TRUE)
names(status) <- list.of.packages
suppressWarnings(if (!all(status)) status[which(status==FALSE)])
rm(list.of.packages,new.packages)












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