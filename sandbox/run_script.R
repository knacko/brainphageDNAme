list.of.packages <- c("DMRcate","minfi","data.table","scMethrix","minfiData","minfiDataEPIC","GEOquery","testthat","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
  BiocManager::install(new.packages)
}
status <- lapply(list.of.packages, require, character.only = TRUE)
names(status) <- list.of.packages
suppressWarnings(if (!all(status)) status[which(status==FALSE)])
rm(list.of.packages,new.packages)

scm <- subset_scMethrix(scm.big,regions = glm_probes_mani)
scm1 <- bin_scMethrix(scm,regions = glm_probes_mani, bin_size = NULL, batch_size = 200)

scm <- subset_scMethrix(scm.big,regions = TCGA_HG38_GR)
scm2 <- bin_scMethrix(scm,regions = TCGA_HG38_GR, bin_size = NULL, batch_size = 200)


