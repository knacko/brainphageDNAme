#' Setups up the environment.
#' 
#' All necessary packages are loaded, then scripts are sourced. Chains are downloaded, probe lists are generated, and cell types are grouped. Two global variables will be created: a list for 'probes', 'cell_types', and 'chains'. Each contains named lists that can be called with other functions. 
#' @param script_dir string; the file path to the directory containing the R scripts
#' @param sig_dir string; the file path to the directory containing the methylation signatures (e.g., Illumina probes, CpG sites)
#' @return Nothing
#' @export
setupEnv <- function(script_dir, sig_dir) {

  list.of.packages <- c("DMRcate","minfi","data.table","scMethrix","minfiData","minfiDataEPIC","GEOquery","testthat","stringr","IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylation27kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19","openxlsx","AnnotationHub","future","ComplexHeatmap","ggplot2","ggforce","filesstrings","tibble","e1071","parallel","preprocessCore","ggpubr","cluster","TCGAbiolinks")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
    BiocManager::install(new.packages)
  }
  status <- lapply(list.of.packages, require, character.only = TRUE)
  names(status) <- list.of.packages
  suppressWarnings(if (!all(status)) status[which(status==FALSE)])
  rm(list.of.packages,new.packages,status)
  
  plan(sequential)
  
  source(paste0(script_dir,"/zzz.R"))
  source(paste0(script_dir,"/accessory_funcs.R"))
  source(paste0(script_dir,"/validate_inputs.R"))
  source(paste0(script_dir,"/raw_data_to_scMethrix.R"))
  source(paste0(script_dir,"/feature_selection.R"))
  source("D:/Git/monobrainDNAme/R/load_data_sets.R")
  
  # Get global variables ---------------------------------------------------------------------------------------------
  chains <<- get_chains()
  probes <<- get_probes()
  cell_types <<- get_cell_types()

}

#---- get_chains -------------------------------------------------------------------------------------------------
#' Retrieves a chain from AnnotationHub for use in LiftOver
#' @return chain; the conversion matrix for CpGs in different genome builds 
#' @export
#'
#' @examples
get_chains <- function() { #chain = c("hg19ToHg38","hg38ToHg19")) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("AnnotationHub")
  
  # chain <- .validateArg(chain,get_chains)
  
  if (!exists("ah")) ah <- AnnotationHub()
  if (!exists("chains")) {
    chains <- query(ah , c("hg19","hg38", "chainfile"))
    chains <- list("hg19ToHg38" = chains[['AH14150']], "hg38ToHg19" = chains[['AH14108']])
    #assign("chains", chains, envir = .GlobalEnv)
  }
  
  return(chains)
  
  # if (is.null(chain)) return(invisible(TRUE))
  # 
  # if (chain == "hg19ToHg38") return(chains[["hg19ToHg38"]])
  # if (chain == "hg38ToHg19") return(chains[["hg38ToHg19"]])
}

#---- get_probes -------------------------------------------------------------------------------------------------
#' Retrieves the master list of probes
#' @return list; a list of [GRanges()] corresponding to sites of interest
#' @export
get_probes <- function() {
  
  if (!exists("probes")) {
    
    file <- "https://api.gdc.cancer.gov/data/d9027b0c-8d24-47ff-98fb-4066852e3ab3"
    filename <- paste0(sig_dir,"/PanGlioma_MethylationSignatures.xlsx")
    if(!file.exists(filename)) downloader::download(file,filename)
    
    glm.probelist <- openxlsx::read.xlsx(filename,sheet=1,startRow = 2,
                                         colNames = FALSE)[,1]
    glm.probes.hg19 <- makeGRfromArrayProbes(IlluminaHumanMethylation450kanno.ilmn12.hg19, probes = glm.probelist) 
    glm.probes.hg38 <- unlist(rtracklayer::liftOver(glm.probes.hg19,chains[["hg19ToHg38"]]))
    glm.probes.hg19.win <- disjointWindow(glm.probes.hg19,window=1000)
    glm.probes.hg38.win <- disjointWindow(glm.probes.hg38,window=1000)
    
    probes <- list("glm.hg19" = glm.probes.hg19, "glm.hg38" = glm.probes.hg38,
                   "glm.win.hg19" = glm.probes.hg19.win, "glm.win.hg38" = glm.probes.hg38.win)
    
    rm(file,filename,glm.probelist,glm.probes.hg19,glm.probes.hg38,glm.probes.hg19.win,glm.probes.hg38.win)
    
    # Generate the probe bedgraphs ------------------------------------------------------------------------------------
    file.450k <- paste0(sig_dir,"/probes.450k.hg19.tsv")
    file.27k <- paste0(sig_dir,"/probes.27k.hg19.tsv")
    file.EPIC <- paste0(sig_dir,"/probes.EPIC.hg19.tsv")
    
    func.19to38 <- function(str) stringr::str_replace(str,"hg19","hg38")
    
    if (!all(file.exists(file.450k,file.27k,file.EPIC))) {
      exportAnnoToBed("IlluminaHumanMethylation450kanno.ilmn12.hg19",file.450k)
      exportAnnoToBed("IlluminaHumanMethylation27kanno.ilmn12.hg19",file.27k)
      exportAnnoToBed("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",file.EPIC)
      liftover_beds(c(file.450k,file.27k,file.EPIC),chain = chains[["hg19ToHg38"]], func = func.19to38)
    }
    
    readAndRemove <- function(file,array = "450k") {makeGRangesFromDataFrame(remove_bad_probes(
      fread(file,sep="\t",header=T),array),keep.extra.columns = T)}
    
    readNotRemove <- function(file,array = "450k") {makeGRangesFromDataFrame(
      fread(file,sep="\t",header=T),keep.extra.columns = T)}
    
    if (!exists("probes.ill")) {
      probes.ill <- list("ill.27k.hg19" = readAndRemove(file.27k),
                         "ill.450k.hg19" = readAndRemove(file.450k),
                         "ill.450k.raw.hg19" = readNotRemove(file.450k),
                         "ill.EPIC.hg19" = readAndRemove(file.EPIC,"EPIC"),
                         "ill.27k.hg38" = readAndRemove(func.19to38(file.27k)),
                         "ill.450k.hg38" = readAndRemove(func.19to38(file.450k)),
                         "ill.EPIC.hg38" = readAndRemove(func.19to38(file.EPIC)),"EPIC")
      
      
      probes.ill[["ill.450k.win.hg38"]] <- disjointWindow(probes.ill[["ill.450k.hg38"]],window=1000)
      probes.ill[["ill.450k.win.red.hg38"]] <- reduceWithMcols(probes.ill[["ill.450k.win.hg38"]])
    }
    
    probes <- append(probes,probes.ill)
    rm(probes.ill)
    
    #---- Get raw CpGs -----------------------------------------------------------------------------------------------------
    file.cpgs.hg38 <- paste0(sig_dir,"/cpgs.hg38.tsv")
    # file.450k <- "D:/Git/thesis_data/methSignatures/bin.cpgs.450k.hg38.tsv"
    
    if (file.exists(file.cpgs.hg38)) {
      cpgs.hg38 <- fread(file.cpgs.hg38)
    } else {
      cpgs.hg38 <- extract_CpGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
      #cpgs.hg38 <- makeGRangesFromDataFrame(ref.gen.hg38)
      fwrite(cpgs.hg38, file = file.cpgs.hg38)
      #ref.gen.hg38.450k <- subsetByOverlaps(gr.hg38,probes.ill[["i450k.hg38.win.red"]])
      #ref.gen.hg38.450k <- as.data.table(ref.gen.hg38)[,-c("width","strand")]
      #setnames(ref.gen.hg38.450k, "seqnames", "chr")
    }
    
    #---- Read the promoters ----------------------------------------------------------------------------------------------
    
    proms.hg38 <- list(hg38 = fread(paste0(sig_dir,"\\epd_promoters.bed"), header = F, col.names = c("chr","start","end"), select  = c(1:3)))
    
    proms.hg38 <- makeGRangesFromDataFrame(proms.hg38)
    proms.hg19 <- rtracklayer::liftOver(proms.hg38,chains[["hg38ToHg19"]])
    proms.hg38 <- promoters(proms.hg38, upstream=2000, downstream=200)
    proms.hg19 <- promoters(proms.hg38, upstream=2000, downstream=200)
    
    probes.prom <- list(proms.hg38 = proms.hg38, proms.hg19 = proms.hg19)
    
    probes <- append(probes,probes.prom)
    
  }
  
  return (probes)
}

#---- get_cell_types -------------------------------------------------------------------------------------------------
#' Retrieves the master list of cell types
#' @return list; a list of strings corresponding to cell types of interest
#' @export
get_cell_types <- function() {
  
  if (!exists("get_cell_types")) {
    cell_types <- list(all = c("NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","WholeBlood","Granulocyte","Endothelial","Immune","CMP","GMP","cMOP","Ly6C","HSCb","HSCm","MPPb","MPPm","Microglia","Inf.microglia","Inf.macrophage","Treg","ImmMix","Ini.Glioma","Glioma","Neuron","Glia","GBM-IDH","GBM-WT","GBM-imm","CLP","Dendritic"),
                       immune = c("NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","Granulocyte","Treg","Dendritic"),
                       brain = c("Microglia","Inf.microglia","Inf.macrophage","Ini.Glioma","Glioma","Neuron","Glia","GBM-IDH","GBM-WT","GBM-imm"),
                       progenitor = c("CMP","GMP","cMOP","HSCb","HSCm","MPPb","MPPm"),
                       basic = c("CD4Tcell","CD8Tcell","Treg","NKcell","Bcell","Monocyte","Granulocyte","Neutrophil","Eosinophils","Neuron","Glia","Endothelial","Glioma","WholeBlood"))
  }
  
  return (cell_types)
}

