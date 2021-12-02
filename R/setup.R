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

source("D:/Git/monobrainDNAme/R/zzz.R")
source("D:/Git/monobrainDNAme/R/accessory_funcs.R")
source("D:/Git/monobrainDNAme/R/validate_inputs.R")
source("D:/Git/monobrainDNAme/R/convert_to_bed.R")
source("D:/Git/monobrainDNAme/R/feature.select.new.R")
source("D:/Git/monobrainDNAme/R/file_input.R")
source("D:/Git/monobrainDNAme/R/CIBERSORT.R")
#source("D:/Git/monobrainDNAme/R/load_data_sets.R")


# Get liftover chains ---------------------------------------------------------------------------------------------
if (!exists("ah")) ah <- AnnotationHub()
if (!exists("chains")) {
  chains <- query(ah , c("hg19","hg38", "chainfile"))
  chains <- list("hg19ToHg38" = chains[['AH14150']], "hg38ToHg19" = chains[['AH14108']])
}

# Get pan-glioma probe list (hg19) --------------------------------------------------------------------------------
file <- "https://api.gdc.cancer.gov/data/d9027b0c-8d24-47ff-98fb-4066852e3ab3"
filename <- "D:/Git/thesis_data/methSignatures/PanGlioma_MethylationSignatures.xlsx"
if(!file.exists(filename)) downloader::download(file,filename)

glm.probelist <- openxlsx::read.xlsx(filename,sheet=1,startRow = 2,
                                     colNames = FALSE)[,1]
glm.probes.hg19 <- makeGRfromArrayProbes(IlluminaHumanMethylation450kanno.ilmn12.hg19, probes = glm.probelist) 
glm.probes.hg38 <- unlist(rtracklayer::liftOver(glm.probes.hg19,chains[["hg19ToHg38"]]))
glm.probes.hg19.win <- disjointWindow(glm.probes.hg19,window=1000)
glm.probes.hg38.win <- disjointWindow(glm.probes.hg38,window=1000)

probes.glm <- list("raw" = glm.probelist, "hg19" = glm.probes.hg19, "hg38" = glm.probes.hg38,
                   "hg19.win" = glm.probes.hg19.win, "hg38" = glm.probes.hg38.win)

rm(file,filename,glm.probelist,glm.probes.hg19,glm.probes.hg38,glm.probes.hg19.win,glm.probes.hg38.win)

# Generate the probe bedgraphs ------------------------------------------------------------------------------------
file.450k <- "D:/Git/thesis_data/methSignatures/probes.450k.hg19.tsv"
file.27k <- "D:/Git/thesis_data/methSignatures/probes.27k.hg19.tsv"
file.EPIC <- "D:/Git/thesis_data/methSignatures/probes.EPIC.hg19.tsv"

func.19to38 <- function(str) stringr::str_replace(str,"hg19","hg38")

if (!all(file.exists(file.450k,file.27k,file.EPIC))) {
  exportAnnoToBed("IlluminaHumanMethylation450kanno.ilmn12.hg19",file.450k)
  exportAnnoToBed("IlluminaHumanMethylation27kanno.ilmn12.hg19",file.27k)
  exportAnnoToBed("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",file.EPIC)
  liftover_beds(c(file.450k,file.27k,file.EPIC),chain = chains[["hg19ToHg38"]], func = func.19to38)
}

readAndRemove <- function(file,array = "450k") {makeGRangesFromDataFrame(remove_bad_probes(
                                 fread(file,sep="\t",header=T),array),keep.extra.columns = T)}

if (!exists("probes.ill")) {
  probes.ill <- list("i27k.hg19" = readAndRemove(file.27k),
                 "i450k.hg19" = readAndRemove(file.450k) ,
                 "iEPIC.hg19" = readAndRemove(file.EPIC,"EPIC"),
                 "i27k.hg38" = readAndRemove(func.19to38(file.27k)),
                 "i450k.hg38" = readAndRemove(func.19to38(file.450k)),
                 "iEPIC.hg38" = readAndRemove(func.19to38(file.EPIC)),"EPIC")
  
  
  probes.ill[["i450k.hg38.win"]] <- disjointWindow(probes.ill[["i450k.hg38"]],window=1000)
  probes.ill[["i450k.hg38.win.red"]] <- reduceWithMcols(probes.ill[["i450k.hg38.win"]])
}

file.cpgs.hg38 <- "D:/Git/thesis_data/methSignatures/cpgs.hg38.tsv"
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

rm(file.450k, file.27k, file.EPIC, file.cpgs.hg38)

# Read the promoters ----------------------------------------------------------------------------------------------

proms.hg38 <- list(hg38 = fread("D:\\Git\\thesis_data\\methSignatures\\epd_promoters.bed", header = F, col.names = c("chr","start","end"), select  = c(1:3)))

proms.hg38 <- makeGRangesFromDataFrame(proms.hg38)
proms.hg38 <- promoters(proms.hg38, upstream=2000, downstream=200)

probes.prom <- list(proms.hg38 = proms.hg38)

# Setup the cell type lists ---------------------------------------------------------------------------------------
cell_types <- list(all = c("NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","WholeBlood","Granulocyte","Endothelial","Immune","CMP","GMP","cMOP","Ly6C","HSCb","HSCm","MPPb","MPPm","Microglia","Inf.microglia","Inf.macrophage","Treg","ImmMix","Ini.Glioma","Glioma","Neuron","Glia","GBM-IDH","GBM-WT","GBM-imm","CLP","Dendritic"),
                   immune = c("NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","Granulocyte","Treg","Dendritic"),
                   brain = c("Microglia","Inf.microglia","Inf.macrophage","Ini.Glioma","Glioma","Neuron","Glia","GBM-IDH","GBM-WT","GBM-imm"),
                   progenitor = c("CMP","GMP","cMOP","HSCb","HSCm","MPPb","MPPm"),
                   basic = c("CD4Tcell","CD8Tcell","Treg","NKcell","Bcell","Monocyte","Granulocyte","Neutrophil","Eosinophils","Neuron","Glia","Endothelial","Glioma","WholeBlood"))

# Load files that need signatures ---------------------------------------------------------------------------------

source("D:/Git/monobrainDNAme/R/load_data_sets.R")
# Load up various data --------------------------------------------------------------------------------------------
#scm.big <- scMethrix::load_scMethrix("D:/Git/thesis_data/GSE151506/exp/")