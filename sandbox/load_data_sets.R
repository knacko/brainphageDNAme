home_dir <- "D:/Git/thesis_data/"
cpg_dir <- "D:/Git/sampleData/ref_cpgs/"
chain_dir <- "D:/Git/thesis_data/chains/"
probe.set <- probes.ill[["i450k.hg38.win.red"]]

cell_types <- c("Neutrophil","NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","WholeBlood","Granulocyte","Eosinophil","Endothelial","Immune","CMP","GMP","cMOP","Ly6C","HSCb","HSCm","MPPb","MPPm","Microglia","Inf.microglia","Inf.macrophage","Treg","ImmMix","Ini.Glioma","Glioma","Neuron","Glia","GBM-IDH","GBM-WT","GBM-imm")

#------------------------------------------------------------------------------------------------------------
# Types: 
# Paper: 
# GEO link: 
# Citation: 
# Genome: hg19
# Platform: 
GEO = ""
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (F) if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CONTROL")] <- "Neuron"
  cell[str_detect(cell, "myeloid")] <- "CMP"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(all(cell %in% cell_types))
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}

#--- GSE35069 -----------------------------------------------------------------------------------------------
# GEO: GSE35069
# Types: CD4+ T cells, CD8+ T cells, CD56+ NK cells, CD19+ B cells, CD14+ monocytes, neutrophils, and eosinophils
# Paper: https://pubmed.ncbi.nlm.nih.gov/22848472/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069
# Citation: Reinius LE, Acevedo N, Joerink M, Pershagen G et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PLoS One 2012;7(7):e41361. PMID: 22848472
# Genome: hg19
# Platform: Illumina 450k
GEO = "GSE35069"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "blood")] <- "WholeBlood"
  cell[str_detect(cell, "Gran")] <- "Granulocyte"
  cell[str_detect(cell, "CD4")] <- "CD4Tcell"
  cell[str_detect(cell, "CD8")] <- "CD8Tcell"
  cell[str_detect(cell, "CD14")] <- "Monocyte"
  cell[str_detect(cell, "CD19")] <- "Bcell"
  cell[str_detect(cell, "CD56")] <- "NKcell"
  cell[str_detect(cell, "Neu")] <- "Neutrophil"
  cell[str_detect(cell, "Eos")] <- "Eosinophil"
  remove_idx <- which(str_detect(cell, "PBMC"))
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE41826 -----------------------------------------------------------------------------------------------
# GEO: GSE41826
# Types: Neuron, glia
# Paper: https://pubmed.ncbi.nlm.nih.gov/23426267/
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41826
# Citation: Guintivano J, Aryee MJ, Kaminsky ZA. A cell epigenotype specific model for the correction of brain cellular heterogeneity bias and its application to age, brain region and major depression. Epigenetics 2013 Mar;8(3):290-302. PMID: 23426267
# Genome: hg19
# Platform: 450k
GEO = "GSE41826"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "G")] <- "Glia"
  cell[str_detect(cell, "N")] <- "Neuron"
  status <- sapply(names(soft@gsms), function(gsm)   soft@gsms[[gsm]]@header$characteristics_ch1[[2]])
  cell[!str_detect(status, "Control")] <- "###-"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE49618 -----------------------------------------------------------------------------------------------
# GEO: GSE49618
# Types: Monocyte, Bcell, HSC
# Paper: https://pubmed.ncbi.nlm.nih.gov/23426267/
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49618
# Citation: Guintivano J, Aryee MJ, Kaminsky ZA. A cell epigenotype specific model for the correction of brain cellular heterogeneity bias and its application to age, brain region and major depression. Epigenetics 2013 Mar;8(3):290-302. PMID: 23426267
# Genome: hg19
# Platform: 450k
GEO = "GSE49618"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CD34")] <- "HSCm"
  cell[str_detect(cell, "PMN")] <- "Granulocyte"
  cell[str_detect(cell, "MONO")] <- "Monocyte"
  cell[str_detect(cell, "CD19")] <- "Bcell"
  cell[str_detect(cell, "PROS")] <- "GMP"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE49667 -----------------------------------------------------------------------------------------------
# GEO: GSE49667
# Types: Treg, CD4Tcell
# Paper: https://pubmed.ncbi.nlm.nih.gov/23974203/
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49667
# Citation: Zhang Y, Maksimovic J, Naselli G, Qian J et al. Genome-wide DNA methylation analysis identifies hypomethylated genes regulated by FOXP3 in human regulatory T cells. Blood 2013 Oct 17;122(16):2823-36. PMID: 23974203
# Genome: hg19
# Platform: 450k
GEO = "GSE49667"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "Treg")] <- "Treg"
  cell[!str_detect(cell, "Treg")] <- "CD4Tcell"
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE50798 -----------------------------------------------------------------------------------------------
# GEO: GSE50798
# Types: Neuron, glia
# Paper: 
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50798
# Citation: 
# Genome: hg19
# Platform: 
GEO = "GSE50798"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "neur")] <- "Neuron"
  cell[str_detect(cell, "nonneu")] <- "Glia"
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)

  # Get data
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE61195 -----------------------------------------------------------------------------------------------
# GEO: GSE61195
# Types: CD4Tcell, CD8Tcell
# Paper: 
# GEO link: 
# Citation: 	Renauer PA, Coit P, Sawalha AH. The DNA methylation signature of human TCRαβ+CD4-CD8- double negative T cells reveals CG demethylation and a unique epigenetic architecture permissive to a broad stimulatory immune response. Clin Immunol 2015 Jan;156(1):19-27. PMID: 25451162
# Genome: hg19
# Platform: 450k
GEO = "GSE61195"
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CD4")] <- "CD4Tcell"
  cell[str_detect(cell, "CD8")] <- "CD8Tcell"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE63409 -----------------------------------------------------------------------------------------------
# GEO: GSE63409 
# Types: CMP, GMP, HSC, MPP
# Paper: 
# GEO link: 
# Citation: 
# Genome: hg19
# Platform: 
GEO = "GSE63409"
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "AML")] <- "No"
  cell <- paste0("###--",cell)
  cell[str_detect(cell, "HSC")] <- "HSCm"
  cell[str_detect(cell, "CMP")] <- "CMP"
  cell[str_detect(cell, "MPP")] <- "MPPm"
  cell[str_detect(cell, "CMP")] <- "CMP"
  cell[str_detect(cell, "GMP")] <- "GMP"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE66351 -----------------------------------------------------------------------------------------------
# GEO: GSE66351 
# Types: Glia, Neuron
# Paper: 
# GEO link: 
# Citation: 
# Genome: hg19
# Platform: 
GEO = "GSE66351"
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell[!str_detect(cell, "CTRL")] <- "No"
  cell <- paste0("###--",cell)
  cell[str_detect(cell, "Glia")] <- "Glia"
  cell[str_detect(cell, "Neuron")] <- "Neuron"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE87196 -----------------------------------------------------------------------------------------------
# Types: Bcell, CD4Tcell, CD8Tcell, CMP, GMP, HSCb, HSCm, Monocyte, MPPb, MPPm, Neutrophil, NKcell
# Paper: https://www.ncbi.nlm.nih.gov/pubmed/27867036
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87196
# Citation: Farlik M, Halbritter F, Müller F, Choudry FA et al. DNA Methylation Dynamics of Human Hematopoietic Stem Cell Differentiation. Cell Stem Cell 2016 Dec 1;19(6):808-822. PMID: 27867036
# Genome: hg38
# Platform: uWGBS
GEO = "GSE87196"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")

if (FALSE) {
  if (.validateExp(exp_path,throws=F)) {
    
    raw_dir = paste0(base_dir,"raw/")
    mkdirs(base_dir,raw_dir,exp_dir)
    
    soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
    soft <- GEOquery::getGEO(filename=soft)
    
    # Get colData
    cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
    cell <- paste0("###-",cell)
    cell[str_detect(cell, "B cell")] <- "Bcell"
    cell[str_detect(cell, "CD4")] <- "CD4Tcell"
    cell[str_detect(cell, "CD8")] <- "CD8Tcell"
    cell[str_detect(cell, "myeloid")] <- "CMP"
    cell[str_detect(cell, "Mono")] <- "Monocyte"
    cell[str_detect(cell, "Neu")] <- "Neutrophil"
    cell[str_detect(cell, "Hematopoietic stem cell from bone marrow")] <- "HSCm"
    cell[str_detect(cell, "Hematopoietic stem cell from peripheral blood")] <- "HSCb"
    cell[str_detect(cell, "Multipotent progenitor from bone marrow")] <- "MPPm"
    cell[str_detect(cell, "Multipotent progenitor from peripheral blood")] <- "MPPb"
    cell[str_detect(cell, "killer")] <- "NKcell"
    cell[str_detect(cell, "Granulocyte")] <- "GMP"
    remove_idx <- which(str_detect(cell, "###-"))
    exclude_id <- names(cell[remove_idx])
    colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
    colData <- colData[-remove_idx,,drop=FALSE]
    stopifnot(length(setdiff(colData$Cell,cell_types))==0)
    
    # Get data
    
    
    
    
    
  
    scm <- liftover_scMethrix(scm,chain,target_genome = "hg38")
    scm <- bin_scMethrix(scm,regions = probe.set)
    assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
  }
  assign(exp_name,scMethrix::load_scMethrix(exp_path))
}

#--- GSE88824 -----------------------------------------------------------------------------------------------
# GEO: GSE88824
# Types: neutrophils, CD4+ T cells, CD8+ T cells, NK cells, B cells and monocytes
# Paper: https://pubmed.ncbi.nlm.nih.gov/30571772/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88824
# Citation: Kennedy DW, White NM, Benton MC, Fox A et al. Critical evaluation of linear regression models for cell-subtype specific methylation signal from mixed blood cell DNA. PLoS One 2018;13(12):e0208915. PMID: 30571772
# Genome: hg19
# Platform: Illumina 450k
GEO = "GSE88824"
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  remove_idx <- which(str_detect(cell, "Case-"))
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "Neu")] <- "Neutrophil"
  cell[str_detect(cell, "Blood")] <- "WholeBlood"
  cell[str_detect(cell, "Mono")] <- "Monocyte"
  cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
  cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
  cell[str_detect(cell, "CD19")] <- "Bcell"
  cell[str_detect(cell, "NKcell")] <- "NKcell"
  remove_idx <- union( which(str_detect(cell, "###-")),remove_idx)
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = paste0(exp_dir,exp_name,".rds")))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE96612 -----------------------------------------------------------------------------------------------
# GEO: GSE96612
# Types: Neuron, glia
# Paper: https://pubmed.ncbi.nlm.nih.gov/30643296/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96612
# Citation: Rizzardi LF, Hickey PF, Rodriguez DiBlasi V, Tryggvadóttir R et al. Neuronal brain-region-specific DNA methylation and chromatin accessibility are associated with neuropsychiatric trait heritability. Nat Neurosci 2019 Feb;22(2):307-316. PMID: 30643296
# Genome: hg19
# Platform: Bulk WGBS
GEO = "GSE110554"
chain = chains[["hg19ToHg38"]]

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  soft <- GEOquery::getGEOfile(GEO)
  soft <- GEOquery::getGEO(filename=soft)
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  colData <- data.table(Sample = names(soft@gsms), Cell = cell, ID = cell)
  remove_idx <- str_detect(colData$Cell,".*_neg .*|.*_pos .*",negate = TRUE)
  #remove_id <- colData$Sample[remove_idx]
  colData <- colData[!remove_idx,]
  colData$Cell <- str_replace(colData$Cell,".*_neg .*","Glia")
  colData$Cell <- str_replace(colData$Cell,".*_pos .*","Neuron")
  colData$ID <- str_remove(colData$ID," \\(bisulfite-Seq\\)")
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "NK")] <- "NKcell"
  cell[str_detect(cell, "WB")] <- "Neutrophil"
  cell[str_detect(cell, "Th")] <- "CD4Tcell"
  cell[str_detect(cell, "Mono")] <- "MonocyteNKcell"
  cell[str_detect(cell, "Tc")] <- "CD8Tcell"
  cell[str_detect(cell, "NK")] <- "NKcell"
  cell[str_detect(cell, "NK")] <- "NKcell"
  cell[str_detect(cell, "NK")] <- "NKcell"
  remove_idx <- which(str_detect(cell, "pulp"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  
  
  
  
  
  
  
  
  
  
  base_dir = paste0(home_dir,GEO,"/")
  raw_dir = paste0(base_dir,"raw/")
  bed_dir = paste0(base_dir,"bed/")
  exp_dir = paste0(base_dir,"exp/")
  mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
  colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
  chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")
  
  # Get colData
  
  data.table::fwrite(colData, file=colData_file, quote=FALSE, sep='\t', row.names = FALSE)
  
  # Get data
  supp_files <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*sorted.CpG_unstranded.txt.gz")
  sapply(supp_files, function(file) GEOquery::gunzip(file, remove=TRUE))
  file.remove(supp_files)
  files = list.files(raw_dir, full.names=TRUE)
  meth_mtx <- raw_files[files %like% "M_matrix"]
  cov_mtx <- raw_files[files %like% "Cov_matrix"]
  
  rrng <- fread(meth_mtx,select=c(1:3))
  
  for (i in 4+which(!remove_idx)) {
    
    message("Parsing ",i)
    
    meth <- fread(meth_mtx,select=i)
    cov <- fread(cov_mtx,select=i)
    beta <- meth/cov
    name <- colData$Sample[which(colData$ID %like% names(beta))]
    
    data.table::fwrite(cbind(rrng,beta), paste0(bed_dir,name,".bedgraph"), append = FALSE, sep = "\t", row.names = FALSE,
                       col.names = FALSE, quote = FALSE)
  }
  
  liftover_beds(files = files,chain = chain_file)
  
  # Compare colData and data
  expect_true(setequal(get_sample_name(files),colData$Sample))
  colData = read.table(file = colData_file, sep = '\t', header = TRUE)
  expect_true(all(colData$Cell %in% cell_types))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))
  
#--- GSE98203 -----------------------------------------------------------------------------------------------
# GEO: GSE98203
# Types: Neuron
# Paper: https://pubmed.ncbi.nlm.nih.gov/28556790/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98203
# Citation: Kozlenkov A, Jaffe AE, Timashpolsky A, Apontes P et al. DNA Methylation Profiling of Human Prefrontal Cortex Neurons in Heroin Users Shows Significant Difference between Genomic Contexts of Hyper- and Hypomethylation and a Younger Epigenetic Age. Genes (Basel) 2017 May 30;8(6). PMID: 28556790
# Genome: hg19
# Platform: Illumina 450k
GEO = "GSE98203"
chain = chains[["hg19ToHg38"]]
array = "450k"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CONTROL")] <- "Neuron"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE110554 ----------------------------------------------------------------------------------------------
# GEO: GSE110554
# Types: neutrophils, monocytes, B-lymphocytes, natural killer (NK) cells, CD4+ T-cells, and CD8+ T-cells
# Paper: https://pubmed.ncbi.nlm.nih.gov/29843789/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554
# Citation: Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
# Genome: hg19
# Platform: Illumina EPIC
GEO = "GSE110554"
chain <- chains[["hg19ToHg38"]]
array = "EPIC"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  # Get colData
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)

  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description[1])
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "NK")] <- "NKcell"
  cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
  cell[str_detect(cell, "Mono")] <- "Monocyte"
  cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
  cell[str_detect(cell, "mix")] <- "ImmMix"
  cell[str_detect(cell, "Bcell")] <- "Bcell"
  cell[str_detect(cell, "Neu")] <- "Neutrophil"
  remove_idx <- which(str_detect(cell, "###-"))
  colData <- colData[-remove_idx,,drop=FALSE]
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  stopifnot(length(setdiff(colData(scm)$Cell,cell_types))==0)
  stopifnot(nrow(scm) == length(probe.set))
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE112618 --------------------------------------------------------------------------------------------
# Types: Mixed Immune
# Paper: https://www.ncbi.nlm.nih.gov/pubmed/29843789
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112618
# Citation: Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
# Genome: hg19
# Platform: EPIC
GEO = "GSE112618"
chain <- chains[["hg19ToHg38"]]
array <- "EPIC"

base_dir <- paste0(home_dir,GEO,"/")
exp_dir <- paste0(base_dir,"exp/")
exp_name <- paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  # Get colData
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  colData <- data.frame(row.names = names(soft@gsms), Cell = rep("ImmMix",6))
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))
  
#--- GSE121483 ----------------------------------------------------------------------------------------------
# GEO: GSE121483
# Types: Microglia-like macrophages
# Paper: https://pubmed.ncbi.nlm.nih.gov/30451869/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121483
# Citation: Lund H, Pieber M, Parsa R, Han J et al. Competitive repopulation of an empty microglial niche yields functionally distinct subsets of microglia-like cells. Nat Commun 2018 Nov 19;9(1):4845. PMID: 30451869
# Genome: hg19
# Platform: Illumina Epic
GEO <- "GSE121483"
chain <- chains[["hg19ToHg38"]]
array <- "EPIC"

base_dir <- paste0(home_dir,GEO,"/")
exp_dir <- paste0(base_dir,"exp/")
exp_name <- paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "myeloid")] <- "CMP"
  cell[str_detect(cell, "Granulocyte")] <- "GMP"
  cell[str_detect(cell, "Ly6C")] <- "Ly6C"
  cell[str_detect(cell, "Naive")] <- "Microglia"
  cell[str_detect(cell, "proliferating")] <- "Inf.microglia"
  cell[str_detect(cell, "CNS")] <- "Inf.macrophage"
  cell[str_detect(cell, "monocyte")] <- "cMOP"
  remove_idx <- which(str_detect(cell, "pulp"))
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE128654 --------------------------------------------------------------------------------------------
# Types: Glioma, Ini.Glioma
# Paper: 
# GEO link: 
# Citation: 
# Genome: hg19
# Platform: 
GEO = "GSE128654"
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"
annotation = "ilmn12.hg19"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "tumor")] <- "Glioma"
  cell[str_detect(cell, "cell")] <- "Ini.Glioma"
  remove_idx <- which(str_detect(cell, "###-"))
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))



#------------------------------------------------------------------------------------------------------------
# Types: 
# Paper: 
# GEO link: 
# Citation: 
# Genome: hg19
# Platform: 
GEO = "GSE137845"
chain = chains[["hg19ToHg38"]]
array = "IlluminaHumanMethylation450k"
annotation = "ilmn12.hg19"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  break()
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CONTROL")] <- "Neuron"
  cell[str_detect(cell, "myeloid")] <- "CMP"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(all(cell %in% cell_types))
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, exclude_id = exclude_id,array = array, annotation = annotation)
  scm <- liftover_scMethrix(scm,chain,target_genome = "hg38")
  scm <- bin_scMethrix(scm,regions = probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE144804 ----------------------------------------------------------------------------------------------
# GEO: GSE144804
# Types: Endothelial
# Paper: https://pubmed.ncbi.nlm.nih.gov/32231389/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144804
# Citation: Rhead B, Shao X, Quach H, Ghai P et al. Global expression and CpG methylation analysis of primary endothelial cells before and after TNFa stimulation reveals gene modules enriched in inflammatory and infectious diseases and associated DMRs. PLoS One 2020;15(3):e0230884. PMID: 32231389
# Genome: hg19
# Platform: Illumina EPIC

GEO = "GSE144804"
chain = chains[["hg19ToHg38"]]
array = "EPIC"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  remove_idx <- which(str_detect(cell, "_TNF"))
  cell[str_detect(cell, "Huvec")] <- "Neuron"
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE151506 ----------------------------------------------------------------------------------------------
# GEO: GSE151506
# Types: GBM-IDH,GBM-WT,GBM-imm
# Paper: https://pubmed.ncbi.nlm.nih.gov/34594037/
# GEO: Chaligne, Ronan, Federico Gaiti, Dana Silverbush, Joshua S. Schiffman, Hannah R. Weisman, Lloyd Kluegel, Simon Gritsch et al. "Epigenetic encoding, heritability and plasticity of glioma transcriptional cell states." Nature Genetics 53, no. 10 (2021): 1469-1479.
# Genome: hg38
# Platform: RRBS
GEO = "GSE151506"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "CONTROL")] <- "Neuron"
  cell[str_detect(cell, "myeloid")] <- "CMP"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(all(cell %in% cell_types))
  
  # Get data
  scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))
  

#--- GSE164149 --------------------------------------------------------------------------------------------
# Types: Treg
# Paper: https://pubmed.ncbi.nlm.nih.gov/34348163/
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164149
# Citation: Lam, Avery J., et al. "Optimized CRISPR-mediated gene knockin reveals FOXP3-independent maintenance of human Treg identity." Cell Reports 36.5 (2021): 109494.
# Genome: hg19
# Platform: EPIC
GEO = "GSE164149"
array = "IlluminaHumanMethylationEPIC"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
  cell <- paste0("###-",cell)
  cell[str_detect(cell, "Cas9")] <- "Treg"
  remove_idx <- which(str_detect(cell, "###-"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(length(setdiff(colData$Cell,cell_types))==0)
  
  # Get data
  scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array=array)
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))

#--- GSE166844 ----------------------------------------------------------------------------------------------
# GEO: GSE166844
# Types: CD4 T-cells, CD8 T-cells, Monocytes, Granulocyte, B-cells
# Paper: https://pubmed.ncbi.nlm.nih.gov/33739972/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166844
# Citation: Hannon E, Mansell G, Walker E, Nabais MF et al. Assessing the co-variability of DNA methylation across peripheral cells and tissues: Implications for the interpretation of findings in epigenetic epidemiology. PLoS Genet 2021 Mar;17(3):e1009443. PMID: 33739972
# Genome: hg19
# Platform: Illumina EPIC
GEO = "GSE166844"
array = "EPIC"

base_dir = paste0(home_dir,GEO,"/")
exp_dir = paste0(base_dir,"exp/")
exp_name = paste0("scm.", GEO)
exp_path = paste0(exp_dir,exp_name,".rds")

if (!.validateType(exp_path,"file",throws=F)) {
  
  raw_dir = paste0(base_dir,"raw/")
  mkdirs(base_dir,raw_dir,exp_dir)
  
  soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
  soft <- GEOquery::getGEO(filename=soft)
  
  # Get colData
  ids <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
  cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
  cell[str_detect(cell, "Monocyte")] <- "Monocyte"
  cell[str_detect(cell, "B-cells")] <- "Bcell"
  cell[str_detect(cell, "CD4")] <- "CD4Tcell"
  cell[str_detect(cell, "CD8")] <- "CD8Tcell"
  cell[str_detect(cell, "blood")] <- "WholeBlood"
  cell[str_detect(cell, "Granulocyte")] <- "Granulocyte"
  cell[str_detect(cell, "CD8")] <- "CD8Tcell"
  remove_idx <- which(str_detect(cell, "Nasal|Buccal"))
  exclude_id <- names(cell[remove_idx])
  colData <- data.frame(row.names = names(soft@gsms), Cell = cell, ID = ids)
  colData <- colData[-remove_idx,,drop=FALSE]
  stopifnot(all(colData$cell %in% cell_types))

  scm <- var.proc.to.scMethrix(GEO = GEO, raw_dir = raw_dir,colData = colData, array = array, proc_file = paste0(raw_dir, "GSE166844_Variance_processed_signals.csv.gz"))
  scm <- standardize.scMethrix(scm, chain, probe.set)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
}
assign(exp_name,scMethrix::load_scMethrix(exp_path))
