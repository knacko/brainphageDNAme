# delete_data_set <- function(...) {
#   Exps <- unlist(...)
#   sapply(Exps, function (Exp) {
#     file.remove(paste0(dataset_dir, Exp, "/", "exp/", "scm.", Exp, ".rds"))
#   })
# }

get_data_set <- function(dataset_dir, GEO = names(type), regions = "all", genome = "hg38", assign.it = F, merge = F, cells = NULL) {

  soft = TRUE
  
  type <- list(singh =    c(!soft, NA,   "hg38", "EPIC"),
               GSE35069 =  c(soft, "soft", "hg19", "n450k"),
               GSE41826 =  c(soft, "soft", "hg19", "n450k"),
               GSE49618 =  c(soft, "idat", "hg19", "n450k"),
               GSE49667 =  c(soft, "soft", "hg19", "n450k"),
               GSE50798 =  c(soft, "soft", "hg19", "n450k"),
               GSE61195 =  c(soft, "soft", "hg19", "n450k"),
               GSE63409 =  c(soft, "idat", "hg19", "n450k"),
               GSE66351 =  c(soft, "idat", "hg19", "n450k"),
               GSE83458 =  c(soft, "soft", "hg19", "n450k"),
               GSE87196 =  c(soft, "bed",  "hg38",  NA),
               GSE88824 =  c(soft, "idat", "hg19", "n450k"),
               #GSE96612 =  c(soft, "bed", "hg19",  NULL),
               GSE98203 =  c(soft, "idat", "hg19", "n450k"),
               GSE103211 = c(soft, "soft", "hg19", "n450k"),
               GSE103659 = c(soft, "soft", "hg19", "n450k"),
               GSE104293 = c(soft, "idat", "hg19", "n450k"),
               GSE110554 = c(soft, "idat", "hg19", "EPIC"),
               GSE112618 = c(soft, "idat", "hg19", "EPIC"),
               GSE121483 = c(soft, "idat", "hg19", "EPIC"),
               GSE128654 = c(soft, "idat", "hg19", "n450k"),
               GSE144804 = c(soft, "idat", "hg19", "EPIC"),
               #GSE151506 = c(soft, "bed",  NULL),
               #GSE151506.mgsig = c(!soft, "bed",  "hg38", NA),
               GSE164149 = c(soft, "idat", "hg19", "EPIC"),
               GSE166844 = c(soft, "var", "hg19",  "EPIC")
               #test = c()
  )
  
  if (length(GEO) > 1) {
    
    exps <- lapply(GEO,get_data_set,dataset_dir = dataset_dir, genome = genome, cells = cells, regions = regions)
    names(exps) <- GEO

    if (assign.it || !merge) return(exps)

    if (merge == T) {
      scm <- exps[[1]]
      if (length(exps) > 1) {
        for (i in 2:length(exps)) {
          message("Adding ",names(exps)[i])
          scm <- merge_scMethrix(scm, exps[[i]], by = "col",verbose=F)
        }
      }
    }
    
    return(scm)
    
  } else {

    has_soft = type[[GEO]][1]
    import_type = type[[GEO]][2]
    src_genome = type[[GEO]][3]
    array = type[[GEO]][4]
    remove_idx <- cd <- NULL
    
    message("Getting ",GEO)

    base_dir <- paste0(dataset_dir, GEO, "/")
    exp_dir <- paste0(base_dir, "exp/")
    exp_name <- paste0(c("scm",GEO,genome,regions), collapse=".")
    exp_all <- paste0(exp_dir,"scm.",GEO,".",src_genome,".all.rds")
    exp_path <- paste0(exp_dir,exp_name,".rds")

    if (!.validateExp(exp_path, throws = F)) {
      
      if (!.validateExp(exp_all, throws = F)) {
        raw_dir <- paste0(base_dir, "raw/")
        mkdirs(base_dir, raw_dir, exp_dir)
        
        if (has_soft) {
          #if(!exists("soft") && !identical(soft@header$geo_accession, GEO)) {
          soft <<- GEOquery::getGEOfile(GEO, destdir = raw_dir)
          soft <<- GEOquery::getGEO(filename = soft)
          #}
        }

        if (GEO == "singh") {
          
# singh -----------------------------------------------------------------------------------------------------------
# Types: Microglia
# Paper: https://actaneurocomms.biomedcentral.com/articles/10.1186/s40478-021-01249-9#Sec10
# GEO link:
# Citation: Singh, O., Pratt, D., & Aldape, K. (2021). Immune cell deconvolution of bulk DNA methylation data reveals an association with methylation class, key somatic alterations, and cell state in glial/glioneuronal tumors. Acta neuropathologica communications, 9(1), 1-17.
          colData <- data.frame(row.names = "Singh", Cell = "Microglia")
          scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)

        } else if (GEO == "GSE151506.mgsig"){
          
          files <- list.files(raw_dir, full.names = TRUE, pattern = ".*bedgraph$", ignore.case = T)
          
          rnames <- names <- gsub("^.*(P105.*)\\.M.*","\\1",get_sample_name(files))
          rnames[duplicated(rnames)] <- paste0(names[duplicated(rnames)],"_2")
          cell <- gsub("^.*?\\.(.*?)\\..*","\\1",get_sample_name(files))
          
          colData <- data.frame(row.names = get_sample_name(files), Sample = names, Cell = cell)
          
          scm <- read_beds(files,colData=colData, chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4,
                             cov_idx = 5,keep_cov=F, genome_name="hg38")
          
          scm <- collapse_samples(scm,colname="Sample")
          
          colData <- colData[which(!duplicated(colData$Sample)),]
          
          row.names(colData(scm)) <- colData$Sample
          colData(scm)$Cell <- colData$Cell
          colData(scm)$Samples <- NULL
          colData(scm)$n_Samples <- NULL
          
        } else {
          
          if (GEO == "GSE35069") { 
# GSE35069 -------------------------------------------------------------------------------------------------------- 
# Types: CD4+ T cells, CD8+ T cells, CD56+ NK cells, CD19+ B cells, CD14+ monocytes, granulocytes
# Paper: https://pubmed.ncbi.nlm.nih.gov/22848472/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069
    # Citation: Reinius LE, Acevedo N, Joerink M, Pershagen G et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PLoS One 2012;7(7):e41361. PMID: 22848472            
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell[str_detect(cell, "blood")] <- "WholeBlood"
            cell[str_detect(cell, "Gran")] <- "Granulocyte"
            cell[str_detect(cell, "CD4")] <- "CD4Tcell"
            cell[str_detect(cell, "CD8")] <- "CD8Tcell"
            cell[str_detect(cell, "CD14")] <- "Monocyte"
            cell[str_detect(cell, "CD19")] <- "Bcell"
            cell[str_detect(cell, "CD56")] <- "NKcell"
            cell[str_detect(cell, "Neu")] <- "Granulocyte"
            cell[str_detect(cell, "Eos")] <- "Granulocyte"
            remove_idx <- which(str_detect(cell, "PBMC"))
          } else if (GEO =="GSE41826") {
# GSE41826 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "G")] <- "Glia"
            cell[str_detect(cell, "N")] <- "Neuron"
            status <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$characteristics_ch1[[2]])
            cell[!str_detect(status, "Control")] <- "###-"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE49618") {
# GSE49618 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "CD34")] <- "HSCm"
            cell[str_detect(cell, "PMN")] <- "Granulocyte"
            cell[str_detect(cell, "MONO")] <- "Monocyte"
            cell[str_detect(cell, "CD19")] <- "Bcell"
            cell[str_detect(cell, "PROS")] <- "GMP"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE49667") {
# GSE49667 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell[str_detect(cell, "Treg")] <- "Treg"
            cell[!str_detect(cell, "Treg")] <- "CD4Tcell"
          } else if (GEO =="GSE50798") {
# GSE50798 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell[str_detect(cell, "neur")] <- "Neuron"
            cell[str_detect(cell, "nonneu")] <- "Glia"
          } else if (GEO =="GSE61195") {
# GSE61195 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "CD4")] <- "CD4Tcell"
            cell[str_detect(cell, "CD8")] <- "CD8Tcell"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE63409") {
# GSE63409 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell[str_detect(cell, "AML")] <- "No"
            cell <- paste0("###--", cell)
            cell[str_detect(cell, "HSC")] <- "HSCm"
            cell[str_detect(cell, "CMP")] <- "CMP"
            cell[str_detect(cell, "MPP")] <- "MPPm"
            cell[str_detect(cell, "CMP")] <- "CMP"
            cell[str_detect(cell, "GMP")] <- "GMP"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE66351") {
# GSE66351 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell[!str_detect(cell, "CTRL")] <- "No"
            cell <- paste0("###--", cell)
            cell[str_detect(cell, "Glia")] <- "Glia"
            cell[str_detect(cell, "Neuron")] <- "Neuron"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE83458") {
# GSE83458 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "CD14_monocyte")] <- "Monocyte"
            cell[str_detect(cell, "DC_donor")] <- "Dendritic"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE87196") {
# GSE87196 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "B cell")] <- "Bcell"
            cell[str_detect(cell, "CD4")] <- "CD4Tcell"
            cell[str_detect(cell, "CD8")] <- "CD8Tcell"
            cell[str_detect(cell, "myeloid")] <- "CMP"
            cell[str_detect(cell, "Mono")] <- "Monocyte"
            cell[str_detect(cell, "Neu")] <- "Granulocyte"
            cell[str_detect(cell, "Hematopoietic stem cell from bone marrow")] <- "HSCm"
            cell[str_detect(cell, "Hematopoietic stem cell from peripheral blood")] <- "HSCb"
            cell[str_detect(cell, "Multipotent progenitor from bone marrow")] <- "MPPm"
            cell[str_detect(cell, "Multipotent progenitor from peripheral blood")] <- "MPPb"
            cell[str_detect(cell, "killer")] <- "NKcell"
            cell[str_detect(cell, "Granulocyte")] <- "GMP"
            cell[str_detect(cell, "Common lymphoid progenitor")] <- "CLP"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE88824") {
# GSE88824 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            remove_idx <- which(str_detect(cell, "Case-"))
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "Neu")] <- "Granulocyte"
            cell[str_detect(cell, "Blood")] <- "WholeBlood"
            cell[str_detect(cell, "Mono")] <- "Monocyte"
            cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
            cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
            cell[str_detect(cell, "CD19")] <- "Bcell"
            cell[str_detect(cell, "NKcell")] <- "NKcell"
            remove_idx <- union(which(str_detect(cell, "###-")), remove_idx)
         # } else if (GEO =="GSE96612") {
            
          } else if (GEO =="GSE98203") {
# GSE98203 --------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "CONTROL")] <- "Neuron"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE103211") {
# GSE103211 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "DC")] <- "Dendritic"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE103659") {
# GSE103659 -------------------------------------------------------------------------------------------------------
            cell = rep("GBM", length(names(soft@gsms)))
          } else if (GEO =="GSE104293") {
# GSE104293 -------------------------------------------------------------------------------------------------------
            cell = rep("Glioma", length(names(soft@gsms)))
          } else if (GEO =="GSE110554") {
# GSE110554 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description[1])
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "NK")] <- "NKcell"
            cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
            cell[str_detect(cell, "Mono")] <- "Monocyte"
            cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
            cell[str_detect(cell, "mix")] <- "ImmMix"
            cell[str_detect(cell, "Bcell")] <- "Bcell"
            cell[str_detect(cell, "Neu")] <- "Granulocyte"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE112618") {
# GSE112618 -------------------------------------------------------------------------------------------------------
            cell = rep("ImmMix", length(names(soft@gsms)))
          } else if (GEO =="GSE121483") {
# GSE121483 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
            cell[str_detect(cell, "myeloid")] <- "CMP"
            cell[str_detect(cell, "Granulocyte")] <- "GMP"
            cell[str_detect(cell, "Ly6C")] <- "Ly6C"
            cell[str_detect(cell, "Naive")] <- "Microglia"
            cell[str_detect(cell, "proliferating")] <- "Inf.microglia"
            cell[str_detect(cell, "CNS")] <- "Inf.macrophage"
            cell[str_detect(cell, "monocyte")] <- "cMOP"
            remove_idx <- which(str_detect(cell, "pulp"))
          } else if (GEO =="GSE122914") {
            
          } else if (GEO =="GSE128654") {
# GSE128654 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "tumor")] <- "Glioma"
            cell[str_detect(cell, "cell")] <- "Ini.Glioma"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE144804") {
# GSE144804 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            remove_idx <- which(str_detect(cell, "_TNF"))
            cell[str_detect(cell, "Huvec")] <- "Endothelial"
          #} else if (GEO =="GSE151506") {
            #structure(list(Labels=c("MG","MG","MG","M<U+03D5>","MG","MG","M<U+03D5>","M<U+03D5>","MG","MG","M<U+03D5>","MG","M<U+03D5>","M<U+03D5>","M<U+03D5>","MG","MG","M<U+03D5>","M<U+03D5>","M<U+03D5>")),class="data.frame",row.names=c("MGH105C_P8_A10_","MGH105C_P8_B1_","MGH105C_P8_B5_","MGH105C_P8_B7_","MGH105C_P8_A2_","MGH105C_P8_C4_","MGH105C_P8_A3_","MGH105C_P8_C10_","MGH105C_P8_D11_","MGH105C_P8_E1_","MGH105C_P8_E2_","MGH105C_P8_E5_","MGH105C_P8_E6_","MGH105C_P8_A6_","MGH105C_P8_E12_","MGH105C_P8_F4_","MGH105C_P8_F7_","MGH105C_P8_A8_","MGH105C_P8_G11_","MGH105C_P8_H11_"))
            scm <- load_scMethrix("D:\\Git\\thesis_data\\GSE151506\\raw") 
            
            
          } else if (GEO =="GSE164149") {
# GSE164149 -------------------------------------------------------------------------------------------------------
            cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
            cell <- paste0("###-", cell)
            cell[str_detect(cell, "Cas9")] <- "Treg"
            remove_idx <- which(str_detect(cell, "###-"))
          } else if (GEO =="GSE166844") {
# GSE166844 -------------------------------------------------------------------------------------------------------
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
            cd <- data.frame(row.names = names(soft@gsms), ID = ids)
            id_col = "ID"
          } 
          
          colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
          if (!is.null(cd)) colData <- transform(merge(colData,cd,by="row.names"), row.names=Row.names, Row.names=NULL)
          if (length(remove_idx) != 0) colData <- colData[-remove_idx, , drop = FALSE]
          stopifnot(all(colData$cell %in% get_cell_types()$all))

          if (import_type == "idat") {
            scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
          } else if (import_type == "var") {
            scm <- var.proc.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array, 
                                         proc_file = proc_file, id_col = id_col)
          } else if (import_type == "soft") {
            scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
          } else if (import_type == "bed") {
            scm <- tar.bed.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData)
          }
        }
        
        scm <- standardize.scMethrix(scm, GEO)
        scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
      } else {
        scm <- scMethrix::load_scMethrix(dest = exp_all)
      }

      if (regions != "all" || genome != src_genome) {
        
        scm <- standardize.scMethrix(scm, GEO = GEO, src_genome = src_genome, out_genome = genome,  regions = get_probes()[[regions]])
        scm <- scMethrix::save_scMethrix(scm, dest = exp_path)
      }
    } else {
      scm <- scMethrix::load_scMethrix(exp_path)
    }
    
    if (!is.null(cells)) scm <- scm[,(left_match(colData(scm)$Cell,cells))]
    
    if (assign.it) {
      assign(exp_name, scm, envir = .GlobalEnv)
      return(TRUE)
    }
  }
  
  return(scm)
}

#         
#     
#   # #------------------------------------------------------------------------------------------------------------
#   # # Types:
#   # # Paper:
#   # # GEO link:
#   # # Citation:
#   # # Genome: hg19
#   # # Platform:
#   # GEO = ""
#   # chain = chains[["hg19ToHg38"]]
#   # array = "IlluminaHumanMethylation450k"
#   #
#   # base_dir = paste0(home_dir,GEO,"/")
#   # exp_dir = paste0(base_dir,"exp/")
#   # exp_name = paste0("scm.", GEO)
#   # exp_path = paste0(exp_dir,exp_name,".rds")
#   #
#   # if (!.validateType(exp_path,"file",throws=F)) {
#   #
#   #   raw_dir = paste0(base_dir,"raw/")
#   #   mkdirs(base_dir,raw_dir,exp_dir)
#   #
#   #   soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
#   #   soft <- GEOquery::getGEO(filename=soft)
#   #
#   #   # Get colData
#   #   cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#   #   cell <- paste0("###-",cell)
#   #   cell[str_detect(cell, "CONTROL")] <- "Neuron"
#   #   cell[str_detect(cell, "myeloid")] <- "CMP"
#   #   remove_idx <- which(str_detect(cell, "###-"))
#   #   exclude_id <- names(cell[remove_idx])
#   #   colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#   #   colData <- colData[-remove_idx,,drop=FALSE]
#   #   stopifnot(length(setdiff(colData$Cell,cell_types[["all"]]))==0)
#   #   
#   #
#   #   # Get data
#   #   scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#   #   scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#   #   scm <- var.proc.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array, proc_file = paste0(raw_dir, ""))
#   #
#   #   scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#   #   assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path), envir = .GlobalEnv)
#   # } else {assign(exp_name,scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)}
# 
#   #--- Singh --------------------------------------------------------------------------------------------------
#   # Types: Microglia
#   # Paper: https://actaneurocomms.biomedcentral.com/articles/10.1186/s40478-021-01249-9#Sec10
#   # GEO link:
#   # Citation: Singh, O., Pratt, D., & Aldape, K. (2021). Immune cell deconvolution of bulk DNA methylation data reveals an association with methylation class, key somatic alterations, and cell state in glial/glioneuronal tumors. Acta neuropathologica communications, 9(1), 1-17.
#   # Genome: hg19
#   # Platform:
#   if (GEO == "singh") {
#     array <- "IlluminaHumanMethylation450k"
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       # Get colData
#       colData <- data.frame(row.names = "Singh", Cell = "Microglia")
#       stopifnot(length(setdiff(colData$Cell, cell_types[["all"]])) == 0)
#       
#       files <- list.files(raw_dir, full.names = TRUE, pattern = ".*idat$", ignore.case = T)
#       RGset <- minfi::read.metharray.exp(raw_dir, force = TRUE)
#       Mset <- preprocessNoob(RGset)
#       colnames(Mset) <- "Singh"
#       Rset <- minfi::ratioConvert(Mset)
#       GRset <- minfi::mapToGenome(Rset)
#       GRset <- convertArray(GRset, outType = "IlluminaHumanMethylation450k")
# 
#       scm <- scMethrix::as.scMethrix.GRset(GRset = GRset, colData = colData)
#       colnames(colData(scm)) <- "Cell" # TODO: Not sure why this is necessary
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path, envir = .GlobalEnv))
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE35069 -----------------------------------------------------------------------------------------------
#   # GEO: GSE35069
#   # Types: CD4+ T cells, CD8+ T cells, CD56+ NK cells, CD19+ B cells, CD14+ monocytes, granulocytes
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/22848472/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069
#   # Citation: Reinius LE, Acevedo N, Joerink M, Pershagen G et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PLoS One 2012;7(7):e41361. PMID: 22848472
#   # Genome: hg19
#   # Platform: Illumina 450k
# 
#   if ("GSE35069" %in% GEOs) {
#     GEO <- "GSE35069"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
# 
#         # Get colData
#           cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#           cell[str_detect(cell, "blood")] <- "WholeBlood"
#           cell[str_detect(cell, "Gran")] <- "Granulocyte"
#           cell[str_detect(cell, "CD4")] <- "CD4Tcell"
#           cell[str_detect(cell, "CD8")] <- "CD8Tcell"
#           cell[str_detect(cell, "CD14")] <- "Monocyte"
#           cell[str_detect(cell, "CD19")] <- "Bcell"
#           cell[str_detect(cell, "CD56")] <- "NKcell"
#           cell[str_detect(cell, "Neu")] <- "Granulocyte"
#           cell[str_detect(cell, "Eos")] <- "Granulocyte"
#           remove_idx <- which(str_detect(cell, "PBMC"))
#           colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#           colData <- colData[-remove_idx, , drop = FALSE]
#           stopifnot(colData$Cell %allin% cell.list)
#           
#           # Get data
#           scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
# 
#         
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE41826 -----------------------------------------------------------------------------------------------
#   # GEO: GSE41826
#   # Types: Neuron, glia
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/23426267/
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41826
#   # Citation: Guintivano J, Aryee MJ, Kaminsky ZA. A cell epigenotype specific model for the correction of brain cellular heterogeneity bias and its application to age, brain region and major depression. Epigenetics 2013 Mar;8(3):290-302. PMID: 23426267
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE41826" %in% GEOs) {
#     GEO <- "GSE41826"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "G")] <- "Glia"
#         cell[str_detect(cell, "N")] <- "Neuron"
#         status <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$characteristics_ch1[[2]])
#         cell[!str_detect(status, "Control")] <- "###-"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE49618 -----------------------------------------------------------------------------------------------
#   # GEO: GSE49618
#   # Types: Monocyte, Bcell, HSC
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/23426267/
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49618
#   # Citation: Guintivano J, Aryee MJ, Kaminsky ZA. A cell epigenotype specific model for the correction of brain cellular heterogeneity bias and its application to age, brain region and major depression. Epigenetics 2013 Mar;8(3):290-302. PMID: 23426267
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE49618" %in% GEOs) {
#     GEO <- "GSE49618"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "CD34")] <- "HSCm"
#         cell[str_detect(cell, "PMN")] <- "Granulocyte"
#         cell[str_detect(cell, "MONO")] <- "Monocyte"
#         cell[str_detect(cell, "CD19")] <- "Bcell"
#         cell[str_detect(cell, "PROS")] <- "GMP"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(colData$Cell %allin% cell_types[["all"]])
#         
#         # Get data
#         scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#       
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE49667 -----------------------------------------------------------------------------------------------
#   # GEO: GSE49667
#   # Types: Treg, CD4Tcell
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/23974203/
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49667
#   # Citation: Zhang Y, Maksimovic J, Naselli G, Qian J et al. Genome-wide DNA methylation analysis identifies hypomethylated genes regulated by FOXP3 in human regulatory T cells. Blood 2013 Oct 17;122(16):2823-36. PMID: 23974203
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE49667" %in% GEOs) {
#     GEO <- "GSE49667"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#         soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#         soft <- GEOquery::getGEO(filename = soft)
#       }
#       
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       cell[str_detect(cell, "Treg")] <- "Treg"
#       cell[!str_detect(cell, "Treg")] <- "CD4Tcell"
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       stopifnot(all(colData$cell %in% cell_types))
#       
#       # Get data
#       scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE50798 -----------------------------------------------------------------------------------------------
#   # GEO: GSE50798
#   # Types: Neuron, glia
#   # Paper:
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50798
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE50798" %in% GEOs) {
#     GEO <- "GSE50798"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       cell[str_detect(cell, "neur")] <- "Neuron"
#       cell[str_detect(cell, "nonneu")] <- "Glia"
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
# 
# 
#   #--- GSE60274 --------------------------------------------------------------------------------------------------------
# 
# 
# 
#   #--- GSE61195 -----------------------------------------------------------------------------------------------
#   # GEO: GSE61195
#   # Types: CD4Tcell, CD8Tcell
#   # Paper:
#   # GEO link:
#   # Citation: 	Renauer PA, Coit P, Sawalha AH. The DNA methylation signature of human TCRaß+CD4-CD8- double negative T cells reveals CG demethylation and a unique epigenetic architecture permissive to a broad stimulatory immune response. Clin Immunol 2015 Jan;156(1):19-27. PMID: 25451162
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE61195" %in% GEOs) {
#     GEO <- "GSE61195"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "CD4")] <- "CD4Tcell"
#       cell[str_detect(cell, "CD8")] <- "CD8Tcell"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE63409 -----------------------------------------------------------------------------------------------
#   # GEO: GSE63409
#   # Types: CMP, GMP, HSC, MPP
#   # Paper:
#   # GEO link:
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE63409" %in% GEOs) {
#     GEO <- "GSE63409"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       cell[str_detect(cell, "AML")] <- "No"
#       cell <- paste0("###--", cell)
#       cell[str_detect(cell, "HSC")] <- "HSCm"
#       cell[str_detect(cell, "CMP")] <- "CMP"
#       cell[str_detect(cell, "MPP")] <- "MPPm"
#       cell[str_detect(cell, "CMP")] <- "CMP"
#       cell[str_detect(cell, "GMP")] <- "GMP"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE66351 -----------------------------------------------------------------------------------------------
#   # GEO: GSE66351
#   # Types: Glia, Neuron
#   # Paper:
#   # GEO link:
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE66351" %in% GEOs) {
#     GEO <- "GSE66351"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       cell[!str_detect(cell, "CTRL")] <- "No"
#       cell <- paste0("###--", cell)
#       cell[str_detect(cell, "Glia")] <- "Glia"
#       cell[str_detect(cell, "Neuron")] <- "Neuron"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   
#   #--- GSE83458 -----------------------------------------------------------------------------------------------
#   # Types: Dendritic Cells
#   # Paper: https://www.sciencedirect.com/science/article/pii/S2211124717312792?via%3Dihub
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83458
#   # Citation: Rodríguez-Ubreva J, Català-Moll F, Obermajer N, Álvarez-Errico D et al. Prostaglandin E2 Leads to the Acquisition of DNMT3A-Dependent Tolerogenic Functions in Human Myeloid-Derived Suppressor Cells. Cell Rep 2017 Oct 3;21(1):154-167. PMID: 28978469
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE83458" %in% GEOs) {
#     GEO = "GSE83458"
#     array = "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
#   
#     base_dir = paste0(home_dir,GEO,"/")
#     exp_dir = paste0(base_dir,"exp/")
#     exp_name = paste0("scm.", GEO)
#     exp_path = paste0(exp_dir,exp_name,".rds")
#   
#     if (!.validateType(exp_path,"file",throws=F)) {
#   
#       raw_dir = paste0(base_dir,"raw/")
#       mkdirs(base_dir,raw_dir,exp_dir)
#   
#       soft <- GEOquery::getGEOfile(GEO,destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename=soft)
#   
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "CD14_monocyte")] <- "Monocyte"
#       cell[str_detect(cell, "DC_donor")] <- "Dendritic"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
#       
#       # Get data
#       scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE87196 -----------------------------------------------------------------------------------------------
#   # Types: Bcell, CD4Tcell, CD8Tcell, CMP, GMP, HSCb, HSCm, Monocyte, MPPb, MPPm, Granulocyte, NKcell
#   # Paper: https://www.ncbi.nlm.nih.gov/pubmed/27867036
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87196
#   # Citation: Farlik M, Halbritter F, Müller F, Choudry FA et al. DNA Methylation Dynamics of Human Hematopoietic Stem Cell Differentiation. Cell Stem Cell 2016 Dec 1;19(6):808-822. PMID: 27867036
#   # Genome: hg38
#   # Platform: uWGBS
#   if ("GSE87196" %in% GEOs) {
#     GEO <- "GSE87196"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       h5_dir <- paste0(base_dir, "h5/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "B cell")] <- "Bcell"
#       cell[str_detect(cell, "CD4")] <- "CD4Tcell"
#       cell[str_detect(cell, "CD8")] <- "CD8Tcell"
#       cell[str_detect(cell, "myeloid")] <- "CMP"
#       cell[str_detect(cell, "Mono")] <- "Monocyte"
#       cell[str_detect(cell, "Neu")] <- "Granulocyte"
#       cell[str_detect(cell, "Hematopoietic stem cell from bone marrow")] <- "HSCm"
#       cell[str_detect(cell, "Hematopoietic stem cell from peripheral blood")] <- "HSCb"
#       cell[str_detect(cell, "Multipotent progenitor from bone marrow")] <- "MPPm"
#       cell[str_detect(cell, "Multipotent progenitor from peripheral blood")] <- "MPPb"
#       cell[str_detect(cell, "killer")] <- "NKcell"
#       cell[str_detect(cell, "Granulocyte")] <- "GMP"
#       cell[str_detect(cell, "Common lymphoid progenitor")] <- "CLP"
#       remove_idx <- which(str_detect(cell, "###-"))
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- tar.bed.to.scMethrix(GEO = GEO, raw_dir = raw_dir, exp_dir = exp_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
# 
#   #--- GSE88824 -----------------------------------------------------------------------------------------------
#   # GEO: GSE88824
#   # Types: Granulocyte, CD4+ T cells, CD8+ T cells, NK cells, B cells and monocytes
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/30571772/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88824
#   # Citation: Kennedy DW, White NM, Benton MC, Fox A et al. Critical evaluation of linear regression models for cell-subtype specific methylation signal from mixed blood cell DNA. PLoS One 2018;13(12):e0208915. PMID: 30571772
#   # Genome: hg19
#   # Platform: Illumina 450k
#   if ("GSE88824" %in% GEOs) {
#     GEO <- "GSE88824"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       remove_idx <- which(str_detect(cell, "Case-"))
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "Neu")] <- "Granulocyte"
#       cell[str_detect(cell, "Blood")] <- "WholeBlood"
#       cell[str_detect(cell, "Mono")] <- "Monocyte"
#       cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
#       cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
#       cell[str_detect(cell, "CD19")] <- "Bcell"
#       cell[str_detect(cell, "NKcell")] <- "NKcell"
#       remove_idx <- union(which(str_detect(cell, "###-")), remove_idx)
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #XXX--- GSE96612 -----------------------------------------------------------------------------------------------
#   # GEO: GSE96612
#   # Types: Neuron, glia
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/30643296/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96612
#   # Citation: Rizzardi LF, Hickey PF, Rodriguez DiBlasi V, Tryggvadóttir R et al. Neuronal brain-region-specific DNA methylation and chromatin accessibility are associated with neuropsychiatric trait heritability. Nat Neurosci 2019 Feb;22(2):307-316. PMID: 30643296
#   # Genome: hg19
#   # Platform: Bulk WGBS
#   if ("GSE96612" %in% GEOs) {
#     GEO <- "GSE96612"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (F) {
#       if (!.validateType(exp_path, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
# 
#         soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#         soft <- GEOquery::getGEO(filename = soft)
# 
#         # Get colData
#         soft <- GEOquery::getGEOfile(GEO, dest_dir = raw_dir)
#         soft <- GEOquery::getGEO(filename = soft)
#         
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#         colData <- data.table(Sample = names(soft@gsms), Cell = cell, ID = cell)
#         remove_idx <- str_detect(colData$Cell, ".*_neg .*|.*_pos .*", negate = TRUE)
#         # remove_id <- colData$Sample[remove_idx]
#         colData <- colData[!remove_idx, ]
#         colData$Cell <- str_replace(colData$Cell, ".*_neg .*", "Glia")
#         colData$Cell <- str_replace(colData$Cell, ".*_pos .*", "Neuron")
#         colData$ID <- str_remove(colData$ID, " \\(bisulfite-Seq\\)")
#         stopifnot(length(setdiff(colData$Cell, cell_types[["all"]])) == 0)
#         
# 
# 
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#         cell[str_detect(cell, "NK")] <- "NKcell"
#         cell[str_detect(cell, "WB")] <- "Granulocyte"
#         cell[str_detect(cell, "Th")] <- "CD4Tcell"
#         cell[str_detect(cell, "Mono")] <- "MonocyteNKcell"
#         cell[str_detect(cell, "Tc")] <- "CD8Tcell"
#         cell[str_detect(cell, "NK")] <- "NKcell"
#         cell[str_detect(cell, "NK")] <- "NKcell"
#         cell[str_detect(cell, "NK")] <- "NKcell"
#         remove_idx <- which(str_detect(cell, "pulp"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#         base_dir <- paste0(home_dir, GEO, "/")
#         raw_dir <- paste0(base_dir, "raw/")
#         bed_dir <- paste0(base_dir, "bed/")
#         exp_dir <- paste0(base_dir, "exp/")
#         mkdirs(base_dir, raw_dir, bed_dir, exp_dir)
#         colData_file <- paste0(base_dir, "colData/", GEO, "_colData.tsv")
#         chain_file <- paste0(chain_dir, "hg19ToHg38.over.chain")
# 
#         # Get colData
# 
#         data.table::fwrite(colData, file = colData_file, quote = FALSE, sep = "\t", row.names = FALSE)
# 
#         # Get data
#         supp_files <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir, 1, nchar(raw_dir) - 1), filter_regex = ".*sorted.CpG_unstranded.txt.gz")
#         sapply(supp_files, function(file) GEOquery::gunzip(file, remove = TRUE))
#         file.remove(supp_files)
#         files <- list.files(raw_dir, full.names = TRUE)
#         meth_mtx <- raw_files[files %like% "M_matrix"]
#         cov_mtx <- raw_files[files %like% "Cov_matrix"]
# 
#         rrng <- fread(meth_mtx, select = c(1:3))
# 
#         for (i in 4 + which(!remove_idx)) {
#           message("Parsing ", i)
# 
#           meth <- fread(meth_mtx, select = i)
#           cov <- fread(cov_mtx, select = i)
#           beta <- meth / cov
#           name <- colData$Sample[which(colData$ID %like% names(beta))]
# 
#           data.table::fwrite(cbind(rrng, beta), paste0(bed_dir, name, ".bedgraph"),
#             append = FALSE, sep = "\t", row.names = FALSE,
#             col.names = FALSE, quote = FALSE
#           )
#         }
# 
#         liftover_beds(files = files, chain = chain_file)
# 
#         # Compare colData and data
#         expect_true(setequal(get_sample_name(files), colData$Sample))
#         colData <- read.table(file = colData_file, sep = "\t", header = TRUE)
#         expect_true(all(colData$Cell %in% cell_types))
#       } else {
#         assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#       }
#     }
#   }
# 
#   #--- GSE98203 -----------------------------------------------------------------------------------------------
#   # GEO: GSE98203
#   # Types: Neuron
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/28556790/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98203
#   # Citation: Kozlenkov A, Jaffe AE, Timashpolsky A, Apontes P et al. DNA Methylation Profiling of Human Prefrontal Cortex Neurons in Heroin Users Shows Significant Difference between Genomic Contexts of Hyper- and Hypomethylation and a Younger Epigenetic Age. Genes (Basel) 2017 May 30;8(6). PMID: 28556790
#   # Genome: hg19
#   # Platform: Illumina 450k
#   if ("GSE98203" %in% GEOs) {
#     GEO <- "GSE98203"
#     array <- "450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "CONTROL")] <- "Neuron"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
# 
# 
#   #--- GSE103211 -------------------------------------------------------------------------------------------------------
#   # Types: Dendritic
#   # Paper: https://www.sciencedirect.com/science/article/pii/S2211124717312792?via%3Dihub
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103211
#   # Citation: Rodríguez-Ubreva J, Català-Moll F, Obermajer N, Álvarez-Errico D et al. Prostaglandin E2 Leads to the Acquisition of DNMT3A-Dependent Tolerogenic Functions in Human Myeloid-Derived Suppressor Cells. Cell Rep 2017 Oct 3;21(1):154-167. PMID: 28978469
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE103211" %in% GEOs) {
#     GEO <- "GSE103211"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "DC")] <- "Dendritic"
#       remove_idx <- which(str_detect(cell, "###-"))
#       exclude_id <- names(cell[remove_idx])
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
#       
# 
#       # Get data
#       scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE103659 -------------------------------------------------------------------------------------------------------
#   # Types: GBM
#   # Paper:
#   # GEO link:
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE103659" %in% GEOs) {
#     GEO <- "GSE103659"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       colData <- data.frame(row.names = names(soft@gsms), Cell = rep("GBM", length(names(soft@gsms))))
#       stopifnot(all(colData$cell %in% cell_types))
#       
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       assign(paste0(exp_name,".std"), scMethrix::save_scMethrix(scm, dest = paste0(exp_dir,exp_name,".std.rds")), envir = .GlobalEnv)
#       if (standardize) scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE104293 ----------------------------------------------------------------------------------------------
#   # Types: Glioma
#   # Paper: https://www.ncbi.nlm.nih.gov/pubmed/29368212
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#   # Citation: Bady P, Kurscheid S, Delorenzi M, Gorlia T et al. The DNA methylome of DDR genes and benefit from RT or TMZ in IDH mutant low-grade glioma treated in EORTC 22033. Acta Neuropathol 2018 Apr;135(4):601-615. PMID: 29368212
#   # Genome: hg19
#   # Platform: 450k
#   if ("GSE104293" %in% GEOs) {
#     GEO <- "GSE104293"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         colData <- data.frame(row.names = names(soft@gsms), Cell = rep("Glioma", length(names(soft@gsms))))
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
#   
#   #--- GSE110554 ----------------------------------------------------------------------------------------------
#   # GEO: GSE110554
#   # Types: Granulocyte, monocytes, B-lymphocytes, natural killer (NK) cells, CD4+ T-cells, and CD8+ T-cells
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/29843789/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554
#   # Citation: Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
#   # Genome: hg19
#   # Platform: Illumina EPIC
#   if ("GSE110554" %in% GEOs) {
#     GEO <- "GSE110554"
#     array <- "EPIC"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       # Get colData
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description[1])
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "NK")] <- "NKcell"
#       cell[str_detect(cell, "CD4T")] <- "CD4Tcell"
#       cell[str_detect(cell, "Mono")] <- "Monocyte"
#       cell[str_detect(cell, "CD8T")] <- "CD8Tcell"
#       cell[str_detect(cell, "mix")] <- "ImmMix"
#       cell[str_detect(cell, "Bcell")] <- "Bcell"
#       cell[str_detect(cell, "Neu")] <- "Granulocyte"
#       remove_idx <- which(str_detect(cell, "###-"))
#       colData <- colData[-remove_idx, , drop = FALSE]
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE112618 --------------------------------------------------------------------------------------------
#   # Types: Mixed Immune
#   # Paper: https://www.ncbi.nlm.nih.gov/pubmed/29843789
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112618
#   # Citation: Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
#   # Genome: hg19
#   # Platform: EPIC
#   if ("GSE112618" %in% GEOs) {
#     GEO <- "GSE112618"
#     array <- "EPIC"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       # Get colData
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       colData <- data.frame(row.names = names(soft@gsms), Cell = rep("ImmMix", length(names(soft@gsms))))
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE121483 ----------------------------------------------------------------------------------------------
#   # GEO: GSE121483
#   # Types: Microglia-like macrophages
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/30451869/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121483
#   # Citation: Lund H, Pieber M, Parsa R, Han J et al. Competitive repopulation of an empty microglial niche yields functionally distinct subsets of microglia-like cells. Nat Commun 2018 Nov 19;9(1):4845. PMID: 30451869
#   # Genome: hg19
#   # Platform: Illumina Epic
#   if ("GSE121483" %in% GEOs) {
#     GEO <- "GSE121483"
#     array <- "EPIC"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#       cell[str_detect(cell, "myeloid")] <- "CMP"
#       cell[str_detect(cell, "Granulocyte")] <- "GMP"
#       cell[str_detect(cell, "Ly6C")] <- "Ly6C"
#       cell[str_detect(cell, "Naive")] <- "Microglia"
#       cell[str_detect(cell, "proliferating")] <- "Inf.microglia"
#       cell[str_detect(cell, "CNS")] <- "Inf.macrophage"
#       cell[str_detect(cell, "monocyte")] <- "cMOP"
#       remove_idx <- which(str_detect(cell, "pulp"))
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
# 
#   #XXX--- GSE122994 --------------------------------------------------------------------------------------------
#   # Types: Glioblastoma
#   # Paper: https://www.ncbi.nlm.nih.gov/pubmed/32064499
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122994
#   # Citation: Wick A, Kessler T, Platten M, Meisner C et al. Superiority of temozolomide over radiotherapy for elderly patients with RTK II methylation class, MGMT promoter methylated malignant astrocytoma. Neuro Oncol 2020 Aug 17;22(8):1162-1172. PMID: 32064499
#   # Genome: hg19
#   # Platform: 450k/EPIC
#   if ("GSE122994" %in% GEOs) {
#     GEO <- "GSE122994"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (F) {
#       if (!.validateType(exp_path, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
# 
#         soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#         soft <- GEOquery::getGEO(filename = soft)
# 
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "CONTROL")] <- "Neuron"
#         cell[str_detect(cell, "myeloid")] <- "CMP"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
# 
#         # Get data
#         scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#         scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#         scm <- var.proc.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array, proc_file = paste0(raw_dir, ""))
# 
#         scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       } else {
#         assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#       }
#     }
#   }
# 
#   #--- GSE128654 --------------------------------------------------------------------------------------------
#   # Types: Glioma, Ini.Glioma
#   # Paper:
#   # GEO link:
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE128654" %in% GEOs) {
#     GEO <- "GSE128654"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       cell <- paste0("###-", cell)
#       cell[str_detect(cell, "tumor")] <- "Glioma"
#       cell[str_detect(cell, "cell")] <- "Ini.Glioma"
#       remove_idx <- which(str_detect(cell, "###-"))
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
# 
# 
#   #XXX--- GSE137845 ----------------------------------------------------------------------------------------------
#   # Types:
#   # Paper:
#   # GEO link:
#   # Citation:
#   # Genome: hg19
#   # Platform:
#   if ("GSE137845" %in% GEOs) {
#     GEO <- "GSE137845"
#     array <- "IlluminaHumanMethylation450k"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (F) {
#       if (!.validateType(exp_path, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
# 
#         soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#         soft <- GEOquery::getGEO(filename = soft)
# 
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "CONTROL")] <- "Neuron"
#         cell[str_detect(cell, "myeloid")] <- "CMP"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
# 
#         # Get data
#         scm <- soft.to.scMethrix(soft = soft, colData = colData, exclude_id = exclude_id, array = array, annotation = annotation)
#         scm <- liftover_scMethrix(scm, chain, target_genome = "hg38")
#         scm <- bin_scMethrix(scm, regions = probe.set)
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       } else {
#         assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#       }
#     }
#   }
# 
#   #--- GSE144804 ----------------------------------------------------------------------------------------------
#   # GEO: GSE144804
#   # Types: Neurons
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/32231389/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144804
#   # Citation: Rhead B, Shao X, Quach H, Ghai P et al. Global expression and CpG methylation analysis of primary endothelial cells before and after TNFa stimulation reveals gene modules enriched in inflammatory and infectious diseases and associated DMRs. PLoS One 2020;15(3):e0230884. PMID: 32231389
#   # Genome: hg19
#   # Platform: Illumina EPIC
#   if ("GSE144804" %in% GEOs) {
#     GEO <- "GSE144804"
#     array <- "EPIC"
#     message("Getting ",GEO)
# 
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0("scm.", GEO)
#     exp_path <- paste0(exp_dir, exp_name, ".rds")
# 
#     if (!.validateType(exp_path, "file", throws = F)) {
#       raw_dir <- paste0(base_dir, "raw/")
#       mkdirs(base_dir, raw_dir, exp_dir)
# 
#       soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#       soft <- GEOquery::getGEO(filename = soft)
# 
#       # Get colData
#       cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#       remove_idx <- which(str_detect(cell, "_TNF"))
#       cell[str_detect(cell, "Huvec")] <- "Endothelial"
#       colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#       colData <- colData[-remove_idx, , drop = FALSE]
#       stopifnot(all(colData$cell %in% cell_types))
#       
# 
#       # Get data
#       scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#       scm <- standardize.scMethrix(scm, GEO, chain, probe.set)
#       assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#     
#     
#     
#     array <- "EPIC"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#         remove_idx <- which(str_detect(cell, "_TNF"))
#         cell[str_detect(cell, "Huvec")] <- "Endothelial"
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#     
#     
#   }
# 
#   #XXX--- GSE151506 ----------------------------------------------------------------------------------------------
#   # GEO: GSE151506
#   # Types: GBM-IDH,GBM-WT,GBM-imm
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/34594037/
#   # GEO: Chaligne, Ronan, Federico Gaiti, Dana Silverbush, Joshua S. Schiffman, Hannah R. Weisman, Lloyd Kluegel, Simon Gritsch et al. "Epigenetic encoding, heritability and plasticity of glioma transcriptional cell states." Nature Genetics 53, no. 10 (2021): 1469-1479.
#   # Genome: hg38
#   # Platform: RRBS
# if(F)  if ("GSE151506" %in% GEOs) {
#     GEO <- "GSE151506"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "CONTROL")] <- "Neuron"
#         cell[str_detect(cell, "myeloid")] <- "CMP"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- soft.to.scMethrix(soft = soft, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE164149 --------------------------------------------------------------------------------------------
#   # Types: Treg
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/34348163/
#   # GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164149
#   # Citation: Lam, Avery J., et al. "Optimized CRISPR-mediated gene knockin reveals FOXP3-independent maintenance of human Treg identity." Cell Reports 36.5 (2021): 109494.
#   # Genome: hg19
#   # Platform: EPIC
#   if ("GSE164149" %in% GEOs) {
#     GEO <- "GSE164149"
#     array <- "EPIC"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
#         cell <- paste0("###-", cell)
#         cell[str_detect(cell, "Cas9")] <- "Treg"
#         remove_idx <- which(str_detect(cell, "###-"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- raw.idat.to.scMethrix(GEO = GEO, raw_dir = raw_dir, colData = colData, array = array)
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# 
#   #--- GSE166844 ----------------------------------------------------------------------------------------------
#   # GEO: GSE166844
#   # Types: CD4 T-cells, CD8 T-cells, Monocytes, Granulocyte, B-cells, Whole blood
#   # Paper: https://pubmed.ncbi.nlm.nih.gov/33739972/
#   # GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166844
#   # Citation: Hannon E, Mansell G, Walker E, Nabais MF et al. Assessing the co-variability of DNA methylation across peripheral cells and tissues: Implications for the interpretation of findings in epigenetic epidemiology. PLoS Genet 2021 Mar;17(3):e1009443. PMID: 33739972
#   # Genome: hg19
#   # Platform: Illumina EPIC
#   if ("GSE166844" %in% GEOs) {
#     GEO <- "GSE166844"
#     array <- "EPIC"
#     message("Getting ",GEO)
#     
#     base_dir <- paste0(home_dir, GEO, "/")
#     exp_dir <- paste0(base_dir, "exp/")
#     exp_name <- paste0(c("scm",GEO,region), collapse=".")
#     exp_all <- paste0(exp_dir,"scm.",GEO,".all.rds")
#     exp_path <- paste0(exp_dir,exp_name,".rds")
#     
#     if (!.validateType(exp_path, "file", throws = F)) {
#       
#       if (!.validateType(exp_all, "file", throws = F)) {
#         raw_dir <- paste0(base_dir, "raw/")
#         mkdirs(base_dir, raw_dir, exp_dir)
#         
#         if (exists("soft") && !identical(soft@header$geo_accession, GEO)) {
#           soft <- GEOquery::getGEOfile(GEO, destdir = raw_dir)
#           soft <- GEOquery::getGEO(filename = soft)
#         }
#         
#         # Get colData
#         ids <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
#         cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
#         cell[str_detect(cell, "Monocyte")] <- "Monocyte"
#         cell[str_detect(cell, "B-cells")] <- "Bcell"
#         cell[str_detect(cell, "CD4")] <- "CD4Tcell"
#         cell[str_detect(cell, "CD8")] <- "CD8Tcell"
#         cell[str_detect(cell, "blood")] <- "WholeBlood"
#         cell[str_detect(cell, "Granulocyte")] <- "Granulocyte"
#         cell[str_detect(cell, "CD8")] <- "CD8Tcell"
#         remove_idx <- which(str_detect(cell, "Nasal|Buccal"))
#         exclude_id <- names(cell[remove_idx])
#         colData <- data.frame(row.names = names(soft@gsms), Cell = cell, ID = ids)
#         colData <- colData[-remove_idx, , drop = FALSE]
#         stopifnot(all(colData$cell %in% cell_types))
#         
#         # Get data
#         scm <- var.proc.to.scMethrix(
#           GEO = GEO, raw_dir = raw_dir, colData = colData, array = array,
#           proc_file = "GSE166844_Variance_processed_signals.csv.gz"
#         )
#         scm <- standardize.scMethrix(scm, GEO, chain, region = NULL)
#         scm <- scMethrix::save_scMethrix(scm, dest = exp_all)
#       } else {
#         scm <- scMethrix::load_scMethrix(dest = exp_all)
#       }
#       
#       if (region == "all") {
#         assign(exp_name, scm, envir = .GlobalEnv)
#       } else { 
#         scm <- standardize.scMethrix(scm, GEO = GEO, chain = chain, bin_region = probe.set[[region]])
#         assign(exp_name, scMethrix::save_scMethrix(scm, dest = exp_path), envir = .GlobalEnv)
#       }
#     } else {
#       assign(exp_name, scMethrix::load_scMethrix(exp_path), envir = .GlobalEnv)
#     }
#   }
# }
