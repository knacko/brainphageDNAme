home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"lineage/")
exp_name <- "scm.lineage.base"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)

cells <- c("Bcell", "CD4Tcell", "CD8Tcell", "Treg", "NKcell","Granulocyte",
           "Microglia", "Dendritic", "Monocyte", "CMP","GMP","cMOP","HSCb","HSCm","MPPb","MPPm")

if (!.validateType(exp_path,"file",throws=F)) {
  
  exps <- c("singh","GSE35069","GSE49667","GSE121483","GSE83458","GSE88824","GSE103211","GSE87196")
  exps <- get_data_set(exps,"proms.hg38")
  scms <- lapply(exps,function(exp) {
    m <- get(paste0("scm.",exp))
    names(rowRanges(m)) <- paste0("rid_",1:nrow(scm))
    m
  })
  
  scm <- scms[[1]]
  for (i in 2:length(exps)) {
    message("Adding ",exps[i])
    scm <- merge_scMethrix(scm, scms[[i]], by = "col",verbose=F)
  }
  
  remove(list = paste0("scm.",exps))
  remove(scms,exps)
  
  col_idx <- sapply(colData(scm)$Cell, `%in%`, cells)
  scm <- scm[,col_idx]
  
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
} else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}

scm <- get(exp_name)
beep()

# Do feature finder -----------------------------------------------------------------------------------------------

exp_name <- "scm.lineage.topfeat"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)
setwd(exp_dir)

prog <-  c("CMP","GMP","cMOP","HSCb","HSCm","MPPb","MPPm","CLP")
idx <- (colData(scm)$Cell %in% prog)

scm.prog <- scm[,which(idx)]

if (!.validateType(exp_path,"file",throws=F)) {
  
  set.seed("123")
  
  # k = round(sqrt(ncol(scm.prog))/2)
  # if((k %% 2) == 0) k = k+1
  # message("k = ",k)
  # 
  # scm.prog <- impute_regions(scm.prog,k=k)
  names(rowRanges(scm.prog)) <- paste0("rid_",1:nrow(scm.prog))
  scm.prog <- do_methylcibersort(scm.prog, assay="score",MaxDMPs = 500, deltaBeta = 0.25)
  
  assign(exp_name, scMethrix::save_scMethrix(scm.prog,dest = exp_path))
} else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}

#scm <- get(exp_name)

scm.prog <- subset_scMethrix(scm,regions = rowRanges(scm.prog))

k = round(sqrt(ncol(scm.prog))/2)
if((k %% 2) == 0) k = k+1
scm.prog <- impute_regions(scm.prog,k=k)


file.sigmtx <- list.files(getwd(), full.names = TRUE, pattern = ".*Signature.txt$", ignore.case = T)
feats.top <- rowRanges(scm)
beep()

# Graph the features ----------------------------------------------------------------------------------------------

hmap <- generate_heatmap(scm.prog,assay="impute")
hmap

scm.prog <- dim_red_scMethrix(scm.prog,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
plot_dim_red(scm.prog,"tSNE",color_anno="Cell")


# Find the lineage ------------------------------------------------------------------------------------------------




red.dim <- reducedDim(scm,"tSNE")
lineages <- slingshot(data=red.dim,start.clus=c("HSCm","HSCb"), clusterLabels = colData(scm)$Cell)



