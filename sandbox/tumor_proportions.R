home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"tumorprop/")
exp_name <- "scm.tumorprop.base"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)

cells <- c("GBM-WT","GBM-IDH","Glioma")
# 
# if (!.validateType(exp_path,"file",throws=F)) {
#   
#   exps <- c("GSE104293","GSE103659")
#   get_data_set(exps)
#   scms <- lapply(exps,function(exp) get(paste0("scm.",exp)))
#   
#   scm <- scms[[1]]
#   for (i in 2:length(exps)) {
#     message("Adding ",exps[i])
#     scm <- merge_scMethrix(scm, scms[[i]], by = "col",verbose=F)
#   }
#   
#   remove(list = paste0("scm.",exps))
#   remove(scms,exps)
#   
#   col_idx <- sapply(colData(scm)$Cell, `%in%`, cells)
#   scm <- scm[,col_idx]
#   
#   assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
# } else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}


# Classify --------------------------------------------------------------------------------------------------------

query <- GDCquery(project= "TCGA-GBM",
                  data.category = "DNA methylation",
                  barcode = c("TCGA-06-0122","TCGA-14-1456"),
                  platform = "Illumina Human Methylation 27",
                  legacy = TRUE)
GDCdownload(query)
data.hg19 <- GDCprepare(query)
classification <- gliomaClassifier(scm.tumors)

# Sort ------------------------------------------------------------------------------------------------------------

scm.tcga.probe <- subset_scMethrix(scm.tumors,regions = probes.glm[["hg38"]])
scm.tcga.probe <- dim_red_scMethrix(scm.tcga.probe,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
print(plot1 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Cell"))
print(plot2 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Model"))
print(plot3 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Class"))

grid.arrange(plot1, plot2, plot3, ncol=3)

hmap <- generate_heatmap(scm.tcga.probe,assay="impute")
hmap

scm.dimred <- dim_red_scMethrix(scm,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Cell"))
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Model"))
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Class"))
grid.arrange(plot1, plot2, plot3, ncol=3)

scm <- dim_red_scMethrix(scm,type="UMAP",assay="impute",perplexity=85,max_iter=2000)
plot <- plot_dim_red(scm,"UMAP",color_anno="Cell") 
plot

scm <- scm.tumors

names(rowRanges(scm)) <- paste0("rid_",1:nrow(scm))
scm.feat <- do_methylcibersort(scm, assay="impute",MaxDMPs = 200, deltaBeta = 0.25,col="Model")

scm.feat <- dim_red_scMethrix(scm.feat,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
print(plot <- plot_dim_red(scm.feat,"tSNE",color_anno="Cell"))
print(plot <- plot_dim_red(scm.feat,"tSNE",color_anno="Model"))
print(plot <- plot_dim_red(scm.feat,"tSNE",color_anno="Class"))
grid.arrange(plot1, plot2, plot3, ncol=3)

scm.feat <- dim_red_scMethrix(scm.feat,type="PCA",assay="impute",perplexity=85,max_iter=2000)
print(plot <- plot_dim_red(scm.feat,"PCA",color_anno="Cell"))
print(plot <- plot_dim_red(scm.feat,"PCA",color_anno="Model"))
print(plot <- plot_dim_red(scm.feat,"PCA",color_anno="Class"))
grid.arrange(plot1, plot2, plot3, ncol=3)

# Check immune cell proportions -----------------------------------------------------------------------------------

sig.mtx <- read.delim("D:\\Git\\thesis_data\\deconvolution\\methylCibersort_0.25_200_Signature.txt", row.names = 1, header = TRUE)

scm.bin <- standardize.scMethrix(scm.tumors,GEO="",chain=chain,probe.set = probe.set)

betas <- get_matrix(scm.bin,assay="impute")

betas <- betas[rownames(sig.mtx), ]*100

results <- CIBERSORT(sig_matrix = sig.mtx,
                     mixture_file = as.data.frame(betas),
                     perm = 100,
                     QN = FALSE,
                     absolute = FALSE,
                     abs_method = 'sig.score')


# graph the immune portion --------------------------------------------------------------------------------------------


melt <- setNames(reshape2::melt(results), c('GSM', 'CellFound', 'Value'))
cd <- colData(scm.tumors)
cd$GSM <- row.names(cd)

cellprops <- as.data.frame(merge(melt,cd, by="GSM"))

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Model)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Class)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

