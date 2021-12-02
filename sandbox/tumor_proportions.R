home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"tumorprop/")
exp_name <- "scm.tumorprop.base"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)

cells <- c("GBM","GBM","Glioma")

# if (!.validateType(exp_path,"file",throws=F)) {

scm.tcga.probe <- get_data_set(c("GSE104293","GSE103659"),merge=T, cells = cells)

scm.tcga.probe <- impute_regions(scm.tcga.probe)


 
# 
#   assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
# } else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}



# Classify --------------------------------------------------------------------------------------------------------

query <- TCGAbiolinks::GDCquery(project= "TCGA-GBM",
                  data.category = "DNA methylation",
                  barcode = c("TCGA-06-0122","TCGA-14-1456"),
                  platform = "Illumina Human Methylation 27",
                  legacy = TRUE)
GDCdownload(query)
data.hg19 <- GDCprepare(query)
classification <- gliomaClassifier(data.hg19)

# Sort ------------------------------------------------------------------------------------------------------------

scm.tcga.probe <- subset_scMethrix(scm.tcga.probe,regions = probes.ill[["i450k.hg19"]])                 

row_idx <- match(start(probes.ill[["i450k.hg19"]]), start(rowRanges(scm.tcga.probe)))

scm.tcga.probe <- scm.tcga.probe[na.omit(row_idx),]

rowRanges(scm.tcga.probe)$ID <- NA

row.names(scm.tcga.probe) <- 

IDs <- sapply(start(rowRanges(scm)), function (st) {
  
  probes.ill[["i450k.hg19"]][match(start(rowRanges(scm)[i]),start(probes.ill[["i450k.hg19"]])) ]$ID 
  
  # rowRanges(scm)[i]$ID <- probes.ill[["i450k.hg19"]][i]$ID

})


names(rowRanges(scm)) <- rowRanges(scm)$ID
  
rowRanges(scm)$ID <- probes.ill[["i450k.hg19"]][IDs]$ID


match(start(rowRanges(scm)),start(probes.ill[["i450k.hg19"]]))
  
  
  
xx <- match(start(probes.ill[["i450k.hg19"]]), start(rowRanges(scm.tcga.probe)))

match(start(rowRanges(scm.tcga.probe)),start(probes.ill[["i450k.hg19"]]))

rowRanges(scm.tcga.probe)$ID[row_idx] <-  probes.ill[["i450k.hg19"]]$ID




rowRanges(scm.tcga.probe)$ID <- probes.ill[["i450k.hg19"]][-remov]$ID






s.s <- start(rowRanges(scm.tcga.probe))
p.s <- start(probes.ill[["i450k.hg19"]])

s.g <- intersect(s.s,p.s)

setdiff(s.s,s.g)
setdiff(p.s,s.g)







names(rowRanges(scm.tcga.probe)) <- probes.glm[["hg38"]]$ID

classification <- gliomaClassifier(get_matrix(scm.tcga.probe))

hmap <- generate_heatmap(scm.tcga.probe,assay="impute")
hmap
scm.tcga.probe <- subset_scMethrix(scm.tumors,regions = probes.glm[["hg38"]])

scm.tcga.probe <- dim_red_scMethrix(scm.tcga.probe,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
print(plot1 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Cell"))
print(plot2 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Model"))
print(plot3 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Class"))

grid.arrange(plot1, plot2, plot3, ncol=3)

scm.dimred <- dim_red_scMethrix(scm,type="tSNE",assay="impute",perplexity=85,max_iter=2000)
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Cell"))
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Model"))
print(plot <- plot_dim_red(scm.dimred,"tSNE",color_anno="Class"))
grid.arrange(plot1, plot2, plot3, ncol=3)

scm <- dim_red_scMethrix(scm,type="UMAP",assay="impute",perplexity=85,max_iter=2000)
plot <- plot_dim_red(scm,"UMAP",color_anno="Cell") 
plot

# Feature select  ------------------------------------------------------------------------------------------------------------

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

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Cell)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Class)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

