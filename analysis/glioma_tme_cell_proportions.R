home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"tumorprop/")
exp_name <- "scm.tumorprop.base"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)

cells <- c("GBM","Glioma")
scm.tumors <- get_data_set(c("GSE104293","GSE103659"), merge=T, genome = "hg19", cells = cells, regions = "ill.450k.raw.hg19")
scm.tumors <- impute_regions(scm.tumors)
mtx <- get_matrix(scm.tumors,"impute")

match_idx <- match(start(rowRanges(scm.tumors)),start(probes[["ill.450k.raw.hg19"]]))

rownames(mtx) <- probes[["ill.450k.hg19"]]$ID
mtx <- as.data.table(mtx,keep.rownames=TRUE)
setorder(mtx, rn)

query <- TCGAbiolinks::GDCquery(project= "TCGA-GBM",
                                data.category = "DNA methylation",
                                barcode = c("TCGA-26-1442"),
                                platform = "Illumina Human Methylation 450",
                                legacy = TRUE)
GDCdownload(query)
hg19.450k <- GDCprepare(query)
hg19.450k <- as.data.table(assays(hg19.450k)[[1]],keep.rownames=TRUE)

mtx2 <- merge(hg19.450k,mtx,by="rn",all.x=TRUE)
mtx2[ , 2 := NULL ]
mtx2 <- as.matrix(mtx2,rownames=1)

# classification <- gliomaClassifier(mtx2)
classification <- load_scMethrix("D:\\Git\\thesis_data\\tumorprop\\classification.rds")
 gsm.class <- data.frame(row.names = classification$final.classification$Sample, Class = classification$final.classification$Final_classification, Model = classification$model.classifications$glioma.idh.model)


scm.tumors <- get_data_set(c("GSE104293","GSE103659"), merge=T, genome = "hg38", cells = cells, regions = "ill.450k.hg38")
scm.tumors <- impute_regions(scm.tumors)

# 
# hmap <- generate_heatmap(scm.tcga.probe,assay="impute")
# hmap
# scm.tcga.probe <- subset_scMethrix(scm.tumors,regions = probes.glm[["hg38"]])

colData(scm.tumors) <- cbind(colData(scm.tumors),gsm.class)

scm.tcga.probe <- subset_scMethrix(scm.tumors,regions = probes[["glm.hg19"]])

scm.tcga.probe <- dim_red_scMethrix(scm.tcga.probe,type="tSNE",assay="impute",perplexity=75,max_iter=2000)
plot1 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Cell")
plot2 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Model")
plot3 <- plot_dim_red(scm.tcga.probe,"tSNE",color_anno="Class")

grid.arrange(plot1, plot2, plot3, ncol=3)

scm.tcga.probe <- dim_red_scMethrix(scm.tcga.probe,type="PCA",assay="impute",perplexity=85,max_iter=2000)
plot1 <- plot_dim_red(scm.tcga.probe,"PCA",color_anno="Cell")
plot2 <- plot_dim_red(scm.tcga.probe,"PCA",color_anno="Model")
plot3 <- plot_dim_red(scm.tcga.probe,"PCA",color_anno="Class")
grid.arrange(plot1, plot2, plot3, ncol=3)

scm.tcga.probe <- dim_red_scMethrix(scm.tcga.probe,type="UMAP",assay="impute",perplexity=85,max_iter=2000)
plot1 <- plot_dim_red(scm.tcga.probe,"UMAP",color_anno="Cell")
plot2 <- plot_dim_red(scm.tcga.probe,"UMAP",color_anno="Model")
plot3 <- plot_dim_red(scm.tcga.probe,"UMAP",color_anno="Class")
grid.arrange(plot1, plot2, plot3, ncol=3)

# Feature select for model  ------------------------------------------------------------------------------------------------------------

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

# Feature select for class ----------------------------------------------------------------------------------------

names(rowRanges(scm)) <- paste0("rid_",1:nrow(scm))
scm.feat <- do_methylcibersort(scm, assay="impute",MaxDMPs = 200, deltaBeta = 0.25,col="Class")

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
feats.top <- readRDS(file.features)
sig.mtx <- read.delim("D:\\Git\\thesis_data\\deconvolution\\methylCibersort_0.25_200_Signature.txt", row.names = 1, header = TRUE)
scm <- bin_scMethrix(scm.tumors,regions = feats.top,fill=T)
scm <- impute_regions(scm,k=7)
names(rowRanges(scm)) <- names(feats.top)
# scm.bin <- standardize.scMethrix(scm.tumors,GEO="",chain=chain,probe.set = probe.set)
betas <- get_matrix(scm,assay="impute")
betas <- betas[rownames(sig.mtx), ]*100
results <- CIBERSORT(sig_matrix = sig.mtx,
                     mixture_file = as.data.frame(betas),
                     perm = 100,
                     QN = TRUE,
                     absolute = FALSE,
                     abs_method = 'sig.score')

# graph the immune portion --------------------------------------------------------------------------------------------
scores <- as.matrix(results)[,1:(ncol(results)-3)]
melt <- setNames(reshape2::melt(scores), c('GSM', 'CellFound', 'Value'))
cd <- colData(scm.tumors)
cd$GSM <- row.names(cd)

cellprops <- as.data.frame(merge(melt,cd, by="GSM"))

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Cell)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Class)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")

ggplot(cellprops, aes(x=CellFound, y=Value, fill=Model)) + 
  geom_boxplot() +
  facet_wrap(~CellFound, scale="free")
