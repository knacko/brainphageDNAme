home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"/deconvolution/")
exp_name <- "scm.deconvolution.rds"
exp_path <- paste0(exp_dir,exp_name)
mkdirs(exp_dir)

if (!.validateType(exp_path,"file",throws=F)) {
  scm.sig <- list(scm.GSE35069,scm.GSE88824,scm.GSE110554,scm.GSE66351,scm.GSE49667,scm.GSE164149,
                    scm.GSE128654,scm.GSE144804,scm.singh,scm.GSE103211)
  scm <- Reduce(function(x,y) {merge_scMethrix(x, y, by = "col",verbose=F)}, scm.sig)
  assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
} else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}


cells <- c("Bcell", "CD4Tcell", "CD8Tcell", "Treg", "NKcell","Eosinophil", "Neutrophil",
             "Granulocyte", "Microglia", "Dendritic", "Monocyte", "Neuron", "Glia", "Glioma",
               "WholeBlood","Endothelial")


  
shapes <- c(21,22,23,24,21,22,23,24,21,22,23,21,22,23,24,25) # Ensure shapes of cell type groups are different
colors <- c('#ffe119', '#4363d8', '#a9a9a9', '#ffffff', '#000000')
colors <- sapply(shapes,function(x) colors[[x-20]])

dim_red_theme <- quote(scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + scale_fill_discrete(cells))

col_idx <- sapply(colData(scm)$Cell, `%in%`, cells)
scm <- scm[,col_idx]

# Do feature finder -----------------------------------------------------------------------------------------------

set.seed("123")

k = round(sqrt(ncol(scm)))
if((k %% 2) == 0) k = k+1

feat <- impute_regions(scm,k=5)
generate_heatmap(feat,assay="impute")

# feat <- transform_assay(feat,assay="impute",new_assay="bin",trans=binarize)
feat <- do_methylcibersort(feat, assay="impute",MaxDMPs = 80)

scale_colour_manual(name = "Cell",values = myColors)

feat <- dim_red_scMethrix(feat,type="UMAP",assay="impute")
plot_dim_red(feat,"UMAP",color_anno="Cell")

feat <- dim_red_scMethrix(feat,type="tSNE",assay="impute",perplexity=85,max_iter=5000)
print(plot_dim_red(feat,"tSNE",color_anno="Cell",axis_labels=list(Y="",X=paste0("Perplexity=",i))) + geom_mark_hull(expand=0.01,aes(fill=Cell)))

reds <- list()
for (i in 5*1:25) {
  feat <- dim_red_scMethrix(feat,type="tSNE",assay="impute",perplexity=i,max_iter=5000)
  reds[[paste0("P",i)]] <- reducedDim(feat,"tSNE")
  print(plot_dim_red(feat,"tSNE",color_anno="Cell",axis_labels=list(Y="",X=paste0("Perplexity=",i))))
}

feat <- dim_red_scMethrix(feat,type="PCA",assay="impute")
plot_dim_red(feat,"PCA",color_anno="Cell")











# Make a synthetic mix --------------------------------------------------------------------------------------------

#source("./synth_mix_2019.R")
require(FlowSorted.Blood.450k)
idx <- which(FlowSorted.Blood.450k$CellType %in% c("Bcell","CD4T","CD8T","NK"))
bvals <- getBeta(FlowSorted.Blood.450k[, idx])
b.sig <- read.delim("./methylCibersort_0.2_100_Signature.txt", row.names = 1, header = TRUE)
bvals <- bvals[rownames(b.sig), ]
all(rownames(bvals)==rownames(b.sig))
a <- Sys.time()
make.synth.mix(input.data = bvals,
               pop.factor = as.factor(FlowSorted.Blood.450k$CellType[idx]),
               pop.rows = 10,
               n.cores = 1,
               output.name = "test")
b <- Sys.time()
b-a
## Time difference of 1.76159 secs




# Test the signatures ---------------------------------------------------------------------------------------------
results <- CIBERSORT(sig_matrix = "./methylCibersort_0.2_100_Signature.txt",
                     mixture_file = "./test_synth_mix.txt",
                     perm = 1000,
                     QN = FALSE,
                     absolute = FALSE,
                     abs_method = 'sig.score')