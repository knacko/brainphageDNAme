home_dir <- "D:/Git/thesis_data/"
mkdirs(home_dir)
proms <- probes.ill

# Create merged object to generate the signature matrix -----------------------------------------------------------
exp_dir <- paste0(home_dir,"deconvolution/")
exp_name <- "scm.deconvolution.base"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)

cells <- c("Bcell", "CD4Tcell", "CD8Tcell", "Treg", "NKcell",
           "Microglia", "Dendritic", "Monocyte","Granulocyte","WholeBlood", "Neuron", "Glia", "Glioma",
           "Endothelial")
# 
# cells <- c("WholeBlood", "Granulocyte", "CD4Tcell", "CD8Tcell", "Monocyte", 
#            "Bcell", "NKcell")

#if (!.validateType(exp_path,"file",throws=F)) {
  
  exps <- c("singh","GSE35069","GSE49667","GSE66351",
            "GSE83458","GSE88824","GSE103211","GSE128654")#,"GSE151506.mgsig")
  # exps <- get_data_set(exps, regions = "proms.hg38", assign = F)
  # #scms <- lapply(exps,function(exp) get(paste0("scm.",exp)))
  # 
  scm <- get_data_set(exps, regions = "proms.hg38", merge=T, cells = cells)
  # 
  # if (length(exps) > 1) {
  #   for (i in 2:length(exps)) {
  #     message("Adding ",names(exps[i]))
  #     scm <- merge_scMethrix(scm, exps[[i]], by = "col",verbose=F)
  #   }
  # }
  
  # remove(list = paste0("scm.",exps))
  # remove(scms,exps)

  # col_idx <- sapply(colData(scm)$Cell, `%in%`, cells)
  # scm <- scm[,col_idx]
  
#   assign(exp_name, scMethrix::save_1scMethrix(scm,dest = exp_path))
# } else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}
# 
# scm <- get(exp_name)
# beep()
# shapes <- c(21,22,23,24,21,22,23,24,21,22,23,21,22,23,24,25) # Ensure shapes of cell type groups are different
# names(shapes) <- cells
# colors <- c('#ffe119', '#4363d8', '#3cb44b', '#e6194B', '#a9a9a9', '#000000')
# colors <- sapply(shapes,function(x) colors[[x-20]])
# names(colors) <- cells
# 
# colData(scm)$Shape <- shapes[colData(scm)$Cell]
# colData(scm)$Color <- colors[colData(scm)$Cell]
# 
# dim_red_theme <- quote(scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + scale_fill_discrete(cells))

# Do feature finder -----------------------------------------------------------------------------------------------
exp_name <- "scm.deconvolution.topfeat"
exp_path <- paste0(exp_dir,exp_name,".rds")
mkdirs(exp_dir)
setwd(exp_dir)

# if (!.validateType(exp_path,"file",throws=F)) {

# cells <- c("WholeBlood", "Granulocyte", "CD4Tcell", "CD8Tcell", "Monocyte", 
#   "Bcell", "NKcell")


  DMPs = 300
  set.seed("123")
 
  k = round(sqrt(ncol(scm))/2)
  if((k %% 2) == 0) k = k+1
  message("k = ",k)

  scm <- impute_regions(scm,k=k)
  names(rowRanges(scm)) <- paste0("rid_",1:nrow(scm))
  scm <- do_methylcibersort(scm, assay="impute",MaxDMPs = DMPs, deltaBeta = 0.25)

  # assign(exp_name, scMethrix::save_scMethrix(scm,dest = exp_path))
# s} else {assign(exp_name,scMethrix::load_scMethrix(exp_path))}

# scm <- get(exp_name)

file.sigmtx <- list.files(getwd(), full.names = TRUE, pattern = paste0(".*",DMPs,"_Signature.txt$"), ignore.case = T)
feats.top <- rowRanges(scm)
# beep()
# Graph the features ----------------------------------------------------------------------------------------------

# ord <- unlist(sapply(cells,function(c)which(cd(scm)$Cell == c)))
# names(ord) <- NULL
# scm <- scm[,ord]

# generate_heatmap <- function(scm,assay = "score",grouping = NULL, n_cpg = NULL,...) {
#   
#   if (!is.null(n_cpg)) scm <- reduce_scMethrix(scm,assay=assay,n_cpg = n_cpg,var="rand")
#   
#   mat <- as.matrix(get_matrix(scm,assay))
#   type <- colData(scm)$Cell
#   # ha = HeatmapAnnotation(
#   #   df = data.frame(type = type),
#   #   foo = anno_block(gp = gpar(fill = get_palette(length(unique(colData(scm)$Cell))))),
#   #   annotation_height = unit(8, "mm"),gp = gpar(col = "black")
#   # )
#   
#   Heatmap(mat, name = "β-value", km = 5, column_title_rot = 90, border = TRUE,#top_annotation = ha, 
#           column_split = colData(scm)$Cell,
#           #column_names_gp = 14,#cluster_columns = agnes(t(mat)), 
#           show_row_names = FALSE, row_title = NULL, column_dend_side = "bottom", column_names_side = "top", show_column_names = FALSE,...) 
#   
# }
# 
# 
# hmap <- generate_heatmap(scm,assay="impute")
# draw(hmap)
# 
# Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\sigmtx_hmp.png",
#              type="png",
#              units="px", 
#              width=500, 
#              height=500, 
#              pointsize=20, 
#              dpi="auto")
# 
# draw(hmap, padding = unit(c(2, 2, 10, 2), "mm"))
# 
# dev.off()
# 
# 
# legend.text=element_text(size=X)
# 
# 
# 
# 
# scm <- dim_red_scMethrix(scm,type="tSNE",assay="impute",perplexity=65,max_iter=2000)
# dimred <- plot_dim_red(scm,"tSNE",color_anno="Cell",axis_labels=list(X="tSNE1",Y="tSNE2"),base_size=20,legend.title = element_blank(),legend.position = "none",
#              #legend.justification="bottom",
#              legend.margin=margin(5,5,5,5),
#              legend.box.margin=margin(-18,0,0,0),
#              legend.box.background = element_blank(),
#              legend.background = element_blank(),
#              legend.text=element_text(size=14)) 
# 
# 
# Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\sigmtx_tsne.png",
#              type="png",
#              units="px", 
#              width=500, 
#              height=500, 
#              pointsize=24, 
#              dpi="auto")
# 
# dimred
# 
# dev.off()










# 
# 
# 
#   f <- function(x) {deparse(substitute(x))}
#   g <- function(x) {eval(parse(text=x))}
# 
#   
#   hull <- geom_mark_hull(expand=0.02,aes(fill=Color,label=NULL))
#   
#   g(paste0(f(plot), " + ", f(hull)))
#   
# 
# q <- plot + list(geom_polygon(data = Color, alpha = 0.5))
# 
# feat <- dim_red_scMethrix(feat,type="UMAP",assay="impute")
# plot_dim_red(feat,"UMAP",color_anno="Cell",axis_labels=list(Y="",X="PCA")) 
# 
# reds <- list()
# for (i in 5*1:25) {
#   feat <- dim_red_scMethrix(feat,type="tSNE",assay="impute",perplexity=i,max_iter=5000)
#   reds[[paste0("P",i)]] <- reducedDim(feat,"tSNE")
#   print(plot_dim_red(feat,"tSNE",color_anno="Cell",axis_labels=list(Y="",X=paste0("Perplexity=",i))))
# }
# 
# feat <- dim_red_scMethrix(feat,type="PCA",assay="impute")
# plot_dim_red(feat,"PCA",color_anno="Cell",axis_labels=list(Y="",X="PCA"))

## Make a synthetic mix of GSE35069 --------------------------------------------------------------------------------

scm <- get_data_set(c("GSE110554"), cells = cells,merge = T)
scm <- bin_scMethrix(scm,regions = feats.top,fill=T)
scm <- impute_regions(scm,k=7)
names(rowRanges(scm)) <- names(feats.top)
betas <- get_matrix(scm,assay="impute")
sig.mtx <- read.delim(file.sigmtx, row.names = 1, header = TRUE)
#row.names(b.sig) <- paste0("id",1:nrow(b.sig))
betas <- betas[rownames(sig.mtx), ]
make.synth.mix(input.data = betas*100,
               pop.factor = as.factor(colData(scm)$Cell),
               pop.rows = 10,
               n.cores = 1,
               output.name = "GSE35069")

file.synmix <- list.files(getwd(), full.names = TRUE, pattern = ".*GSE35069_synth_mix.txt$", ignore.case = T)
file.synprop <- list.files(getwd(), full.names = TRUE, pattern = ".*GSE35069_synth_prop.txt$", ignore.case = T)

scm.synmix <- scm

## Test the signatures ---------------------------------------------------------------------------------------------
sig.mtx = read.delim(file.sigmtx, row.names = 1, header = TRUE)
syn.mix = read.delim(file.synmix, row.names = 1, header = TRUE, sep="\t")

stopifnot(rownames(syn.mix) == rownames(sig.mtx))

results <- CIBERSORT(sig_matrix = sig.mtx,
                     mixture_file = syn.mix,
                     perm = 200,
                     QN = TRUE,
                     absolute = FALSE,
                     abs_method = 'sig.score')

props <- read.delim(file.synprop, header = TRUE,sep=" ")
row.names(props) <- paste0("synth.",1:nrow(props))

stats <- results[,(ncol(results)-2):ncol(results)]
scores <- as.matrix(results)[,1:(ncol(results)-3)]

scores <- scores[,colnames(props)]

props.melt <- setNames(reshape2::melt(as.matrix(props)), c('Mix', 'Cell', 'X'))
scores.melt <- setNames(reshape2::melt(scores), c('Mix', 'Cell', 'Y'))

syn.stats <- merge(scores.melt,props.melt,by=c("Mix","Cell"))
syn.stats$X <- round(syn.stats$X,3)
syn.stats$Y <- round(syn.stats$Y,3)

corr <- cor.test(c(as.matrix(props)), c(as.matrix(scores)), method="pearson")
corr.lbl <- paste("Spearman R = ",round(corr$estimate,3),"\np << 0.001")


prop.vec <- c(as.matrix(props))
res.vec <- c(as.matrix(scores))

gg <- ggplot(syn.stats, aes(x=X, y=Y,fill = Cell)) + geom_point(shape = 21, size =2, alpha = .7) + 
  geom_abline(intercept = 0, slope = 1) + annotate("text", x=0.35, y=0.75, label=corr.lbl) +xlab("Synthetic Proportion") + ylab("Calculated Proportion") + scMethrix_theme(base_size=20,legend.text=element_text(size=16))


Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\sigmtx_corr.png",
             type="png",
             units="px",
             width=500,
             height=300,
             pointsize=20,
             dpi="auto")

gg

dev.off()





  

# Make synthetic mix of PBMC --------------------------------------------------------------------------------------
# 
# get_data_set(c("GSE110554","GSE112618"))
# scm.syn <- merge_scMethrix(scm.GSE110554,scm.GSE112618,by="col")
# scm <- subset_scMethrix(scm,regions = feats.top)
# scm <- impute_regions(scm,k=7)
# scm <- scm[,which(colData(scm)$Cell == "ImmMix")]
# names(rowRanges(scm)) <- names(feats.top)
# betas <- get_matrix(scm,assay="impute")
# sig.mtx <- read.delim(file.sigmtx, row.names = 1, header = TRUE)
# #row.names(b.sig) <- paste0("id",1:nrow(b.sig))
# betas <- betas[rownames(sig.mtx), ]*100
# 
# results <- CIBERSORT(sig_matrix = sig_matrix,
#                      mixture_file = betas,
#                      perm = 100,
#                      QN = FALSE,
#                      absolute = FALSE,
#                      abs_method = 'sig.score')
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
# 
# make.synth.mix(input.data = betas,
#                pop.factor = as.factor(colData(scm)$Cell),
#                pop.rows = 15,
#                n.cores = 1,
#                output.name = "ImmMix")
# 
# file.synmix <- list.files(getwd(), full.names = TRUE, pattern = ".*ImmMix_synth_mix.txt$", ignore.case = T)
# file.synprop <- list.files(getwd(), full.names = TRUE, pattern = ".*ImmMix_synth_prop.txt$", ignore.case = T)
# 
# scm.synmix <- scm
# 
# ## Test the signatures ---------------------------------------------------------------------------------------------
# sig_matrix <- read.delim("./methylCibersort_0.2_100_Signature.txt", row.names = 1, header = TRUE)
# mixture_file <- read.delim("./test_synth_mix.txt", row.names = 1, header = TRUE)
# 
# 
# 
# 




