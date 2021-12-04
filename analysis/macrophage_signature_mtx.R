
# Import the RNA-seq data -----------------------------------------------------------------------------------------

tpm <- data.table::fread(file="D:\\Git\\thesis_data\\GSE151506\\RNAseq\\GSM4579804_MGH105C_CD45pos.tpm.txt.gz",header = T)
tpm <- data.frame(tpm)
row.names(tpm) <- tpm$gene_name
tpm$gene_name <- NULL

cell_colors <- c("#FFA900","#0056FF")


 mg_genes <- c("P2ry12", "Slc2a5", "Tmem119","Fcrls")#"TMEM119","P2ry13","P2ry12","Gpr34","Slc2a5","Siglec-H","Olfml3","Fcrls")#"P2RY13","GPR34","SLC2A5","OLFML3",,"Siglech","Fcrls"
 mp_genes <- c("Gda" , "Hp", "Sell" , "Emilin2")#"F10","Emilin2","F5","C3","Gda","Mki67","Sell","Hp")#"EMILIN2","C3","GDA","SELL",,"TREM2","HP","F10","F5","Mki67"

# mg_genes <- c("TMEM119","p2ry12")
# mp_genes <- c("Itga4","Fn1")

genes <- c(mp_genes,mg_genes)

counts <- tpm
counts <- tpm[which(tolower(row.names(tpm)) %in% tolower(genes)),,drop=FALSE]
counts <- log(counts + 1) / log(2)
counts <- M3Drop::M3DropConvertData(counts)

# Create an scMethrix expeirment ----------------------------------------------------------------------------------
scm <- create_scMethrix(assay=list(score = counts), colData = colnames(counts), rowRanges = GRanges(seqnames = "chr1", ranges = IRanges(start = 1:nrow(counts)*100, width = 3)))

colData(scm)$mg <- round(colSums(counts[tolower(row.names(counts)) %in% tolower(mg_genes),],na.rm=T),1)
colData(scm)$mp <- round(colSums(counts[tolower(row.names(counts)) %in% tolower(mp_genes),],na.rm=T),1)
colData(scm)$fold <- round(foldchange(colData(scm)$mg, colData(scm)$mp),1)

scm <- remove_col(scm,"MGH105C_P8_F12_")

cd <- colData(scm)
cd$mg[is.infinite(cd$fold) & cd$fold > 0] <- cd$mg[is.infinite(cd$fold) & cd$fold > 0]*1000
cd$mp[is.infinite(cd$fold) & cd$fold < 0] <- cd$mp[is.infinite(cd$fold) & cd$fold < 0]*-1000
colData(scm) <- cd

MGs <- names(sort(colData(scm)$mg,decreasing=T)[1:10])
MPs <- names(sort(colData(scm)$mp)[1:10])

#scm <- scm[,which(colnames(scm) %in% c(MGs,MPs))]

counts <- tpm
#counts <- tpm[which(tolower(row.names(tpm)) %in% tolower(genes)),,drop=FALSE]
counts <- log(counts + 1) / log(2)
scm <- create_scMethrix(assay=list(score = counts), colData = colnames(counts),
                        rowRanges = GRanges(seqnames = "chr1", ranges = IRanges(start = 1:nrow(counts)*100, width = 3)))
scm <- scm[,which(colnames(scm) %in% c(MGs,MPs))]

colData(scm)$Labels <- "Other"
colData(scm)$Labels[match(MGs,colnames(scm))] <- "MG"
colData(scm)$Labels[match(MPs,colnames(scm))] <- "M\u03D5"

dim_red = "UMAP"

lapply(c("tSNE","UMAP","PCA"), function(dim_red) {
  scm <- dim_red_scMethrix(scm,type=dim_red)

  dimred <- plot_dim_red(scm,dim_red,color_anno = "Labels",axis_labels=list(X="UMAP1",Y="UMAP2"),
                         base_size=20,legend.title = element_blank(),legend.position="bottom",
                         #legend.justification="bottom",
                         legend.margin=margin(5,5,5,5),
                         legend.box.margin=margin(-18,0,0,0),
                         legend.box.background = element_blank(),
                         legend.background = element_blank(),
                         legend.text=element_text(size=20))
  print(dimred)
})

  
  
  
  
Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\microglia_selection2.png",
      type="png",
      units="px", 
      width=300, 
      height=300, 
      pointsize=24, 
      dpi="auto")

dimred
  
dev.off()


# Get the experiment ----------------------------------------------------------------------------------------------

scm.GSE151506 <- load_scMethrix("D:\\Git\\thesis_data\\GSE151506\\raw") 
scm.GSE151506 <- subset_scMethrix(scm.GSE151506,samples=c(MGs,MPs))

sapply(c(MGs,MPs),function(x) colData(scm.GSE151506)$Sample %like% x)

files <- list.files("D:\Git\thesis_data\GSE151506.mgsig\\raw", full.names=TRUE, pattern = "",ignore.case = T)

row.names(colData(scm))
names <- str_remove(row.names(colData(scm)),".R1.fastq_bismark_bt2_pe.bismark")

scmmm <- read_beds(files,colData=colData, chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4,
          cov_idx = 5,keep_cov=F)

scmmmm = collapse_samples(scmmm,colname="Collapse")

row.names(x) <- x$OldName
x$OldName <- NULL
colData <- x

colData = colData[!duplicated(colData$Collapse),]

colData(scmmmm)$Cell <- colData$Cell 
colData(scmmmm)$Name <- colData$Name

colData(scmmmm) <- colData
row.names(colData(scmmmm)) <- row.names(colData)




scm <- bin_scMethrix(scm,regions = probes.ill[["i450k.hg38.win"]])

scm.bin <- scm


scm1 <- scm[,which(colData(scm)$Cell == "Microglia")]

scm1 <- impute_regions(scm1,k=3,regions=range(rowRanges(scm1)))

scm2 <- scm[,which(colData(scm)$Cell != "Microglia")]
scm2 <- impute_regions(scm2,k=5)

scm <- merge_scMethrix(scm1,scm2,by="col")



# Get the selected fiels ------------------------------------------------------------------------------------------








# Compare to the reference microglia ------------------------------------------------------------------------------










scm <- get_data_set(c("singh","GSE151506.mgsig"), regions = "i450k.hg38.win", merge=T, cells = "Microglia")
scm1 <- scm[,which(row.names(colData(scm)) != "Singh")]
scm1 <- get_rowdata_stats(scm1)
row_idx <- which(rowData(scm1)$cells != 0)
#scm <- scm[row_idx,]

mtx <- get_matrix(scm)

for (i in 1:9) {
  # row_idx <- which(!is.na(score(scm)[,i]))
  # cor.test(score(scm)[row_idx,10], score(scm)[row_idx,i])
  
  row_idx <- 1:nrow(scm)#which(is.na(mtx[,i]))
  row_idx <- sample(row_idx,length(row_idx)*.90, replace = F)
  
  mtx[row_idx,i] <- mtx[row_idx,10]
  
}

assay(scm, "fill", withDimnames=FALSE) <- mtx

scm <- impute_regions(scm,assay="fill",k=10)

assay(scm, "score", withDimnames=FALSE) <- get_matrix(scm,"impute")






mtx <- get_matrix(scm,"impute")







xx <- get_data_set("GSE151506.mgsig")
all.dens <- plot_density(xx,na.rm=T,base_size = 20)
trimodal <- plot_density(scm,"impute",base_size = 20)
generate_heatmap(scm,assay="impute", n_cpg =  5000)




Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\microglia_sig.png",
             type="png",
             units="px", 
             width=1025, 
             height=300, 
             pointsize=24, 
             dpi="auto")

egg::ggarrange(dimred,all.dens, trimodal, ncol = 3, widths = c(300,300,300),labels = c("  A", "  B", "  C"),label.args = list(gp = grid::gpar(fontface = "bold", fontsize =20)))  

dev.off()



# Heatmap ---------------------------------------------------------------------------------------------------------





generate_heatmap <- function(scm,assay = "score", type_anno="Sample", n_cpg = NULL, grouping = NULL, ...) {
  
  if (!is.null(n_cpg)) scm <- reduce_scMethrix(scm,assay=assay,n_cpg = n_cpg,var="rand")
  
  mat <- as.matrix(get_matrix(scm,assay))
  type <- colData(scm)[[type_anno]]
  # ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=c("#FF0000","#0000FF")),
  #                                         labels = c("Serum", "2i"),labels_gp = gpar(col = "white", fontsize = 20))#,
  #                        # df = data.frame(type = type)
  #                        #annotation_height = unit(8, "mm")
  # )
  
  Heatmap(mat, name = "β-value", km = 5, #top_annotation = ha,# cluster_columns = agnes(t(mat)), 
          show_row_names = FALSE, column_dend_side = "top", column_names_side = "top", show_column_names = FALSE,
          column_title = NULL, row_title = NULL, column_split = colData(scm)$Medium) 
}

hm <- generate_heatmap(scm, n_cpg =  5000)
hm


Cairo::Cairo(file="D:\\Documents\\School\\Thesis\\Report\\microglia_heatmap.png",
             type="png",
             units="px", 
             width=300, 
             height=300, 
             pointsize=24, 
             dpi="auto")

hm

dev.off()




# Concordence -----------------------------------------------------------------------------------------------------

for (i in 1:9) {
  print(cor.test(get_matrix(scm,"impute")[,i],get_matrix(scm,"impute")[,10]))
}




reduce(scm,n_=10000)

mtx <- get_matrix(reduce_scMethrix(scm,assay="impute",n_cpg = 25000),"impute")







