gbm.big <- scMethrix::load_HDF5_scMethrix("D:/Git/thesis_data/GSE151506/exp/")
gbm.collapsed <- scMethrix::load_HDF5_scMethrix("D:/Git/thesis_data/GSE151506/gbm.collapsed/")
gbm.collapsed <- subset_scMethrix(gbm.collapsed,regions = disjointWindow(probes.ill[["i450k.hg38"]],10))
#quickResaveHDF5SummarizedExperiment(scm.big, verbose=FALSE)


# Renaming rows to [Patient]_[Cell number] ------------------------------------------------------------------------
# 
# cd <- colData(scm.big.immune)
# cd$ID <- cd$Patient
# 
# for (patient in unique(cd$Patient)) {
#   row_idx <- which(cd$Patient == patient)
#   ids <- paste0(cd$Patient[row_idx],"_C",1:length(row_idx))
#   ids <- str_replace(ids,"MGH","P")
#   cd$ID[row_idx] <- ids
# }

# rownames(colData(scm.big.immune)) <- cd$ID

# Collapse non-immune samples -------------------------------------------------------------------------------------

# scm.big.cancer <- scm.big[,colData(scm.big)$Cell != "Immune"]
# scm.big.cancer <- collapse_samples(scm.big.cancer,colname="Patient",batch_size = 1000000)
# scm.big.cancer <- scMethrix::save_HDF5_scMethrix(scm.big.cancer,"D:/Git/thesis_data/GSE151506/gbm.cancer.collapsed")
# scm.big.immune <- scm.big[,colData(scm.big)$Cell == "Immune"]
# scm.collapsed <- merge_scMethrix(scm.big.cancer,scm.big.immune,by="col")
# scm.collapsed <- scMethrix::save_HDF5_scMethrix(scm.collapsed,"D:/Git/thesis_data/GSE151506/gbm.collapsed")



# Dim red the collapsed -------------------------------------------------------------------------------------------

gbm <- bin_scMethrix(gbm.collapsed,regions = disjointWindow(probes.ill[["i450k.hg38"]],1000))


#---- Replicate Gaiti data

scm <- readRDS("D:/Git/thesis_data/GSE151506/bin/scm_binned_to_1300_glioma_probes_my_code.rds")
scm <- readRDS("D:/Git/thesis_data/GSE151506/bin/scm_binned_to_1300_glioma_probes_gaiti_code.rds")

scm.1300 <- subset_scMethrix(scm.big,regions = glm_probes_mani)
summ <- get_region_summary(scm.1300,regions = glm_probes_mani)

regions = TCGA_HG38_GR
regions = glm_probes_mani

scm <- subset_scMethrix(scm.big, regions = regions)
scm <- bin_scMethrix(scm,regions = regions, bin_size = NULL, batch_size = 100)
beep()

scm <- subset_scMethrix(scm,regions = glm_probes_mani)

scm <- scm.1056
#scm <- scm2

scm <- convert_HDF5_scMethrix(scm)

scm <- impute_regions(scm,assay="score",type="kNN",k=5,regions=glm_probes_mani)

scm <- impute_regions(scm,assay="score",type="kNN",k=5)
scm <- scm[,-which(is.na(get_matrix(scm,"impute")[1,]))]
#scm <- mask_by_variance(scm,assay="impute",)
#scm <- remove_uncovered(scm)
scm <- transform_assay(scm,assay="impute", new_assay="binarize",trans = binarize)

#assay="binarize"
assay="impute"

scm <- dim_red_scMethrix(scm,assay=assay,type="UMAP")
plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cell",)

scm <- dim_red_scMethrix(scm,assay=assay,type="PCA")
plot_dim_red(scm,dim_red = "PCA",color_anno = "Cell")


plots <- lapply(c(10,25,50,75,100,125,150,200), function(perp) {
scm <- dim_red_scMethrix(scm,assay=assay,type="tsne",perplexity=perp)
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Cell")
})

scm <- cluster_scMethrix(scm,assay="impute",type="part",n_clusters=3)

plot_dim_red(scm,dim_red = "PCA",color_anno = "Patient", shape_anno="Cell")
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Patient", shape_anno="Cell")
plot_dim_red(scm,dim_red = "UMAP",color_anno = "Patient", shape_anno="Cell")

plot_dim_red(scm,dim_red = "PCA",color_anno = "Cluster")
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Cluster")
plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cluster")

# 10 bp window

scm <- subset_scMethrix(scm.big,regions = glm_probes)
scm <- bin_scMethrix(scm,regions = glm_probes, bin_size = NULL, batch_size = 100)

scm <- convert_HDF5_scMethrix(scm)

scm <- impute_regions(scm,type="knn")

scm <- dim_red_scMethrix(scm,assay="impute",type="UMAP")
scm <- dim_red_scMethrix(scm,assay="impute",type="PCA")
scm <- dim_red_scMethrix(scm,assay="impute",type="tSNE")

plot_dim_red(scm,dim_red = "PCA",color_anno = "Cell")
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Cell")
plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cell")





#---- 

scm <- subset_scMethrix(scm.big,regions = anno450k)
scm <- bin_scMethrix(scm,regions = hg38_proms)
scm <- get_rowdata_stats(scm)
cpg_idx <- which(rowData(scm)$n_cpgs <= 5)
scm <- scm[cpg_idx,]
scm <- mask_by_stat(scm,threshold = 10, by="row", stat="count",op ="<")
scm <- mask_by_stat(scm,threshold = 0.05, by = "row", stat="var", op = "<")
scm <- convert_HDF5_scMethrix(scm)
scm <- remove_uncovered(scm)

scm <- impute_regions(scm)

scm <- dim_red_scMethrix(scm,assay="impute",type="UMAP")
scm <- dim_red_scMethrix(scm,assay="impute",type="PCA")
scm <- dim_red_scMethrix(scm,assay="impute",type="tSNE")

plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cell")
plot_dim_red(scm,dim_red = "PCA",color_anno = "Cell")
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Cell")




scm <- convert_HDF5_scMethrix(scm)
beep()

scm <- mask_by_variance(scm,low_threshold = 0.05)
saveRDS(scm,"D:/Git/thesis_data/GSE151506/bin/scm_450k_hg38proms.rds")
scm <- readRDS("D:/Git/thesis_data/GSE151506/bin/scm_450k_hg38proms.rds")















scm <- bin_scMethrix(scm,regions = hg38_proms)

scm <- mask_by_sample(scm,low_threshold = 20)
scm <- remove_uncovered(scm)

scm <- impute_regions(scm)

scm <- dim_red_scMethrix(scm,assay="impute",type="UMAP")
scm <- dim_red_scMethrix(scm,assay="impute",type="PCA")
scm <- dim_red_scMethrix(scm,assay="impute",type="tSNE")

plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cell")
plot_dim_red(scm,dim_red = "PCA",color_anno = "Cell")
plot_dim_red(scm,dim_red = "tSNE",color_anno = "Cell")


scm2 <- readRDS("D:/Git/thesis_data/GSE151506/bin/scm_bin_prom.rds")
proms <- trim(promoters(genes(txdb), upstream = 1000, downstream = 1000))
scm <- bin_scMethrix(scm,regions = proms)

scm <- mask_by_variance(scm,low_threshold = 0.1)
scm <- mask_by_sample(scm,low_threshold = 20)
scm <- remove_uncovered(scm)
scm <- get_coldata_stats(scm)
scm <- get_rowdata_stats(scm)


scm <- impute_regions(scm)
scm <- saveRDS(scm, "D:/Git/thesis_data/GSE151506/bin/scm_imp.rds")
scm <- readRDS("D:/Git/thesis_data/GSE151506/bin/scm_imp.rds")

scm <- reduce_scMethrix(scm)

scm <- dim_red_scMethrix(scm,assay="impute",type="PCA")
plot_dim_red(scm,dim_red = "UMAP",color_anno = "Cell")




reducedDim(scm,"UMAP")





require(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
proms <- promoters(genes(txdb), upstream = 1000, downstream = 1000)




####### Taken calculate_windows_around_450K_array_CpGs.Rmd

resources.dir = "D:/Git/monobrainDNAme/"

library(readr)
library(readxl)
library(dplyr)
#library(randomForest)
library(doMC)
#library(e1071)

# Control random number generation
set.seed(210) # set a seed to RNG

# register number of cores to be used for parallel evaluation
registerDoMC(cores = parallel::detectCores())

# DNA methylation matrix
file <- paste0(resources.dir, "GSE151506/LGG.GBM.meth.txt") # This file is provided by the study "Molecular profiling
# file <- paste0(resources.dir, "DNAmetBulkSignatures/LGG.GBM.meth.txt") # This file is provided by the study "Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma", Ceccarelli et al., accessible via TCGA database
LGG.GBM <- as.data.frame(readr::read_tsv(file))
rownames(LGG.GBM) <- LGG.GBM$Composite.Element.REF
idx <- grep("TCGA",colnames(LGG.GBM))
colnames(LGG.GBM)[idx] <- substr(colnames(LGG.GBM)[idx], 1, 12) # reduce complete barcode to sample identifier (first 12 characters) 

# metadata with samples molecular subtypes
library(DT)
file <- paste0(resources.dir, "GSE151506/mmc2.xlsx") # This file is provided by the study "Molecular profiling reveals
# file <- paste0(resources.dir, "DNAmetBulkSignatures/mmc2.xlsx") # This file is provided by the study "Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma"
metadata <-  readxl::read_excel(file, sheet = "S1A. TCGA discovery dataset", skip = 1)
DT::datatable(metadata[,c("Case",
                          "Pan-Glioma DNA Methylation Cluster",
                          "Supervised DNA Methylation Cluster",
                          "IDH-specific DNA Methylation Cluster")])

#Probes metadata information are downloaded from http://zwdzwd.io/InfiniumAnnotation This will be used to remove probes that should be masked from the training.
# file <- "http://zwdzwd.io/InfiniumAnnotation/20170313/EPIC/EPIC.manifest.hg38.rda"
# if(!file.exists(basename(file))) downloader::download(file,basename(file))
# load(basename(file))
load(paste0(resources.dir,"GSE151506/EPIC.manifest.hg38.rda"))
library(ChAMPdata)
data(EPIC.manifest.hg38)

# load signatures
signatures.file <- paste0(resources.dir, "GSE151506/PanGlioma_MethylationSignatures.xlsx") # This file is provided by the
# signatures.file <- paste0(resources.dir, "DNAmetBulkSignatures/PanGlioma_MethylationSignatures.xlsx") # This file is provided by the study "Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma"
trainingcol <- 'IDH status'


# prepare training data
sheet <- "1,300 pan-glioma tumor specific"
trainingset <- grep("Mutant|WT|NA",unique(metadata$`IDH status`),value = T)
#trainingset <- grep("LGm1|LGm2|LGm3|LGm4|LGm5|LGm6",unique(metadata$'Pan-Glioma DNA Methylation Cluster'),value = T)

# The DNA methylation matrix will be subset to the DNA methylation signatures and samples with classification.
plat <- "EPIC"
signature.probes <-  read_excel(signatures.file,  sheet = sheet)  %>% pull(1) 
#samples <- dplyr::filter(metadata, 'IDH-specific DNA Methylation Cluster' %in% trainingset)
samples <- as.data.frame(metadata[metadata$`IDH status` %in% trainingset,])
RFtrain <- LGG.GBM[signature.probes, colnames(LGG.GBM) %in% as.character(samples$Case)] %>% na.omit 

# Filter signatures with EPIC?
# consider.probs <- intersect(signature.probes, names(EPIC.manifest.hg38))
# signature.probes <- consider.probs[!EPIC.manifest.hg38[consider.probs]$MASK.general]
# RFtrain <- RFtrain[signature.probes,]

#merge the samples with their classification. In the end, we will have samples in the row, and probes and classification as columns.
trainingdata <- t(RFtrain)
trainingdata <- merge(trainingdata, metadata[,c("Case", trainingcol)], by.x=0,by.y="Case", all.x=T)
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata$Row.names <- NULL

library(gplots)
pdf(paste0(resources.dir,"GSE151506/TCGA_1300_heatmap.pdf"))
heatmap.2(as.matrix(RFtrain), trace="none") # Shortcut to final result
dev.off()

LGG.GBM$loci <- paste0(LGG.GBM$Chromosome,":",LGG.GBM$Genomic_Coordinate)
bulk.loci.signature <- unique(LGG.GBM$loci)

######################################################################################3

library(GenomicRanges)
library(gUtils)

WINDOW = 500 # This is half the window size (left margin and right margin, the full window in this case will be 1kb)
TCGA_HG38_loci = read.table(paste0(resources.dir, "GSE151506/TCGA_HG38.csv"), sep = ",", header = 1)
#patients_names <- colnames(TCGA_HG38_loci)[grep("TCGA",colnames(TCGA_HG38_loci))]
TCGA_HG38_loci$loci <- paste0(TCGA_HG38_loci$Chromosome,":",TCGA_HG38_loci$Genomic_Coordinate)

LGm_probes = read.table(paste0(resources.dir, "GSE151506/1300_pan_glioma_tumor_specific.csv"), sep = ",", header = 1)

# In case you want just the 1300 CpGs use this
TCGA_HG38_loci = TCGA_HG38_loci[TCGA_HG38_loci$Composite.Element.REF %in% LGm_probes$probeID,]

# In case you want all the CpGs you will need alot of memory or to filter based on another strategy Here filtering based on variability across samples:
# rowSTD <- apply(dplyr::select(TCGA_HG38_loci,contains("TCGA")),1,sd,na.rm = TRUE)
# hist(rowSTD)
# keepCpG <- (rowSTD > 0.05)
# keepCpG[is.na(keepCpG)] <- FALSE
# TCGA_HG38_loci <- TCGA_HG38_loci[keepCpG,]

TCGA_HG38_loci$left_coordinate_window <- unlist(lapply(TCGA_HG38_loci$Genomic_Coordinate - WINDOW, function(x) {max(x, 0)}))
TCGA_HG38_loci$right_coordinate_window <- unlist(lapply(TCGA_HG38_loci$Genomic_Coordinate + WINDOW, function(x) {max(x, 0)}))

# Check that all windows managed to be WINDOW size
window_size <- unlist(TCGA_HG38_loci$right_coordinate_window) - unlist(TCGA_HG38_loci$left_coordinate_window)

TCGA_HG38_loci.raw <- TCGA_HG38_loci

# Check if windows are overlapping:

# Create ranges
TCGA_HG38_loci_for_mapping <- as.data.frame((TCGA_HG38_loci))
TCGA_HG38_loci_for_mapping <- TCGA_HG38_loci_for_mapping[!is.na(TCGA_HG38_loci_for_mapping[,'left_coordinate_window']),]
TCGA_HG38_loci_for_mapping <- TCGA_HG38_loci_for_mapping[order(TCGA_HG38_loci_for_mapping[,'Chromosome'],TCGA_HG38_loci_for_mapping[,'left_coordinate_window']),]
rownames(TCGA_HG38_loci_for_mapping) <- 1:nrow(TCGA_HG38_loci_for_mapping)
TCGA_HG38_loci_for_mapping <- droplevels(TCGA_HG38_loci_for_mapping)

TCGA_HG38_loci_for_mapping$center = ((TCGA_HG38_loci_for_mapping$left_coordinate_window + TCGA_HG38_loci_for_mapping$right_coordinate_window)/2)

left_coordinate_window <- GRanges((TCGA_HG38_loci_for_mapping$Chromosome), IRanges(start=TCGA_HG38_loci_for_mapping$center, end=TCGA_HG38_loci_for_mapping$center))
right_coordinate_window <- GRanges((TCGA_HG38_loci_for_mapping$Chromosome), IRanges(start=TCGA_HG38_loci_for_mapping$center, end=TCGA_HG38_loci_for_mapping$center))


# calculate overlaps
library(gUtils)
df.distances<-as.data.frame(gr.dist(left_coordinate_window,right_coordinate_window))
df.distances[df.distances==0]<-10000000
all_overlaps <- which(df.distances < (2*WINDOW), arr.ind = TRUE)
all_overlaps<-as.data.frame(all_overlaps)

# What is the maximun distnace?
all_overlaps$gc1 <-  TCGA_HG38_loci_for_mapping[all_overlaps[,1],"Genomic_Coordinate"]
all_overlaps$gc2 <-  TCGA_HG38_loci_for_mapping[all_overlaps[,2],"Genomic_Coordinate"]
all_overlaps$gc1_Chromosome <-  TCGA_HG38_loci_for_mapping[all_overlaps[,1],"Chromosome"]
all_overlaps$gc2_Chromosome <-  TCGA_HG38_loci_for_mapping[all_overlaps[,2],"Chromosome"]
all_overlaps$gc1_left <-  TCGA_HG38_loci_for_mapping[all_overlaps[,1],"left_coordinate_window"]
all_overlaps$gc1_right <-  TCGA_HG38_loci_for_mapping[all_overlaps[,1],"right_coordinate_window"]
all_overlaps$gc2_left <-  TCGA_HG38_loci_for_mapping[all_overlaps[,2],"left_coordinate_window"]
all_overlaps$gc2_right <-  TCGA_HG38_loci_for_mapping[all_overlaps[,2],"right_coordinate_window"]

isProblem = function(x){
  if (x["gc1"] < x["gc2"]){
    is.it <- (x["gc2_left"] < x["gc1_right"])
  } else {
    is.it <- (x["gc1_left"] < x["gc2_right"])
  }
  return(is.it)
}
all_overlaps$problem <-apply(all_overlaps, 1, function(x) isProblem(x))
all_overlaps$distance <- (all_overlaps$gc1 - all_overlaps$gc2) 
all_overlaps <- all_overlaps[all_overlaps$problem,]
all_overlaps <- all_overlaps[all_overlaps$distance > 0,]

# Create a graph from the problematic CpGs (these that are too close to one another)
GCs <- unique(c(paste0(all_overlaps[,c("gc1_Chromosome")],":",all_overlaps[,c("gc1")]),paste0(all_overlaps[,c("gc2_Chromosome")],":",all_overlaps[,c("gc2")])))
number_gc <- c(1:length(GCs))
names(number_gc) <- GCs

a<-lapply(paste0(all_overlaps[,c("gc1_Chromosome")],":",all_overlaps[,c("gc1")]), function(x) number_gc[as.character(x)])
b<-lapply(paste0(all_overlaps[,c("gc2_Chromosome")],":",all_overlaps[,c("gc2")]), function(x) number_gc[as.character(x)])

library(igraph)
all_overlaps$gc1_num <- a
all_overlaps$gc2_num <- b
g <- graph_from_data_frame((all_overlaps[,c("gc1_num","gc2_num")]), directed = FALSE)
g <- set.vertex.attribute(g, "name", value=GCs)
g <- set.edge.attribute(g, "weight", value=all_overlaps[,c("distance")])

# find the CpG and update its left and right
setNewTile = function(Chromosome, Genomic_Coordinate, loci, left, right) {

  row.num <- match(TRUE, TCGA_HG38_loci$loci == loci)
  a = TCGA_HG38_loci[row.num,]
  #test
  if ((a$Chromosome != Chromosome) | (a$Genomic_Coordinate != Genomic_Coordinate)) {
    print ("Something is wrong in setNewTile")
  }
  else {
    rowname <- rownames(TCGA_HG38_loci[TCGA_HG38_loci$loci == loci,])
    TCGA_HG38_loci[rowname, "left_coordinate_window"] = left
    TCGA_HG38_loci[rowname,"right_coordinate_window"] = right
    print(TCGA_HG38_loci[TCGA_HG38_loci$loci == loci,"right_coordinate_window"])
  }

  assign('TCGA_HG38_loci',TCGA_HG38_loci,envir=.GlobalEnv)
}


findNewCenteredTiles = function(x) {
  print (x)
  print(V(g)$name[cl$membership %in% x])
  cluster <- c(V(g)$name[cl$membership %in% x])
  
  # Get the genomic coordinate
  cluster_Genomic_Coordinates <- unlist(lapply(cluster,function (x) {as.numeric(strsplit(x, ":")[[1]][2])}))
  space <- max(cluster_Genomic_Coordinates) - min(cluster_Genomic_Coordinates)
  
  # If the cluster of close CpGs is less than the WINDOW then find the new center
  if (space < (WINDOW)*2) {
    new_center <- ((max(cluster_Genomic_Coordinates) + min(cluster_Genomic_Coordinates))/2)
    left <- (new_center - WINDOW)
    right <- (new_center + WINDOW)
    
    for (CpG in cluster)  {
      chromosome = (strsplit(CpG,":")[[1]][1])
      genomic.coordinate = (strsplit(CpG,":")[[1]][2])
      setNewTile(chromosome, genomic.coordinate, CpG, floor(as.numeric(left)), floor(as.numeric(right)))
    }
    
  }
  # If the cluster is more than the WINDOW then tile the space to minimum number of non overlapping tiles
  else {

    # extra space needed form each side 
    extra_space = (2*WINDOW - (space - WINDOW * floor(space/WINDOW)))-1
    
    # The whole new space left coordinate and right coordinate
    left <- (min(cluster_Genomic_Coordinates) - (extra_space/2))  
    right <- (max(cluster_Genomic_Coordinates) + (extra_space/2)) 
    
    # Create a range from the whole new space begining to end
    new.range <-  IRanges(left,right)
    
    # Create a range from the whole new space beginnign to end
    cluster_Chromosome <- unlist(lapply(cluster,function (x) {(strsplit(x, ":")[[1]][1])}))
    
    # the min and max should be he same chromosome
    #new.range <- GRanges(c(min(cluster_Chromosome)),c(left,right))
    
    # Break it to tiles
    tiles <- tile(new.range, width = 2*WINDOW)
    target_tiles <- GRanges(min(cluster_Chromosome) ,tiles[[1]])
    
    # Match each point to a tile and set its new window
    d <- GRanges(min(cluster_Chromosome), IRanges(cluster_Genomic_Coordinates, width=1))
    OL <- findOverlaps(d, target_tiles)
    hits <- target_tiles[as.data.frame(OL)[,2],]
    df.d <- as.data.frame(d)
    df.hits <- as.data.frame(hits)
    for (hit in as.data.frame(OL)$queryHits) {
      setNewTile(df.d[hit,]$seqnames, df.d[hit,]$start, cluster[hit],df.hits[1,]$start, df.hits[1,]$end)
    }
    
  }
}

cl <- clusters(g)
lapply(seq_along(cl$csize)[cl$csize > 1], function(x) findNewCenteredTiles(x) )

# If you run overlaps now it should come out empty

# Create ranges
TCGA_HG38_loci_for_mapping <- as.data.frame((TCGA_HG38_loci))
TCGA_HG38_loci_for_mapping <- TCGA_HG38_loci_for_mapping[order(TCGA_HG38_loci_for_mapping[,'Chromosome'],TCGA_HG38_loci_for_mapping[,'left_coordinate_window']),]
rownames(TCGA_HG38_loci_for_mapping) <- 1:nrow(TCGA_HG38_loci_for_mapping)
TCGA_HG38_loci_for_mapping <- droplevels(TCGA_HG38_loci_for_mapping)

TCGA_HG38_loci_for_mapping$center = ((TCGA_HG38_loci_for_mapping$left_coordinate_window + TCGA_HG38_loci_for_mapping$right_coordinate_window)/2)

left_coordinate_window <- GRanges((TCGA_HG38_loci_for_mapping$Chromosome), IRanges(start=TCGA_HG38_loci_for_mapping$center, end=TCGA_HG38_loci_for_mapping$center))
right_coordinate_window <- GRanges((TCGA_HG38_loci_for_mapping$Chromosome), IRanges(start=TCGA_HG38_loci_for_mapping$center, end=TCGA_HG38_loci_for_mapping$center))


# calculate overlaps
distances<-gr.dist(left_coordinate_window,right_coordinate_window)
df.distances <- as.data.frame(distances)
df.distances[df.distances==0]<-10000000
all_overlaps <- which(df.distances < (2*WINDOW), arr.ind = TRUE)


TCGA_HG38_loci$rrbs_group <- TCGA_HG38_loci$left_coordinate_window
TCGA_HG38_loci_for_mapping <- as.data.frame((TCGA_HG38_loci))
TCGA_HG38_loci_for_mapping <- TCGA_HG38_loci_for_mapping[order(TCGA_HG38_loci_for_mapping[,'Chromosome'],TCGA_HG38_loci_for_mapping[,'left_coordinate_window']),]
rownames(TCGA_HG38_loci_for_mapping) <- 1:nrow(TCGA_HG38_loci_for_mapping)
TCGA_HG38_loci_for_mapping <- droplevels(TCGA_HG38_loci_for_mapping)
target_range <- GRanges((TCGA_HG38_loci_for_mapping$Chromosome), IRanges(start=TCGA_HG38_loci_for_mapping$left_coordinate_window, end=TCGA_HG38_loci_for_mapping$right_coordinate_window),ID=TCGA_HG38_loci_for_mapping$Composite.Element.REF)

TCGA_HG38_GR <- TCGA_HG38_loci_for_mapping[c("Chromosome","left_coordinate_window","right_coordinate_window", "Composite.Element.REF")]
colnames(TCGA_HG38_GR) <- c("chr","start","end","ID")

TCGA_HG38_GR <- makeGRangesFromDataFrame(TCGA_HG38_GR, keep.extra.columns=T)


# Save bulk windows
#write.table(TCGA_HG38_loci, paste0(resources.dir,"DNAmetBulkSignatures/all_varried_cpgs/bulk_centered_Tiles",2*WINDOW,"_HG38.csv"), quote = FALSE, sep = "\t")

#adjust sc data to windows:

# Calc 1kb window around the bulk CpGs
mapSCToRRBSTiles = function(file) {
  
  dataset <- read.table(paste0(coverage.dir,file), header=FALSE, sep="\t")
  dataset$filename <- basename(file)
  
  # Create loci coulmn chromosome:location
  dataset$V2 <- sub("^[0]+", "", dataset$V2) 
  dataset$V1 <- sub("chr", "", dataset$V1) 
  
  colnames(dataset) <- c("Chromosome","Genomic_Coordinate","Genomic_coordinate_end","percentage","methylated","unmethylated","filename")
  
  # Filter to loci of signatre
  dataset$loci <- paste0(dataset$V1,":",dataset$V2)
  
  dataset <- dataset[order(dataset[,'Chromosome'],dataset[,'Genomic_Coordinate']),]
  rownames(dataset) <- 1:nrow(dataset)
  d <- GRanges((dataset$Chromosome), IRanges(as.numeric(dataset$Genomic_Coordinate), width=1))
  
  OL <- findOverlaps(d, target_range)
  hits <- TCGA_HG38_loci_for_mapping[as.data.frame(OL)[,2],]
  hits_in_dataset <- dataset[as.data.frame(OL)[,1],]
  hits_in_dataset$rrbs_group <- hits$rrbs_group
  return(hits_in_dataset)
}


f = function(x, chromosome_field = 'Chromosome', genomic_Coordinate_field = 'Genomic_Coordinate') {
  
  chromosome <- (as.integer(x[chromosome_field]))
  Genomic_Coordinate <- as.integer(x[genomic_Coordinate_field])
  
  c <- b[(b$chr == chromosome)]
  rrbs.group.id <- c[(c$start <= Genomic_Coordinate) & (c$end >= Genomic_Coordinate)]$withinGroupID.V1
  return (rrbs.group.id)
}

#load sc
coverage.dir<-paste0(scRRBS.dir,"/COV/")
file.list <- list.files(path = coverage.dir, recursive = TRUE, pattern = "cov$")

rm(dataset)
i<-0
for (file in file.list){
  
  print(file)
  print (i)
  i <- i +1
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- mapSCToRRBSTiles(file)
    
    
    
  } else if (exists("dataset")){
    
    temp_dataset <-  mapSCToRRBSTiles(file)
    dataset <- rbind(dataset, temp_dataset)
  }
}


# Save sc windows
write.table(dataset, paste0(results.dir,"/all_var_sc_tiles_L",2*WINDOW,"_mean.csv"), quote = FALSE, sep = "\t")
```



