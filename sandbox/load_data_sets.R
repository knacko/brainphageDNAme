home_dir <- "D:/Git/thesis_data/"
cpg_dir <- "D:/Git/sampleData/ref_cpgs/"
chain_dir <- "D:/Git/thesis_data/chains/"

cell_types <- c("Neutrophil","NKcell","Bcell","CD4Tcell","CD8Tcell","Monocyte","WholeBlood","Endothelial","Immune")

#------------------------------------------------------------------------------------------------------------
# GEO: 
# Types: 
# Paper: 
# GEO: 
# Citation: 
# Genome: hg19
# Platform: 

GEO = ""
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData

# Get data

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))

#------------------------------------------------------------------------------------------------------------
# GEO: GSE121483
# Types: Microglia-like macrophages
# Paper: https://pubmed.ncbi.nlm.nih.gov/30451869/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121483
# Citation: Lund H, Pieber M, Parsa R, Han J et al. Competitive repopulation of an empty microglial niche yields functionally distinct subsets of microglia-like cells. Nat Commun 2018 Nov 19;9(1):4845. PMID: 30451869
# Genome: hg19
# Platform: Illumina Epic

GEO = "GSE121483"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))

# GSE121483 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
#           chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = length(files))
# 
# GSE121483 <- scMethrix::liftover_CpGs(scm = GSE121483, chain = chain_file,target_genome="hg38")
# 
# save_HDF5_scMethrix(GSE121483, h5_dir=exp_dir,replace=TRUE)
# GSE121483 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE35069
# Types: CD4+ T cells, CD8+ T cells, CD56+ NK cells, CD19+ B cells, CD14+ monocytes, neutrophils, and eosinophils
# Paper: https://pubmed.ncbi.nlm.nih.gov/22848472/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069
# Citation: Reinius LE, Acevedo N, Joerink M, Pershagen G et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PLoS One 2012;7(7):e41361. PMID: 22848472
# Genome: hg19
# Platform: Illumina 450k

# Setup dirs
GEO = "GSE35069"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData



# Get data
soft.to.bed(in_file,bed_dir,"(.*)")#_.*_.*")
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(setequal(get_sample_name(files),colData$Sample))
expect_true(all(colData$Cell %in% cell_types))

# 
# 
# 
# GSE35069 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
#                                  chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 30)
# 
# GSE35069 <- scMethrix::liftover_CpGs(scm = GSE35069, chain = chain_file,target_genome="hg38")
# 
# save_HDF5_scMethrix(GSE35069, h5_dir=exp_dir,replace=TRUE)
# GSE35069 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)


#------------------------------------------------------------------------------------------------------------
# GEO: GSE88824
# Types: neutrophils, CD4+ T cells, CD8+ T cells, NK cells, B cells and monocytes
# Paper: https://pubmed.ncbi.nlm.nih.gov/30571772/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88824
# Citation: Kennedy DW, White NM, Benton MC, Fox A et al. Critical evaluation of linear regression models for cell-subtype specific methylation signal from mixed blood cell DNA. PLoS One 2018;13(12):e0208915. PMID: 30571772
# Genome: hg19
# Platform: Illumina 450k

# Setup dirs
GEO = "GSE88824"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)
cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)
colData <- data.table(Sample = names(soft@gsms), Cell = cell)
remove_idx <- str_detect(colData$Cell,"^Case.*|.*WBC$")
remove_id <- colData$Sample[remove_idx]
colData <- colData[!remove_idx,]
colData$Cell <- str_remove(colData$Cell,"Control-")
colData$Cell <- str_replace(colData$Cell,"CD19B","Bcell")
colData$Cell <- str_replace(colData$Cell,"CD4T","CD4Tcell")
colData$Cell <- str_replace(colData$Cell,"CD8T","CD8Tcell")
data.table::fwrite(colData, file=colData_file, quote=FALSE, sep='\t', row.names = FALSE)

# Get data
idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")
files = list.files(bed_dir, full.names=TRUE)
file.remove(files[get_sample_name(files) %in% remove_id])
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
files = list.files(bed_dir, full.names=TRUE)
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(setequal(get_sample_name(files),colData$Sample))
expect_true(all(colData$Cell %in% cell_types))




# 
# 
# 
# 
# GSE88824 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
#                                   chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 20)
# 
# GSE88824 <- scMethrix::liftover_CpGs(scm = GSE88824, chain = chain_file,target_genome="hg38")
# 
# save_HDF5_scMethrix(GSE88824, h5_dir=exp_dir,replace=TRUE)
# GSE88824 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE166844
# Types: CD4 T-cells, CD8 T-cells, Monocytes, Granulocyte, B-cells
# Paper: https://pubmed.ncbi.nlm.nih.gov/33739972/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166844
# Citation: Hannon E, Mansell G, Walker E, Nabais MF et al. Assessing the co-variability of DNA methylation across peripheral cells and tissues: Implications for the interpretation of findings in epigenetic epidemiology. PLoS Genet 2021 Mar;17(3):e1009443. PMID: 33739972
# Genome: hg19
# Platform: Illumina EPIC

base_dir = paste0(home_dir,"GSE166844/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE166844_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get data and colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)

ids <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$description)
cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$source_name_ch1)

colData <- data.table(Sample = names(soft@gsms), ID = ids, Cell = cell)
data.table::fwrite(colData, file=colData_file, quote=FALSE, sep='\t', row.names = FALSE)

header <- fread(file=paste0(raw_dir,"GSE166844_Variance_processed_signals.csv"),nrows=1,header = F)
cols <- which(apply(header, 2, function(x) !any(grepl("Pval", x))))

meths <- as.data.frame(fread(file = paste0(raw_dir,"GSE166844_Variance_processed_signals.csv"),select=as.integer(cols)))
rownames(meths) <- meths[,1]
meths$V1 <- NULL

colMatch <- match(colnames(meths),colData$ID)

stopifnot(identical(colnames(meths)[1],colData$ID[colMatch[1]]))

colnames(meths) <- colData$Sample[colMatch]
colData$ID <- NULL

remove_id <- which(colData$Cell %in% c("Buccal","Nasal"))
colData <- colData[-remove_id,]
meths <- meths[ , order(names(meths))]
meths <- meths[,-remove_id] 

RSet = RatioSet(Beta = meths)
annotation(RSet) = annotation(minfiDataEPIC::RGsetEPIC)

ratioset.to.bed(RSet, out_dir = bed_dir, regex = "(.*)")
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
files = list.files(bed_dir, full.names=TRUE)
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(setequal(get_sample_name(files),colData$Sample))
expect_true(all(colData$Cell %in% cell_types))

# files = list.files(bed_dir, full.names=TRUE)
# GSE166844 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
#                                   chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 30)
# 
# GSE166844 <- scMethrix::liftover_CpGs(scm = GSE166844, chain = chain_file,target_genome="hg38")
# 
# save_HDF5_scMethrix(GSE166844, h5_dir=exp_dir,replace=TRUE)
# GSE166844 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE110554
# Types: neutrophils, monocytes, B-lymphocytes, natural killer (NK) cells, CD4+ T-cells, and CD8+ T-cells
# Paper: https://pubmed.ncbi.nlm.nih.gov/29843789/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554
# Citation: Salas LA, Koestler DC, Butler RA, Hansen HM et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 2018 May 29;19(1):64. PMID: 29843789
# Genome: hg19
# Platform: Illumina EPIC

GEO = "GSE110554"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)

colData <- sapply(names(soft@gsms), function(gsm) {
  gsm <- soft@gsms[[gsm]]@header$characteristics_ch1
  gsm <- gsm[which(gsm %like% "cell type:")]
  gsm <- gsub(".*: (.*)","\\1",gsm)
  return(gsm)
})

colData <- data.frame(Sample = names(colData), Cell = colData,row.names = NULL)
remove_mix <- colData$Sample[which(colData$Cell == "MIX")]
colData <- colData[-(which(colData$Cell == "MIX")),]
data.table::fwrite(colData, file=colData_file, quote=FALSE, sep='\t', row.names = FALSE) 
###### ^ This needs manual editting for cell names

# Get data
supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
supp_file <- rownames(supp_file)
supp_files <- untar(tarfile = supp_file, list=TRUE)
supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
supp_files <- untar(tarfile = supp_file, exdir = raw_dir, files = supp_files)
file.remove(supp_file)
sapply(supp_files, function(file) GEOquery::gunzip(file, remove=TRUE))
mix_files <- list.files(raw_dir, full.names=TRUE)
file.remove(mix_files[which(rowSums(sapply(remove_mix, like, vector = mix_files)) == 1)])
idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))

#------------------------------------------------------------------------------------------------------------
# GEO: GSE96612
# Types: Neuron, glia
# Paper: https://pubmed.ncbi.nlm.nih.gov/30643296/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96612
# Citation: Rizzardi LF, Hickey PF, Rodriguez DiBlasi V, Tryggvadóttir R et al. Neuronal brain-region-specific DNA methylation and chromatin accessibility are associated with neuropsychiatric trait heritability. Nat Neurosci 2019 Feb;22(2):307-316. PMID: 30643296
# Genome: hg19
# Platform: Bulk WGBS

GEO = "GSE96612"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

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

#------------------------------------------------------------------------------------------------------------
# GEO: GSE151506
# Types: tumor and immune
# Paper: 
# GEO: 
# Citation: 
# Genome: hg19
# Platform: 

GEO = "GSE151506"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)
cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$type)
colData <- data.table(Sample = names(soft@gsms), Cell = cell, ID = cell)



# Get data

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))




#------------------------------------------------------------------------------------------------------------
# GEO: GSE144804
# Types: Endothelial
# Paper: https://pubmed.ncbi.nlm.nih.gov/32231389/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144804
# Citation: Rhead B, Shao X, Quach H, Ghai P et al. Global expression and CpG methylation analysis of primary endothelial cells before and after TNFa stimulation reveals gene modules enriched in inflammatory and infectious diseases and associated DMRs. PLoS One 2020;15(3):e0230884. PMID: 32231389
# Genome: hg19
# Platform: Illumina EPIC

GEO = "GSE144804"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)
cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
colData <- data.table(Sample = names(soft@gsms), Cell = cell)
remove_idx <- str_detect(colData$Cell,".*_TNF$")
remove_gsm <- colData$Sample[remove_idx]
colData <- colData[!remove_idx,]

# Get data
supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
supp_file <- rownames(supp_file)
supp_files <- untar(tarfile = supp_file, list=TRUE)
supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
supp_files <- untar(tarfile = supp_file, exdir = raw_dir, files = supp_files)
file.remove(supp_file)
sapply(supp_files, function(file) GEOquery::gunzip(file, remove=TRUE))
files <- list.files(raw_dir, full.names=TRUE)
file.remove(files[which(rowSums(sapply(remove_gsm, like, vector = files)) == 1)])
idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))

#------------------------------------------------------------------------------------------------------------
# GEO: GSE98203
# Types: Neuron
# Paper: https://pubmed.ncbi.nlm.nih.gov/28556790/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98203
# Citation: Kozlenkov A, Jaffe AE, Timashpolsky A, Apontes P et al. DNA Methylation Profiling of Human Prefrontal Cortex Neurons in Heroin Users Shows Significant Difference between Genomic Contexts of Hyper- and Hypomethylation and a Younger Epigenetic Age. Genes (Basel) 2017 May 30;8(6). PMID: 28556790
# Genome: hg19
# Platform: Illumina 450k

GEO = "GSE98203"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
mkdirs(base_dir,raw_dir,bed_dir,exp_dir)
colData_file = paste0(base_dir,"colData/",GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

# Get colData
soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)
cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)
colData <- data.table(Sample = names(soft@gsms), Cell = cell)
remove_idx <- str_detect(colData$Cell,".*_TNF$")
remove_gsm <- colData$Sample[remove_idx]
colData <- colData[!remove_idx,]

# Get data
supp_file <- GEOquery::getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = substr(raw_dir,1,nchar(raw_dir)-1), filter_regex = ".*RAW.tar")
supp_file <- rownames(supp_file)
supp_files <- untar(tarfile = supp_file, list=TRUE)
supp_files <- supp_files[grepl(".*idat.gz$", supp_files,ignore.case = TRUE)]
supp_files <- untar(tarfile = supp_file, exdir = raw_dir, files = supp_files)
file.remove(supp_file)
sapply(supp_files, function(file) GEOquery::gunzip(file, remove=TRUE))
files <- list.files(raw_dir, full.names=TRUE)
file.remove(files[which(rowSums(sapply(remove_gsm, like, vector = files)) == 1)])
idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")
files = list.files(bed_dir, full.names=TRUE)
liftover_beds(files = files,chain = chain_file)

# Compare colData and data
expect_true(setequal(get_sample_name(files),colData$Sample))
colData = read.table(file = colData_file, sep = '\t', header = TRUE)
expect_true(all(colData$Cell %in% cell_types))


#------------------------------------------------------------------------------------------------------------
# GEO: GSE151506
# Types: 
# Paper: 
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151506
# Citation: 
# Genome: hg19
# Platform: 

GEO = "GSE151506"
base_dir = paste0(home_dir,GEO,"/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_dir = paste0(base_dir,"colData/")
mkdirs(c(base_dir,raw_dir,bed_dir,exp_dir,colData_dir))
colData_file = paste0(colData_dir,GEO,"_colData.tsv")
chain_file = paste0(chain_dir,"hg19ToHg38.over.chain")

soft <- GEOquery::getGEOfile(GEO)
soft <- GEOquery::getGEO(filename=soft)

cell <- sapply(names(soft@gsms), function(gsm) soft@gsms[[gsm]]@header$title)





gsm2file <- data.table::fread("D:/Git/thesis_data/GSE151506/colData/GSE151506_SRRs_and_original_rawfile_names.txt",
                                blank.lines.skip = TRUE, sep = "|", header = FALSE)
row_idx <-  which(apply(gsm2file, 1, function(x) str_detect(x,"_r1:")))
gsm <- unlist(gsm2file[row_idx,])
gsm <- str_remove(gsm,"SRR.*; ")
gsm <- str_remove(gsm,"_r1:")
file.name <- unlist(gsm2file[row_idx+1,])
file.name <- str_remove(file.name,".fastq.*")
gsm2file = data.table(GSM = gsm,file = file.name)    




