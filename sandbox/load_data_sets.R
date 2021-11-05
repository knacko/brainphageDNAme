home_dir <- "D:/Git/thesis_data/"
cpg_dir <- "D:/Git/sampleData/ref_cpgs/"
chain_dir <- "D:/Git/thesis_data/chains/"

#------------------------------------------------------------------------------------------------------------
# GEO: GSE121483
# Types: Microglia-like macrophages
# Paper: https://pubmed.ncbi.nlm.nih.gov/30451869/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121483
# Citation: Lund H, Pieber M, Parsa R, Han J et al. Competitive repopulation of an empty microglial niche yields functionally distinct subsets of microglia-like cells. Nat Commun 2018 Nov 19;9(1):4845. PMID: 30451869
# Genome: hg19
# Platform: Illumina Epic

base_dir = paste0(home_dir,"GSE121483/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE121483_colData.tsv")

idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = cpg_dir, genome = "BSgenome.Hsapiens.UCSC.hg19")
colData = read.table(file = colData_file, sep = '\t', header = TRUE)

GSE121483 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
          chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = length(files))

GSE121483 <- scMethrix::liftover_CpGs(scm = GSE121483, chain = "D:/Git/thesis_data/chains/hg19ToHg38.over.chain",target_genome="hg38")

save_HDF5_scMethrix(GSE121483, h5_dir=exp_dir,replace=TRUE)
GSE121483 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE35069
# Types: CD4+ T cells, CD8+ T cells, CD56+ NK cells, CD19+ B cells, CD14+ monocytes, neutrophils, and eosinophils
# Paper: https://pubmed.ncbi.nlm.nih.gov/22848472/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069
# Citation: Reinius LE, Acevedo N, Joerink M, Pershagen G et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PLoS One 2012;7(7):e41361. PMID: 22848472
# Genome: hg19
# Platform: Illumina 450k

base_dir = paste0(home_dir,"GSE35069/")
raw_dir = paste0(base_dir,"raw/")
in_file = list.files(raw_dir, full.names=TRUE)
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE35069_colData.tsv")

soft.to.bed(in_file,bed_dir,"(.*)")#_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = cpg_dir, genome = "BSgenome.Hsapiens.UCSC.hg19")
colData = read.table(file = colData_file, sep = '\t', header = TRUE)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE88824
# Types: neutrophils, CD4+ T cells, CD8+ T cells, NK cellsa, B cells and monocytes
# Paper: https://pubmed.ncbi.nlm.nih.gov/30571772/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88824
# Citation: Kennedy DW, White NM, Benton MC, Fox A et al. Critical evaluation of linear regression models for cell-subtype specific methylation signal from mixed blood cell DNA. PLoS One 2018;13(12):e0208915. PMID: 30571772
# Genome: hg19
# Platform: Illumina 450k

base_dir = paste0(home_dir,"GSE88824/")
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE88824_colData.tsv")

idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = cpg_dir, genome = "BSgenome.Hsapiens.UCSC.hg19")
colData = read.table(file = colData_file, sep = '\t', header = TRUE)

GSE88824 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
                                  chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 20)

GSE88824 <- scMethrix::liftover_CpGs(scm = GSE88824, chain = chain_dir,target_genome="hg38")

save_HDF5_scMethrix(GSE88824, h5_dir=exp_dir,replace=TRUE)
GSE121483 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)

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

ref_cpgs = load_ref_cpgs(dir = cpg_dir, genome = "BSgenome.Hsapiens.UCSC.hg19")

soft <- GEOquery::getGEOfile("GSE166844")
soft <- GEOquery::getGEO(filename=colData)

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
GSE166844 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
                                  chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 30)

GSE166844 <- scMethrix::liftover_CpGs(scm = GSE166844, chain = "D:/Git/thesis_data/chains/hg19ToHg38.over.chain",target_genome="hg38")

save_HDF5_scMethrix(GSE166844, h5_dir=exp_dir,replace=TRUE)
GSE166844 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)
