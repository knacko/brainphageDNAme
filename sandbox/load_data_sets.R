#------------------------------------------------------------------------------------------------------------
# GEO: GSE121483
# Types: Microglia-like macrophages
# Paper: https://pubmed.ncbi.nlm.nih.gov/30451869/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121483
# Citation: Lund H, Pieber M, Parsa R, Han J et al. Competitive repopulation of an empty microglial niche yields functionally distinct subsets of microglia-like cells. Nat Commun 2018 Nov 19;9(1):4845. PMID: 30451869
# Genome: hg19

base_dir = "D:/Git/thesis_data/GSE121483/"
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE121483_colData.tsv")

idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = "D:/Git/sampleData/ref_cpgs/", genome = "BSgenome.Hsapiens.UCSC.hg19")
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

base_dir = "D:/Git/thesis_data/GSE35069/"
raw_dir = paste0(base_dir,"raw/")
in_file = list.files(raw_dir, full.names=TRUE)
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE35069_colData.tsv")

soft.to.bed(in_file,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = "D:/Git/sampleData/ref_cpgs/", genome = "BSgenome.Hsapiens.UCSC.hg19")
colData = read.table(file = colData_file, sep = '\t', header = TRUE)

#------------------------------------------------------------------------------------------------------------
# GEO: GSE88824
# Types: neutrophils, CD4+ T cells, CD8+ T cells, NK cellsa, B cells and monocytes
# Paper: https://pubmed.ncbi.nlm.nih.gov/30571772/
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88824
# Citation: Kennedy DW, White NM, Benton MC, Fox A et al. Critical evaluation of linear regression models for cell-subtype specific methylation signal from mixed blood cell DNA. PLoS One 2018;13(12):e0208915. PMID: 30571772
# Genome: hg19

base_dir = "D:/Git/thesis_data/GSE88824/"
raw_dir = paste0(base_dir,"raw/")
bed_dir = paste0(base_dir,"bed/")
exp_dir = paste0(base_dir,"exp/")
colData_file = paste0(base_dir,"colData/GSE88824_colData.tsv")

idat.to.bed(raw_dir,bed_dir,"(.*)_.*_.*")

files = list.files(bed_dir, full.names=TRUE)
ref_cpgs = load_ref_cpgs(dir = "D:/Git/sampleData/ref_cpgs/", genome = "BSgenome.Hsapiens.UCSC.hg19")
colData = read.table(file = colData_file, sep = '\t', header = TRUE)

GSE88824 <- scMethrix::read_beds(files=files, h5=TRUE, ref_cpgs = ref_cpgs,
                                  chr_idx=1, start_idx=2, end_idx=3, beta_idx=5, colData = colData, n_threads=0, batch_size = 20)

GSE88824 <- scMethrix::liftover_CpGs(scm = GSE88824, chain = "D:/Git/thesis_data/chains/hg19ToHg38.over.chain",target_genome="hg38")

save_HDF5_scMethrix(GSE88824, h5_dir=exp_dir,replace=TRUE)
GSE121483 <- scMethrix::load_HDF5_scMethrix(dir=exp_dir)