






library(AnnotationHub)
ah = AnnotationHub()
hg38_genome <- ah[["AH97949"]]
hg38_genes <- genes(hg38_genome)
hg38_tss <- promoters(hg38_genes,upstream=1000,downstream=1000)

library(annotatr)
build_annotations(genome, annotations)
builtin_annotations()

#--- Get promoter regions 
library(AnnotationHub)
ah = AnnotationHub()
hg38_genome <- ah[["AH97949"]]
hg38_genes <- genes(hg38_genome)
hg38_proms = promoters(hg38_genes, upstream=1000, downstream=1000)









