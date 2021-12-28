library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Randomly select 1000 CpGs to be significantly differentially methylated
sigcpgs <- sample(rownames(ann),1000,replace=FALSE)

# All CpG sites tested
allcpgs <- rownames(ann)

# GO testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO", plot.bias = TRUE, prior.prob = TRUE)

# Total number of GO categories significant at 5\% FDR
table(gst$FDR<0.05)

# Table of top GO results
topGO(gst)

# GO testing ignoring bias
gst.bias <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "GO", prior.prob=FALSE)

# Total number of GO categories significant at 5\% FDR ignoring bias
table(gst.bias$FDR<0.05)

# Table of top GO results ignoring bias
topGO(gst.bias)

# KEGG testing
kegg <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = "KEGG", prior.prob=TRUE)

# Table of top KEGG results
topKEGG(kegg)












gometh <- function(sig.cpg, all.cpg=NULL, collection=c("GO","KEGG"), 
                   array.type = c("450K","EPIC"), plot.bias=FALSE, 
                   prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE, 
                   fract.counts = TRUE, 
                   genomic.features = c("ALL", "TSS200","TSS1500","Body",
                                        "1stExon","3'UTR","5'UTR","ExonBnd"),
                   sig.genes = FALSE)
  # Gene ontology testing or KEGG pathway analysis for Illumina methylation 
  # arrays based on goseq
  # Takes into account probability of differential methylation based on
  # numbers of probes on array per gene
  # Belinda Phipson
  # 28 January 2015. Last updated 1 September 2020.
  # EPIC functionality contributed by Andrew Y.F. Li Yim
{
  
  browser()
  
  array.type <- match.arg(toupper(array.type), c("450K","EPIC"))    
  collection <- match.arg(toupper(collection), c("GO","KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", "TSS200","TSS1500",
                                                    "Body", "1stExon","3'UTR",
                                                    "5'UTR","ExonBnd"), 
                                several.ok = TRUE)
  
  if(length(genomic.features) > 1 & any(grepl("ALL", genomic.features))){
    message("All input CpGs are used for testing.") 
    genomic.features <- "ALL"
  } 
  
  if(array.type == "450K" & any(grepl("ExonBnd", genomic.features))){
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n
           please remove it from your genomic.feature parameter\n
           specification.") 
  }
  
  if(collection == "GO"){
    go <- .getGO()
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
    rownames(result) <- result$GOID
    
  } else if(collection == "KEGG"){
    kegg <- .getKEGG()
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=kegg$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(kegg$idTable,result,by.x="PathwayID",by.y="row.names")
    rownames(result) <- result$PathwayID
  }
  
  result[,-1]
}  

.getGO <- function(){
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db package required but not installed.")
  egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
  d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
  GeneID.PathID <- GeneID.PathID[d, ]
  GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                      keys=unique(GeneID.PathID$go_id), 
                                                      columns=c("GOID","ONTOLOGY","TERM"), 
                                                      keytype="GOID"))
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
  
  list(idList=go, idTable=GOID.TERM)
}

.getKEGG <- function(){
  GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", convert = TRUE)
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "Hsa", 
                                                remove.qualifier = TRUE)
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by="PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, list)
  
  list(idList = kegg, idTable = PathID.PathName)
}  

.plotBias <- function(D,bias)
  # Plotting function to show gene level CpG density bias
  # Belinda Phipson
  # 5 March 2015
{
  o <- order(bias)
  splitf <- rep(1:100,each=200)[1:length(bias)]
  avgbias <- tapply(bias[o],factor(splitf),mean)
  sumDM <- tapply(D[o],factor(splitf),sum)
  propDM <- sumDM/table(splitf)
  graphics::par(mar=c(5,5,2,2))
  graphics::plot(avgbias,as.vector(propDM),
                 xlab="Number of CpGs per gene (binned)",
                 ylab="Proportion Differential Methylation",cex.lab=1.5,
                 cex.axis=1.2)
  graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}

.estimatePWF <- function(D,bias)
  # An alternative to goseq function nullp, which is transformation invariant
  # Belinda Phipson and Gordon Smyth
  # 6 March 2015
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                               columns=c("ENTREZID","SYMBOL"), 
                                               keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}

# .reduceMultiMap <- function(flat){
#   mm <- table(flat$cpg)
#   mm <- names(mm[mm > 2])
#   red <- tapply(flat$entrezid[flat$cpg %in% mm],
#                 flat$cpg[flat$cpg %in% mm], sample, 1)
#   key <- paste(flat$cpg,flat$entrezid,sep=".")
#   rkey <- paste(names(red),red,sep = ".")
#   w <- which(key %in% rkey)
#   
#   newmm <- flat[w,] 
#   newflat <- flat[-which(flat$cpg %in% mm),]
#   newflat <- rbind(newflat, newmm)
#   newflat
# }

