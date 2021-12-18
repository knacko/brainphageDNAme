#---- do_methylcibersort -------------------------------------------------------------------------------------
#' @param scm scMethrix; The scMethrix object
#' @param assay string; The assay to use
#' @param cell_types string; Cells to ignore
#' @param cell_col string; The column in colData(scm) that contains the cell annotation
#' @param ... Arguments passed to feature.select.new
#'
#' @return scMethrix object containing only the selected features
#' @export
#'
#' @examples
do_methylcibersort <- function(scm, assay= "score", ignore_cells = NULL, cell_col = "Cell", ...) {
  
  .validateExp(scm)
  .validateAssay(scm,assay)
  .validateType(ignore_cells,c("string","null"))
  .validateType(cell_col,"string")
  
  #Remove cells, if specified
  if (is.null(ignore_cells)) {
    col_idx <- 1:ncol(scm)
  } else {
    col_idx <- which(!is.na(match(colData(scm)[[cell_col]],ignore_cells)))
  }

  beta <- get_matrix(scm,assay)[,col_idx]
  rownames(beta) <- names(rowRanges(scm))
  
  cell_types <- as.factor(colData(scm)[[cell_col]][col_idx])
  
  fet <- feature.select.new(Stroma.matrix = beta, Phenotype.stroma = cell_types,
                            export = TRUE, silent = FALSE,...)
  
  return(scm[row.names(fet$SignatureMatrix),col_idx])
}


#---- feature.select.new -------------------------------------------------------------------------------------

## [ Edited Feature Selection Function from MethylCibersort 0.2.1 - AL (2021) ] ////////////////////

## Changes made to function:
## added argument to specify output directory instead of defaulting to getwd()


## [ Edited Feature Selection Function from MethylCibersort 0.2.1 - YG, DW (2018) ] ////////////////////

## Changes made to function:
## added argument & functionality for internal conversion for M-values
## added arguments & functionality for exporting of various intermediate objects
## Fixed the way DMPs are filtered to prevent the factor order / naming order from influencing the CpG result

## Base working of the function is retained and the function takes and outputs the same objects, as before.

## FUNCTION START //////////////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////////////////////////////

feature.select.new <- function(MaxDMPs = 100,                       ## maximum differentially methylated probes to use, takes n/2 from top & n/2 from bottom
                               deltaBeta = 0.2,                     ## cutoff for the minimum difference between pairwise groups by delta-beta
                               useM = FALSE,                        ## option to use M-values not beta-values, largely not useful
                               CellLines.matrix = NULL,             ## input matrix for cell line 'cancer' data
                               export = TRUE,                       ## save a table of signature results
                               export.fit = TRUE,                   ## export the limma result
                               export.cpg = TRUE,                   ## export the CpGs selected
                               output.dir = getwd(),
                               sigName = "methylCibersort",         ## name appended to start of filename
                               Stroma.matrix = NULL,                ## matrix of betas for populations
                               Phenotype.stroma = NULL,             ## pheno that corresponds to Stroma.matrix
                               FDR = 0.01,                          ## FDR cutoff
                               silent = TRUE){                      ## run function without returning output 
  
  ### import packages
  require(magrittr)
  require(matrixStats)
  require(BiocGenerics)
  require(MethylCIBERSORT)
  
  message(paste0("methylCibersort ", packageVersion("MethylCIBERSORT")))
  message(paste0("deltaBeta: ", deltaBeta, " MaxDMPs: ", MaxDMPs))
  
  ## process pheno information
  if (!is.null(ncol(CellLines.matrix))) {
    Pheno1 <- c(rep("Cancer", ncol(CellLines.matrix)))
    Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma))
  } else { 
    Pheno2 <- as.character(Phenotype.stroma)
  } # end if else 
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  if (useM){ 
    Mat3 <- minfi::logit2(Mat2) 
  } else {
    Mat3 <- Mat2
  } # end if
  
  message("Setting up for pairwise feature selection")
  ContrastMatrix <- design.pairs(levels(factor(Pheno2)))
  
  ## set up limma comparisons
  Des <- model.matrix(~0 + Pheno2)
  colnames(Des) <- rownames(ContrastMatrix)
  
  ## do limma comparisons
  Fit <- limma::lmFit(Mat3, Des) %>% 
    limma::contrasts.fit(ContrastMatrix) %>% 
    limma::eBayes()
  
  FitList <- list()
  
  ## get top results for each comparison
  for (i in 1:ncol(ContrastMatrix)) {
    tmp.name <- gsub("\\-", "\\.", colnames(ContrastMatrix)[i])
    cat(tmp.name)
    FitList[[i]] <- limma::topTable(Fit, coef = i, number = nrow(Mat2)) %>% 
      dplyr::mutate(ID = rownames(.)) %>% 
      dplyr::filter(adj.P.Val < FDR)
    cat(" ... done", "\n")
    names(FitList)[i] <- tmp.name
  } # end i
  
  if (all(c(export,export.fit))){
    message("Saving limma fit results")
    fN <- paste(sigName, deltaBeta, MaxDMPs, "Fit_list.rds", sep = "_")
    saveRDS(FitList, paste0(output.dir,"/",fN))
  } # end if
  
  message("Calculating population medians")
  ## slow
  #Split <- apply(Mat2, 1, function(x){tapply(x, as.factor(Pheno2), median)}) 
  
  Transformed <- data.frame(t(Mat2))
  
  ## get means by population
  ## faster
  Split <- split(Transformed, Pheno2) 
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split) 
  rownames(Split) <- rownames(Mat2)
  
  dbList <- list()
  message("Getting Delta Beta estimates")
  
  ## get DB estimates
  for (i in 1:ncol(ContrastMatrix)) {
    tmp.name <- gsub("\\-", "\\.", colnames(ContrastMatrix)[i])
    cat(tmp.name)
    dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
    dB <- data.frame(dB = dB, ID = rownames(Split))
    dbList[[i]] <- dB
    cat(" ... done", "\n")
    names(dbList)[i] <- tmp.name
  } # end i
  
  message("Filtering by Delta Beta")
  ## filter by DB threshold
  dbList <- lapply(dbList, function(x) dplyr::filter(x, abs(dB) > deltaBeta))
  
  message("Filtering by max DMP number")
  for (i in 1:length(FitList)) {
    A1 <- FitList[[i]]
    A1 <- dplyr::filter(A1, ID %in% dbList[[i]]$ID)
    A1 <- A1 %>% .[rev(order(.$t)), ]
    if (nrow(A1) > MaxDMPs) {
      A2 <- head(A1[, ], MaxDMPs/2) ## ///////////// Edited to correct issues with ordering of factors
      A3 <- tail(A1[, ], MaxDMPs/2) ## ///////////// affecting resulting CpG selection & erroneous 
      A1 <- rbind(A2, A3) ## /////////////////////// correlation effects, might need to be better adressed in future
    }
    FitList[[i]] <- A1
  } # end i
  
  message("Bulding unique CpG list and annotating probes")
  Nonzeros <- lapply(FitList, function(x) dplyr::select(x, ID))
  
  Nonzeros <- do.call(rbind, Nonzeros)
  #Nonzeros <- dplyr::filter(Nonzeros, !duplicated(ID))
  if (all(c(export,export.cpg))){
    message("Saving complete list of CpGs before filtering")
    fN <- paste(sigName, deltaBeta, MaxDMPs, "all_Nonzeros.txt", sep = "_")
    write.table(data.frame(ID = Nonzeros,
                           PopID = rownames(Nonzeros)), 
                file = paste0(output.dir,"/",fN), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  } # end if

  Nonzeros <- data.frame(ID = Nonzeros[-which(duplicated(Nonzeros$ID)), ], 
                         row.names = rownames(Nonzeros)[-which(duplicated(Nonzeros$ID))])
  
  Nonzeros <- data.frame(ID = Nonzeros$ID, PopID = gsub("\\.[0-9]+", "", rownames(Nonzeros)))
  Nonzeros.anno <- minfi::getAnnotation(minfiData::MsetEx, lociNames = Nonzeros$ID)
  tmp <- Nonzeros.anno[as.character(Nonzeros$ID), ]
  rownames(tmp) <- 1:nrow(tmp)
  Nonzeros.anno <- cbind(Nonzeros, tmp)

  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros.anno$ID, ]
  nrow(Mat3)
  Mat3 <- 100 * Mat3
  
  if (export) {
    
    message("Writing text files")
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), 
                                  stringsAsFactors = F), Collapsed)
    
    fN <- paste(sigName, deltaBeta, MaxDMPs, "Signature.txt", sep = "_")
    write.table(Collapsed, 
                file = paste0(output.dir,"/",fN), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    
    fN <- paste(sigName, deltaBeta, MaxDMPs, "CpG_Annotation.txt", sep = "_")
    write.table(Nonzeros.anno, 
                file = paste0(output.dir,"/",fN), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  } # end export
  
  if(!silent) return(list(SignatureMatrix = Mat3))
} # end function

## FUNCTION END ////////////////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////////////////////////////

#' 
#' CIBERSORT <- function(response,features, transform, usegenes, norm=T, nu= c(0.25,0.5,0.75), optim.nu =F, mc.cores=3, ...) {
#'   
#'   if (norm) features <- apply(features, 2, function(x) x / mean(x) * mean(response))
#'   
#'   response <- transform(response) -> allresponse
#'   response <- response[usegenes]
#'   allresponse <- allresponse[intersect(rownames(features),names(allresponse))]
#'   
#'   features <- transform(features[usegenes, ])
#'   #
#'   #   features <- apply(features,2,scale)
#'   #   response <- scale(response)
#'   
#'   #proper way to choose nu
#'   #set up CV scheme
#'   if (optim.nu & length(nu)>1) {
#'     cvorder <- sample(1:nrow(features), nrow(features))
#'     cvset <- rep(c(1:5), length.out = nrow(features))
#'     nuerror <- sapply(nu, function(n) {
#'       tuning <- do.call(rbind,parallel::mclapply(1:5, function(i) {
#'         train <- cvorder[cvset != i]
#'         test <- cvorder[cvset == i]
#'         svr <- e1071::svm(features[train,], y = response[train], type="nu-regression", kernel="linear", nu = n,  ...)
#'         out <- data.frame(real = response[test], predicted = e1071:::predict.svm(svr, newdata = features[test,]))
#'       }, mc.cores=mc.cores))
#'       sqrt(mean((tuning$predicted - tuning$real) ^2))
#'     })
#'     nu <- nu[which.min(nuerror)]
#'   }
#'   
#'   SVR <- parallel::mclapply(nu, function(n) e1071::svm(features, y = response, type="nu-regression", kernel="linear", nu = n,  ...))
#'   
#'   get_RMSE <- function(SVR, truey) {
#'     predicted <-   e1071:::predict.svm(SVR)
#'     sqrt(mean((predicted - truey) ^2))
#'   }
#'   RMSE <- sapply(SVR, get_RMSE, truey = response)
#'   SVR <- SVR[[which.min(RMSE)]]
#'   
#'   coef <- t(SVR$coefs) %*% SVR$SV
#'   coef[coef<0] <- 0
#'   
#'   
#'   #P value estimation.
#'   predicted <- e1071:::predict.svm(SVR)
#'   test_statistics <- cor(predicted , response)
#'   
#'   
#'   
#'   #CIBERSORT is then run on m*i to produce a vector of estimated cellular fractions, f*i. CIBERSORT determines the correlation coefficient R*i between the random mixture m*i and the reconstituted mixture, f*i ?? B. This process is repeated for I iterations (I = 500 in this work) to produce R*.
#'   
#'   out <- t(coef/sum(coef))
#'   attr(out, "p") <-  cor.test(predicted , response)$p.value
#'   attr(out, "SV") <- SVR$SV
#'   
#'   out
#' }
#' 
#' 
#' 
#' #'Runs CIBERSORT for decomposing bulk RNA-seq samples using a single cell reference
#' #'
#' #' This is a custom implementation of the algorithm described by Newman et al (Nautre Methods 12:453-457). CIBERSORT is an algorithm for estimating the cell type composition of a bulk sample, given a gene expression profile of the sample and a known gene expression profile for each cell type potentially contributing to the sample.
#' #'@param exprs A data frame or matrix of raw read counts of \epmh{bulk} RNA-seq samples. Column names correspond to sample names, row names to genes.
#' #'@param base A matrix of read counts representing the gene expression profiles of all cell types that might contribute to the bulk samples. See examples below for how to generate this object from an object of class  \code{seurat}.
#' #'@param design A named vector assigning sample names to sample class, see examples below.
#' #'@param markergenes A vector of genes to be included in the analysis, defaults to \code{intersect( rownames(mean_by_cluster),  rownames(exprs) )}
#' #'@param transform A function to be applied to columns of \code{exprs} and \{base} following normalization. Defaults to no transformation since bulk RNA-seq profiles are generated by pooling up RNA from constituent cell types. In the original CIBERSORT paper, a logarithmic transform was used.
#' #'@param nu Different values of nu to evaluate support vector regression at, see \code{\link{[e1071]svm}}. Nu defines how many support vectors (i.e. genes) to use in regression.
#' #'@param optim.nu In the original CIBERSORT implementation, SVR is evaluated at several values of nu and the value with the best RSME is chosen. This can lead to overfitting. If \code{optim.nu} is set to \code{TRUE}, the value for nu is chosen by cross validation, which leads to longer runtimes.
#' #'@param mc.cores Number of cores used, e.g. for the parallel evaluation at different balues of nu.
#' #'@param ... Parameters passed to \code{\link[e1071]svm}
#' #'@return A data frame in long format suitable for plotting with ggplot2.
#' #'@examples
#' #'\dontrun{
#' #'#See also package vignette CIBERORT.Rmd
#' #'#1. identify marker genes from a seurat object
#' #'NicheMarkers10x <- FindAllMarkers(NicheData10x, test.use = "roc")
#' #'usegenes <- unique(NicheMarkers10x$gene[(NicheMarkers10x$myAUC > 0.8 |NicheMarkers10x$myAUC < 0.2) ])
#' #'
#' #'#2. compute mean expression per cell type
#' #'mean_by_cluster <- do.call(cbind, lapply(unique(NicheData10x@ident), function(x) {
#' #'apply(NicheData10x@raw.data[usegenes,NicheData10x@cell.names][,NicheData10x@ident == x], 1,mean )
#' #'}))
#' #'colnames(mean_by_cluster) <- unique(NicheData10x@ident)
#' #'
#' #'#3. Create a vector that maps samples to biological class
#' #'LCM_design <- NicheMetaDataLCM$biological.class
#' #'names(LCM_design) <- NicheMetaDataLCM$id
#' #'
#' #'4. Run CIBERSORT
#' #'CIBER <- runCIBERSORT(NicheDataLCM, mean_by_cluster,usegenes, LCM_design, mc.cores=3)
#' #'}
#' #'@export
#' runCIBERSORT <- function(exprs, base,design, markergenes = intersect( rownames(mean_by_cluster),  rownames(exprs) ),transform=function(x) x,nu = c(0.25,0.5,0.75), optim.nu = F, mc.cores= 3, ...) {
#'   
#'   res <- list()
#'   for (i in 1:ncol(exprs)) {
#'     x <- exprs[,i]
#'     names(x) <- rownames(exprs)
#'     res[[i]] <- CIBERSORT(x, features=base, transform=transform, usegenes = intersect(markergenes, rownames(exprs)), nu=nu, optim.nu = optim.nu, mc.cores = mc.cores, ...)
#'   }
#'   #out <- apply(exprs,2, CIBERSORT, features=mean_by_cluster, kernel=kernel, cost =cost, method = method,alpha=alpha, gamma=gamma, transform=transform, usegenes = intersect(markergenes, rownames(exprs)), norm=norm, nu=nu)
#'   pvals <- data.frame(
#'     pvals = sapply(res, attr, "p"),
#'     samples = colnames(exprs)
#'   )
#'   
#'   svs <- lapply(res, attr, "SV")
#'   
#'   
#'   out <- do.call(cbind,res)
#'   colnames(out) <- colnames(exprs)
#'   rownames(out) <- colnames(mean_by_cluster)
#'   out <- reshape2::melt(out)
#'   out$experiment <- design[out$Var2]
#'   colnames(out) <- c("CellType","SampleID","Fraction","SampleClass")
#'   out
#' }
#' 
#' plotCIBER <- function(ciber,nrow=NULL) {
#'   ggplot(aes(x = Var1, y = value, fill=Var1), data=ciber) + geom_bar(stat="summary", fun.y=mean) + scale_x_discrete(labels = annos) + scale_fill_manual(values=colors, labels=annos, guide=F) + facet_wrap(~experiment,nrow=nrow) + geom_errorbar(stat="summary",fun.ymin = function(x) mean(x)-sd(x), fun.ymax = function(x) mean(x)+sd(x), width=0.2) + theme(axis.text.x =  element_text(angle=90, size=8)) + ylab("% of population") + xlab("")
#' }
#' 
#' 
#' get_sample <- function(runid, usegenes,npop=10,ncell=1000, replicates = 1, kernel ="radial",cost = 1, transform="log2", method="SVR",alpha=0.5, gamma=NULL, nu = c(0.25,0.5,0.75)){
#'   raw <- seurat@raw.data[,seurat@cell.names]
#'   
#'   #if (transform == "log2") transform <- log2p1 else transform <- lin
#'   
#'   what.to.sample <- sample(as.character(unique(seurat@ident)), npop)
#'   fractions <- rdirichlet(1, rep(1,npop))
#'   fractions <- round(ncell*rdirichlet(1, rep(1,npop)))
#'   
#'   out <- replicate(replicates, {
#'     gex <- lapply(1:length(what.to.sample), function(x) {
#'       use <- sample(seurat@cell.names[seurat@ident == what.to.sample[x]], fractions[x], replace=T)
#'       if (length(use)==1) raw[,use] else apply(raw[,use],1,sum)
#'     })
#'     gex <- do.call(cbind,gex)
#'     gex <- apply(gex,1,sum)
#'     
#'     
#'     # all <- rep(0, length(unique(seurat@ident)) )
#'     # names(all ) <- unique(seurat@ident)
#'     # all[what.to.sample] <- fractions
#'     #
#'     # gex <- mean_by_cluster %*% all
#'     # gex <- gex[,1]
#'     
#'     
#'     usegenes <- intersect(usegenes, names(gex))#,size = ngene)
#'     
#'     
#'     
#'     #features <- transform(mean_by_cluster[usegenes, ])
#'     test_ciber <- CIBERSORT(gex, mean_by_cluster, kernel, cost, transform, usegenes,method = method,alpha=alpha, gamma=gamma, nu=nu)
#'   })
#'   
#'   all <- rep(0, length(unique(seurat@ident)) )
#'   names(all ) <- unique(seurat@ident)
#'   all[what.to.sample] <- fractions/ncell
#'   
#'   
#'   test_ciber <- apply(out,1,mean)
#'   if (method == "SVR") merged <- data.frame(truth = all, result = test_ciber, cluster = names(all),runid = runid) else merged <- data.frame(truth = all, result = test_ciber, cluster = names(all),runid = runid)
#'   merged
#'   
#' }
#' 
#' runSampling <- function(usegenes, nsamples=15,npop=10,replicates=1,ncell=1000, kernel ="radial",cost = 1, transform="log2", method = "SVR",alpha=0.5, gamma = NULL, nu=c(0.25,0.5,0.75)) {
#'   samples <- lapply(1:nsamples, get_sample, usegenes=usegenes, npop=npop,ncell=ncell, replicates=replicates,kernel =kernel,cost = cost, transform=transform, method=method, alpha=alpha, gamma=gamma, nu=nu)
#'   rsqrs <- sapply(samples, function(x) cor(x$truth,x$result)^2)
#'   rho <- sapply(samples, function(x) cor(x$truth,x$result,method="spearman"))
#'   out <- do.call(rbind, samples)
#'   out$gamma <- gamma; out$npop <- npop
#'   list(rsqr = mean(rsqrs), rho = mean(rho), samples = out, allrho = rho, allrsqrs = rsqrs)
#' }
#' 
#' 
#' 
#' 
#' get_sample_fixed <- function(runid, populations,frequencies,usegenes,ncell=1000, kernel ="radial",cost = 1, transform="log2", method="SVR",alpha=0.5, gamma=NULL){
#'   raw <- as.matrix(raw)
#'   # if (transform == "log2") transform <- log2p1 else transform <- lin
#'   
#'   what.to.sample <- populations
#'   fractions <- round(ncell*frequencies)
#'   gex <- lapply(1:length(what.to.sample), function(x) {
#'     use <- sample(cell.names[ident == what.to.sample[x]], fractions[x], replace=T)
#'     if (length(use)==1) raw[,use] else apply(raw[,use],1,sum)
#'   })
#'   gex <- do.call(cbind,gex)
#'   gex <- apply(gex,1,sum)
#'   all <- rep(0, length(unique(ident)) )
#'   names(all ) <- unique(ident)
#'   all[what.to.sample] <- fractions/ncell
#'   
#'   # all <- rep(0, length(unique(seurat@ident)) )
#'   # names(all ) <- unique(seurat@ident)
#'   # all[what.to.sample] <- fractions
#'   #
#'   # gex <- mean_by_cluster %*% all
#'   # gex <- gex[,1]
#'   
#'   
#'   usegenes <- intersect(usegenes, names(gex))#,size = ngene)
#'   # mean_by_cluster <- apply(mean_by_cluster, 2, function(x) x / mean(x) * mean(gex))
#'   
#'   
#'   
#'   # features <- transform(mean_by_cluster[usegenes, ])
#'   # test_ciber <- CIBERSORT(transform(gex[usegenes]), features, kernel, cost, method = method,alpha=alpha)
#'   test_ciber <- CIBERSORT(gex, mean_by_cluster, kernel, cost, transform = transform, usegenes = usegenes, method = method,alpha=alpha, gamma = gamma)
#'   if (method == "SVR") merged <- data.frame(truth = all, result = test_ciber[,1], cluster = names(all),runid = runid, stringsAsFactors = F) else merged <- data.frame(truth = all, result = test_ciber, cluster = names(all),runid = runid, stringsAsFactors = F)
#'   merged
#'   
#' }

## generate random proportions for  cell types
generate.random.prop <- function(nr = 100, ## number of rows
                                 nc = 12) ## number of columns
{
  ## create empty matrix
  prop.mat <- matrix(data = 0, nrow = nr, ncol = nc)
  dimnames(prop.mat) <- list(1:nrow(prop.mat), NULL)
  ## generate a set sequence of proportions for the population being generated
  ctoi <- seq(0, 1, 1/nr)
  ## fill in the rest with values from uniform distribution 
  for (i in 1:nrow(prop.mat)){
    tmp.set <- ctoi[i] ## set your static value
    tmp.max <- 1-tmp.set ## maximum possible value given 1-static value
    tmp.val <- rep(0, nc-1) ## populate the slots with 0
    for (j in sample(1:(nc-2))){
      tmp.val[j] <- runif(1, min = 0, max = tmp.max-sum(tmp.val)) ## generate random values
    } # end j
    tmp.val[nc-1] <- tmp.max-sum(tmp.val) ## generate the final value for the row to add up to 1
    prop.mat[i, ] <- c(sample(tmp.val, nc-1, replace = FALSE), tmp.set) ## sample without replacement then append the set value
  } # end i
  return(prop.mat)
} # end function

## generate synthetic proportions for use in benchmarking
make.synth.prop <- function(synth.colnames = NULL,      ## factor names, 1 per population
                            pop.rows = 100)             ## number of rows in matrix to populate per population
{ 
  synth.list <- list()
  ## run through the list of factors, for each factor generate mixes where 
  ## for n pop.rows the proportion is a set sequence and the rest are random
  for (i in 1:length(synth.colnames)){
    idx <- 1:(length(synth.colnames)-1)
    idx <- R.utils::insert(idx, ats = i, length(synth.colnames))
    synth.list[[synth.colnames[i]]] <- generate.random.prop(nr = pop.rows, 
                                                            nc = length(synth.colnames))[, idx]
    colnames(synth.list[[synth.colnames[i]]]) <- synth.colnames
  } # end i
  synth.prop.mat <- do.call("rbind", synth.list)
  rownames(synth.prop.mat) <- paste0("synth.", 1:nrow(synth.prop.mat))
  return(synth.prop.mat)
} # end function

make.synth.mix <- function(input.data = NULL,                                  ## data matrix to be converted to synth mixes
                           pop.factor = NULL,                                  ## factor with levels of the populations of data matrix
                           pop.rows = 100,                                     ## how many rows to generate for each population
                           output.dir = getwd(),                               ## location to save proportions table and resulting mixtures
                           output.name = gsub("-", "_", Sys.Date()),           ## name to append to the filenames
                           n.cores = 1){                                       ## do in parallel for n cores
  synth.prop.mat <- make.synth.prop(levels(pop.factor), pop.rows)
  write.table(synth.prop.mat, 
              paste0(output.dir, "/", output.name, "_synth_prop.txt"), 
              quote = FALSE, 
              row.names = FALSE)
  soi.mean <- t(apply(input.data, 1, function(x){tapply(x, pop.factor, mean)}))
  require(parallel)
  synth.list <- mclapply(1:nrow(synth.prop.mat), 
                         function(x) {apply(soi.mean[, colnames(synth.prop.mat)], 
                                            1, 
                                            function(y) weighted.mean(y, 
                                                                      w = synth.prop.mat[x, ]))},
                         mc.cores = n.cores)
  synth.df <- as.data.frame(do.call(cbind, synth.list))
  synth.df <- cbind(rownames(soi.mean), synth.df)
  colnames(synth.df) <- c("NAMES", rownames(synth.prop.mat))
  write.table(synth.df, 
              paste0(output.dir, "/", output.name, "_synth_mix.txt"), 
              quote = FALSE, 
              row.names = FALSE,
              sep = "\t")
} # end function 

# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
# install.packages('e1071') 
# install.packages('parallel')
# install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
# library(e1071)
# library(parallel)
# library(preprocessCore)

#Core algorithm
#' Title
#'
#' @param X
#' @param y
#' @param absolute
#' @param abs_method
#'
#' @export
#'
CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#do permutations
#' Title
#'
#' @param perm Permutation
#' @param X
#' @param Y
#' @param absolute
#' @param abs_method
#'
#' @export
#'
doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}




#' CIBERSORT is an analytical tool developed by Newman et al. to provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data.

#' @param sig_matrix  Cell type GEP barcode matrix: row 1 = sample labels; column 1 = gene symbols; no missing values; default =LM22.txt download from CIBERSORT (https://cibersort.stanford.edu/runcibersort.php)
#' @param mixture_file  GEP matrix: row 1 = sample labels; column 1 = gene symbols; no missing values
#' @param perm Set permutations for statistical analysis (≥100 permutations recommended).
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param absolute  Run CIBERSORT in absolute mode (default = FALSE)
#' - note that cell subsets will be scaled by their absolute levels and will not be represented as fractions
#' (to derive the default output, normalize absolute levels such that they sum to 1 for each mixture sample)
#' - the sum of all cell subsets in each mixture sample will be added to the ouput ('Absolute score').
#' If LM22 is used, this score will capture total immune content.
#' @param abs_method  if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#' - sig.score = for each mixture sample, define S as the median expression
#' level of all genes in the signature matrix divided by the median expression
#' level of all genes in the mixture. Multiple cell subset fractions by S.
#' - no.sumto1 = remove sum to 1 constraint
#' @author Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' @return cibersrot with immune cell fractions
#' @export
#' @import e1071
#' @import parallel
#' @import preprocessCore
#' @import tidyverse
#' @examples
#' cibersort<-CIBERSORT(sig_matrix = lm22, mixture_file = eset_ec, perm = 1000, QN=TRUE, absolute=FALSE)
#' head(cibersort)

CIBERSORT <- function(sig_matrix = lm22, mixture_file, perm, QN = TRUE, absolute, abs_method='sig.score'){
  
  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0){
    stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score')
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  #read in data
  X<- sig_matrix
  # X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- rownames_to_column(mixture_file,var = "symbol")
  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  dups <- dim(Y)[1] - length(unique(Y[,1]))
  if(dups > 0) {
    warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
    rownames(Y) <- make.names(Y[,1], unique=TRUE)
  }else {rownames(Y) <- Y[,1]}
  Y <- Y[,-1]
  ###################################
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y,na.rm=T) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  # library(preprocessCore)
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Yorig),1)
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    if(absolute && abs_method == 'sig.score') {
      w <- w * median(Y[,itor]) / Ymedian
    }
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  # write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
  obj
}

# design.pairs <-
#   function(levels) {
#     n <- length(levels)
#     design <- matrix(0,n,choose(n,2))
#     rownames(design) <- levels
#     colnames(design) <- 1:choose(n,2)
#     k <- 0
#     for (i in 1:(n-1))
#       for (j in (i+1):n) {
#         k <- k+1
#         design[i,k] <- 1
#         design[j,k] <- -1
#         colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
#       }
#     design
#   }

#' @export
FeatureSelect.V4 <- function( CellLines.matrix = NULL, Heatmap = FALSE, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, deltaBeta, FDR , MaxDMRs = 1000) {
  
  #This function uses pairwise limma
  #This may be rather slow but whatever - it shoul do better on the basis of the nearest cell types
  
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)
  require(limma)
  
  if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) } else { Pheno2 <- as.character(Phenotype.stroma)}
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  
  message("Setting up for pairwise feature selection")
  
  
  #Work out the linear model fit here
  
  ContrastMatrix <- design.pairs(levels(factor(Pheno2)))
  Des <- model.matrix(~0 + Pheno2)
  colnames(Des) <- rownames(ContrastMatrix)
  Fit <- lmFit(Mat2, Des)%>%
    contrasts.fit(., ContrastMatrix)%>%
    eBayes(.)
  
  
  FitList <- list()
  for(i in 1:ncol(ContrastMatrix)) {
    
    FitList[[i]] <- topTable(Fit, coef = i, number = nrow(Mat2))%>%
      mutate(ID = rownames(.))%>%
      filter(adj.P.Val < FDR)
    
    message(paste0(i, " done"))
    
    
  }
  
  
  
  #Here goes the pairwise dB thing
  
  Transformed <- data.frame(t(Mat2))
  Split <- split(Transformed, Pheno2)
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split)
  rownames(Split) <- rownames(Mat2)
  
  #Then I need to actually annotate each one of these comparison topTables
  
  dbList <- list()
  message("Getting Delta Beta estimates")
  for(i in 1:ncol(ContrastMatrix)) {
    
    dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
    dB <- data.frame(dB = dB, ID = rownames(Split))
    dbList[[i]] <- dB
    message(paste0 (i, " done"))
  }
  
  
  #Filter by thresholds
  
  
  dbList <- lapply(dbList, function(x) filter(x, abs(dB) > deltaBeta))
  for(i in 1:length(FitList)) {
    
    A1 <- FitList[[i]]
    A1 <- filter(A1 , ID %in% dbList[[i]]$ID)
    A1 <- A1%>%.[rev(order(.$t)),]
    if(nrow(A1) > MaxDMRs) { A1 <-  A1[1:MaxDMRs,]                   }
    FitList[[i]] <- A1
  }
  
  Nonzeros <- lapply(FitList, function(x) dplyr::select(x,ID))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]
  
  if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }
  
  #Print number of selected probes
  nrow(Mat3)
  Mat3 <- 100 * Mat3
  
  #Then , export object
  if(export) {
    
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), stringsAsFactors = F),
                       Collapsed)
    
    fN <- paste0(sigName, "_Signature.txt")
    write.table(Collapsed, file = fN, sep = "\t", row.names = FALSE, quote = FALSE )
    
  }
  
  return(list(SignatureMatrix = Mat3))
  
}



#This function creates the pairs for the pairwise matrices
design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }






