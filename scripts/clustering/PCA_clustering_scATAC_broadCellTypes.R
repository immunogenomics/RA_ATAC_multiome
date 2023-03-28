print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Matrix.utils)
    library(Seurat)
    library(Signac)
    library(RANN)
    library(parallel)
    library(harmony)
    library(umap)
    library(ggplot2)
    library(viridis)
    library(irlba)
    library(tidyr)
    library(plyr)
    library(dplyr)
    library(tibble)
    library(gridExtra)
    library(lme4)
    library(argparse)
})

#we did not use this option, but keeping for posterity
FindVariableGenesSeurat <- function (data, x.low.cutoff = 0.1, x.high.cutoff = 8,
                                     y.cutoff = 1, y.high.cutoff = Inf, num.bin = 0,
                                     binning.method = "equal_width", sort.results = TRUE,
                                     display.progress = TRUE, ...)
{
    genes.use <- rownames(data)
    if (class(data) != "dgCMatrix") {
        data <- as(as.matrix(data), "dgCMatrix")
    }
    ## (1) get means and variances
    gene.mean <- FastExpMean(data, display.progress)
    names(gene.mean) <- genes.use
    gene.dispersion <- FastLogVMR(data, display.progress)
    names(gene.dispersion) <- genes.use

    gene.dispersion[is.na(x = gene.dispersion)] <- 0
    gene.mean[is.na(x = gene.mean)] <- 0

    mv.df <- data.frame(gene.mean, gene.dispersion)
    rownames(mv.df) <- rownames(data)

    ## (OPTIONAL) do the binning correction
    if (num.bin > 0) {
      if (binning.method == "equal_width") {
          data_x_bin <- cut(x = gene.mean, breaks = num.bin)
      }
      else if (binning.method == "equal_frequency") {
          data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean >
              0], probs = seq(0, 1, length.out = num.bin))))
     }
      else {
          stop(paste0("Invalid selection: '", binning.method,
              "' for 'binning.method'."))
      }
      names(x = data_x_bin) <- names(x = gene.mean)
      mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = mean)
      sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
          FUN = sd)
      gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
      gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0

      mv.df$gene.dispersion.scaled <- gene.dispersion.scaled
    }

    return(mv.df)
}
environment(FindVariableGenesSeurat) <- asNamespace("Seurat")

#we did not use these plots in the final manuscript
plot_overlays <- function(toPlot_df, toPlot_x, toPlot_y, toPlot_covar, toPlot_outPrefix, toRet=FALSE){
    if(toRet & length(toPlot_covar)!=1) stop('Returning plot can only have one covariate to overlay')
    
    outFile <- paste(sep="",toPlot_outPrefix,".png")
    g <- ggplot(toPlot_df, aes_string(x=toPlot_x, y=toPlot_y)) +
        geom_point(size=1,alpha=0.5) +
        theme_classic()
    if(!toRet) ggsave(outFile,plot=g)
    
    for(i in toPlot_covar){
        outFile <- paste(sep="",toPlot_outPrefix,"_col-",i,".png")
        
        if(is.numeric(toPlot_df[,i])){
            if(any(is.na(toPlot_df[,i]))){
                thisPlot_df <- toPlot_df
            } else {
                thisPlot_df <- toPlot_df[order(toPlot_df[,i]),]
            }
            g <- ggplot(thisPlot_df, aes_string(x=toPlot_x, y=toPlot_y,color=i)) +
                geom_point(size=1,alpha=0.5) +
                theme_classic() + 
                scale_color_viridis_c(option = "plasma")
        } else {
            g <- ggplot(toPlot_df, aes_string(x=toPlot_x, y=toPlot_y,color=i)) +
                geom_point(size=1,alpha=0.5) +
                theme_classic()
        }
        
        if(!toRet) ggsave(outFile,plot=g)
    }
    
    if(toRet) return(g)

}

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0) 
{
    my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", 
        eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn_res) <- row.names(data.use)
    colnames(snn_res) <- row.names(data.use)
    return(snn_res)
}
environment(BuildSNNSeurat) <- asNamespace("Seurat") 

getClusters <- function(adt_harmony_res, resolution_list){
    snn_ref <- BuildSNNSeurat(adt_harmony_res, nn.eps = 0)
    
    ids_ref <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
        Seurat:::RunModularityClustering(SNN = snn_ref, modularity = 1,
            resolution = res_use, algorithm = 1, n.start = 20,
            n.iter = 20, random.seed = 100, print.output = FALSE,
            temp.file.location = NULL, edge.file.name = NULL)
    }, mc.preschedule = FALSE, mc.cores = min(20, length(resolution_list))))

    ids_ref <- as.data.frame(ids_ref)
    rm(snn_ref)
    colnames(ids_ref) <- sprintf("res_%.2f", resolution_list)
    rownames(ids_ref) <- rownames(adt_harmony_res)
    
    return(ids_ref)
}

#we did not use these plots in the final manuscript
cluster_UMAP_funct <- function(toDo_mat, res_list, toDo_type, toDo_metadata, toDo_covar, toDo_outPrefix){
    print_time(paste(toDo_type,"Clusters"))
    cluster_res <- getClusters(toDo_mat,res_list)
    outFile <- paste(sep="_",toDo_outPrefix,toDo_type,"Clusters.rds")
    saveRDS(cluster_res,outFile)
    cluster_df <- data.frame(apply(cluster_res, 2, as.character),stringsAsFactors=FALSE)
    
    
    print_time(paste(toDo_type,"UMAP"))
    umap_res <- umap(toDo_mat, n_neighbors = 30, metric = "cosine", min_dist = .3)
    outFile <- paste(sep="_",toDo_outPrefix,toDo_type,"UMAP.rds")
    saveRDS(umap_res,outFile)
    
    
    print_time(paste("Plotting",toDo_type,"UMAP"))
    umap_only_df <- data.frame("UMAP1"=umap_res$layout[, 1],"UMAP2"=umap_res$layout[, 2])
    umap_df <- cbind(umap_only_df,toDo_metadata,cluster_df)
    umap_df_scrambled <- umap_df[sample(nrow(umap_df), nrow(umap_df)), ]
    saveRDS(umap_df_scrambled, file=paste(sep="_",toDo_outPrefix,toDo_type,"UMAP",'scrambledDF.rds'))
    plot_overlays(umap_df_scrambled, "UMAP1", "UMAP2", c(toDo_covar,colnames(cluster_res)),
                  paste(sep="_",toDo_outPrefix,toDo_type,"UMAP"))
    
}


parser <- ArgumentParser(description='Clustering pipeline for raw counts matrix')
parser$add_argument('mat_file', metavar='mat_file', help='Matrix file')
parser$add_argument('cellData_file', metavar='cellData_file', help='Cell metadata file')
parser$add_argument('numPCs', metavar='numPCs', type='integer', help='Number of PCs to calculate')
parser$add_argument('covars', metavar='covars', help='Covariates to Harmonize upon; separated by commas')
parser$add_argument('featSel', metavar='featSel', choices=c('top', 'varNum', 'varDisp', 'perc', 'not'), help='Feature Selection methods: top, varNum, varDisp, perc, not')
parser$add_argument('featArg', metavar='featArg', type='double', help='Feature Argument')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
parser$add_argument('--seed', metavar='seed', type='integer', default=1234567890,help='Seed for randomization')
args <- parser$parse_args()

mat_file <- args$mat_file
cellData_file <- args$cellData_file
numPCs <- args$numPCs
covars <- args$covars
featSel <- args$featSel
featArg <- args$featArg
outDir <- args$outDir
prefix <- args$prefix
seed <- args$seed


print_time("Argument Checking")
if(!all(file.exists(c(mat_file,cellData_file)))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time(paste("Set Seed:",seed))
set.seed(seed)


print_time("Load Matrix")
ori_mat <- readRDS(mat_file)
dim(ori_mat)
ori_mat[1:3,1:3]


print_time("Cell Data Matrix")
cellData_df <- readRDS(cellData_file)
dim(cellData_df)
head(cellData_df)

covariates <- unlist(strsplit(covars,","))
if(!all(covariates %in% colnames(cellData_df))) stop("Covariates not in cell metadata matrix.")


#options: nothing, features in x% of cells, top features, most variable features
print_time("Feature Selection")
if(featSel == "top"){
    if(featArg != as.integer(featArg) | featArg < 0 | featArg > nrow(ori_mat)) stop("Feature Argument for 'top' must be an integer between 0 and the number of features in the matrix.")
    
    featSums <- rowSums(ori_mat)
    topFeat <- names(sort(featSums,decreasing = TRUE)[1:featArg])
    mat <- ori_mat[topFeat,]
} else if(featSel == "varNum" | featSel == "varDisp"){
    if(featSel == "varNum" & (featArg != as.integer(featArg) | featArg < 0 | featArg > nrow(ori_mat))) stop("Feature Argument for 'varNum' must be an integer between 0 and the number of features in the matrix.")
    
    mean_disp_df <- FindVariableGenesSeurat(data = ori_mat)
    dim(mean_disp_df)
    head(mean_disp_df)
    saveRDS(mean_disp_df,paste(sep="",outPrefix,'_Mean-Disp.rds'))
    
    g <- ggplot(mean_disp_df, aes_string(x='gene.dispersion', y="gene.mean")) +
            geom_point(size=1,alpha=0.5) +
            theme_classic()
    ggsave(paste(sep="",outPrefix,"_Mean-Disp.png"),plot=g)
    
    if(featSel == "varNum") mat <- ori_mat[rownames(mean_disp_df[sort(mean_disp_df$gene.dispersion,decreasing = TRUE, index.return=TRUE)$ix[1:featArg],]),]
    if(featSel == "varDisp"){
        mat <- ori_mat[which(mean_disp_df$gene.dispersion > featArg),]
        if(nrow(mat) == 0) stop("Feature Argument for 'varDisp' resulted in an empty matrix.")
    }
} else if(featSel == "perc"){
    if(featArg < 0 | featArg > 1 ) stop("Feature Argument for 'perc' must be between 0 and 1.")
    
    meta = data.frame(
        feat_row = rownames(ori_mat),
        nCells = rep(ncol(ori_mat), nrow(ori_mat))
    )
    meta$nCells_prest <- Matrix::rowSums(ori_mat > 0)
    meta$nCells_prest_perc <- meta$nCells_prest / meta$nCells
    
    feat_use <- meta[which(meta$nCells_prest_perc > featArg), ]$feat_row
    
    mat <- ori_mat[feat_use,]
} else if(featSel == "not"){
    mat <- ori_mat[1:nrow(ori_mat),1:ncol(ori_mat)]
} else {
    stop("Feature Selection choice not valid")
}
dim(mat)
mat[1:3,1:3]
outFile <- paste(sep="",outPrefix,"_initial_featMat.rds")
saveRDS(mat,outFile)


print_time("Cell metadata for subsetted matrix")
if(!all(colnames(mat) %in% rownames(cellData_df))) stop("colnames from matrix not in cell data matrix.")
meta_data <- as.data.frame(cellData_df[colnames(mat),])
dim(meta_data)
head(meta_data)
outFile <- paste(sep="",outPrefix,"_initial_metadata.rds")
saveRDS(meta_data,outFile)


print_time("TF-IDF")
mat_norm <- TF.IDF(data=mat)
outFile <- paste(sep="",outPrefix,"_TF-IDF.rds")
saveRDS(mat_norm,outFile)


print_time("PCA")
pcaRes <- prcomp_irlba(t(mat_norm), n = numPCs, center = TRUE, scale = TRUE)
outFile <- paste(sep="",outPrefix,"_PCA.rds")
saveRDS(pcaRes,outFile)


print_time("Formatting PC Matrix")
pcaMat <- pcaRes$x
identical(rownames(pcaMat),rownames(meta_data))
rownames(pcaMat) <- colnames(mat)
if(!identical(rownames(pcaMat),rownames(meta_data))) stop("Rownames not the same between PCA and meta data matrices.")
if(ncol(pcaMat) < numPCs) stop("Asked for too many PCs for this matrix.")


print_time("Saving matrices")
saveRDS(pcaMat,paste(sep="",outPrefix,"_pcaMat.rds"))
saveRDS(mat,paste(sep="",outPrefix,"_featMat.rds"))
saveRDS(meta_data,paste(sep="",outPrefix,"_metadata.rds"))


print_time("PCA Plots part 1")
pca_df <- data.frame('PC1'=pcaMat[,'PC1'],'PC2'=pcaMat[,'PC2'])
pca_df <- cbind(pca_df,meta_data)
pca_df_scrambled <- pca_df[sample(nrow(pca_df), nrow(pca_df)), ]
plot_overlays(pca_df_scrambled, "PC1", "PC2", colnames(meta_data), paste(sep="",outPrefix,"_PCA"))


print_time("Cluster + UMAP function for PCA matrix")
res_list <- c(0.1, 0.5, 0.9)
cluster_UMAP_funct(pcaMat, res_list, "PCA", meta_data, colnames(meta_data), outPrefix)


print_time("Harmony")
harmony_embed <- HarmonyMatrix(pcaMat,meta_data,covariates,do_pca = FALSE)
outFile <- paste(sep="",outPrefix,"_Harmony.rds")
saveRDS(harmony_embed,outFile)


print_time("Harmony Plots")
harmony_df <- data.frame('PC1'=harmony_embed[,'PC1'],'PC2'=harmony_embed[,'PC2'])
harmony_df <- cbind(harmony_df,meta_data)
harmony_df_scrambled <- harmony_df[sample(nrow(harmony_df), nrow(harmony_df)), ]
plot_overlays(harmony_df_scrambled, "PC1", "PC2", colnames(meta_data), paste(sep="",outPrefix,"_Harmony"))


print_time("Cluster + UMAP function for Harmony matrix")
cluster_UMAP_funct(harmony_embed, res_list, "Harmony", meta_data, colnames(meta_data), outPrefix)


print_time("Done.")

