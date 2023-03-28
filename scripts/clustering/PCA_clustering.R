print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Seurat)
    library(irlba)
    library(ggplot2)
    library(umap)
    library(ggrepel)
    library(parallel)
    library(harmony)
    library(aricode)
    library(gtools)
    library(stringr)
    library(plyr)
    library(argparse)
})

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

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}

plot_order <- function(uniq_val){
    if(all(uniq_val %in% grep('[a-zA-Z]+-[0-9]+:',uniq_val,value=TRUE)))
        ret <- mixedsort_new(uniq_val)
    else if(all(uniq_val %in% grep('[a-zA-z]+_c[0-9]+',uniq_val,value=TRUE)))
        ret <- mixedsort(uniq_val)
    else if(!any(is.na(as.numeric(uniq_val))))
        ret <- as.character(sort(as.numeric(uniq_val)))
    else
        ret <- sort(uniq_val)
    
    return(ret)
}

plot_dot <- function(toPlot, plotX, plotY, plotTit, plotPrefix,
                     plotHeight=6, plotWidth=6, plotMeta=FALSE, plotMetaCols=NA,
                     plotColors=NA){
    g <- ggplot(toPlot,aes_string(x=plotX,y=plotY)) + geom_point(size=1,alpha=0.5) + 
            theme_classic(base_size = 15) + ggtitle(plotTit)
    ggsave(paste(sep="",plotPrefix,".png"),
           plot=g,units='in',height=plotHeight,width=plotWidth)
    
    if(plotMeta){
        for(col in plotMetaCols){
            widthAdd = floor(max(nchar(unique(toPlot[,col])),nchar(col))/7)
            if(length(unique(toPlot[,col]))>17) widthAdd=widthAdd*2
            
            g <- ggplot(toPlot,aes_string(x=plotX,y=plotY,color=col)) + 
                    geom_point(size=1,alpha=0.5) + theme_classic(base_size = 15) + 
                    ggtitle(plotTit)
            
            if(all(unique(toPlot[,col]) %in% names(plotColors))){
                g <- g + scale_color_manual(values=plotColors)
            }
            ggsave(paste(sep="",plotPrefix,"_col-",col,".png"),
                   plot=g,units='in',height=plotHeight,width=plotWidth+widthAdd)
        }
    }
}

#modified from gtools - to allow for the hyphen to be used as a delimiter or a negative number.
mixedorder_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE, 
    numeric.type = c("decimal", "roman"), roman.case = c("upper", 
        "lower", "both"),keepNegative=FALSE) 
{
    numeric.type <- match.arg(numeric.type)
    roman.case <- match.arg(roman.case)
    if (length(x) < 1) 
        return(NULL)
    else if (length(x) == 1) 
        return(1)
    if (!is.character(x)) 
        return(order(x, decreasing = decreasing, na.last = na.last))
    delim = "\\$\\@\\$"
    if (numeric.type == "decimal") {
        if(keepNegative)
            regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
        else
            regex <- "((?:(?i)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)))"
        numeric <- function(x) as.numeric(x)
    }
    else if (numeric.type == "roman") {
        regex <- switch(roman.case, both = "([IVXCLDMivxcldm]+)", 
            upper = "([IVXCLDM]+)", lower = "([ivxcldm]+)")
        numeric <- function(x) roman2int(x)
    }
    else stop("Unknown value for numeric.type: ", numeric.type)
    nonnumeric <- function(x) {
        ifelse(is.na(numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    delimited <- gsub(regex, paste(delim, "\\1", delim, sep = ""), 
        x, perl = TRUE)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) {x[x > ""]})
    suppressWarnings(step1.numeric <- lapply(step1, numeric))
    suppressWarnings(step1.character <- lapply(step1, nonnumeric))
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) {sapply(step1.numeric, 
        function(x) {x[i]})})
    step1.character.t <- lapply(1:maxelem, function(i) {sapply(step1.character, 
        function(x) {x[i]})})
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) {as.numeric(factor(x))})
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric), 
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric, 
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0) 
        if (is.na(na.last)) 
            order.frame[which.nas, ] <- NA
        else if (na.last) 
            order.frame[which.nas, ] <- Inf
        else order.frame[which.nas, ] <- -Inf
    if (length(which.blanks) > 0) 
        if (is.na(blank.last)) 
            order.frame[which.blanks, ] <- NA
        else if (blank.last) 
            order.frame[which.blanks, ] <- 1e+99
        else order.frame[which.blanks, ] <- -1e+99
    order.frame <- as.list(order.frame)
    order.frame$decreasing <- decreasing
    order.frame$na.last <- NA
    retval <- do.call("order", order.frame)
    return(retval)
}

mixedsort_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE, 
    numeric.type = c("decimal", "roman"), roman.case = c("upper", 
        "lower", "both"),keepNegative=FALSE) 
{
    ord <- mixedorder_new(x, decreasing = decreasing, na.last = na.last, 
        blank.last = blank.last, numeric.type = numeric.type, 
        roman.case = roman.case, keepNegative = keepNegative)
    x[ord]
}


parser <- ArgumentParser(description='PCA & Harmony clustering of a feature selected file')
parser$add_argument('featSel_file', metavar='featSel_file', help='Feature selection file (cells x feature) - assumes it is center/scale already')
parser$add_argument('--centerScale', action='store_true', help='If given, center/scale')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='double', help='randomization seed; default 1234567890')
parser$add_argument('plotTitle', metavar='plotTitle', help='title for plots')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('--meta_cols', metavar='meta_cols', help='columns to plot in meta file; comma delimited list; if not given, all columns.')
parser$add_argument('harmony_covars', metavar='harmony_covars', help='meta columns to use in harmony')
parser$add_argument('--color_file', metavar='color_file', help='file for colors')
parser$add_argument('--donor_color_file', metavar='donor_color_file', help='sample colors - vector in RDS file; if not given, just uses default ggplot colors')
parser$add_argument('this_nPCs', metavar='this_nPCs', type='integer', help='number of PCs to calculate')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

featSel_file <- args$featSel_file
centerScale <- args$centerScale
seed <- args$seed
plotTitle <- args$plotTitle
meta_file <- args$meta_file
meta_cols <- args$meta_cols
harmony_covars <- args$harmony_covars
color_file <- args$color_file
donor_color_file <- args$donor_color_file
this_nPCs <- args$this_nPCs
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(featSel_file,meta_file))) stop("Input file(s) don't exist.")

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


print_time("Overall seed")
set.seed(seed)


print_time('Load feature selection file')
featSel_mat <- readRDS(featSel_file)
if(centerScale){
    featSel_mat <- scale(featSel_mat,center=TRUE)
}
featSel_mat[1:3,1:3]
cat(paste("first column: mean:",mean(featSel_mat[,1]),"var: ",var(featSel_mat[,1]),"\n"))
summary(featSel_mat[,1])


print_time('Load meta file')
meta <- readRDS(meta_file)
meta_col_toUse <- parse_columns(meta_cols)
harmony_col_toUse <- parse_columns(harmony_covars)
meta_col_toUse <- unique(c(meta_col_toUse,harmony_col_toUse)) #definitely want to plot the harmony covariates!

use_meta <- FALSE
if(all(rownames(featSel_mat) %in% rownames(meta))){
    meta <- meta[rownames(featSel_mat),]

    if(any(meta_col_toUse %in% colnames(meta))){
        meta_col_toUse <- meta_col_toUse[which(meta_col_toUse %in% colnames(meta))]
        meta <- meta[,meta_col_toUse]
        use_meta <- TRUE
    }
}

color_vec <- NA
if(use_meta & !is.null(color_file)){
    if(file.exists(color_file)){
        print_time('colors file')
        color_vec <- readRDS(color_file)
    }
}

if(use_meta & !is.null(donor_color_file)){
    if(file.exists(donor_color_file)){
        print_time('donor colors file')
        donor_color_vec <- readRDS(donor_color_file)
        
        if(any(is.na(color_vec))){
            cat('donor colors are the only colors\n')
            color_vec <- donor_color_vec
        } else if(any(names(color_vec) %in% names(donor_color_vec))){
            cat('not adding donor colors since there is conflict with the color vec\n')
        } else {
            cat('adding donor colors to existing color vec\n')
            color_vec <- c(color_vec,donor_color_vec)
        }
    }
}


print_time("PCA")
set.seed(seed)
pcaRes <- prcomp_irlba(featSel_mat, n = this_nPCs, center = FALSE, scale = FALSE) #matrix already has center/scale
saveRDS(pcaRes,paste(sep='',outPrefix,'_PC',this_nPCs,'_res.rds'))

pcaMat <- pcaRes$x
rownames(pcaMat) <- rownames(featSel_mat)
dim(pcaMat)
pcaMat[1:3,1:3]
saveRDS(pcaMat,paste(sep='',outPrefix,'_PC',this_nPCs,'_mat.rds'))


print_time('Clustering')
res_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)
pca_cluster_res <- getClusters(pcaMat,res_list)
saveRDS(pca_cluster_res, paste(sep='',outPrefix,'_PC',this_nPCs,'_cluster.rds'))
pca_cluster_df <- data.frame(apply(pca_cluster_res, 2, as.character),stringsAsFactors=FALSE)
cluster_cols <- colnames(pca_cluster_df)


print_time('UMAP')
set.seed(seed)
pca_umap <- umap(pcaMat, n_neighbors = 30, metric = "cosine", min_dist = .3) 
saveRDS(pca_umap,paste(sep='',outPrefix,'_PC',this_nPCs,'_UMAP.rds'))
pca_umap_df <- data.frame("UMAP1"=pca_umap$layout[, 1],"UMAP2"=pca_umap$layout[, 2])


print_time('PCA full DF')
pca_full_df <- cbind(meta,pcaMat,pca_cluster_df,pca_umap_df)
set.seed(seed)
pca_full_df <- pca_full_df[sample(nrow(pca_full_df),nrow(pca_full_df)),]
saveRDS(pca_full_df,paste(sep='',outPrefix,'_PC',this_nPCs,'_fullDF.rds'))


print_time('Plot initial PCs')
plot_dot(pca_full_df,'PC1','PC2',plotTitle, 
         paste(sep="",outPrefix,'_PC',this_nPCs,"_PC1-PC2"), 
         plotMeta=use_meta, plotMetaCols=c(meta_col_toUse,cluster_cols), plotColors=color_vec)


print_time("Plot UMAPs")
plot_dot(pca_full_df,'UMAP1','UMAP2',
         paste(plotTitle,this_nPCs,'PCs'), 
         paste(sep="",outPrefix,'_PC',this_nPCs,'_UMAP'), 
         plotMeta=use_meta, plotMetaCols=c(meta_col_toUse,cluster_cols), plotColors=color_vec)


print_time('Harmony')
set.seed(seed)
harmonyObj <- HarmonyMatrix(pcaMat, meta, harmony_col_toUse, 
                                 do_pca = FALSE, return_object = TRUE,
                                 max.iter.harmony = 30, plot_convergence = FALSE) 

harmony_mat <- as.matrix(t(harmonyObj$Z_corr))
rownames(harmony_mat) <- rownames(pcaMat)
colnames(harmony_mat) <- paste(sep='','h',colnames(pcaMat))
saveRDS(harmony_mat,paste(sep='',outPrefix,'_Harmony',this_nPCs,'_mat.rds'))


print_time('Harmony Clustering')
harmony_cluster_res <- getClusters(harmony_mat,res_list)
colnames(harmony_cluster_res) <- paste(sep='','h',colnames(harmony_cluster_res))
saveRDS(harmony_cluster_res, paste(sep='',outPrefix,'_Harmony',this_nPCs,'_cluster.rds'))
harmony_cluster_df <- data.frame(apply(harmony_cluster_res,2,as.character),
                                         stringsAsFactors=FALSE)
harmony_cluster_cols <- colnames(harmony_cluster_df)


print_time('Harmony UMAP')
set.seed(seed)
harmony_umap <- umap(harmony_mat, n_neighbors = 30, metric = "cosine", min_dist = .3) 
saveRDS(harmony_umap,paste(sep='',outPrefix,'_Harmony',this_nPCs,'_UMAP.rds'))
harmony_umap_df <- data.frame("hUMAP1"=harmony_umap$layout[, 1],"hUMAP2"=harmony_umap$layout[, 2])


print_time('Harmony full DF')
harmony_full_df <- cbind(meta,harmony_mat,harmony_cluster_df,harmony_umap_df)
set.seed(seed)
harmony_full_df <- harmony_full_df[sample(nrow(harmony_full_df),nrow(harmony_full_df)),]
saveRDS(harmony_full_df,paste(sep='',outPrefix,'_Harmony',this_nPCs,'_fullDF.rds'))


print_time('Plot initial hPCs')
plot_dot(harmony_full_df,'hPC1','hPC2',plotTitle, 
         paste(sep="",outPrefix,'_Harmony',this_nPCs,'_hPC1-hPC2'), 
         plotMeta=use_meta, plotMetaCols=c(meta_col_toUse,harmony_cluster_cols), plotColors=color_vec)


print_time("Plot hUMAPs")
plot_dot(harmony_full_df,'hUMAP1','hUMAP2',
         paste(plotTitle,this_nPCs,'hPCs'), 
         paste(sep="",outPrefix,'_Harmony',this_nPCs,'_UMAP'), 
         plotMeta=use_meta, plotMetaCols=c(meta_col_toUse,harmony_cluster_cols), plotColors=color_vec)


print_time('Done.')

