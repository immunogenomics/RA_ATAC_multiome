print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(class)
    library(RANN)
    library(Matrix)
    library(stringr)
    library(plyr)
    library(harmony)
    library(symphony)
    library(ggplot2)
    library(argparse)
})


revKNN <- function(thisRef, thisQuery, toPred_col, thisMeta=NULL, toQuery_cells=NULL, k_numNN=5, thisSeed=0, ref_cell_col='cell',
                   resClust_col='pred_cluster_name', resProb_col='pred_cluster_prob', ref_umap1_col='UMAP1',ref_umap2_col='UMAP2',
                   plotTitle='',plotPrefix='',plotColors=NA){
    
    if(!all(c('meta_data','Z') %in% ls(thisQuery))) stop('query object must have meta_data and Z slots')
    if(!identical(colnames(thisQuery$Z),rownames(thisQuery$meta_data))) stop('query structure off: Z cells do not match metadata cells')
    if(!all(c('Z_corr','meta_data','umap') %in% ls(thisRef))) stop('ref object must have Z_corr, meta_data, and umap slots')
    if(!('embedding' %in% ls(thisRef$umap))) stop('ref object must have thisRef$umap$embedding slot')
    if(!all(c(ref_umap1_col,ref_umap2_col) %in% colnames(thisRef$umap$embedding))) stop('UMAP columns must be in thisRef$umap$embedding')
    
    if(!is.null(toQuery_cells)){
        if(!all(toQuery_cells %in% rownames(thisQuery$meta_data))){
            stop('if toQuery_cells given, they all must exist in the query metadata.')
        } else {
            toPred_cells <- toQuery_cells
        }
    } else {
        toPred_cells <- rownames(thisQuery$meta_data)
    }
    
    #priority given to thisMeta
    if(!is.null(thisMeta)){
        if(!(toPred_col %in% colnames(thisMeta))){
            stop('if thisMeta is given, it is expected to find toPred_col column.')
        } else {
            if(!all(rownames(thisQuery$meta_data) %in% rownames(thisMeta))){
                stop('all query metadata rownames must be in meta rownames')
            } else{
                thisQuery$meta_data[,toPred_col] <- thisMeta[rownames(thisQuery$meta_data),toPred_col]
            }
        }
    } else {
        if(!(toPred_col %in% colnames(thisQuery$meta_data))) stop('toPred_col not in query metadata')
    }
    
    refMeta <- cbind(thisRef$meta_data,thisRef$umap$embedding)
    if(!is.null(ref_cell_col)){
        if(!(ref_cell_col %in% colnames(thisRef$meta_data))) {
            stop('if reference cell name column given, it should exist in thisRef$meta_data')
        } else {
            rownames(refMeta) <- thisRef$meta_data[,ref_cell_col]
        }
    } else {
        print('No reference cell name column given, so not changing.\n')
    }
    
    set.seed(thisSeed)
    knn_pred = class::knn(t(thisQuery$Z[,toPred_cells]),t(thisRef$Z_corr),
                          thisQuery$meta_data[toPred_cells,toPred_col], k = k_numNN, prob = TRUE)
    knn_prob = attributes(knn_pred)$prob
    
    refMeta[,resClust_col] <- as.character(knn_pred)
    refMeta[,resProb_col] <- knn_prob

    g <- ggplot(refMeta,aes_string(x=ref_umap1_col,y=ref_umap2_col,color=resClust_col)) + geom_point(shape='.') +
            theme_bw(base_size=15) +
            labs(color='classified\ncluster',title=plotTitle) +
            theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank())
    if(all(unique(refMeta[,resClust_col]) %in% names(plotColors))){
        g <- g + scale_color_manual(values=plotColors)
    }
    ggsave(file=paste(sep='',plotPrefix,'_oriUMAP_predClusters.png'),plot=g,units='in',height=5,width=8)

    g <- ggplot(refMeta,aes_string(y=resClust_col,x=resProb_col,color=resClust_col)) +
            geom_boxplot() + theme_classic(base_size=15) + 
            theme(legend.position="none") +
            labs(y='classified cluster',x='mapping score',title=plotTitle)
    if(all(unique(refMeta[,resClust_col]) %in% names(plotColors))){
        g <- g + scale_color_manual(values=plotColors)
    }
    ggsave(file=paste(sep='',plotPrefix,'_mappingScore.png'),plot=g,units='in',height=8,width=5)
    
    return(refMeta)
}


avgUMAP <- function(thisRef, refMeta, 
                    thisQuery, query_umap_df, toQuery_cells=NULL, query_umap1_col='UMAP1', query_umap2_col='UMAP2',
                    k_numNN=5, thisSeed=0,
                    ref_cell_col='cell',ref_umap1_col='pred_UMAP1',ref_umap2_col='pred_UMAP2',resClust_col='pred_cluster_name',
                    plotTitle='',plotPrefix='',plotColors=NA){
    
    if(!('Z' %in% ls(thisQuery))) stop('query object must have Z')
    if(!all(colnames(thisQuery$Z) %in% rownames(query_umap_df))) stop('cellnames not the same for thisQuery$Z and query_umap_df')
    if(!all(c('Z_corr','meta_data') %in% ls(thisRef))) stop('ref object must have Z_corr and meta_data slots')

    if(!is.null(toQuery_cells)){
        if(!all(toQuery_cells %in% rownames(thisQuery$meta_data))){
            stop('if toQuery_cells given, they all must exist in the query metadata.')
        } else {
            toPred_cells <- toQuery_cells
        }
    } else {
        toPred_cells <- rownames(thisQuery$meta_data)
    }
    
    set.seed(thisSeed)
    all_nn <- nn2(t(thisQuery$Z[,toPred_cells]), query=t(thisRef$Z_corr), k = k_numNN, eps = 0)

    nn_mat <- matrix(0, nrow = ncol(thisRef$Z_corr), ncol = ncol(thisQuery$Z[,toPred_cells]))
    for(idx in 1:nrow(all_nn$nn.idx)){
        nn_mat[idx,all_nn$nn.idx[idx,]] <- 1
    }
    ll <- rowSums(nn_mat)
    if(!all(ll==k_numNN)) stop('Number of NN not right in nn_mat.')
    
    query_umap_mat <- as.matrix(query_umap_df[colnames(thisQuery$Z[,toPred_cells]),c(query_umap1_col,query_umap2_col)])
    ref_umap_mat <- nn_mat %*% query_umap_mat
    ref_umap_mat <- ref_umap_mat/k_numNN

    if(!is.null(ref_cell_col)){
        if(!(ref_cell_col %in% colnames(thisRef$meta_data))) {
            stop('if reference cell name column given, it should exist in thisRef$meta_data')
        } else {
            rownames(ref_umap_mat) <- thisRef$meta_data[,ref_cell_col]
        }
    } else {
        rownames(ref_umap_mat) <- rownames(thisRef$meta_data)
    }
    if(!identical(rownames(refMeta),rownames(ref_umap_mat))) stop('rownames different between ref_umap_mat and refMeta')
    
    refMeta[,ref_umap1_col] <- ref_umap_mat[,1]
    refMeta[,ref_umap2_col] <- ref_umap_mat[,2]

    g <- ggplot(refMeta, aes_string(x=ref_umap1_col,y=ref_umap2_col,color=resClust_col)) + 
            geom_point(shape='.') + theme_bw(base_size=15) + 
            labs(color='classified\ncluster',title=plotTitle) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank())
    if(all(unique(refMeta[,resClust_col]) %in% names(plotColors))){
        g <- g + scale_color_manual(values=plotColors)
    }
    ggsave(file=paste(sep='',plotPrefix,'_predUMAP_predClusters.png'),plot=g,units='in',height=5,width=8)
    
    return(refMeta)
}


parser <- ArgumentParser(description='symphony reverse KNN prediction and average UMAP coordinates')
parser$add_argument('ref_file', metavar='ref_file', help='symphony reference file')
parser$add_argument('--ref_cellCol', metavar='ref_cellCol', help='symphony reference object cell column - or rownames')
parser$add_argument('query_file', metavar='query_file', help='symphony query file')
parser$add_argument('query_meta_file', metavar='query_meta_file', help='query metadata file')
parser$add_argument('query_predCol', metavar='query_predCol', help='query prediction column')
parser$add_argument('--query_umap1Col', metavar='query_umap1Col', default='UMAP1', help='query UMAP1 column; default=UMAP1')
parser$add_argument('--query_umap2Col', metavar='query_umap2Col', default='UMAP2', help='query UMAP2 column; default=UMAP2')
parser$add_argument('--revKNN_n', metavar='revKNN_n', default=5, type='integer', help='Number of nearest neighbors to use for reverse knnPredict; default=5')
parser$add_argument('--avgUMAP_n', metavar='avgUMAP_n', default=5, type='integer', help='Number of nearest neighbors to use for average UMAP positions; default=5')
parser$add_argument('--color_file', metavar='color_file', help='celltype colors - vector in RDS file; if not given, just uses default ggplot colors')
parser$add_argument('--seed', metavar='seed', default=0, type='double', help='randomization seed; default 0')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

ref_file <- args$ref_file
ref_cellCol <- args$ref_cellCol
query_file <- args$query_file
query_meta_file <- args$query_meta_file
query_predCol <- args$query_predCol
query_umap1Col <- args$query_umap1Col
query_umap2Col <- args$query_umap2Col
revKNN_n <- args$revKNN_n
avgUMAP_n <- args$avgUMAP_n
color_file <- args$color_file
seed <- args$seed
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(ref_file,query_file,query_meta_file))) stop("Input files don't exist.")

opt_files <- c(color_file)
opt_file_check <- function(fl){
    if(length(fl)!=0){
        if(!file.exists(fl)) stop(paste("if given,",fl,"file must exist."))
      }
}
invisible(lapply(opt_files,opt_file_check))

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


print_time('Load files')
ref <- readRDS(ref_file)
query <- readRDS(query_file)
query_meta <- readRDS(query_meta_file)

if(!all(rownames(query$meta_data) %in% rownames(query_meta))) stop('all query$meta_data cells must be in query_meta')
query_meta <- query_meta[rownames(query$meta_data),]
if(!identical(rownames(query_meta),rownames(query$meta_data))) stop('cells not the same between query_meta and query$meta_data')
if(!(query_predCol %in% unique(c(colnames(query_meta),colnames(query$meta_data))))) stop('query_predCol must be in columns of either query_meta and query$meta_data')

color_vec <- NA
if(!is.null(color_file)){
    if(file.exists(color_file)){
        print_time('colors file')
        color_vec <- readRDS(color_file)
    }
}


print_time('reverse knnPredict')
ref_meta <- revKNN(ref, query, query_predCol, thisMeta=query_meta, k_numNN=revKNN_n, thisSeed=seed, ref_cell_col=ref_cellCol,
                   plotPrefix=outPrefix,plotColors=color_vec)


print_time('average UMAP locations')
ref_meta <- avgUMAP(ref, ref_meta, query, query_meta, query_umap1_col=query_umap1Col, query_umap2_col=query_umap2Col,
                    k_numNN=avgUMAP_n, thisSeed=seed,ref_cell_col=ref_cellCol,plotPrefix=outPrefix,plotColors=color_vec)


print_time('Save')
saveRDS(ref_meta,paste(sep='',outPrefix,'.rds'))


print_time('Done.')
