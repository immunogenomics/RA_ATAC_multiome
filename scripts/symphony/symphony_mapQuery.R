print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Seurat)
    library(harmony)
    library(presto)
    library(symphony)

    library(Rcpp)
    library(Matrix)
    library(MASS)
    library(matrixStats)
    library(rstatix)
    library(aricode)
    library(data.table)
    library(irlba)
    library(umap)
    library(uwot)


    library(magrittr)
    library(plyr)
    library(dplyr)
    library(stringr)

    library(reticulate)

    library(ggplot2)
    library(ggthemes)
    library(ggbeeswarm)
    library(ggrepel)
    library(ggrastr)
    library(gplots)
    library(gtools)
    library(grid)
    library(gridExtra)
    library(viridis)
    library(RColorBrewer)
    library(scales)
    library(patchwork)
    library(pheatmap)
    
    library(argparse)
})


#from Seurat
NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
    A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    A@x <- scaling_factor * A@x
    if (do_ftt) {
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else {
        A@x <- log(1 + A@x)
    }
    return(A)
}

#KW adapted from Seurat TF.IDF function - from Symphony paper
TF.IDF_multOut <- function(data,verbose = TRUE, refIDF = NULL){
    if (is.data.frame(x = data)) {
        data <- as.matrix(x = data)
    }
    if (!inherits(x = data, what = "dgCMatrix")) {
        data <- as(object = data, Class = "dgCMatrix")
    }
    if (verbose) {
        message("Performing TF-IDF normalization")
    }
    
    res <- list()
    
    npeaks <- colSums(x = data)
    tf <- t(x = t(x = data)/npeaks)
    
    #KW added option to input IDF vector
    if(is.null(refIDF)){
        if (verbose) {
            message("Using IDF calculated from current data")
        }
        idf <- ncol(x = data)/rowSums(x = data)
    } else {
        if (verbose) {
            message("Using provided IDF")
        }
        idf <- refIDF
    }
    
    norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
    norm.data[is.na(norm.data)] <- 0 #KW changed b/c errored: norm.data[which(x = is.na(x = norm.data))] <- 0
    
    #KW added more outputs
    res$TF <- tf
    res$IDF <- idf
    res$norm <- norm.data
    
    return(res)
}

#added more options
mapQuery3 <- function (exp_query, metadata_query, ref_obj, vars = NULL, verbose = TRUE, 
    do_normalize = FALSE, do_umap = TRUE, sigma = 0.1, do_log_normalize = TRUE) 
{
    if (do_normalize) {
        if (verbose) 
            message("Normalizing")
        exp_query = normalizeData(exp_query, 10000, "log")
    }
    
    #KW changed ref_obj$vargenes to ref_obj$varpeaks and means/sds to TF-IDF - from symphony paper
    if (verbose) 
        message("Normalizing and synchronizing query peak counts")
    idx_shared_peaks = which(names(ref_obj$varpeaks) %in% rownames(exp_query))
    shared_peaks_IDF = ref_obj$varpeaks[idx_shared_peaks]
    if (verbose) 
        message("Found ", length(shared_peaks_IDF), " reference variable peaks in query dataset")
    exp_query_scaled_res = TF.IDF_multOut(exp_query[names(shared_peaks_IDF),],refIDF = shared_peaks_IDF) #KW
    if(do_log_normalize) 
        exp_query_scaled_res$norm <- LogNormalize(exp_query_scaled_res$norm) #KW 
    exp_query_scaled = as.matrix(exp_query_scaled_res$norm) #KW
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$varpeaks), ncol = ncol(exp_query))
    exp_query_scaled_sync[idx_shared_peaks, ] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$varpeaks
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    
    
    if (verbose) 
        message("Project query cells using reference gene loadings")
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    
    
    if (verbose) 
        message("Clustering query cells to reference centroids")
    Z_pca_query_cos = cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, sigma)
    
    
    if (verbose) 
        message("Correcting query batch effects")
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[, vars] %>% as.data.frame()
        onehot = design %>% purrr::map(function(.x) {
            if (length(unique(.x)) == 1) {
                rep(1, length(.x))
            }
            else {
                stats::model.matrix(~0 + .x)
            }
        }) %>% purrr::reduce(cbind)
        Xq = cbind(1, intercept = onehot) %>% t()
    }
    else {
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), 
            sparse = TRUE)
    }
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq), 
        as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
    colnames(Z_pca_query) = row.names(metadata_query)
    rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
    colnames(Zq_corr) = row.names(metadata_query)
    rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
    
    
    umap_query = NULL
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        if (verbose) 
            message("UMAP")
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path, 
            verbose = FALSE)
        umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
        colnames(umap_query) = c("UMAP1", "UMAP2")
    }
    if (verbose) 
        message("All done!")
    return(list(exp = exp_query, meta_data = metadata_query,
                Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq, 
                umap = umap_query))
}
environment(mapQuery3) <- asNamespace("symphony")

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


plotBasic = function(umap_labels,             # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = 'Query',         # Plot title
                     color.by = 'cell_type',  # metadata column name for coloring
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.8, stroke = 0.5, shape = 16)
        if (!is.null(color.mapping) & all(unique(umap_labels[,color.by]) %in% names(color.mapping))) { p = p + scale_color_manual(values = color.mapping) }
    
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = FALSE)

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    return(p)
}


parser <- ArgumentParser(description='symphony mapQuery')
parser$add_argument('modality', metavar='modality', choices=c('RNA','ATAC'), help='modality: either RNA or ATAC')
parser$add_argument('ref_file', metavar='ref_file', help='symphony reference file')
parser$add_argument('--harmClust_file', metavar='harmClust_file', help='harmony cluster file for reference; not needed if cellType_col is in ref$meta_data')
parser$add_argument('cellType_col', metavar='cellType_col', help='cell type column in reference')
parser$add_argument('donor_col', metavar='donor_col', help='donor column in reference')
parser$add_argument('allCell_meta_file', metavar='allCell_meta_file', help='all cell symphony metadata file')
parser$add_argument('--onlyQuery_meta_flag', action='store_true', help='if given, allCell_meta_file should not be subsetted by ref_query column')
parser$add_argument('--meta_cellName_col', metavar='meta_cellName_col', help='if given, allCell_meta_file column containing cell names; if not given, rownames assumed.')
parser$add_argument('--query_cellType_col', metavar='query_cellType_col', help='if given, use as comparison query cell type. if not, default to unknown.')
parser$add_argument('--allCell_cellType_col', metavar='allCell_cellType_col', help='all cell symphony metadata cell type column')
parser$add_argument('--allCell_cellType_values', metavar='allCell_cellType_values', help='all cell symphony metadata cell type(s) to extract - comma-delimited if multiple cell types')
parser$add_argument('allCell_mat_file', metavar='allCell_mat_file', help='all cells features x cells file')
parser$add_argument('--uwot_file', metavar='uwot_file', help='symphony reference uwot file')
parser$add_argument('--logNorm_flag', action='store_true', help='if given, logNormalize the (ATAC) mapQuery TFxIDF')
parser$add_argument('--cellType_string', metavar='cellType_string',default='Cell Type', help='cell type string to use in figure titles')
parser$add_argument('--dataset_string', metavar='dataset_string', default='dataset', help='dataset string to use in figure titles')
parser$add_argument('--ref_string', metavar='ref_string', default='reference', help='reference string to use in figure titles')
parser$add_argument('--color_file', metavar='color_file', help='celltype colors - vector in RDS file; if not given, just uses default ggplot colors')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

modality <- args$modality
ref_file <- args$ref_file
harmClust_file <- args$harmClust_file ##optional; if cellType_col in reference$meta_data
cellType_col <- args$cellType_col
donor_col <- args$donor_col
allCell_meta_file <- args$allCell_meta_file
onlyQuery_meta_flag <- args$onlyQuery_meta_flag ##optional - if not given, subset of ref_query
meta_cellName_col <- args$meta_cellName_col ##optional - if not given, use rownames
query_cellType_col <- args$query_cellType_col ##optional - if not given, use unknown
allCell_cellType_col <- args$allCell_cellType_col ##optional - if not given, use everything
allCell_cellType_values <- args$allCell_cellType_values ##optional - if not given, use everything
allCell_mat_file <- args$allCell_mat_file
uwot_file <- args$uwot_file ##optional; can use the reference's uwot file already there
logNorm_flag <- args$logNorm_flag ##optional; if not given, FALSE (the original version); else TRUE (for the new log(TF.IDF) versions)
cellType_string <- args$cellType_string
dataset_string <- args$dataset_string
ref_string <- args$ref_string
color_file <- args$color_file ##optional
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(ref_file,allCell_meta_file,allCell_mat_file))) stop("Input files don't exist.")

opt_files <- c(harmClust_file,uwot_file,color_file)
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


set.seed(0)

print_time("Load Reference")

reference <- readRDS(ref_file)
ls(reference)
head(reference$meta_data)

if(length(harmClust_file)!=0){
    print_time("Harmony Clusters")
    harmony_clusters <- readRDS(harmClust_file)
    dim(harmony_clusters)
    head(harmony_clusters)
    if(!identical(rownames(reference$meta_data),rownames(harmony_clusters))) stop('Harmony cluster rownames are not identical to the reference metadata rownames')
    
    reference$meta_data <- cbind(reference$meta_data,harmony_clusters)
    head(reference$meta_data)
}

print_time("Verify reference cell type column")
if(!(cellType_col %in% colnames(reference$meta_data))) stop('Cell type column not in reference metadata')
reference$meta_data[,cellType_col] <- as.character(reference$meta_data[,cellType_col])
table(reference$meta_data[,cellType_col],useNA='ifany')

NA_cellType <- 'NL: Not Labeled'
while(NA_cellType %in% unique(reference$meta_data[,cellType_col])) NA_cellType <- paste(sep='',NA_cellType,'1')
if(any(is.na(reference$meta_data[,cellType_col]))){
    cat(paste(sep='','There are NA cell type values. They will be renamed to ',NA_cellType,'.\n'))
    
    reference$meta_data[is.na(reference$meta_data[,cellType_col]),cellType_col] <- NA_cellType
        
    table(reference$meta_data[,cellType_col],useNA='ifany')
}
reference$meta_data[,cellType_col] <- factor(reference$meta_data[,cellType_col],
                                             levels=mixedsort_new(unique(reference$meta_data[,cellType_col])))


print_time("Colors vector")
colors_vec <- hue_pal()(length(unique(reference$meta_data[,cellType_col])))
names(colors_vec) <- sort(unique(reference$meta_data[,cellType_col]))
if(length(color_file)!=0){
    cv <- readRDS(color_file)
    if(!any(unique(reference$meta_data[,cellType_col]) %in% names(cv))){
        cat("No color names overlap with reference cell types, so just using defaults\n")
    } else {
        colors_vec <- cv
        cat("Using colors from input file.\n")
    }
}
if(NA_cellType %in% unique(reference$meta_data[,cellType_col])) colors_vec[NA_cellType] <- 'grey15'

print_time("Adding cellname (i.e. rowname) column")
cellName_col <- 'cellName'
reference$meta_data[,cellName_col] <- rownames(reference$meta_data)

#we did not use these plots in the final manuscript
print_time("Reference UMAPs")
umap_labels = cbind(reference$meta_data[,c(cellName_col,donor_col,cellType_col)], reference$umap$embedding)
saveRDS(umap_labels,paste(sep='',outPrefix,'_ref_UMAP_labels.rds'))

g <- ggplot() +
  geom_point(
    data = umap_labels[sample(nrow(umap_labels)),],
    mapping = aes_string(x = "UMAP1", y = "UMAP2", color = cellType_col),size=0.5,alpha=0.5) +
  theme_bw(base_size = 15) +
  ggtitle(paste(ref_string,cellType_string,'reference')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=colors_vec)
ggsave(paste(sep='',outPrefix,'_ref_UMAP-cellType.png'),plot=g,units='in',height=7,width=12)


print_time("Query metadata")
all_cell_symphony_meta <- readRDS(allCell_meta_file)
head(all_cell_symphony_meta)

if(!(donor_col %in% colnames(all_cell_symphony_meta))) stop(paste('all_cell_symphony_meta requires the donor col:',donor_col))

if(length(meta_cellName_col)==0){
    all_cell_symphony_meta[,cellName_col] <- rownames(all_cell_symphony_meta)
} else {
    if(!(meta_cellName_col %in% colnames(all_cell_symphony_meta))) stop('if given, meta_cellName_col should exist in all_cell_symphony_meta')
    all_cell_symphony_meta[,cellName_col] <- all_cell_symphony_meta[,meta_cellName_col]
}

if(length(query_cellType_col)==0){
    all_cell_symphony_meta[,cellType_col] <- 'unknown'
} else {
    if(!(query_cellType_col %in% colnames(all_cell_symphony_meta))) stop('if given, query_cellType_col should exist in all_cell_symphony_meta')
    all_cell_symphony_meta[,cellType_col] <- all_cell_symphony_meta[,query_cellType_col]
}

if(length(allCell_cellType_col)!=0 & length(allCell_cellType_values)!=0){
    if(!(allCell_cellType_col %in% colnames(all_cell_symphony_meta))) stop('subsetting cell type column in allCell metadata file not found.')
    if(!(all(unlist(strsplit(allCell_cellType_values,",")) %in% unique(all_cell_symphony_meta[,allCell_cellType_col])))) stop('not all given cell types found in allCell metadata file.')
    
    cat('Subsetting query\n')    
    CT_toGet <- unlist(strsplit(allCell_cellType_values,",")) #checked in if statement
    all_cell_symphony_meta <- all_cell_symphony_meta[which(all_cell_symphony_meta[,allCell_cellType_col] %in% CT_toGet),]

}

predCT_col <- 'cell_type_pred_knn'
refQuery_col <- 'ref_query'

if(onlyQuery_meta_flag){
    query_metadata <- all_cell_symphony_meta[,c(cellName_col,donor_col,cellType_col)]
} else {
    if(!refQuery_col %in% colnames(all_cell_symphony_meta)) stop('not all required columns in all_cell_symphony_meta')
    
    paste("query nrow: ",nrow(all_cell_symphony_meta[which(all_cell_symphony_meta[,refQuery_col]=='query'),]))
    
    query_metadata <- all_cell_symphony_meta[which(all_cell_symphony_meta[,refQuery_col]=='query'),
                                             c(cellName_col,donor_col,cellName_col)]
}


print_time("Query features x cells")
allCell_mat <- readRDS(allCell_mat_file)

paste("allCells features x cells matrix dim:",paste(collapse=' ',dim(allCell_mat)))
if(!all(rownames(query_metadata) %in% colnames(allCell_mat))){
    if(length(meta_cellName_col)!=0){
        if(all(query_metadata[,cellName_col] %in% colnames(allCell_mat))){
            cat(paste(sep='','Not all query metadata rownames were in the features x cells colnames, but the meta_cellName_col (',meta_cellName_col,') were, so the query metadata rownames were replaced with the cellName column!\n'))
            
            cat('Affected cells:')
            print(table(rownames(query_metadata)!=query_metadata[,cellName_col]))
            
            rownames(query_metadata) <- query_metadata[,cellName_col]
        }
    }
}
if(!all(rownames(query_metadata) %in% colnames(allCell_mat))) stop("Not all rownames in query metadata in allCells features x cells matrix.")

query_mat <- allCell_mat[,rownames(query_metadata)]
paste("query features x cell matrix",paste(collapse=' ',dim(query_mat)))
saveRDS(query_mat,paste(sep='',outPrefix,'_query_featMat.rds'))


print_time("Save query metadata")
dim(query_metadata)
head(query_metadata)
saveRDS(query_metadata,paste(sep='',outPrefix,'_query_meta_initial.rds'))


if(modality=='RNA'){
    print_time('Normalize query matrix')
    query_exprs_norm <- query_mat %>% NormalizeDataSeurat()
    saveRDS(query_exprs_norm,paste(sep='',outPrefix,'_query_featMat_norm.rds'))   
}


if(length(uwot_file)!=0){
    print_time("Change reference uwot path")
    print(paste(sep='','original reference uwot path: ',reference$save_uwot_path))
    
    reference$save_uwot_path <- uwot_file
    print(paste(sep='','new reference uwot path: ',reference$save_uwot_path))
}


print_time("Map Query")
if(modality=='ATAC'){
    query = mapQuery3(query_mat,                  # query peak counts (peaks x cells)
                      query_metadata,             # query metadata (cells x attributes)
                      reference,                  # Symphony reference object
                      vars = c(donor_col),        # also integrate over query batches
                      do_normalize = FALSE,       # perform log(CP10k) normalization on query - no, normalize with TF-IDF
                      do_log_normalize = logNorm_flag,   #added after symphony paper
                      do_umap = TRUE)             # project query cells into reference UMAP 
} else if(modality=='RNA'){
    query = mapQuery(query_exprs_norm,     # query gene expression (genes x cells)
                    query_metadata,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    vars = c(donor_col),   # also integrate over query batches
                    do_normalize = FALSE,  # perform log(CP10k) normalization on query
                    do_umap = TRUE)        # project query cells into reference UMAP
} else {
    stop('Should not have gotten past the two modality choices for mapping a query')
}

print_time("Predict query cell types using k-NN")
query = knnPredict(query, reference, reference$meta_data[,cellType_col], k = 5,confidence=TRUE)
#saved at the end

head(query$meta_data)
saveRDS(query$meta_data,paste(sep='',outPrefix,'_query_meta_postPredict.rds'))


#we did not use these plots in the final manuscript
print_time("Query freq")
CT_ct_df <- as.data.frame(table(query$meta_data[,predCT_col]))
colnames(CT_ct_df) <- c(predCT_col,'Freq')
head(CT_ct_df)

head(CT_ct_df[order(CT_ct_df$Freq),])

paste("nrow(query$meta_data):",nrow(query$meta_data))
paste("sum(CT_ct_df$Freq):",sum(CT_ct_df$Freq))

CT_ct_df$perc <- CT_ct_df$Freq/nrow(query$meta_data)
saveRDS(CT_ct_df,paste(sep='',outPrefix,'_query_cellType_freq.rds'))

g <- ggplot(CT_ct_df,aes_string(x=predCT_col,y='Freq',fill=predCT_col)) + 
        geom_bar(position="dodge", stat="identity") +
        labs(x="Predicted Cell Type",fill='Predicted\nCell Type',y='Number of query cells') +
        theme_grey(base_size=15) +
        ggtitle(paste(dataset_string,cellType_string,'\nPredicted Cell Type Counts')) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual(values=colors_vec)
ggsave(paste(sep='',outPrefix,'_query_predictedCTct.png'),plot=g,units='in',height=6, width=6)


#we did not use these plots in the final manuscript
print_time("mapping quality")
predProb_col <- paste(sep='',predCT_col,'_prob')
g <- ggplot(query$meta_data, aes_string(x=predCT_col,y=predProb_col,color=predCT_col)) +
        geom_boxplot() + theme_bw(base_size=15) + theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #expand_limits(y=0) +
        labs(x='Predicted Cell Type',y='Proportion of 5 NN\nwith winning vote') +
        ggtitle(paste(dataset_string,cellType_string,'Mapping Quality')) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=colors_vec)
ggsave(paste(sep='',outPrefix,'_query_mappingQuality.png'),plot=g,units='in',height=6, width=6)

table(query$meta_data[,c(predCT_col,predProb_col)])

print_time("Merge reference and query metadata tables")
#Sync the column names for both data frames
reference$meta_data[,predCT_col] = NA
reference$meta_data[,predProb_col] = NA
reference$meta_data[,refQuery_col] = 'reference'
query$meta_data[,refQuery_col] = 'query'

#Add the UMAP coordinates to the metadata
meta_data_combined = rbind(query$meta_data[, c(cellName_col,donor_col,cellType_col,#'cluster_number',
                                               predCT_col,predProb_col, refQuery_col)],
                           reference$meta_data[, c(cellName_col,donor_col,cellType_col,#'cluster_number',
                                                   predCT_col,predProb_col, refQuery_col)])
umap_combined = rbind(query$umap, reference$umap$embedding)
umap_combined_labels = cbind(meta_data_combined, umap_combined)
#saved at the end

#we did not use these plots in the final manuscript
print_time("ref/query UMAPs")
g <- ggplot() +
  geom_point(
    data = umap_combined_labels[rev(order(umap_combined_labels[,refQuery_col])),],
    mapping = aes_string(x = "UMAP1", y = "UMAP2", color = refQuery_col),
    size=0.5,alpha=0.5
  ) +
  theme_bw(base_size = 15) +
  ggtitle(paste(dataset_string,cellType_string,'reference/query')) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste(sep='',outPrefix,'_refQuery_UMAP.png'),plot=g,units='in',height=7,width=10)

#we did not use these plots in the final manuscript
print_time("Query UMAPs")
toPlot <- umap_combined_labels[which(umap_combined_labels[,refQuery_col]=='query'),]
g <- ggplot(toPlot,aes_string(x='UMAP1',y='UMAP2',color=predCT_col)) + geom_point(size=0.5,alpha=0.5) +
        theme_bw(base_size=15) +
        ggtitle(paste(dataset_string,cellType_string,'query')) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=colors_vec)
ggsave(paste(sep='',outPrefix,'_query_UMAP-cellType.png'),plot=g,units='in',height=7,width=12)


#we did not use these plots in the final manuscript
print_time("ref query facet UMAP")
umap_combined_labels[,cellType_col] <- factor(umap_combined_labels[,cellType_col],
                                              levels=mixedsort_new(unique(as.character(umap_combined_labels[,cellType_col]))))

paste("nrow(umap_combined_labels):",nrow(umap_combined_labels))
comboCT = c(as.character(query$meta_data[,predCT_col]), as.character(reference$meta_data[,cellType_col]))
paste("length(comboCT):",length(comboCT))
head(comboCT)

temp_umap_combined_labels <- cbind(umap_combined_labels,comboCT)
temp_umap_combined_labels$comboCT <- factor(temp_umap_combined_labels$comboCT,
                                            levels=mixedsort_new(unique(as.character(temp_umap_combined_labels$comboCT))))
table(temp_umap_combined_labels$comboCT)
rbind(head(temp_umap_combined_labels,n=3),tail(temp_umap_combined_labels,n=3))
saveRDS(temp_umap_combined_labels,paste(sep='',outPrefix,'_refQuery_meta_wComboCT.rds'))

g <- plotBasic(temp_umap_combined_labels,
               title = paste(dataset_string,cellType_string,'reference and query cells'),
               color.by = 'comboCT', facet.by = refQuery_col, color.mapping = colors_vec) + 
        theme_bw(base_size=20)
ggsave(paste(sep='',outPrefix,'_refQuery_UMAP_facet_query-colors.png'),plot=g,units='in',height=6,width=14)


print_time("Saving")
rbind(head(umap_combined_labels,n=3),tail(umap_combined_labels,n=3))
saveRDS(umap_combined_labels,paste(sep='',outPrefix,'_refQuery_metadata.rds'))

saveRDS(reference,paste(sep='',outPrefix,'_reference_postScript.rds'))
saveRDS(query,paste(sep='',outPrefix,'_query_postScript.rds'))

print_time('Done.')
