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
    library(plyr)
    library(dplyr)
    library(stringr)
    library(argparse)
})

aggMat <- function(thisPxC,thisMeta,thisCTcol,fun){

    theseCT <- sort(unique(thisMeta[,thisCTcol]))
    for(idx in 1:length(theseCT)){
        thisCT <- theseCT[idx]
        these_cells <- rownames(thisMeta[which(thisMeta[,thisCTcol]==thisCT),])
        this_subset <- thisPxC[,these_cells]

        if(length(these_cells)==1){
            thisAgg <- this_subset
        } else {
            thisAgg <- as.matrix(apply(this_subset,1,fun))
        }

        if(idx==1){
            catAgg <- thisAgg
        } else {
            catAgg <- cbind(catAgg,thisAgg)
        }
        colnames(catAgg)[ncol(catAgg)] <- thisCT
    }
    
    return(as(catAgg,'dgCMatrix'))

}

agg_merge <- function(half1, half2){
    return(rbind(half1,half2))
}       
    
agg_byHalf <- function(thisPxC,thisMeta,thisCTcol,fun,num_rows_max=133000){
    if(nrow(thisPxC)<=num_rows_max){
        return(aggMat(thisPxC,thisMeta,thisCTcol,fun))
    } else {
        mid_row_idx <- floor(nrow(thisPxC)/2)
        a <- agg_byHalf(thisPxC[1:mid_row_idx,],thisMeta,thisCTcol,fun,num_rows_max=num_rows_max)
        b <- agg_byHalf(thisPxC[(mid_row_idx+1):nrow(thisPxC),],thisMeta,thisCTcol,fun,num_rows_max=num_rows_max)
        agg_merge(a,b)
    }
}

checkMat <- function(thisMat){
    rs <- rowSums(thisMat)
    cat(paste(sep='','Number of all zero peaks: ',length(rs[rs==0]),'\n'))

    cs <- colSums(thisMat)
    cat(paste(sep='','Number of all zero cell (types): ',length(cs[cs==0]),'\n'))
}


parser <- ArgumentParser(description='aggregate normalized pxc by mean CT')
parser$add_argument('pxc_norm_file', metavar='pxc_norm_file', help='peaks x cells normalized file')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file (can be bigger, but needs all the cells)')
parser$add_argument('CT_col', metavar='CT_col', help='cell type column')
parser$add_argument('--simplifyCT', action='store_true', help='if given, change X-#: description to X#')
parser$add_argument('--nrow_forHalf', metavar='nrow_forHalf', type='integer', default=133000, help='max number of rows to give to aggregate matrix')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

pxc_norm_file <- args$pxc_norm_file
meta_file <- args$meta_file
CT_col <- args$CT_col
simplifyCT <- args$simplifyCT
nrow_forHalf <- args$nrow_forHalf
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(pxc_norm_file,meta_file))) stop("Input file(s) don't exist.")

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


print_time("load pxc")
pxc_norm <- readRDS(pxc_norm_file)
dim(pxc_norm)
max(pxc_norm)
pxc_norm[1:3,1:3]
checkMat(pxc_norm)

print_time("load meta")
meta <- readRDS(meta_file)
if(!(CT_col %in% colnames(meta))) stop('given cell type column not in meta_file')
if(!all(colnames(pxc_norm) %in% rownames(meta))) stop('not all pxc cells are in meta_file')
meta <- meta[colnames(pxc_norm),]
if(!identical(colnames(pxc_norm),rownames(meta))) stop('cells not identical between pxc and meta')

if(simplifyCT){
    if(all(grepl('^[A-Za-z]+-[0-9]+: ',unique(meta[,CT_col]),perl=TRUE))){
        newCT_col <- 'cellType'
        while(newCT_col %in% colnames(meta)){
            newCT_col <- paste(newCT_col,1)
        }
        meta[,newCT_col] <- str_replace(str_split_fixed(meta[,CT_col],":",2)[,1],"-","")
        print(table(meta[,c(CT_col,newCT_col)]))
        CT_col <- newCT_col
        
    } else {
        cat('keeping original cell type names since "X-#:"" pattern was not found for all cell types.\n')
    }
}


print_time('Aggregate cell types by mean')
pxCTmean <- agg_byHalf(pxc_norm,meta,CT_col,mean,num_rows_max=nrow_forHalf)
saveRDS(pxCTmean,paste(sep='',outPrefix,'.rds'))
checkMat(pxCTmean)

print_time('Done.')

