print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(presto)
    library(argparse)
})


parser <- ArgumentParser(description='differential genes')
parser$add_argument('data_file', metavar='data_file', help='Matrix file: normalized genes x cells')
parser$add_argument('meta_file', metavar='meta_file', help='Cell metadata file')
parser$add_argument('cellType_col', metavar = 'cellType_col', help='cell type column; default: cellType')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

data_file <- args$data_file
meta_file <- args$meta_file
cellType_col <- args$cellType_col
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(data_file,meta_file))) stop("Input files don't exist.")

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


print_time("Load files")
data <- readRDS(data_file)
meta <- readRDS(meta_file)
if(!(cellType_col %in% colnames(meta))) stop("Cell Type column not in meta file")
if(!all(colnames(data) %in% rownames(meta))) stop("Not all gxc matrix cells in meta file")
if(!identical(rownames(meta),colnames(data))) meta <- meta[colnames(data),]


print_time("Differential Genes")
diffGenes_df <- wilcoxauc(data, meta[,cellType_col])


print_time("Save")
saveRDS(diffGenes_df,paste(sep='',outPrefix,'_diffGenes.rds'))


print_time('Done.')

