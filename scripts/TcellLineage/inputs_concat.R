print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(gtools)
    library(argparse)
})

parse_list <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }

    return(list_toUse)
}



parser <- ArgumentParser(description='concat iterPeaks jobs and optional MHTC')
parser$add_argument('inDir', metavar='inDir', help='Input directory')
parser$add_argument('suffix', metavar='suffix', help='file suffix')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
parser$add_argument('--pval_cols', metavar='pval_cols', default='pval', help='p-value columns; default: pval')
parser$add_argument('--padj', metavar='padj', default='BH',choices=c('holm','hochberg','hommel','bonferroni','BH','BY','fdr'), help='if given, do pvalue-adjust; options: holm, hochberg, hommel, bonferroni, BH, BY, fdr; default=BH')
args <- parser$parse_args()

inDir <- args$inDir
suffix <- args$suffix
outDir <- args$outDir
prefix <- args$prefix
pval_cols <- args$pval_cols
padj <- args$padj


print_time("Argument Checking")
if(!all(file.exists(inDir))) stop("Input directory don't exist.")

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


print_time("Load Files")
input_files <- mixedsort(list.files(inDir,pattern=suffix,full.names=TRUE))
length(input_files)

for(idx in 1:length(input_files)){
    input_file <- input_files[idx]
    
    if(idx==1){
        all_df <- readRDS(input_file)
    } else {
        all_df <- rbind(all_df,readRDS(input_file))
    }
    print(nrow(all_df))
}
head(all_df)


print_time('pvalue columns')
pval_cols <- parse_list(pval_cols)
cat(paste('Given number of columns:',length(pval_cols),'\n'))
pval_cols <- pval_cols[which(pval_cols %in% colnames(all_df))]
cat(paste('Number of columns in df:',length(pval_cols),'\n'))
pval_cols


if(length(padj)!=0){
    print_time("MHTC")
    for(pval_col in pval_cols){
        if(pval_col=='pval'){
            padj_col <- 'padj'
            log10padj_col <- 'log10padj'
        } else {
            padj_col <- paste(sep='',pval_col,'_padj')
            log10padj_col <- paste(sep='',pval_col,'_log10padj')
        }
        
        all_df[,padj_col] <- p.adjust(all_df[,pval_col],method=padj)
        all_df[,log10padj_col] <- -log10(all_df[,padj_col])
    }
}


print_time('Save')
saveRDS(all_df,paste(sep='',outPrefix,'.rds'))


print_time('Done.')

