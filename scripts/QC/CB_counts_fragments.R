library(plyr)
library(argparse)

print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

parser <- ArgumentParser(description='counting cell barcodes from a fragment file')
parser$add_argument('frag_file', metavar='frag_file', help='fragment file with 5 columns: chr, start, stop, cell barcode, duplicate count')
parser$add_argument('sample_ID', metavar='sample_ID', help='sample ID for cell barcode concatenation')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

frag_file <- args$frag_file
sample_ID <- args$sample_ID
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(c(frag_file)))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)


print_time("Loading File")
frag_df <- read.table(frag_file, sep="\t", stringsAsFactors=FALSE)
colnames(frag_df) <- c('chr','start','stop','CB','dupCt')


print_time("Counting Cell Barcodes")
fragCt_byCB <- count(frag_df,'CB')


print_time("Adding Prefix to CB")
fragCt_byCB$CB <- paste(sep="#",sample_ID,fragCt_byCB$CB)
colnames(fragCt_byCB) <- c("CB","fragCt")


print_time("Saving File")
outFile <- paste(sep="", outPrefix, "_fragCt.txt")
write.table(fragCt_byCB,file=outFile,sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)


print_time("Done.")



