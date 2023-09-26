suppressMessages({
    library(GenomicRanges)
    library(Matrix)
    library(Matrix.utils)
    library(gtools)
    library(argparse)
})

gRange_str <- function(gr) {
    return(paste(sep="",seqnames(gr),":",start(gr),"-",end(gr)))
}

substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}


parser <- ArgumentParser(description='Peaks by cells matrix with number of overlapping fragments')
parser$add_argument('peaks_file', metavar='peaks_file', help='peaks file with max 10 columns')
parser$add_argument('fragment_file', metavar='fragment_file', help='fragment file with 5 columns: chr, start, stop, cell barcode, duplicate count')
parser$add_argument('--sample_ID', metavar='sample_ID', help='sample ID for cell barcode concatenation')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

peaks_file <- args$peaks_file
fragment_file <- args$fragment_file
sample_ID <- args$sample_ID
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(c(peaks_file,fragment_file)))) stop("Input file(s) don't exist.")

if(is.null(sample_ID)) sample_ID <- prefix

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}

print_time("Starting fragment GR")
fragment_df <- read.table(fragment_file,sep="\t",stringsAsFactors=FALSE)
colnames(fragment_df)<-c('chr','start','end','cell','dupCt')
fragment_df$donor <- sample_ID
fragment_df$CBD <- paste(sep="#",fragment_df$donor,fragment_df$cell)
head(fragment_df)
fragment_GR <- GRanges(fragment_df)
fragment_GR
length(unique(fragment_GR$cell))


print_time("Staring peaks GR")
peaks_df <- read.table(peaks_file,sep="\t",stringsAsFactors=FALSE)
full_colnames <- c('chr','start','end','name','score','strand','sigVal','pval','qval','summit')
colnames(peaks_df)<-full_colnames[1:ncol(peaks_df)]
head(peaks_df)
peaks_GR <- GRanges(peaks_df)
peaks_GR


print_time("Starting initial findOverlaps for reordering")
#This version of findOverlaps seems to miss the first few overlaps/hits, so I do it twice, reordering in the middle such that the non-hits come first.
fragment_GR_forReorder <- fragment_GR[seqnames(fragment_GR) == as.character(seqnames(fragment_GR)[1])]
peaks_GR_forReorder <- peaks_GR[seqnames(peaks_GR) == as.character(seqnames(fragment_GR)[1])]
hits_init <- findOverlaps(fragment_GR_forReorder,peaks_GR_forReorder)
peaks_GR <- c(peaks_GR[peaks_GR != peaks_GR_forReorder[subjectHits(hits_init)[1]]],peaks_GR[peaks_GR == peaks_GR_forReorder[subjectHits(hits_init)[1]]])
peaks_GR
fragment_GR <- c(fragment_GR[fragment_GR != fragment_GR_forReorder[queryHits(hits_init)[1]]],fragment_GR[fragment_GR == fragment_GR_forReorder[queryHits(hits_init)[1]]])
fragment_GR
rm(fragment_GR_forReorder,peaks_GR_forReorder)


print_time("Starting findOverlaps")
hits <- findOverlaps(fragment_GR,peaks_GR)
hits_SM <- sparseMatrix(i=queryHits(hits),j=subjectHits(hits),dimnames = list(gRange_str(fragment_GR),gRange_str(peaks_GR)))
hits_SM[1:3,1:3]
dim(hits_SM)


print_time("Starting peaks by cells")
rownames(hits_SM)<-fragment_GR$CBD
cells_by_peaks <- aggregate.Matrix(hits_SM, row.names(hits_SM))
peaks_by_cells <- t(cells_by_peaks)
rm(cells_by_peaks)
peaks_by_cells <- peaks_by_cells[mixedsort(rownames(peaks_by_cells)),]
peaks_by_cells[1:3,1:3]
dim(peaks_by_cells)
saveRDS(peaks_by_cells,paste(sep="",outDir,prefix,"_peaksXcells.RDS"))


print_time("Starting all-zero cell stats")
cs <- colSums(peaks_by_cells)
peaks_by_cells_non0 <- peaks_by_cells[,cs!=0]
dim(peaks_by_cells_non0)

print_time("Starting peaks by donors")
rs <- rowSums(peaks_by_cells)
peaksXdonor <- data.frame(V1=ifelse(rs==0,0,1))
peaksXdonor$V1 <- as.integer(peaksXdonor$V1)
colnames(peaksXdonor) <- prefix
table(peaksXdonor)
saveRDS(peaksXdonor,paste(sep="",outDir,prefix,"_peaksXdonor.RDS"))


