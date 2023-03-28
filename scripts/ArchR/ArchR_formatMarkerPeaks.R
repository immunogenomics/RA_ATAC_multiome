print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(argparse)
})


parser <- ArgumentParser(description='extract marker peaks from ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project directory')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--FDR_cutoff', metavar='FDR_cutoff', default=0.1, type='double', help='FDR cutoff used for motif analysis; default=0.1')
parser$add_argument('--Log2FC_cutoff', metavar='Log2FC_cutoff', default=0.5, type='double', help='Log2FC cutoff used for motif analysis; default=0.5')
parser$add_argument('--save_intermediate_mat', action='store_true', help='if given, will also save the FDR and Log2FC matrices')
args <- parser$parse_args()

project_dir <- args$project_dir 
wk_dir <- args$wk_dir
extract_dir <- args$extract_dir
FDR_cutoff <- args$FDR_cutoff
Log2FC_cutoff <- args$Log2FC_cutoff
save_intermediate_mat <- args$save_intermediate_mat

print_time("Argument Checking")
if(!all(file.exists(project_dir))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')

markerPeak_sumExp_file <- paste(sep='',extract_dir,basename(project_dir),'_markerPeaks.rds')
if(!file.exists(markerPeak_sumExp_file)) stop("Marker Peak file does not exist!")

if(substrRight(extract_dir,1)!='/') extract_dir<-paste(sep="",extract_dir,'/')
dir.create(extract_dir,showWarnings=FALSE)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Set wd")
if(length(wk_dir)==0){
    wk_dir <- paste(sep='',dirname(project_dir),'/')
} else {
    if(substrRight(wk_dir,1)!='/') wk_dir<-paste(sep="",wk_dir,'/')
    if(!file.exists(wk_dir)) wk_dir <- paste(sep='',dirname(project_dir),'/')
}
cat(paste('New working dir:',wk_dir,'\n'))
setwd(wk_dir)


print_time("Extract FDR and Log2FC matrices")
markerPeak_sumExp <- readRDS(markerPeak_sumExp_file)

peakNames <- paste(sep='',as.character(rowData(markerPeak_sumExp)$seqnames),':',
                   as.numeric(rowData(markerPeak_sumExp)$start),'-',as.numeric(rowData(markerPeak_sumExp)$end))

FDR_mat <- assays(markerPeak_sumExp)$FDR
if(!identical(rownames(rowData(markerPeak_sumExp)),rownames(FDR_mat))) stop('FDR matrix rownames are unknown.')
rownames(FDR_mat) <- peakNames
if(save_intermediate_mat) {saveRDS(FDR_mat,paste(sep='',extract_dir,basename(project_dir),'_markerPeaks_FDR.rds'))}

Log2FC_mat <- assays(markerPeak_sumExp)$Log2FC
if(!identical(rownames(rowData(markerPeak_sumExp)),rownames(Log2FC_mat))) stop('Log2FC matrix rownames are unknown.')
rownames(Log2FC_mat) <- peakNames
if(save_intermediate_mat) {saveRDS(Log2FC_mat,paste(sep='',extract_dir,basename(project_dir),'_markerPeaks_Log2FC.rds'))}


print_time('Make membership matrix')
if(!identical(colnames(FDR_mat),colnames(Log2FC_mat))) stop('FDR and Log2FC matrices have different colnames.')

mp_df <- as.data.frame(matrix(data = FALSE, nrow = nrow(FDR_mat), ncol = ncol(FDR_mat)))
rownames(mp_df) <- rownames(FDR_mat)
colnames(mp_df) <- colnames(FDR_mat)

for(cs in colnames(FDR_mat)){
    FDR_peaks <- rownames(FDR_mat[which(FDR_mat[,cs]<=FDR_cutoff),])
    Log2FC_peaks <- rownames(Log2FC_mat[which(Log2FC_mat[,cs]>=Log2FC_cutoff),])
    cs_peaks <- FDR_peaks[which(FDR_peaks %in% Log2FC_peaks)]
    mp_df[cs_peaks,cs] <- TRUE
}
head(mp_df[which(mp_df[,1]==TRUE),],n=3)
print(colSums(mp_df))
table(rowSums(mp_df))

saveRDS(mp_df,paste(sep='',extract_dir,basename(project_dir),'_markerPeaks_logicalDF.rds'))


print_time('Done.')

