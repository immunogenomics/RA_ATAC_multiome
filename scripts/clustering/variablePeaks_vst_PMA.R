print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(Seurat)
    library(symphony)
    library(argparse)
})


parser <- ArgumentParser(description='minially accessible, norm, and variable peaks using symphony vargenes_vst')
parser$add_argument('pxc_file', metavar='pxc_file', help='peaks x cells file')
parser$add_argument('--peak_acc_crit', metavar='peak_acc_crit', default=0.005, type='double', help='peak accessible in at least X cells (if integer) or X proportion of cells (if decimal); default 0.005')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file for sample')
parser$add_argument('sample_col', metavar='sample_col', help='sample column in metadata file')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='double', help='randomization seed; default 1234567890')
parser$add_argument('--numPeaks', metavar='numPeaks', default=500, type='double', help='number of variables peaks to ask for (to start); default 500')
parser$add_argument('--step_crit', metavar='step_crit', default=50, type='double', help='step criterion; if integer - number of peaks to decrease; if decimal - percentage of peaks to decrease; default 50')
parser$add_argument('--end_crit', metavar='end_crit', default=0.8, type='double', help='end criterion; variable peaks less than X proportion of cells; default 0.8')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

pxc_file <- args$pxc_file
peak_acc_crit <- args$peak_acc_crit
meta_file <- args$meta_file
sample_col <- args$sample_col
seed <- args$seed
numPeaks <- args$numPeaks
step_crit <- args$step_crit
end_crit <- args$end_crit
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(pxc_file,meta_file))) stop("Input file(s) don't exist.")

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


print_time("Load matrix")
pxc <- readRDS(pxc_file)
dim(pxc)
pxc[1:3,1:3]
ll <- max(pxc)
if(ll!=1){
    cat(paste('Max peak counts',ll,'Binarizing.\n'))
    suppressMessages(library(Signac))
    pxc <- BinarizeCounts(pxc)
    saveRDS(pxc,paste(sep='',outPrefix,'_binary.rds'))
}
nCells_tot <- ncol(pxc)


print_time("Peak Accessibility Cutoff")
pxc_rs <- rowSums(pxc)
length(pxc_rs)
head(pxc_rs)
summary(pxc_rs)

cat(paste('\nPeaks accessibility cutoff:',peak_acc_crit,'\n'))
if(peak_acc_crit<1) nCells_cutoff <- floor(nCells_tot * peak_acc_crit) else nCells_cutoff <- peak_acc_crit
cat(paste('Peaks must be accessible in at least',nCells_cutoff,'cells\n'))

peaks_keeping <- names(pxc_rs[pxc_rs >= nCells_cutoff])
cat(paste('Peaks kept:',length(peaks_keeping),'\n\n'))
saveRDS(peaks_keeping,paste(sep='',outPrefix,'_accPeaks-',nCells_cutoff,'.rds'))

pxc <- pxc[peaks_keeping,]
dim(pxc)
pxc[1:3,1:3]
nPeaks_tot <- nrow(pxc)


print_time("Load meta")
meta <- readRDS(meta_file)
dim(meta)
head(meta)
if(!(sample_col %in% colnames(meta))) stop('sample column not in meta file')
if(!all(colnames(pxc) %in% rownames(meta))) stop('Not all pxc cells in meta cells')
meta <- meta[colnames(pxc),]
if(!identical(colnames(pxc),rownames(meta))) stop('cellnames not the same between pxc and meta files')


print_time("Normalize")
set.seed(seed)
pxc_norm <- TF.IDF(pxc)
pxc_norm <- LogNormalize(pxc_norm)
saveRDS(pxc_norm,paste(sep='',outPrefix,'_normE.rds'))


print_time("Variable")
nPeaks_end <- floor(nCells_tot * end_crit)
if(step_crit<1) nPeaks_step <- floor(nPeaks_tot * step_crit)*-1 else nPeaks_step <- step_crit*-1
cat(paste("Total number of (subsetted) peaks:",nPeaks_tot,"\n"))
cat(paste("Ending peak criteria:",nPeaks_end,"\n"))
cat(paste("Peak Step:",nPeaks_step,"\n"))
cat(paste("Initially asking for:",numPeaks,"\n"))
for(nPeaks_ask in seq(numPeaks,1,nPeaks_step)){
    
    set.seed(seed)
    var_peaks = vargenes_vst(pxc_norm, groups = meta[,sample_col], topn = nPeaks_ask)
    saveRDS(var_peaks,paste(sep='',outPrefix,'_varPeaks.rds'))
    nPeaks_var <- length(var_peaks)
    print_time(paste('Asked for',nPeaks_ask,'Got',nPeaks_var,'\n'))
    
    if(nPeaks_var < nPeaks_end) break
}


print_time("Scale")
pxc_scaled <- t(scale(t(pxc_norm[var_peaks,]))) #want features scaled to mean 0 and var 1

pxc_scaled[1:3,1:3]
cat('First column feature\n')
mean(pxc_scaled[1,]) #feature
var(pxc_scaled[1,]) #feature
cat('First row cell\n')
mean(pxc_scaled[,1]) #cell
var(pxc_scaled[,1]) #cell


print_time("Transpose")
cxp <- t(pxc_scaled)
dim(cxp)
cxp[1:3,1:3]
mean(cxp[,1])
var(cxp[,1])
identical(var_peaks,colnames(cxp))
identical(rownames(cxp),rownames(meta))

saveRDS(cxp,paste(sep='',outPrefix,'_normE_scaled_transpose.rds'))


print_time("Done.")
