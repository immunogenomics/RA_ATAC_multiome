print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(lme4)
    library(stringr)
    library(ggplot2)
    library(ggrepel)
    library(argparse)
})

file_parts <- function(fl) {
    x = unlist(strsplit(fl,"[.]"))
    ret = paste(x[-length(x)],collapse=".")
    return(c(ret,x[length(x)]))
}


load_file <- function(fl) {
    if(tolower(file_parts(fl)[2])=="rds"){
        contents <- readRDS(fl)
    } else if(tolower(file_parts(fl)[2])=="txt") {
        contents <- read.table(fl,sep="\t",stringsAsFactors=FALSE,header=FALSE)
        if(ncol(contents)==1) contents <- as.vector(contents$V1)
    } else {
        stop(paste("File format of input file not recognized:",fl))
    }

    return(contents)
}



parser <- ArgumentParser(description='differential peaks')
parser$add_argument('data_file', metavar='data_file', help='Matrix file: binary peaks x cells')
parser$add_argument('meta_file', metavar='meta_file', help='Cell metadata file')
parser$add_argument('--sample_col', metavar='sample_col', default='Sample', help='sample/donor column in meta file; default: Sample')
parser$add_argument('--frag_col', metavar='frag_col', default='nFrags', help='number of fragments column in meta file; default: nFrags')
parser$add_argument('cellType_file', metavar='cellType_file', help='file with cell types; can be the same as the meta_file')
parser$add_argument('--cellType_col', metavar = 'cellType_col', default='cellType', help='cell type column; default: cellType')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
parser$add_argument('--peak2gene_file', metavar='peak2gene_file', help='peak to gene file')
parser$add_argument('--subsetPeak_file', metavar='subsetPeak_file', help='subsetPeak_file')
args <- parser$parse_args()

data_file <- args$data_file
meta_file <- args$meta_file
sample_col <- args$sample_col
frag_col <- args$frag_col
cellType_file <- args$cellType_file
cellType_col <- args$cellType_col
outDir <- args$outDir
prefix <- args$prefix
peak2gene_file <- args$peak2gene_file
subsetPeak_file <- args$subsetPeak_file


if(is.null(peak2gene_file)){
    cat('peak 2 gene file not provided. using peak names as gene names.\n')
}
if(is.null(subsetPeak_file)){
    cat('subset peak file not provided. using all peak in data mat.\n')
}

print_time("Argument Checking")
if(!all(file.exists(data_file,meta_file,cellType_file))) stop("Input files don't exist.")

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
if(!all(c(sample_col,frag_col) %in% colnames(meta))) stop("Not all required columns in meta file")
cellType <- readRDS(cellType_file)
if(!(cellType_col %in% colnames(cellType))) stop("Required column not in cell type file")
if(!all(rownames(cellType) %in% colnames(data))) stop("Not all cell type file cells in pxc matrix")
if(!all(rownames(cellType) %in% rownames(meta))) stop("Not all cell type file cells in meta df")
dim(data)
dim(meta)
dim(cellType)


print_time("subset by cell type cells")
data <- data[,rownames(cellType)]
meta <- meta[rownames(cellType),]
dim(data)
dim(meta)
dim(cellType)
if(!(identical(rownames(cellType),colnames(data)) & identical(rownames(cellType),rownames(meta)))) stop('cells not identical between data structures')


print_time('checking subset peak file')
if(!is.null(subsetPeak_file)){
    peaksToRun <- load_file(subsetPeak_file)
    if(!is.vector(peaksToRun)) stop('peaks to run should be a vector')
    if(length(unique(peaksToRun))!=length(peaksToRun)) stop('peaks should not be duplicated')
    if(!all(peaksToRun %in% rownames(data))) stop("Peak(s) you asked for is not in the data matrix. Please try again.")
} else {
    peaksToRun <- rownames(data)
}
cat(paste(sep='','Initial number of peaksToRun: ',length(peaksToRun),'\n'))
#don't need to run all zero or all one peaks! Binomial model should have both 0s and 1s.
rs <- rowSums(data[peaksToRun,])
cat(paste(sep='','Number of all 0 peaks: ',length(rs[rs==0]),'\n','Number of NOT all 0 peaks: ',length(rs[rs!=0]),'\n'))
cat(paste(sep='','Number of all 1 peaks: ',length(rs[rs==ncol(data)]),'\n','Number of NOT all 1 peaks: ',length(rs[rs!=ncol(data)]),'\n'))
peaksToRun <- names(rs[which(rs!=0 & rs!=ncol(data))])
cat(paste(sep='','Number of peaksToRun: ',length(peaksToRun),'\n'))
head(peaksToRun)


print_time('checking peak2gene file')
if(!is.null(peak2gene_file)){
    peak2gene <- readRDS(peak2gene_file)
    if(!identical(colnames(peak2gene),c('peak','gene'))) stop('peak2gene colnames should be peak and gene')
    if(length(unique(peak2gene$peak))!=nrow(peak2gene)) stop('peak can map to at most 1 gene')
    if(!all(rownames(data) %in% peak2gene$peak)){
        cat('Not all peaks in peak2gene\n')
        print(table(rownames(data) %in% peak2gene$peak))
        
        missingPeaks <- rownames(data)[which(!(rownames(data) %in% peak2gene$peak))]
        print(length(missingPeaks))
        
        peak2gene <- rbind(peak2gene, data.frame('peak'=missingPeaks,
                                                 'gene'=rep(NA,length(missingPeaks)),
                                                 stringsAsFactors=FALSE))
        
        print(table(rownames(data) %in% peak2gene$peak))
    }
} else {
    peak2gene <- data.frame('peak'=rownames(data),'gene'=rep(NA,nrow(data)),stringsAsFactors=FALSE)
}
rbind(head(peak2gene,n=3),tail(peak2gene,n=3))


print_time('for loops over genes and cell types')
donor_col <- meta[,sample_col]
fragments_scale_log10_col <- scale(log10(meta[,frag_col]))
pval_df <- data.frame('peak'=character(),'gene'=character(),'cellType'=character(),
                      'log2FC'=numeric(),'absLog2FC'=numeric(),'pval'=numeric(),
                      stringsAsFactors=FALSE)

for(peak in peaksToRun){
    gene <- peak2gene[which(peak2gene$peak==peak),'gene']
    
    data_vec <- as.vector(t(data[peak,]))
    
    for(CT in sort(unique(cellType[,cellType_col]))){
        CT_col <- ifelse(cellType[,cellType_col]==CT,1,0)
        
        resGlmer <- tryCatch({
            m1 <- glmer(data_vec ~ CT_col + (1|donor_col) + fragments_scale_log10_col,family='binomial',
                        nAGQ=0, control = glmerControl(optimizer = "nloptwrap"))
            m2 <- glmer(data_vec ~ (1|donor_col) + fragments_scale_log10_col,family='binomial',
                        nAGQ=0, control = glmerControl(optimizer = "nloptwrap"))
            res <- lmtest::lrtest(m1,m2)
        }, error = function(err) {
            cat(paste(sep="\n",paste("ERROR: glmer",peak,CT),err))
        })
        if(!is.null(resGlmer)){ 
            
            result <- tryCatch({
                log2FC <- summary(m1)$coefficient['CT_col','Estimate']
            }, warning = function(war) {
                cat(paste(sep="\n",paste("WARNING: log2FC",peak,CT),war))
            }, error = function(err) {
                cat(paste(sep="\n",paste("ERROR: log2FC",peak,CT),err))
            })
            if(is.null(result)) log2FC <- NA
            
            pval_df <- rbind(pval_df,
                             data.frame('peak'=peak,'gene'=gene,'cellType'=CT,'log2FC'=log2FC,'absLog2FC'=abs(log2FC),
                                        'pval'=res[nrow(res),ncol(res)],
                                        stringsAsFactors=FALSE))
        } else {
            cat(paste('SKIPPING:',peak,CT,'\n'))
        }
        
    }   
}


print_time('BH MHTC')
pval_df$padj <- p.adjust(pval_df$pval,method='BH')
pval_df$log10padj <- -log10(pval_df$padj)


print_time("Save")
saveRDS(pval_df,paste(sep='',outPrefix,'_pval.rds'))


print_time('Done.')

