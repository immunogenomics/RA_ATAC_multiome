print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(MASS)
    library(stringr)
    library(gtools)
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

createDonorMat <- function(this_meta,this_sample_col){
    donors <- unique(this_meta[,this_sample_col])
    for(idx in 1:length(donors)){
        thisDonor <- donors[idx]
        thisDonorCol <- as.integer(ifelse(this_meta[,this_sample_col]==thisDonor,1,0))

        if(idx==1){
            this_donor_df <- data.frame('donor'=thisDonorCol)
            rownames(this_donor_df) <- rownames(this_meta)
            colnames(this_donor_df) <- thisDonor
        } else {
            this_donor_df[,thisDonor] <- thisDonorCol
        }   
    }
    ll <- colSums(this_donor_df)
    rr <- table(this_meta[,this_sample_col])
    ss <- as.numeric(rr)
    names(ss) <- names(rr)
    theseNames <- sort(names(rr))
    if(!identical(ll[theseNames],ss[theseNames])){
        print(ll[theseNames])
        print(ss[theseNames])
        stop("donor issues")
    }
    
    return(this_donor_df)
}



parser <- ArgumentParser(description='differential (pseudobulk) peaks via anova')
parser$add_argument('data_file', metavar='data_file', help='Matrix file: peaks x pseudobulk')
parser$add_argument('meta_file', metavar='meta_file', help='pseudobulk metadata file')
parser$add_argument('--sample_col', metavar='sample_col', default='sample', help='sample/donor column in meta file; default: sample')
parser$add_argument('--frag_col', metavar='frag_col', default='sample_CITE_nFrag', help='number of fragments column in meta file; default: sample_CITE_nFrag')
parser$add_argument('--cellType_col', metavar = 'cellType_col', default='CITE_abbr', help='cell type column; default: CITE_abbr')
parser$add_argument('--simplifyCT', action='store_true', help='if given, change X-#: description to X#')
parser$add_argument('--cellType_ref', metavar='cellType_ref', help='reference cell type to NOT include in the 1-hot cell type matrix')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
parser$add_argument('--peak2gene_file', metavar='peak2gene_file', help='peak to gene file')
parser$add_argument('--subsetPeak_file', metavar='subsetPeak_file', help='subsetPeak_file')
args <- parser$parse_args()

data_file <- args$data_file
meta_file <- args$meta_file
sample_col <- args$sample_col
frag_col <- args$frag_col
cellType_col <- args$cellType_col
simplifyCT <- args$simplifyCT
cellType_ref <- args$cellType_ref
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
if(!all(c(sample_col,frag_col,cellType_col) %in% colnames(meta))) stop("Not all required columns in meta file")
if(!identical(colnames(data),rownames(meta))) stop("Pseudobulks not the same")
dim(data)
dim(meta)


print_time('checking subset peak file')
if(!is.null(subsetPeak_file)){
    peaksToRun <- load_file(subsetPeak_file)
    if(!is.vector(peaksToRun)) stop('peaks to run should be a vector')
    if(length(unique(peaksToRun))!=length(peaksToRun)) stop('peaks should not be duplicated')
    if(!all(peaksToRun %in% rownames(data))){
        cat(paste("Peak(s) you asked for is not in the data matrix. Subsetting from",length(peaksToRun)))
        peaksToRun <- peaksToRun[which(peaksToRun %in% rownames(data))]
        cat(paste(sep=''," to ",length(peaksToRun),".\n"))
    }
} else {
    peaksToRun <- rownames(data)
}
cat(paste(sep='','Initial number of peaksToRun: ',length(peaksToRun),'\n'))
rs <- rowSums(data[peaksToRun,]) 
cat(paste(sep='','Number of all 0 peaks: ',length(rs[rs==0]),'\n','Number of NOT all 0 peaks: ',length(rs[rs!=0]),'\n'))
peaksToRun <- names(rs[which(rs!=0)])
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


print_time("cell type design matrix")
CT_binDF <- createDonorMat(meta,cellType_col)
if(simplifyCT & all(grepl('^[A-Za-z]+-[0-9]+: ',colnames(CT_binDF),perl=TRUE))) colnames(CT_binDF) <- str_replace(str_split_fixed(colnames(CT_binDF),":",2)[,1],"-","")
CT_binDF <- CT_binDF[,mixedsort(colnames(CT_binDF))]
if(!is.null(cellType_ref)){
    if(cellType_ref %in% colnames(CT_binDF)){
        CT_binMat <- as.matrix(CT_binDF[,which(colnames(CT_binDF)!=cellType_ref)]) #need to exclude one cell type to not overdetermine
        cat(paste('Removing requested cell type reference: ',cellType_ref,'\n'))
    } else {
        cat(paste(sep='','WARNING: The requested cell type reference (',cellType_ref,
                  ') does not exist in this dataset, so removing the last alphabetical cell type (',
                  colnames(CT_binDF)[ncol(CT_binDF)],')\n'))
        CT_binMat <- as.matrix(CT_binDF[,-ncol(CT_binDF)]) #need to exclude one cell type to not overdetermine
    }
} else {
    CT_binMat <- as.matrix(CT_binDF[,-ncol(CT_binDF)]) #need to exclude one cell type to not overdetermine
}
dim(CT_binMat)
head(CT_binMat)


print_time("sample design matrix")
sample_binDF <- createDonorMat(meta,sample_col)
sample_binDF <- sample_binDF[,sort(colnames(sample_binDF))]
if(paste(sep='',sample_col,'_nFrag') %in% colnames(meta)){
    cat(paste('Excluding lowest sample count\n'))
    ll <- unique(meta[,c(sample_col,paste(sep='',sample_col,'_nFrag'))])
    ll <- ll[order(ll[,paste(sep='',sample_col,'_nFrag')]),]
    sample_binMat <- as.matrix(sample_binDF[,which(colnames(sample_binDF)!=ll[1,sample_col])]) #need to exclude one cell type to not overdetermine
} else {
    cat(paste('Excluding last name sorted sample\n'))
    sample_binMat <- as.matrix(sample_binDF[,-ncol(sample_binDF)]) #need to exclude one cell type to not overdetermine
}
colnames(sample_binMat) <- str_replace(colnames(sample_binMat),"-","")
dim(sample_binMat)
head(sample_binMat)


print_time('for loop over genes')
fragments_scale_log10_col <- scale(log10(meta[,frag_col]))
if(nrow(CT_binMat)!=nrow(sample_binMat) | nrow(CT_binMat)!=length(fragments_scale_log10_col) | nrow(CT_binMat)!=ncol(data)) stop('dimensions issues')

pval_df <- data.frame('peak'=character(),'gene'=character(),'theta'=numeric(),'ResidDF'=numeric(),'2xLL'=numeric(),
                      'df'=numeric(),'LRstat'=numeric(),'pval'=numeric(),
                      stringsAsFactors=FALSE) 
stats_df <- data.frame('peak'=character(),'gene'=character(),'feature'=character(),'value'=numeric(),stringsAsFactors=FALSE)

for(peak in peaksToRun){
    gene <- peak2gene[which(peak2gene$peak==peak),'gene']
    
    data_vec <- as.vector(t(data[peak,]))
        
    resGlm <- tryCatch(withCallingHandlers({
        m1 <- glm.nb(data_vec ~ CT_binMat + sample_binMat + fragments_scale_log10_col)
        m2 <- glm.nb(data_vec ~ sample_binMat + fragments_scale_log10_col)
        
        res <- anova(m1,m2)
        
        glm_coef <- as.data.frame(summary(m1)$coefficient)
    }, error = function(err) {
        cat(paste(sep="\n",paste("ERROR: glm",peak),err))
    }),
    error = function(err){
        doNothing = 1
    })
    if(exists("res")){ 
        pval_df <- rbind(pval_df, 
                         data.frame('peak'=peak,'gene'=gene,'theta'=res[nrow(res),'theta'],'ResidDF'=res[nrow(res),'Resid. df'],
                                    '2xLL'=res[nrow(res),'   2 x log-lik.'],'df'=res[nrow(res),'   df'],
                                    'LRstat'=res[nrow(res),'LR stat.'],'pval'=res[nrow(res),'Pr(Chi)'],
                                    stringsAsFactors=FALSE))
        
        glm_type <- 'glm'
        numFeatures <- nrow(glm_coef)*4
        this_df <- data.frame('peak'=rep(peak,numFeatures),
                              'gene'=rep(gene,numFeatures),
                              'feature'=c(paste(sep='_',glm_type,rownames(glm_coef),'est'), 
                                          paste(sep='_',glm_type,rownames(glm_coef),'stdErr'),
                                          paste(sep='_',glm_type,rownames(glm_coef),'Z'),
                                          paste(sep='_',glm_type,rownames(glm_coef),'PrZ')),
                              'value'=c(glm_coef[,'Estimate'],glm_coef[,'Std. Error'],
                                        glm_coef[,'z value'],glm_coef[,'Pr(>|z|)']),
                              stringsAsFactors=FALSE)
        stats_df <- rbind(stats_df,this_df)
        
        rm(res)
    } else {
        cat(paste('SKIPPING:',peak,'\n'))
    }
}


print_time('BH MHTC')
pval_df$padj <- p.adjust(pval_df$pval,method='BH')
pval_df$log10padj <- -log10(pval_df$padj)


print_time("save")
saveRDS(pval_df,paste(sep='',outPrefix,'_pval.rds'))
saveRDS(stats_df,paste(sep='',outPrefix,'_stats.rds'))


print_time('Done.')

