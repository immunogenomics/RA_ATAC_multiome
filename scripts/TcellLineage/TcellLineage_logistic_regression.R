print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(lme4)
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



parser <- ArgumentParser(description='identity versus function linear model')
parser$add_argument('data_file', metavar='data_file', help='Matrix file: binary peaks x cells')
parser$add_argument('meta_file', metavar='meta_file', help='Cell metadata file')
parser$add_argument('--sample_col', metavar='sample_col', default='sample', help='sample/donor column in meta file; default: sample')
parser$add_argument('--frag_col', metavar='frag_col', default='nFrags', help='number of fragments column in meta file; default: nFrags')
parser$add_argument('--cellType_col', metavar = 'cellType_col', default='cluster_name', help='cell type column; default: cluster_name')
parser$add_argument('--simplifyCT', action='store_true', help='if given, change X-#: description to X#')
parser$add_argument('lineage_file', metavar='lineage_file', help='cells x lineage file: 1 column: 1 if CD4+ CD8A-; 0 if double pos/neg; -1 if CD4- CD8A+')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
parser$add_argument('--peak2gene_file', metavar='peak2gene_file', help='peak to gene file')
parser$add_argument('--subsetPeak_file', metavar='subsetPeak_file', help='subsetPeak_file')
parser$add_argument('--seed', metavar='seed', type='integer', default=1234567890,help='Seed for randomization')
parser$add_argument('--padj', metavar='padj', default='BH',choices=c('holm','hochberg','hommel','bonferroni','BH','BY','fdr'), help='if given, do pvalue-adjust; options: holm, hochberg, hommel, bonferroni, BH, BY, fdr; default=BH')
args <- parser$parse_args()

data_file <- args$data_file
meta_file <- args$meta_file
sample_col <- args$sample_col
frag_col <- args$frag_col
lineage_file <- args$lineage_file
cellType_col <- args$cellType_col
simplifyCT <- args$simplifyCT
outDir <- args$outDir
prefix <- args$prefix
peak2gene_file <- args$peak2gene_file
subsetPeak_file <- args$subsetPeak_file
seed <- args$seed
padj <- args$padj


if(is.null(peak2gene_file)){
    cat('peak 2 gene file not provided. using peak names as gene names.\n')
}
if(is.null(subsetPeak_file)){
    cat('subset peak file not provided. using all peak in data mat.\n')
}

print_time("Argument Checking")
if(!all(file.exists(data_file,meta_file,lineage_file))) stop("Input files don't exist.")

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
lineage <- load_file(lineage_file)
if(!identical(c(-1,0,1),sort(unique(lineage)))) stop("lineage should only be -1,0,1")

dim(data)
dim(meta)
length(lineage)

if(!identical(sort(names(lineage)),sort(colnames(data)))) stop("Not all lineage file cells in pxc matrix")
if(!identical(sort(names(lineage)),sort(rownames(meta)))) stop("Not all lineage file cells in meta df")

data <- data[,rownames(meta)]
lineage <- lineage[rownames(meta)]
if(!identical(names(lineage),colnames(data))) stop("Not all lineage file cells in pxc matrix")
if(!identical(names(lineage),rownames(meta))) stop("Not all lineage file cells in meta df")


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
#don't need to run all zero or all one peaks! note, this method only works with binary peaks x cells
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


print_time("cell type design matrix")
CT_binDF <- createDonorMat(meta,cellType_col)
if(simplifyCT & all(grepl('^[A-Za-z]+-[0-9]+: ',colnames(CT_binDF),perl=TRUE))) colnames(CT_binDF) <- str_replace(str_split_fixed(colnames(CT_binDF),":",2)[,1],"-","")
CT_binDF <- CT_binDF[,mixedsort(colnames(CT_binDF))]
CT_binMat <- as.matrix(CT_binDF[,-ncol(CT_binDF)]) #need to exclude one cell type to not overdetermine
dim(CT_binMat)
head(CT_binMat)


print_time('for loop over genes')
donor_col <- meta[,sample_col]
fragments_scale_log10_col <- scale(log10(meta[,frag_col]))
lineage_col <- lineage[1:length(lineage)]
pval_df <- data.frame('peak'=character(),'gene'=character(),
                      'LogLik'=numeric(),'Chisq'=numeric(),'pval'=numeric(),
                      'lineage_beta'=numeric(),'lineage_pval'=numeric(),
                      stringsAsFactors=FALSE) 
stats_df <- data.frame('peak'=character(),'gene'=character(),'feature'=character(),'value'=numeric(),stringsAsFactors=FALSE)

for(peak in peaksToRun){
    gene <- peak2gene[which(peak2gene$peak==peak),'gene']
    
    data_vec <- as.vector(t(data[peak,]))
    
    resGlmer <- tryCatch({
        set.seed(seed)
        
        m1 <- glmer(data_vec ~ lineage_col + CT_binMat + (1|donor_col) + fragments_scale_log10_col,family='binomial',
                    nAGQ=0, control = glmerControl(optimizer = "nloptwrap"))
        m2 <- glmer(data_vec ~ CT_binMat + (1|donor_col) + fragments_scale_log10_col,family='binomial',
                    nAGQ=0, control = glmerControl(optimizer = "nloptwrap"))
        
        res <- lmtest::lrtest(m1,m2)
        
        glmer_coef <- as.data.frame(summary(m1)$coefficient)
    }, error = function(err) {
        cat(paste(sep="\n",paste("ERROR: glmer",peak),err))
    })
    if(!is.null(resGlmer)){ 
        pval_df <- rbind(pval_df, 
                         data.frame('peak'=peak,'gene'=gene,'LogLik'=res[nrow(res),'LogLik'],
                                    'Chisq'=res[nrow(res),'Chisq'],'pval'=res[nrow(res),'Pr(>Chisq)'],
                                    'lineage_beta'=glmer_coef['lineage_col','Estimate'],
                                    'lineage_pval'=glmer_coef['lineage_col','Pr(>|z|)'],
                                    stringsAsFactors=FALSE))  
        
        glm_type <- 'glmer'
        numFeatures <- nrow(glmer_coef)*4
        this_df <- data.frame('peak'=rep(peak,numFeatures),
                              'gene'=rep(gene,numFeatures),
                              'feature'=c(paste(sep='_',glm_type,rownames(glmer_coef),'est'), 
                                          paste(sep='_',glm_type,rownames(glmer_coef),'stdErr'),
                                          paste(sep='_',glm_type,rownames(glmer_coef),'Z'),
                                          paste(sep='_',glm_type,rownames(glmer_coef),'PrZ')),
                              'value'=c(glmer_coef[,'Estimate'],glmer_coef[,'Std. Error'],
                                        glmer_coef[,'z value'],glmer_coef[,'Pr(>|z|)']),
                              stringsAsFactors=FALSE)
        stats_df <- rbind(stats_df,this_df)
    } else {
        cat(paste('SKIPPING:',peak,'\n'))
    }
}


print_time('MHTC')
pval_df$padj <- p.adjust(pval_df$pval,method=padj)
pval_df$log10padj <- -log10(pval_df$padj)


print_time("Save")
saveRDS(pval_df,paste(sep='',outPrefix,'_pval.rds'))
saveRDS(stats_df,paste(sep='',outPrefix,'_stats.rds'))


print_time('Done.')

