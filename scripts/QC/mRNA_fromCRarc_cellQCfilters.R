print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(magrittr)
    library(stringr)
    library(ggplot2)
    library(viridis)
    library(argparse)

})


cellQC_barPlots <- function(thisMeta,nGene_TH,pMito_CO,tle,plotPrefix,height=4,width=12){
    g <- (ggplot(thisMeta) +
        geom_histogram(aes(x = nUMI))) +
            theme_bw(base_size = 20) +
            ggtitle(tle) +
            theme(text = element_text(size = 20),
                legend.position = "none",
                panel.grid = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1)
                  )+ 
    (ggplot(thisMeta) +
        geom_histogram(aes(x = nGene)) +
        geom_vline(xintercept = nGene_TH, col = "darkred")) +
            theme_bw(base_size = 20) +
            theme(text = element_text(size = 20),
                legend.position = "none",
                panel.grid = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1)
                  ) +
    (ggplot(thisMeta) +
        geom_histogram(aes(x = percent_mito)) +
        geom_vline(xintercept = pMito_CO, col = "darkred")) +
            theme_bw(base_size = 20) +
            theme(text = element_text(size = 20),
                legend.position = "none",
                panel.grid = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1)
                  )
    ggsave(paste(sep='',plotPrefix,'_barPlot.png'),plot=g,units='in',height=height,width=width)
    
    return(g)
}

get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix,iy)
    return(dens$z[ii])
}

cellQC_densityPlots <- function(thisMeta,nGene_TH,pMito_CO,tle,plotPrefix,height=5,width=7){
    plot_dens <- get_density(thisMeta$percent_mito, thisMeta$nGene, n = 100)
    
    g <- ggplot(thisMeta, aes(x = percent_mito, y = nGene, color = plot_dens)) +
        geom_point(size = .3) +
        scale_x_continuous(trans = "log2") + 
        scale_y_continuous(trans = "log2") + 
        geom_hline(yintercept = nGene_TH, color = "darkred") +
        geom_vline(xintercept = pMito_CO, color = "darkred") +
        scale_color_viridis() +
        labs(color = "density") +
        theme_classic(base_size = 20) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(color="black", size=22, face = "italic")
         ) + 
        ggtitle(tle)
    ggsave(paste(sep='',plotPrefix,'_nGxpM_log2density.png'),plot=g,units='in',height=height,width=width)

    g <- ggplot(meta, aes(x = percent_mito, y = nGene, color = plot_dens)) +
        geom_point(size = .3) +
        geom_hline(yintercept = nGene_TH, color = "darkred") +
        geom_vline(xintercept = pMito_CO, color = "darkred") +
        scale_color_viridis() +
        labs(color = "density") +
        theme_classic(base_size = 20) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(color="black", size=22, face = "italic")
         ) + 
        ggtitle(tle)
    ggsave(paste(sep='',plotPrefix,'_nGxpM_density.png'),plot=g,units='in',height=height,width=width)
}


parser <- ArgumentParser(description='mRNA CR-arc matrix cell QC filters')
parser$add_argument('matrix_dir', metavar='matrix_dir', help='CR-arc matrix output directory; e.g., filtered_feature_bc_matrix')
parser$add_argument('sample_ID', metavar='sample_ID', help='sample ID used in output files')
parser$add_argument('-nGene_thresh', metavar='nGene_thresh', default=500, type='integer', help='cells must have at least this number of unique genes; default: 500')
parser$add_argument('-pMito_cutoff', metavar='pMito_cutoff', default=0.20, type='double', help='cells must have less than this percent of mitochondrial genes; default: 0.20')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
args <- parser$parse_args()

matrix_dir <- args$matrix_dir
sample_ID <- args$sample_ID
nGene_thresh <- args$nGene_thresh
pMito_cutoff <- args$pMito_cutoff
outDir <- args$outDir

print_time("Argument Checking")
if(!all(file.exists(matrix_dir))) stop("Input files don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)

cat("Arguments\n")
outFile <- paste(sep="",outDir,sample_ID,"_mRNA_fromCRarc_cellQCfilters_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load 10x genes x cells matrix and rename cells")

dat <- Seurat::Read10X(matrix_dir)
ls(dat)

mRNA_exprs <- dat$`Gene Expression` %>% as("dgCMatrix")
dim(mRNA_exprs)

mRNA_cellnames <- str_split_fixed(colnames(mRNA_exprs),'-',2)
new_mRNA_colnames <- paste(sep='',sample_ID,"_",mRNA_cellnames[,1])
head(cbind(colnames(mRNA_exprs),new_mRNA_colnames))
colnames(mRNA_exprs) <- new_mRNA_colnames
mRNA_exprs[1:3,1:3]

saveRDS(mRNA_exprs,paste(sep="",outDir,sample_ID,'_CR-arc_genesXcells.rds'))


print_time("Start Meta")

cellnames_split <- str_split_fixed(colnames(mRNA_exprs),'_',2)
meta <- as.data.frame(cellnames_split,stringsAsFactors=FALSE)
colnames(meta) <- c('sample','CB')
rownames(meta) <- colnames(mRNA_exprs)
meta$cell <- colnames(mRNA_exprs)
dim(meta)
head(meta)


print_time("Pre-Cell QC stats")

meta$nUMI <- Matrix::colSums(mRNA_exprs)
meta$nGene <- Matrix::colSums(mRNA_exprs > 0)
mito_genes <- grep("^MT-", rownames(mRNA_exprs), value = TRUE, ignore.case = TRUE)
meta$percent_mito <- Matrix::colSums(mRNA_exprs[mito_genes, ])/Matrix::colSums(mRNA_exprs)
length(mito_genes)
head(meta)

cat(paste(sep='',"Inital number of cells: ",ncol(mRNA_exprs),"\n"))

saveRDS(meta,paste(sep="",outDir,sample_ID,'_CR-arc_meta.rds'))

cellQC_barPlots(meta,nGene_thresh,pMito_cutoff,sample_ID,paste(sep='',outDir,sample_ID,'_cellQC_pre'))
cellQC_densityPlots(meta,nGene_thresh,pMito_cutoff,sample_ID,paste(sep='',outDir,sample_ID,'_cellQC_pre'))


print_time("Filter nUMI/nGene")
meta <- subset(meta, nGene > nGene_thresh & percent_mito < pMito_cutoff)
dim(meta)
mRNA_exprs <- mRNA_exprs[, meta$cell]

cat(paste(sep='',"post cellQC number of cells: ",ncol(mRNA_exprs),"\n"))

saveRDS(mRNA_exprs,paste(sep="",outDir,sample_ID,'_CR-arc_postQC_genesXcells.rds'))
saveRDS(meta,paste(sep="",outDir,sample_ID,'_CR-arc_postQC_mRNA_meta.rds'))

cellQC_barPlots(meta,nGene_thresh,pMito_cutoff,sample_ID,paste(sep='',outDir,sample_ID,'_cellQC_post'))
cellQC_densityPlots(meta,nGene_thresh,pMito_cutoff,sample_ID,paste(sep='',outDir,sample_ID,'_cellQC_post'))

g <- ggplot(meta,aes_string(x='nGene',y='nUMI')) + geom_point(size=1,alpha=0.5) +
        theme_bw(base_size=15) + ggtitle(sample_ID) +
        scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
ggsave(paste(sep='',outDir,sample_ID,'_cellQC_post_nGxnU_dotPlot.png'),plot=g,units='in',height=5,width=6)


print_time("Scrublet input")
writeMM(mRNA_exprs, paste(sep='',outDir,sample_ID,'_CR-arc_postQC_genesXcells.mm'))


print_time('Done.')
