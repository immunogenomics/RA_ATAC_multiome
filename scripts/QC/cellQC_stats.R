suppressMessages({
    library(GenomicRanges)
    library(Matrix)
    library(Matrix.utils)
    library(MASS)
    library(ggplot2)
    library(viridis)
    library(argparse)
})


gRange_str <- function(gr) {
    return(paste(sep="",seqnames(gr),":",start(gr),"-",end(gr)))
}

substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}

make_GR <- function(fl,colnm){
    df <- read.table(fl,sep="\t",stringsAsFactors=FALSE)
    cat(paste(sep="",dim(df),"\n"))
    colnames(df)<-colnm
    GR <- GRanges(df)
    return(GR)
}

fragment_to_cell_overlap <- function(fragments,criteria,binarize){
    fragments_oL_crit <- countOverlaps(fragments,criteria)
    if(binarize) fragments_oL_crit <- ifelse(fragments_oL_crit==0,0,1)
    names(fragments_oL_crit)<-fragments$cell
    cells_oL_crit <- tapply(fragments_oL_crit,names(fragments_oL_crit),sum)
    
    return(cells_oL_crit)
}

pad_vector <- function(curr, goal){
    cells_to_add<-names(goal)[!(names(goal) %in% names(curr))]
    cells_to_rm<-names(curr)[!(names(curr) %in% names(goal))]
    
    cat(sep="","Number of cells added: ",length(cells_to_add),"\n")
    cat(sep="","Number of cells removed: ",length(cells_to_rm),"\n")
    
    if(length(curr)+length(cells_to_add)-length(cells_to_rm)!=length(goal)) stop("ERROR: lengths differ.")
    
    CBct_to_add <- integer(length(cells_to_add))
    names(CBct_to_add)<-cells_to_add
    
    adjusted <- c(curr[!(names(curr) %in% cells_to_rm)],CBct_to_add)
    adjusted <- adjusted[order(names(adjusted))]
    
    if((length(adjusted)!=length(goal)) || (sum(adjusted)!=sum(curr[!(names(curr) %in% cells_to_rm)]))){
        stop("ERROR: length and/or sum not the same!")
    }
    
    return(adjusted)
}

perc_CB <- function(top,bottom,ratio=FALSE){
    if(length(top)!=length(bottom)) top <- pad_vector(top,bottom)
    
    top <- top[order(names(top))]
    bottom <- bottom[order(names(bottom))]
    
    if(!identical(names(top),names(bottom))) stop("cell names not identical")
    
    percCB <- top/bottom
    
    if(!ratio & !all(min(percCB)>=0,max(percCB)<=1)) stop("bounds not between 0 and 1.")
    
    return(percCB)
    
}

make_CB_vec <- function(fl){
    df <- read.table(fl,sep="\t",stringsAsFactors=FALSE)
    colnames(df)<-c('CB','CBct')
    CBct <- as.vector(df$CBct)
    names(CBct)<-df$CB
    
    return(CBct)
}

readCt_from_frag_GR <- function(GR){
    vec <- GR$dupCt*2
    names(vec) <- GR$cell
    vec_byCell <- tapply(vec,names(vec),sum)
    
    return(vec_byCell)
}

get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix,iy)
    return(dens$z[ii])
}

doEverything <- function(x,y,NODP,WGAS,RDFG,cellThresh,crit,cutoff,oPrefix, histXmax=1, maxGenome=1000000, maxCrit=200000, above=TRUE){
    
    print_time("Getting Percent")
    x_pad <- pad_vector(x,y)
    x_perc <- perc_CB(x_pad,y)
    
    print_time("Histograms")
    tle <- paste(sep="",prefix," ",NODP," % ",WGAS," ",RDFG," in ",crit," minC ",cellThresh)
    png(paste(sep="",paste(sep="_",oPrefix,"QC","hist",crit,NODP,WGAS,RDFG,paste(sep="","cellTh",cellThresh)),".png"))
    hist(x_perc,main=tle,xlab=paste(sep="","% ",WGAS," ",RDFG," in ",crit),breaks=100,xlim=c(0,histXmax))
    dev.off()
    
    
    print_time("Counts")
    if(above) cat(paste("With",crit,"cutoff",cutoff,"keeping: ",length(x_perc[x_perc>cutoff]),"\n")) else cat(paste("With",crit,"cutoff",cutoff,"keeping: ",length(x_perc[x_perc<cutoff]),"\n"))
    
    
    print_time("Plot crit vs reads")
    x_df <- data.frame(
        total <- y,
        criteria <- x_pad
    )
    saveRDS(x_df,paste(sep="",paste(sep="_",oPrefix,"QC","vs",crit,NODP,WGAS,RDFG,paste(sep="","cellTh",cellThresh)),".RDS"))
    
    png(paste(sep="",paste(sep="_",oPrefix,"QC","vs",crit,NODP,WGAS,RDFG,paste(sep="","cellTh",cellThresh)),".png"))
    plot(x_df,xlab=paste(WGAS,NODP),ylab=crit,main=paste(prefix,"minC",cellThresh))
    dev.off()
    
    png(paste(sep="",paste(sep="_",oPrefix,"QC","vs",crit,NODP,WGAS,RDFG,paste(sep="","cellTh",cellThresh)),"_scale.png"))
    plot(x_df,xlab=paste(WGAS,NODP),ylab=crit,main=paste(prefix,"minC",cellThresh),xlim=c(0,maxGenome),ylim=c(0,maxCrit))
    dev.off()
    
    
    print_time("Density Plots")
    result <- tryCatch({
        plot_dens <- get_density(x_df$total,x_df$criteria, n = 100)
        options(repr.plot.height = 6, repr.plot.width = 6)
        ggplot(x_df, aes(x = total, y = criteria, color = plot_dens)) +
            geom_point(size = .3) +
            scale_color_viridis() +
            labs(title=paste(prefix,"minC",cellThresh), x=paste(WGAS,NODP), y=crit, color = "density") +
            theme_bw(base_size = 25) +
            theme(text = element_text(size = 25),
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 60, hjust = 1)
            )
        ggsave(paste(sep="",paste(sep="_",oPrefix,"QC","density",crit,NODP,WGAS,RDFG,paste(sep="","cellTh",cellThresh)),".png"))
    
    }, warning = function(war) {
        
        cat(paste(sep="\n","WARNING: Density Plot",war))
    
    }, error = function(err) {
        
        cat(paste(sep="\n","ERROR: Density Plot not generated",err))

    }, finally = {
        
        print_time("Done.")
        
    })
}


parser <- ArgumentParser(description='cellQC statistics - no actual filtering here - based on peaks, promoters, chrM, blacklist read counts')
parser$add_argument('fragment_file', metavar='fragment_file', help='chr1-22XY fragment file with 5 columns: chr, start, stop, cell barcode, duplicate count')
parser$add_argument('otherChr_fragment_file', metavar='otherChr_fragment_file', help='not chr1-22XY fragment file with 5 columns: chr, start, stop, cell barcode, duplicate count')
parser$add_argument('WG_CB_ct_file', metavar='WG_CB_ct_file', help='read counts by cell barcode for deduplicated chr1-22XY bam file')
parser$add_argument('blacklist_CB_ct_file', metavar='blacklist_CB_ct_file', help='read counts by cell barcode for reads in the blacklist')
parser$add_argument('peaks_CB_ct_file', metavar='peaks_CB_ct_file', help='read counts by cell barcode for reads in peak neighborhoods')
parser$add_argument('promoters_CB_ct_file', metavar='promoters_CB_ct_file', help='read counts by cell barcode for reads in promoters')
parser$add_argument('WG_noDedup_readCtByCell_RDS', metavar='WG_noDedup_readCtByCell_RDS', help='read counts by cell barcode for non-deduplicated genome as R vector; if it does not already exist provide an output filename where it will be saved within this script')
parser$add_argument('WG_dedup_readCtByCell_RDS', metavar='WG_dedup_readCtByCell_RDS', help='read counts by cell barcode for deduplicated chr1-22XY bam file as R vector; if it does not already exist provide an output filename where it will be saved within this script')
parser$add_argument('threshold', metavar='threshold', type='integer', help='Minimal number of reads')
parser$add_argument('--stats_table_file', metavar='stats_table_file', help='for peaks, promoters, chrM, blacklist, read proportion cutoff, plotting total read maximum, plotting criteria read maximum')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

fragment_file <- args$fragment_file
otherChr_fragment_file <- args$otherChr_fragment_file
WG_CB_ct_file <- args$WG_CB_ct_file
blacklist_CB_ct_file <- args$blacklist_CB_ct_file
peaks_CB_ct_file <- args$peaks_CB_ct_file
promoters_CB_ct_file <- args$promoters_CB_ct_file
WG_noDedup_readCtByCell_RDS <- args$WG_noDedup_readCtByCell_RDS
WG_dedup_readCtByCell_RDS <- args$WG_dedup_readCtByCell_RDS
threshold <- args$threshold
stats_table_file <- args$stats_table_file
outDir <- args$outDir
prefix <- args$prefix


print_time("Argument Checking")
if(!all(file.exists(c(fragment_file,otherChr_fragment_file,WG_CB_ct_file,blacklist_CB_ct_file,peaks_CB_ct_file,promoters_CB_ct_file)))) stop("Input file(s) don't exist.")

#if these do not exist; they will be made later.
#WG_noDedup_readCtByCell_RDS
#WG_dedup_readCtByCell_RDS

if(!is.null(stats_table_file)){
    if(!file.exists(stats_table_file)) stop("If given, stats table file must exist.")
    
    statsTable <- read.table(stats_table_file,sep='\t',header=TRUE,row.names=1,stringsAsFactors=FALSE)
    
    if(!identical(sort(rownames(statsTable)),sort(c('peaks','promoters','chrM','blacklist')))){
        stop("stats table requires rownames: peaks, promoters, chrM, and blacklist")
    }
    if(!identical(sort(colnames(statsTable)),sort(c('cutoff','maxGenome','maxCrit')))){
        stop("stats table requires colnames: cutoff, maxGenome, and maxCrit")
    }
    
    if(!all(is.numeric(statsTable$cutoff),is.integer(statsTable$maxGenome),is.integer(statsTable$maxCrit),
            statsTable$cutoff<=1,statsTable$cutoff>=0)){
        stop("stats table needs cutoffs between 0 and 1 and integer maxGenome and maxCrit values")
    }
} else {
    statsTable <- data.frame(cutoff=c(0.5,0.1,0.1,0.1),
                             maxGenome=c(1000000,1000000,1000000,1000000),
                             maxCrit=c(320000,120000,40000,120000))
    rownames(statsTable) <- c('peaks','promoters','chrM','blacklist')
}

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Genome-wide structures")

if(file.exists(WG_dedup_readCtByCell_RDS)){
    WG_dedup_readCtByCell <- readRDS(WG_dedup_readCtByCell_RDS)
    
    OC_dedup_fragment_GR <- make_GR(otherChr_fragment_file,c('chr','start','end','cell','dupCt'))
} else {
    AS_dedup_fragment_GR <- make_GR(fragment_file,c('chr','start','end','cell','dupCt'))
    AS_dedup_readCtByCell <- readCt_from_frag_GR(AS_dedup_fragment_GR)
    
    OC_dedup_fragment_GR <- make_GR(otherChr_fragment_file,c('chr','start','end','cell','dupCt'))
    WG_dedup_fragment_GR <- c(AS_dedup_fragment_GR,OC_dedup_fragment_GR)
    
    WG_dedup_readCtByCell <- readCt_from_frag_GR(WG_dedup_fragment_GR)
    WG_dedup_readCtByCell <- WG_dedup_readCtByCell[order(names(WG_dedup_readCtByCell))]
    
    saveRDS(WG_dedup_readCtByCell,WG_dedup_readCtByCell_RDS)
}

if(file.exists(WG_noDedup_readCtByCell_RDS)){
    WG_noDedup_readCtByCell <- readRDS(WG_noDedup_readCtByCell_RDS)
    
} else {
    WG_noDedup_readCtByCell <- make_CB_vec(WG_CB_ct_file)
    WG_noDedup_readCtByCell <- WG_noDedup_readCtByCell[order(names(WG_noDedup_readCtByCell))]
    
    saveRDS(WG_noDedup_readCtByCell,WG_noDedup_readCtByCell_RDS)
}


print_time(paste("Taking out cells with less than",threshold,"reads"))
cat(paste(sep="",prefix,": Original number of cells in WG dedup: ",length(WG_dedup_readCtByCell),"\n"))
cat(paste(sep="",prefix,": Min number of cells in WG dedup: ",min(WG_dedup_readCtByCell),"\n"))
WG_dedup_readCtByCell <- WG_dedup_readCtByCell[WG_dedup_readCtByCell>=threshold]
cat(paste(sep="",prefix,": New number of cells in WG dedup: ",length(WG_dedup_readCtByCell),"\n"))
cat(paste(sep="",prefix,": Min number of cells in WG dedup: ",min(WG_dedup_readCtByCell),"\n"))

cat(paste(sep="",prefix,": Original number of cells in WG NO dedup: ",length(WG_noDedup_readCtByCell),"\n"))
cat(paste(sep="",prefix,": Min number of cells in WG dedup: ",min(WG_noDedup_readCtByCell),"\n"))
WG_noDedup_readCtByCell <- WG_noDedup_readCtByCell[WG_noDedup_readCtByCell>=threshold]
cat(paste(sep="",prefix,": New number of cells in WG NO dedup: ",length(WG_noDedup_readCtByCell),"\n"))
cat(paste(sep="",prefix,": Min number of cells in WG dedup: ",min(WG_noDedup_readCtByCell),"\n"))


print_time("Mitochondrial")

chrM_dedup_readCtByCell <- readCt_from_frag_GR(OC_dedup_fragment_GR[seqnames(OC_dedup_fragment_GR)=='chrM'])
criterion='chrM'
doEverything(chrM_dedup_readCtByCell,WG_dedup_readCtByCell,"dedup","WG","reads",threshold,"chrM",statsTable[criterion,'cutoff'],outPrefix, maxGenome=statsTable[criterion,'maxGenome'], maxCrit=statsTable[criterion,'maxCrit'], above=FALSE)


print_time("Blacklist")

blacklist_noDedup_readCtByCell <- make_CB_vec(blacklist_CB_ct_file)
criterion='blacklist'
doEverything(blacklist_noDedup_readCtByCell,WG_noDedup_readCtByCell,"NOdedup","WG","reads",threshold,"blacklist",statsTable[criterion,'cutoff'],outPrefix, maxGenome=statsTable[criterion,'maxGenome'], maxCrit=statsTable[criterion,'maxCrit'], above=FALSE)


print_time("Peaks")

peaks_dedup_readCtByCell <- make_CB_vec(peaks_CB_ct_file)
criterion='peaks'
doEverything(peaks_dedup_readCtByCell,WG_dedup_readCtByCell,"dedup","WG","reads",threshold,"peaks",statsTable[criterion,'cutoff'],outPrefix, maxGenome=statsTable[criterion,'maxGenome'], maxCrit=statsTable[criterion,'maxCrit'])


print_time("Promoters")

promoters_dedup_readCtByCell <- make_CB_vec(promoters_CB_ct_file)
criterion='promoters'
doEverything(promoters_dedup_readCtByCell,WG_dedup_readCtByCell,"dedup","WG","reads",threshold,"promoters",statsTable[criterion,'cutoff'],outPrefix, maxGenome=statsTable[criterion,'maxGenome'], maxCrit=statsTable[criterion,'maxCrit'])



